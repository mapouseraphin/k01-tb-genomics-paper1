# =============================================================================
# build_lca_dataset.R
#
# Builds the variant-call-level LCA-ready dataset for Paper 1.
#
# Each row = one (pair_id x position x alt) entry. Columns include:
#   W1, W2          - rule-pass indicators in the two replicates of the pair
#   MAF1, MAF2      - observed minor allele frequencies (NA if not called)
#   DP1, DP2        - observed depths (NA if not called)
#   AD1_1, AD1_2    - observed alt-allele depths (NA if not called)
#   MAF_max_obs     - max(MAF1, MAF2, na.rm)  -- stratification variable
#   DP_min_obs      - min(DP1, DP2, na.rm)
#   maf_bin         - factor: 4 bins under MAF_BINNING="default";
#                            3 bins under MAF_BINNING="merged_low"
#   region          - factor "PE_PPE" / "non_PE_PPE"
#   HIV             - factor "HIV-" / "HIV+" (NA if missing)
#   pair_type       - factor over slot-pair types (EM_S0.1, EM_S0.2, S0.1_S0.2)
#   pair_id         - identifier for the M0 within-visit pair
#   patientId       - patient identifier
#
# Locked design decisions (see session log):
#   - Universe (canonical, UNIVERSE="C"): union of (CHROM,POS,REF,ALT) keys
#       with at least one PASS GATK call across the pair.
#   - Universe (sensitivity, UNIVERSE="A"): union of keys with at least one
#       GATK call (PASS or fail) across the pair. Broader universe; sens
#       under A characterizes joint GATK + rule behavior.
#   - SNV-only: indels dropped (under both UNIVERSE values).
#   - PE/PPE: retained in universe, used as stratifier (NOT excluded). PE/PPE
#       exclusion is operated upstream by the apply_thresholds step when
#       DROP_PPE_APPLY=TRUE; this script honors whatever the variant table
#       provides.
#   - Asymmetric calls: replicate without a call -> W=0, NA covariates.
#   - Pair-type: kept as a column; the primary LCA model collapses across
#       pair-type. Pair-type is a covariate in fit_lca_models.R if needed.
#
# Inputs:
#   PATHS$variants/gh_variants.rds
#   PATHS$meta/gh_pairs.rds
#   PATHS$meta/gh_meta.rds
#   PATHS$calibration/thr_primary.rds   # list(DP_min, AD1_min, MAF_min, MAF_max)
#
# Outputs:
#   PATHS$lca/lca_dataset.rds          # main analytic dataset
#   PATHS$lca/lca_dataset.csv          # human-readable copy
#   PATHS$lca/lca_build_log.txt        # provenance log
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr); library(tibble)
})

banner()

# ---- Validate config --------------------------------------------------------
stopifnot(UNIVERSE %in% c("C", "A"))
stopifnot(MAF_BINNING %in% c("default", "merged_low"))

# ---- Paths ------------------------------------------------------------------
# Use PATHS$lca (defined in 00_pipeline_config.R). This now varies with
# UNIVERSE and MAF_BINNING via TAG$lca_tag.
LCA_PATH <- PATHS$lca
dir.create(LCA_PATH, recursive = TRUE, showWarnings = FALSE)

LOG <- character()
log_msg <- function(...) {
  s <- sprintf(...)
  message(s); LOG <<- c(LOG, s); invisible(s)
}

log_msg("LCA spec: UNIVERSE=%s, MAF_BINNING=%s, lca_tag=%s",
        UNIVERSE, MAF_BINNING, TAG$lca_tag)
log_msg("LCA outputs -> %s/", LCA_PATH)

# ---- Load inputs ------------------------------------------------------------
gh_v    <- readRDS(file.path(PATHS$variants,    "gh_variants.rds"))
gh_pairs<- readRDS(file.path(PATHS$meta,        "gh_pairs.rds"))
gh_meta <- readRDS(file.path(PATHS$meta,        "gh_meta.rds"))
thr     <- readRDS(file.path(PATHS$calibration, "thr_primary.rds"))

stopifnot(!is.null(thr), all(c("DP_min","AD1_min","MAF_min","MAF_max") %in% names(thr)))

log_msg("Inputs loaded:")
log_msg("  gh_variants:  %s rows (%d samples)",
        format(nrow(gh_v), big.mark = ","), n_distinct(gh_v$Sample))
log_msg("  gh_pairs:     %d M0 within-visit pairs", nrow(gh_pairs))
log_msg("  gh_meta:      %d samples (depth_med >= %d)", nrow(gh_meta), MIN_COV_INCLUDE)
log_msg("  thr_primary:  DP>=%g, AD1>=%g, MAF in [%g, %g]",
        thr$DP_min, thr$AD1_min, thr$MAF_min, thr$MAF_max)

# ---- Universe construction --------------------------------------------------
# UNIVERSE = "C" (canonical): enforce FILTER == PASS pre-screen. This is
#                              defensive: gh_variants.rds is typically PASS-only
#                              upstream, but we re-enforce here so the universe
#                              is unambiguous.
# UNIVERSE = "A" (sensitivity): skip the FILTER == PASS gate; all GATK calls
#                                (PASS or fail) enter the universe. Sens/spec
#                                under A characterize joint GATK + rule
#                                behavior; under C they characterize rule-
#                                only behavior given GATK pre-screening.
n0 <- nrow(gh_v)
if (UNIVERSE == "C") {
  if ("FILTER" %in% names(gh_v)) {
    gh_v <- gh_v %>% filter(FILTER == "PASS")
    log_msg("Universe step 1 (UNIVERSE=C) -- FILTER==PASS enforced: %s -> %s",
            format(n0, big.mark = ","), format(nrow(gh_v), big.mark = ","))
  } else {
    log_msg("Universe step 1 (UNIVERSE=C) -- no FILTER column; assuming all PASS upstream.")
  }
} else if (UNIVERSE == "A") {
  filter_summary <- if ("FILTER" %in% names(gh_v)) {
    paste(capture.output(print(table(gh_v$FILTER, useNA = "ifany"))), collapse = " | ")
  } else "(no FILTER column)"
  log_msg("Universe step 1 (UNIVERSE=A) -- FILTER gate SKIPPED. Pre-gate FILTER table: %s",
          filter_summary)
  log_msg("  Retained all %s rows for universe construction.", format(nrow(gh_v), big.mark = ","))
}

n0 <- nrow(gh_v)
gh_v <- filter_snv_only(gh_v)
log_msg("Universe step 2 -- SNV-only: %s -> %s (%.1f%% retained)",
        format(n0, big.mark = ","), format(nrow(gh_v), big.mark = ","),
        100 * nrow(gh_v) / n0)

gh_v <- ensure_maf(gh_v)

# PE/PPE flag (unchanged from canonical)
get_ppe_flag <- function(v) {
  if ("ppe_flag" %in% names(v)) {
    pp <- v$ppe_flag
    if (is.logical(pp))   return(as.integer(pp))
    if (is.numeric(pp))   return(as.integer(pp != 0))
    if (is.character(pp) || is.factor(pp)) {
      drop_levels <- c("PPE","pe/ppe","PE/PPE","TRUE","true","1")
      return(as.integer(as.character(pp) %in% drop_levels))
    }
  }
  if ("PPE" %in% names(v)) {
    pp <- v$PPE
    if (is.logical(pp))   return(as.integer(pp))
    if (is.numeric(pp))   return(as.integer(pp != 0))
    if (is.character(pp)) return(as.integer(pp %in% c("1","TRUE","true","PPE","pe/ppe")))
  }
  rep(NA_integer_, nrow(v))
}
gh_v$ppe_int <- get_ppe_flag(gh_v)
n_ppe_na <- sum(is.na(gh_v$ppe_int))
if (n_ppe_na > 0) log_msg("WARNING: %d variant calls lack a PE/PPE annotation.", n_ppe_na)

# ---- Helper: per-sample call lookup -----------------------------------------
samples_in_pairs <- unique(c(gh_pairs$sampleA, gh_pairs$sampleB))
v_by_sample <- gh_v %>%
  filter(Sample %in% samples_in_pairs) %>%
  mutate(site_key = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
  group_by(Sample) %>% group_split()
names(v_by_sample) <- vapply(v_by_sample, function(d) d$Sample[1], character(1))

empty_v <- gh_v[0, ] %>% mutate(site_key = character())
get_v <- function(s) if (!is.null(v_by_sample[[s]])) v_by_sample[[s]] else empty_v

# ---- Apply calibrated rule (vectorised) -------------------------------------
# Returns logical vector indicating pass/fail per row of the per-sample table.
# Note: rule_pass operates only on quality thresholds (DP/AD1/MAF). It does
# NOT re-check FILTER == PASS; under UNIVERSE=A, non-PASS calls can therefore
# still satisfy the rule if their DP/AD1/MAF meet thresholds. This is the
# intended semantic of UNIVERSE=A as a sensitivity.
rule_pass <- function(v) {
  v$DP  >= thr$DP_min  &
  v$AD1 >= thr$AD1_min &
  v$MAF >= thr$MAF_min &
  v$MAF <= thr$MAF_max
}

# ---- Pair-type derivation ---------------------------------------------------
slot_pair_type <- function(slot_a, slot_b) {
  s <- sort(c(slot_a, slot_b))
  paste(s[1], s[2], sep = "_")
}

# ---- Build per-pair, per-site rows ------------------------------------------
log_msg("\nIterating %d pairs to build call-level dataset...", nrow(gh_pairs))

build_pair_rows <- function(pair_row, pair_id) {
  sA <- pair_row$sampleA; sB <- pair_row$sampleB
  vA <- get_v(sA); vB <- get_v(sB)

  keys <- union(vA$site_key, vB$site_key)
  if (length(keys) == 0) return(NULL)

  iA <- match(keys, vA$site_key)
  iB <- match(keys, vB$site_key)

  pull_col <- function(v, idx, col, type = "num") {
    out <- if (type == "num") rep(NA_real_, length(idx)) else
                              rep(NA_integer_, length(idx))
    has <- !is.na(idx)
    if (any(has)) out[has] <- v[[col]][idx[has]]
    out
  }

  W1 <- integer(length(keys)); W2 <- integer(length(keys))
  for (i in seq_along(keys)) {
    if (!is.na(iA[i])) {
      r <- vA[iA[i], , drop = FALSE]
      W1[i] <- as.integer(rule_pass(r))
    } else W1[i] <- 0L
    if (!is.na(iB[i])) {
      r <- vB[iB[i], , drop = FALSE]
      W2[i] <- as.integer(rule_pass(r))
    } else W2[i] <- 0L
  }

  parse_key <- function(k) {
    parts <- str_split_fixed(k, ":", 4)
    list(CHROM = parts[,1],
         POS   = suppressWarnings(as.integer(parts[,2])),
         REF   = parts[,3],
         ALT   = parts[,4])
  }
  pk <- parse_key(keys)

  ppe <- rep(NA_integer_, length(keys))
  for (i in seq_along(keys)) {
    if (!is.na(iA[i]))      ppe[i] <- vA$ppe_int[iA[i]]
    else if (!is.na(iB[i])) ppe[i] <- vB$ppe_int[iB[i]]
  }

  tibble(
    pair_id     = pair_id,
    sampleA     = sA,
    sampleB     = sB,
    CHROM       = pk$CHROM,
    POS         = pk$POS,
    REF         = pk$REF,
    ALT         = pk$ALT,
    site_key    = keys,
    W1          = W1,
    W2          = W2,
    MAF1        = pull_col(vA, iA, "MAF"),
    MAF2        = pull_col(vB, iB, "MAF"),
    DP1         = pull_col(vA, iA, "DP",  type = "int"),
    DP2         = pull_col(vB, iB, "DP",  type = "int"),
    AD1_1       = pull_col(vA, iA, "AD1", type = "int"),
    AD1_2       = pull_col(vB, iB, "AD1", type = "int"),
    ppe_int     = ppe
  )
}

rows_list <- vector("list", nrow(gh_pairs))
for (p in seq_len(nrow(gh_pairs))) {
  rows_list[[p]] <- build_pair_rows(gh_pairs[p, ], pair_id = p)
}
lca <- bind_rows(rows_list)

log_msg("Pre-annotation row count: %s", format(nrow(lca), big.mark = ","))

# ---- Annotate stratification variables --------------------------------------
# MAF binning branches on MAF_BINNING. The "default" 4-bin scheme is the
# canonical Paper 1 specification. The "merged_low" scheme collapses the two
# thinnest low-MAF bins as a sensitivity to test whether stratum-specific
# posteriors are bin-resolution-driven.
if (MAF_BINNING == "default") {
  maf_breaks <- c(0.02, 0.05, 0.10, 0.25, thr$MAF_max + 1e-9)
  maf_labels <- c("[0.02,0.05)", "[0.05,0.10)", "[0.10,0.25)",
                  paste0("[0.25,", thr$MAF_max, "]"))
} else if (MAF_BINNING == "merged_low") {
  maf_breaks <- c(0.02, 0.10, 0.25, thr$MAF_max + 1e-9)
  maf_labels <- c("[0.02,0.10)", "[0.10,0.25)",
                  paste0("[0.25,", thr$MAF_max, "]"))
}
log_msg("MAF binning (%s): %d bins -- %s", MAF_BINNING, length(maf_labels),
        paste(maf_labels, collapse = " | "))

lca <- lca %>%
  mutate(
    MAF_max_obs = pmax(MAF1, MAF2, na.rm = TRUE),
    DP_min_obs  = pmin(DP1,  DP2,  na.rm = TRUE),
    region      = factor(if_else(ppe_int == 1L, "PE_PPE", "non_PE_PPE"),
                         levels = c("non_PE_PPE", "PE_PPE")),
    maf_bin     = cut(MAF_max_obs, breaks = maf_breaks,
                      labels = maf_labels, right = FALSE,
                      include.lowest = TRUE)
  )

# Slot/pair-type and patient/HIV joins (unchanged from canonical)
lca <- lca %>%
  mutate(slot_a = parse_slot(sampleA),
         slot_b = parse_slot(sampleB),
         pair_type = mapply(slot_pair_type, slot_a, slot_b)) %>%
  left_join(gh_meta %>% select(Sample, patientId, HIV) %>%
              rename(sampleA = Sample, patientId_a = patientId, HIV_a = HIV),
            by = "sampleA") %>%
  left_join(gh_meta %>% select(Sample, HIV) %>%
              rename(sampleB = Sample, HIV_b = HIV),
            by = "sampleB") %>%
  mutate(
    patientId = patientId_a,
    HIV = case_when(
      !is.na(HIV_a) & !is.na(HIV_b) & HIV_a != HIV_b ~ NA_character_,
      !is.na(HIV_a) ~ as.character(HIV_a),
      !is.na(HIV_b) ~ as.character(HIV_b),
      TRUE ~ NA_character_
    ),
    HIV = factor(HIV, levels = c("HIV-", "HIV+"))
  ) %>%
  select(-patientId_a, -HIV_a, -HIV_b)

# Sanity checks
n_w_na <- sum(is.na(lca$W1) | is.na(lca$W2))
if (n_w_na > 0) stop("build_lca_dataset: ", n_w_na, " rows with NA in W1/W2")

if (nrow(lca) < 100) {
  log_msg("WARNING: LCA dataset has only %d rows -- talk to PI before proceeding.",
          nrow(lca))
}

log_msg("\nFinal dataset: %s rows across %d pairs, %d patients",
        format(nrow(lca), big.mark = ","),
        n_distinct(lca$pair_id),
        n_distinct(lca$patientId))

# ---- Cell tabulation snapshot (for the build log) ---------------------------
cells <- lca %>% count(W1, W2) %>% mutate(prop = round(n / sum(n), 3))
log_msg("\nJoint (W1, W2) cells (overall):")
log_msg(paste(capture.output(print(as.data.frame(cells))), collapse = "\n"))

# ---- Save -------------------------------------------------------------------
saveRDS(lca, file.path(LCA_PATH, "lca_dataset.rds"))
write_csv(lca, file.path(LCA_PATH, "lca_dataset.csv"))
writeLines(LOG, file.path(LCA_PATH, "lca_build_log.txt"))

message("\n[build_lca_dataset] Done. Artifacts written to ", LCA_PATH, "/")
