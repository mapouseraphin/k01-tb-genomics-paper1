# =============================================================================
# 02b_calibration_bootstrap.R
#
# Stage 2 of the two-stage data-driven calibration: bootstrap-stabilized
# primary selection.
#
# WORKFLOW
# --------
#   Stage 1 (02_calibrate_hybrid_pareto.R): Pareto selection on observed data
#                                           writes thr_primary_observed.rds
#                                           and thr_primary.rds (initially equal).
#   Stage 2 (this script):                  Pareto selection on B bootstrap
#                                           resamples, take the modal cell as
#                                           the stabilized primary, OVERWRITE
#                                           thr_primary.rds with the modal.
#                                           Save thr_primary_bootstrap.rds.
#
# RATIONALE
# ---------
#   Lex-selection bootstrap revealed AD1 instability (canonical AD1=1 in 48.5%
#   of reps, AD1=6 in 50.4%, joint canonical 19.2%). Switching to Pareto with
#   closest-to-ideal removes the lex-induced low-AD1 bias. Bootstrapping the
#   Pareto procedure provides an honest measure of selection stability under
#   the new rule. Taking the modal cell across reps as the final primary
#   makes the selection completely data-driven: not "the cell the procedure
#   picks on this dataset" but "the cell the procedure picks consistently
#   under resampling."
#
# BOOTSTRAP UNIT
# --------------
#   M0 within-visit replicate pair. Each rep resamples the full set of pairs
#   with replacement (size = n_pairs), recomputes per-grid-point metrics with
#   multiplicity preserved, applies the same guards as Stage 1, and runs the
#   same Pareto + closest-to-ideal selection.
#
# STABILIZATION RULE
# ------------------
#   Final primary = modal cell across OK reps (cells from failed reps are
#   excluded). Ties broken by:
#     (1) higher mean Jaccard across reps where that cell was selected
#     (2) higher mean Confirm
#     (3) higher stringency rank (DP, AD1, MAF)
#
# OUTPUTS (PATHS$bootstrap)
# -------------------------
#   bootstrap_picks.rds         # B-row tibble of primary cells per rep
#   bootstrap_stability.rds     # summary list (joint and marginal stability)
#   bootstrap_summary_table.csv # human-readable summary for supplement
#   bootstrap_log.txt           # provenance log
#
# OUTPUTS (PATHS$calibration; OVERWRITES from Stage 1)
# ----------------------------------------------------
#   thr_primary.rds             # bootstrap-stabilized (modal) primary
#   thr_primary_bootstrap.rds   # backup of the stabilized primary
#   (thr_primary_observed.rds remains untouched -- Stage 1 backup)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr); library(readr); library(purrr)
})

banner()
message("\n[02b_calibration_bootstrap] Stage 2: Pareto bootstrap stabilization")

dir.create(PATHS$bootstrap, recursive = TRUE, showWarnings = FALSE)

LOG <- character()
log_msg <- function(...) {
  args <- list(...)
  s <- if (length(args) == 1) as.character(args[[1]]) else do.call(sprintf, args)
  message(s); LOG <<- c(LOG, s); invisible(s)
}

set.seed(SEED)

# ---- Load -------------------------------------------------------------------
cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))

req_cols <- c("DP_min", "AD1_min", "MAF_min",
              "sampleA", "sampleB",
              "nA", "nB", "n_union", "n_inter",
              "jaccard", "overlap", "confirm_both", "maf_cor")
miss <- setdiff(req_cols, names(cal_pairs))
if (length(miss)) stop("cal_pairs missing required columns: ",
                       paste(miss, collapse = ", "))

# Stage 1 observed primary, for comparison against the bootstrap modal
obs_path <- file.path(PATHS$calibration, "thr_primary_observed.rds")
if (!file.exists(obs_path)) {
  stop("thr_primary_observed.rds not found. Run 02_calibrate_hybrid_pareto.R ",
       "before 02b.")
}
observed_primary <- readRDS(obs_path)
log_msg("Observed-data Pareto primary (Stage 1): DP=%g, AD1=%g, MAF=%g",
        observed_primary$DP_min, observed_primary$AD1_min, observed_primary$MAF_min)

# ---- Pair index -------------------------------------------------------------
all_pairs <- cal_pairs %>%
  distinct(sampleA, sampleB) %>%
  mutate(pair_id = row_number())
n_pairs <- nrow(all_pairs)
log_msg("Total unique M0 within-visit pairs: %d", n_pairs)
log_msg("Bootstrap reps: B = %d", N_BOOTSTRAP_LEX_SEL)

cal_pairs <- cal_pairs %>%
  inner_join(all_pairs, by = c("sampleA", "sampleB"))

# Pre-split cal_pairs by pair_id for fast subsetting
cal_pairs_by_pair <- split(cal_pairs, cal_pairs$pair_id)

# ---- Pareto selection on a resampled cal_pairs subset ----------------------
# Mirrors 02_calibrate_hybrid_pareto.R but compressed to return one cell.
# Returns NULL if no grid points pass guards or no Pareto frontier exists.
is_pareto <- function(J, C) {
  ok <- !is.na(J) & !is.na(C)
  keep <- rep(FALSE, length(J))
  for (i in which(ok)) {
    dom <- (J >= J[i] & C >= C[i]) & (J > J[i] | C > C[i]) & ok
    if (!any(dom)) keep[i] <- TRUE
  }
  keep
}

aggregate_and_pick <- function(cp_subset) {
  cal_grid <- cp_subset %>%
    group_by(DP_min, AD1_min, MAF_min) %>%
    summarise(
      pairs_total      = n(),
      pairs_both0      = sum(nA == 0 & nB == 0),
      pairs_evaluable  = sum(nA > 0 & nB > 0),
      prop_both0       = pairs_both0 / pairs_total,
      jacc_med_eval    = median(jaccard     [n_union > 0 & nA > 0 & nB > 0], na.rm = TRUE),
      conf_med_eval    = median(confirm_both[n_union > 0 & nA > 0 & nB > 0], na.rm = TRUE),
      r2_med_eval      = median((maf_cor^2) [n_inter >= 2 & nA > 0 & nB > 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(pass_guards = (prop_both0 <= MAX_PROP_BOTH0) &
                         (pairs_evaluable >= MIN_EVAL_PAIRS))

  passers <- cal_grid %>% filter(pass_guards)
  if (nrow(passers) == 0) return(NULL)

  passers$is_frontier <- is_pareto(passers$jacc_med_eval,
                                    passers$conf_med_eval)
  frontier <- passers %>% filter(is_frontier)
  if (nrow(frontier) == 0) return(NULL)

  frontier %>%
    mutate(dist_to_ideal = sqrt((1 - jacc_med_eval)^2 + (1 - conf_med_eval)^2),
           S = min_rank(DP_min) + min_rank(AD1_min) + min_rank(MAF_min)) %>%
    arrange(dist_to_ideal,
            desc(r2_med_eval),
            desc(S),
            desc(DP_min), desc(AD1_min), desc(MAF_min)) %>%
    slice(1) %>%
    select(DP_min, AD1_min, MAF_min,
           jacc_med_eval, conf_med_eval, dist_to_ideal,
           pairs_evaluable)
}

# ---- Bootstrap loop ---------------------------------------------------------
log_msg("\nRunning Pareto bootstrap...")
t0 <- Sys.time()

picks <- vector("list", N_BOOTSTRAP_LEX_SEL)
n_failed <- 0L

report_at <- unique(c(1L, 10L, 100L,
                      seq(N_BOOTSTRAP_LEX_SEL %/% 10, N_BOOTSTRAP_LEX_SEL,
                          by = max(1L, N_BOOTSTRAP_LEX_SEL %/% 10))))

for (b in seq_len(N_BOOTSTRAP_LEX_SEL)) {
  resampled_ids <- sample(seq_len(n_pairs), size = n_pairs, replace = TRUE)
  cp_b <- bind_rows(cal_pairs_by_pair[as.character(resampled_ids)])
  pick <- aggregate_and_pick(cp_b)
  if (is.null(pick)) {
    n_failed <- n_failed + 1L
    picks[[b]] <- tibble(rep = b, DP_min = NA_real_, AD1_min = NA_real_,
                         MAF_min = NA_real_, jacc_med_eval = NA_real_,
                         conf_med_eval = NA_real_, dist_to_ideal = NA_real_,
                         pairs_evaluable = NA_integer_, failed = TRUE)
  } else {
    picks[[b]] <- pick %>% mutate(rep = b, failed = FALSE) %>%
      select(rep, everything())
  }
  if (b %in% report_at) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    log_msg("  rep %d / %d  [%.1fs elapsed, %d failed]",
            b, N_BOOTSTRAP_LEX_SEL, elapsed, n_failed)
  }
}

bootstrap_picks <- bind_rows(picks)
log_msg("\nBootstrap complete in %.1f seconds. Failed reps: %d / %d (%.2f%%)",
        as.numeric(difftime(Sys.time(), t0, units = "secs")),
        n_failed, N_BOOTSTRAP_LEX_SEL,
        100 * n_failed / N_BOOTSTRAP_LEX_SEL)

saveRDS(bootstrap_picks, file.path(PATHS$bootstrap, "bootstrap_picks.rds"))

# ---- Stabilized primary = modal cell across OK reps ------------------------
ok <- bootstrap_picks %>% filter(!failed)
n_ok <- nrow(ok)

if (n_ok == 0) {
  stop("All bootstrap reps failed. Cannot stabilize. ",
       "Inspect cal_pairs and guards.")
}

# Cell frequency table with mean (J, C) over reps where each cell was picked
cell_summary <- ok %>%
  group_by(DP_min, AD1_min, MAF_min) %>%
  summarise(n_reps        = n(),
            mean_J_pick   = mean(jacc_med_eval, na.rm = TRUE),
            mean_C_pick   = mean(conf_med_eval, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(prop = n_reps / n_ok,
         S    = min_rank(DP_min) + min_rank(AD1_min) + min_rank(MAF_min)) %>%
  arrange(desc(n_reps),
          desc(mean_J_pick),
          desc(mean_C_pick),
          desc(S))

modal_cell <- cell_summary %>% slice(1)

log_msg("\nBootstrap-stabilized primary (modal across reps):")
log_msg("  DP=%g, AD1=%g, MAF=%g   (selected in %d / %d reps = %.1f%%)",
        modal_cell$DP_min, modal_cell$AD1_min, modal_cell$MAF_min,
        modal_cell$n_reps, n_ok, 100 * modal_cell$prop)

prop_modal <- modal_cell$prop[1]

prop_observed <- ok %>%
  filter(DP_min  == observed_primary$DP_min,
         AD1_min == observed_primary$AD1_min,
         MAF_min == observed_primary$MAF_min) %>%
  nrow() / n_ok

log_msg("\nJoint stability comparison:")
log_msg("  Observed-data primary (DP=%g, AD1=%g, MAF=%g): %.1f%% of reps",
        observed_primary$DP_min, observed_primary$AD1_min, observed_primary$MAF_min,
        100 * prop_observed)
log_msg("  Modal cell             (DP=%g, AD1=%g, MAF=%g): %.1f%% of reps",
        modal_cell$DP_min, modal_cell$AD1_min, modal_cell$MAF_min,
        100 * prop_modal)

modal_eq_observed <- (modal_cell$DP_min  == observed_primary$DP_min) &
                     (modal_cell$AD1_min == observed_primary$AD1_min) &
                     (modal_cell$MAF_min == observed_primary$MAF_min)
if (modal_eq_observed) {
  log_msg("  [CONSISTENT] Modal cell matches observed-data primary; no override.")
} else {
  log_msg("  [DISCREPANT] Modal cell differs from observed-data primary;")
  log_msg("               stabilized primary will OVERRIDE observed-data primary.")
}

# Marginal stability per axis (relative to modal)
marginal <- tibble(
  axis = c("DP_min", "AD1_min", "MAF_min"),
  modal_value = c(modal_cell$DP_min, modal_cell$AD1_min, modal_cell$MAF_min),
  prop_at_modal = c(
    mean(ok$DP_min  == modal_cell$DP_min),
    mean(ok$AD1_min == modal_cell$AD1_min),
    mean(ok$MAF_min == modal_cell$MAF_min)
  )
)

axis_dist <- bind_rows(
  ok %>% count(value = DP_min,  name = "n") %>% mutate(axis = "DP_min"),
  ok %>% count(value = AD1_min, name = "n") %>% mutate(axis = "AD1_min"),
  ok %>% count(value = MAF_min, name = "n") %>% mutate(axis = "MAF_min")
) %>% mutate(prop = n / n_ok) %>% select(axis, value, n, prop)

stability <- list(
  selection_method     = "hybrid_pareto_closest_to_ideal",
  n_reps               = N_BOOTSTRAP_LEX_SEL,
  n_failed             = n_failed,
  n_ok                 = n_ok,
  observed_primary     = observed_primary,
  modal_cell           = list(DP_min = modal_cell$DP_min,
                              AD1_min = modal_cell$AD1_min,
                              MAF_min = modal_cell$MAF_min,
                              MAF_max = MAF_MAX),
  prop_observed_primary = prop_observed,
  prop_modal_joint     = prop_modal,
  modal_eq_observed    = modal_eq_observed,
  cell_summary         = cell_summary,
  marginal_stability   = marginal,
  axis_distribution    = axis_dist
)
saveRDS(stability, file.path(PATHS$bootstrap, "bootstrap_stability.rds"))

# ---- Write stabilized primary as the new thr_primary.rds -------------------
thr_primary_bootstrap <- list(
  DP_min  = modal_cell$DP_min,
  AD1_min = modal_cell$AD1_min,
  MAF_min = modal_cell$MAF_min,
  MAF_max = MAF_MAX
)

saveRDS(thr_primary_bootstrap,
        file.path(PATHS$calibration, "thr_primary_bootstrap.rds"))
saveRDS(thr_primary_bootstrap,
        file.path(PATHS$calibration, "thr_primary.rds"))

log_msg("\n[02b] thr_primary.rds OVERWRITTEN with bootstrap-stabilized cell.")
log_msg("  Backups available in PATHS$calibration:")
log_msg("    thr_primary_observed.rds  = Stage 1 observed-data Pareto primary")
log_msg("    thr_primary_bootstrap.rds = Stage 2 bootstrap-stabilized primary (now active)")

# ---- Sync consolidated thresholds.rds with stabilized primary --------------
# 02_calibrate_hybrid_pareto wrote thresholds.rds with the Stage 1 primary.
# After bootstrap stabilization, the canonical primary is the modal cell.
# Update thresholds.rds so that downstream consumers (03_apply_thresholds.R,
# make_paper1_tables_1to3.R) see the active Stage 2 primary.
thresholds_path <- file.path(PATHS$calibration, "thresholds.rds")
if (file.exists(thresholds_path)) {
  thresholds <- readRDS(thresholds_path)
  thresholds$primary <- thr_primary_bootstrap
  saveRDS(thresholds, thresholds_path)
  log_msg("[02b] thresholds.rds$primary synced to bootstrap-stabilized cell.")
} else {
  warning("[02b] thresholds.rds not found; skipping sync.")
}

# ---- Human-readable summary table ------------------------------------------
summary_tbl <- tibble(
  metric = c(
    "Selection method",
    "Bootstrap reps (B)",
    "Failed reps (no Pareto frontier)",
    "Reps yielding a primary cell (n_ok)",
    "Observed-data primary (Stage 1)",
    "P(rep selects observed-data primary)",
    "Modal cell across reps (Stage 2 stabilized)",
    "P(modal cell)",
    "Modal == observed?",
    "Marginal P(DP_min = modal)",
    "Marginal P(AD1_min = modal)",
    "Marginal P(MAF_min = modal)"
  ),
  value = c(
    "hybrid Pareto, closest-to-ideal",
    sprintf("%d", N_BOOTSTRAP_LEX_SEL),
    sprintf("%d (%.2f%%)", n_failed, 100 * n_failed / N_BOOTSTRAP_LEX_SEL),
    sprintf("%d", n_ok),
    sprintf("(DP=%g, AD1=%g, MAF=%g)",
            observed_primary$DP_min, observed_primary$AD1_min, observed_primary$MAF_min),
    sprintf("%.1f%%", 100 * prop_observed),
    sprintf("(DP=%g, AD1=%g, MAF=%g)",
            modal_cell$DP_min, modal_cell$AD1_min, modal_cell$MAF_min),
    sprintf("%.1f%%", 100 * prop_modal),
    if (modal_eq_observed) "yes" else "NO (override applied)",
    sprintf("%.1f%%", 100 * marginal$prop_at_modal[marginal$axis == "DP_min"]),
    sprintf("%.1f%%", 100 * marginal$prop_at_modal[marginal$axis == "AD1_min"]),
    sprintf("%.1f%%", 100 * marginal$prop_at_modal[marginal$axis == "MAF_min"])
  )
)

write_csv(summary_tbl, file.path(PATHS$bootstrap, "bootstrap_summary_table.csv"))

log_msg("\n=== Bootstrap stability summary ===")
log_msg(paste(capture.output(print(as.data.frame(summary_tbl))), collapse = "\n"))
log_msg("\nTop-5 cells by frequency across reps:")
log_msg(paste(capture.output(print(as.data.frame(head(cell_summary, 5)))),
              collapse = "\n"))

writeLines(LOG, file.path(PATHS$bootstrap, "bootstrap_log.txt"))

message("\n[02b_calibration_bootstrap] Done. Outputs in ", PATHS$bootstrap, "/")
message("  Stabilized primary written to ", PATHS$calibration, "/thr_primary.rds")
