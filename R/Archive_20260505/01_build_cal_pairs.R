# =============================================================================
# 01_build_cal_pairs.R
#
# Builds `cal_pairs` -- the per-pair concordance metric table required by
# 02_calibrate_lexicographic.R.
#
# Canonical specification:
#   - SNV_ONLY = TRUE: indels excluded before grid evaluation
#   - DROP_PPE_CAL: PE/PPE handling for calibration step (canonical: FALSE)
#   - require FILTER == "PASS"
#
# For each (DP, AD1, MAF) grid point and each pair (sampleA, sampleB):
#   nA, nB                    -- variant counts at threshold
#   n_union, n_inter          -- set algebra
#   jaccard, overlap          -- concordance metrics
#   confirm_both              -- fraction of union sites in both
#   maf_cor                   -- MAF correlation across union
#
# Inputs:
#   PATHS$variants/gh_variants.rds  (built in 00_prep_metadata.R)
#   PATHS$meta/gh_pairs.rds         (M0-only canonical)
#
# Outputs (to PATHS$calibration):
#   cal_pairs.rds, cal_pairs.csv
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

# ---- Load inputs ------------------------------------------------------------
gh_v_raw <- readRDS(file.path(PATHS$variants, "gh_variants.rds"))
gh_pairs <- readRDS(file.path(PATHS$meta,     "gh_pairs.rds"))

stopifnot(nrow(gh_v_raw) > 0)
stopifnot(nrow(gh_pairs) > 0)

req_v <- c("Sample", "CHROM", "POS", "REF", "ALT", "DP", "AD1")
miss_v <- setdiff(req_v, names(gh_v_raw))
if (length(miss_v))
  stop("gh_v_raw missing required columns: ", paste(miss_v, collapse = ", "))

req_p <- c("sampleA", "sampleB")
miss_p <- setdiff(req_p, names(gh_pairs))
if (length(miss_p))
  stop("gh_pairs missing required columns: ", paste(miss_p, collapse = ", "))

if (!"pair_class" %in% names(gh_pairs)) {
  gh_pairs$pair_class <- NA_character_
}

message(sprintf("gh_v_raw: %s rows, %d samples.",
                format(nrow(gh_v_raw), big.mark = ","),
                length(unique(gh_v_raw$Sample))))
message(sprintf("gh_pairs: %d replicate pairs.", nrow(gh_pairs)))

# ---- Pre-grid filters -------------------------------------------------------
# Apply variant-class restriction (SNV-only) and PE/PPE handling BEFORE the
# grid is evaluated, so the calibration sees the same variant pool that
# downstream application will use.
n_pre <- nrow(gh_v_raw)
gh_v_filtered <- gh_v_raw

if ("FILTER" %in% names(gh_v_filtered)) {
  gh_v_filtered <- gh_v_filtered %>% filter(FILTER == "PASS")
  message(sprintf("After FILTER == PASS: %s rows",
                  format(nrow(gh_v_filtered), big.mark = ",")))
}

if (DROP_PPE_CAL) {
  gh_v_filtered <- exclude_ppe(gh_v_filtered)
  message(sprintf("After PE/PPE exclusion: %s rows",
                  format(nrow(gh_v_filtered), big.mark = ",")))
}

if (SNV_ONLY) {
  gh_v_filtered <- filter_snv_only(gh_v_filtered)
  message(sprintf("After SNV-only filter: %s rows (%.1f%% retained from PASS)",
                  format(nrow(gh_v_filtered), big.mark = ","),
                  100 * nrow(gh_v_filtered) / n_pre))
}

gh_v_filtered <- ensure_maf(gh_v_filtered)

# ---- Grid -------------------------------------------------------------------
grid <- tidyr::expand_grid(DP_min  = DP_GRID,
                           AD1_min = AD1_GRID,
                           MAF_min = MAF_GRID)

message(sprintf("\nCalibration grid: %d points (DP x AD1 x MAF = %d x %d x %d).",
                nrow(grid), length(DP_GRID), length(AD1_GRID), length(MAF_GRID)))

# ---- Per-sample variant lookup ----------------------------------------------
samples_in_pairs <- unique(c(gh_pairs$sampleA, gh_pairs$sampleB))
v_by_sample <- gh_v_filtered %>%
  filter(Sample %in% samples_in_pairs) %>%
  split(.$Sample)

empty_var_tbl <- gh_v_filtered[0, ]
get_v <- function(s) v_by_sample[[s]] %||% empty_var_tbl

# ---- Threshold filter (per-pair, used inside the grid loop) -----------------
apply_thr_inline <- function(v, DP_min, AD1_min, MAF_min) {
  v %>% filter(DP >= DP_min, AD1 >= AD1_min,
               MAF >= MAF_min, MAF <= MAF_MAX)
}

# ---- Compute per-pair concordance at each grid point ------------------------
message("\nIterating grid x pairs...")
n_grid <- nrow(grid)
n_pair <- nrow(gh_pairs)

results_list <- vector("list", n_grid * n_pair)
idx <- 1L
t0 <- Sys.time()

for (g in seq_len(n_grid)) {
  DP_min  <- grid$DP_min[g]
  AD1_min <- grid$AD1_min[g]
  MAF_min <- grid$MAF_min[g]

  for (p in seq_len(n_pair)) {
    sA <- gh_pairs$sampleA[p]
    sB <- gh_pairs$sampleB[p]

    a_tbl <- apply_thr_inline(get_v(sA), DP_min, AD1_min, MAF_min)
    b_tbl <- apply_thr_inline(get_v(sB), DP_min, AD1_min, MAF_min)

    a_sites <- site_key(a_tbl)
    b_sites <- site_key(b_tbl)

    nA <- length(a_sites)
    nB <- length(b_sites)
    n_union <- length(union(a_sites, b_sites))
    n_inter <- length(intersect(a_sites, b_sites))

    jacc <- if (n_union == 0) 0 else n_inter / n_union
    over <- if (min(nA, nB) == 0) 0 else n_inter / min(nA, nB)
    cr   <- confirm_rate_pair(a_tbl, b_tbl)
    mc   <- if (n_inter >= 2) maf_cor_pair(a_tbl, b_tbl) else NA_real_

    results_list[[idx]] <- tibble(
      DP_min          = DP_min,
      AD1_min         = AD1_min,
      MAF_min         = MAF_min,
      sampleA         = sA,
      sampleB         = sB,
      pair_class      = gh_pairs$pair_class[p],
      nA              = nA,
      nB              = nB,
      n_union         = n_union,
      n_inter         = n_inter,
      jaccard         = jacc,
      overlap         = over,
      confirm_both    = cr["confirm_both"],
      confirm_a_to_b  = cr["confirm_a_to_b"],
      confirm_b_to_a  = cr["confirm_b_to_a"],
      maf_cor         = mc
    )
    idx <- idx + 1L
  }

  if (g %% 10 == 0 || g == n_grid) {
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf("  grid point %d/%d (%.1f sec elapsed)",
                    g, n_grid, elapsed))
  }
}

cal_pairs <- bind_rows(results_list)

message(sprintf("\ncal_pairs: %s rows", format(nrow(cal_pairs), big.mark = ",")))

# ---- Write artifacts --------------------------------------------------------
saveRDS(cal_pairs,    file.path(PATHS$calibration, "cal_pairs.rds"))
readr::write_csv(cal_pairs,
                 file.path(PATHS$calibration, "cal_pairs.csv"))

message("\n[01_build_cal_pairs] Done. Artifacts written to ",
        PATHS$calibration, "/")
