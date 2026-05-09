# =============================================================================
# 04_external_validity.R
#
# External validity / transportability analysis. Implements the full
# sensitivity grid for Paper 1.
#
# Estimands (per cohort x cluster definition):
#   - Sensitivity (pair-level): P(shared CHROM:POS:REF:ALT iSNV | same cluster)
#   - Specificity (pair-level): P(no shared iSNV | different cluster)
#   - PPV: P(same cluster | shared iSNV)
#   - NPV: P(different cluster | no shared iSNV)
#   - PPV enrichment = PPV / within-cluster prevalence
#       (>1: test enriches above baseline -> transmission signal recovery
#         ~1: negative-control behavior, test no more informative than chance)
#   - PR: same_cluster vs different, modified Poisson with log link
#   - Mean Jaccard within vs between cluster
#
# Grid axes (full factorial, see 00_pipeline_config.R):
#   - Threshold tier:      loose | primary | tighter
#   - Cluster definition:  GH primary <=12 SNP-only; sensitivity <=5 SNP-only,
#                          <=12+spatial, <=5+spatial; FL surveillance
#   - Zero-iSNV handling:  all samples (canonical) | restrict to >=1 iSNV
#   - PE/PPE handling:     governed by which TAG of 03 outputs we read in
#
# Standard errors:
#   - Florida primary: dyadic-clustered (Cameron-Gelbach-Miller 2011)
#   - Bootstrap: resample sample IDs -> reconstruct pairs -> recompute,
#     N_BOOTSTRAP reps. Returns CIs for PR, PPV, and enrichment jointly.
#
# Inputs:
#   PATHS$thresholded/{fl,gh}_variants_thr_{loose,primary,tighter}.rds
#   PATHS$thresholded/{fl,gh}_calibrated_isnv_count.rds
#   PATHS$meta/{fl,gh}_meta.rds
#   data_raw/gis-coordinates.rds
#   data_raw/snp_distance_primary.rds
#
# Outputs (PATHS$external):
#   external_validity_grid.csv  -- full sensitivity grid results (incl. PPV)
#   external_validity_main.csv  -- canonical primary row + key sensitivity rows
#   bootstrap_summary.csv       -- bootstrap CIs for PR / PPV / enrichment
#   {fl,gh}_pairs.rds           -- pair tables (for 06 diagnostics)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(sandwich)
  library(lmtest)
})

set.seed(SEED)

# ---- Paths ------------------------------------------------------------------
fl_meta_path  <- file.path(PATHS$meta,        "fl_meta.rds")
gh_meta_path  <- file.path(PATHS$meta,        "gh_meta.rds")
fl_count_path <- file.path(PATHS$thresholded, "fl_calibrated_isnv_count.rds")
gh_count_path <- file.path(PATHS$thresholded, "gh_calibrated_isnv_count.rds")

gh_gis_path      <- file.path(PATHS$data_raw, "gis-coordinates.rds")
gh_snp_dist_path <- file.path(PATHS$data_raw, "snp_distance_primary.rds")

fl_var_paths <- list(
  loose   = file.path(PATHS$thresholded, "fl_variants_thr_loose.rds"),
  primary = file.path(PATHS$thresholded, "fl_variants_thr_primary.rds"),
  tighter = file.path(PATHS$thresholded, "fl_variants_thr_tighter.rds")
)
gh_var_paths <- list(
  loose   = file.path(PATHS$thresholded, "gh_variants_thr_loose.rds"),
  primary = file.path(PATHS$thresholded, "gh_variants_thr_primary.rds"),
  tighter = file.path(PATHS$thresholded, "gh_variants_thr_tighter.rds")
)

# ---- Load --------------------------------------------------------------------
fl_meta    <- readRDS(fl_meta_path)
gh_meta    <- readRDS(gh_meta_path)
fl_isnv    <- readRDS(fl_count_path)
gh_isnv    <- readRDS(gh_count_path)
gh_gis     <- readRDS(gh_gis_path)
gh_snp_mat <- readRDS(gh_snp_dist_path)

fl_meta <- fl_meta %>% inner_join(fl_isnv, by = "Sample")
gh_meta <- gh_meta %>% inner_join(gh_isnv, by = "Sample")

# ---- Pair-level utilities --------------------------------------------------
make_pairs <- function(samples) {
  ids <- combn(samples, 2)
  tibble(sample_a = ids[1, ], sample_b = ids[2, ])
}

combo_num <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")
combo_chr <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "_")

lookup_pair <- function(mat, a, b) {
  out <- vapply(seq_along(a), function(i) {
    if (a[i] %in% rownames(mat) && b[i] %in% colnames(mat)) mat[a[i], b[i]]
    else NA_real_
  }, numeric(1))
  out
}

build_variant_sets <- function(variant_tbl, Samples) {
  v <- variant_tbl %>%
    filter(Sample %in% Samples) %>%
    mutate(varkey = paste(CHROM, POS, REF, ALT, sep = ":"))
  vsets <- split(v$varkey, v$Sample)
  miss <- setdiff(Samples, names(vsets))
  for (id in miss) vsets[[id]] <- character(0)
  vsets
}

pair_isnv_metrics <- function(pairs, vsets) {
  pairs %>%
    rowwise() %>%
    mutate(
      a_set      = list(vsets[[sample_a]] %||% character(0)),
      b_set      = list(vsets[[sample_b]] %||% character(0)),
      n_a        = length(a_set),
      n_b        = length(b_set),
      n_inter    = length(intersect(a_set, b_set)),
      n_union    = length(union(a_set, b_set)),
      shared_iSNV = as.integer(n_inter >= 1),
      jaccard    = if (n_union == 0) 0 else n_inter / n_union
    ) %>%
    ungroup() %>%
    select(-a_set, -b_set)
}

# ---- Pair-level 2x2 metrics (raw counts) -----------------------------------
# Compute predictive values from raw counts in the pair table. PPV/NPV are
# defined empirically (no model). PPV enrichment = PPV / within-cluster
# prevalence is the principal transportability summary at the pair level.
compute_2x2_metrics <- function(pairs, same_col = "same_cluster") {
  pairs <- pairs %>% mutate(.same = .data[[same_col]])

  TP <- sum(pairs$.same == 1L & pairs$shared_iSNV == 1L, na.rm = TRUE)
  FP <- sum(pairs$.same == 0L & pairs$shared_iSNV == 1L, na.rm = TRUE)
  FN <- sum(pairs$.same == 1L & pairs$shared_iSNV == 0L, na.rm = TRUE)
  TN <- sum(pairs$.same == 0L & pairs$shared_iSNV == 0L, na.rm = TRUE)

  total <- TP + FP + FN + TN
  prev  <- if (total > 0)         (TP + FN) / total else NA_real_
  ppv   <- if ((TP + FP) > 0)     TP / (TP + FP)   else NA_real_
  npv   <- if ((TN + FN) > 0)     TN / (TN + FN)   else NA_real_
  enrich <- if (!is.na(prev) && prev > 0) ppv / prev else NA_real_

  list(
    TP = TP, FP = FP, FN = FN, TN = TN,
    prevalence     = prev,
    ppv            = ppv,
    npv            = npv,
    ppv_enrichment = enrich,
    sens_raw       = if ((TP + FN) > 0) TP / (TP + FN) else NA_real_,
    spec_raw       = if ((TN + FP) > 0) TN / (TN + FP) else NA_real_
  )
}

# ---- Florida: spatial / cluster info already in fl_meta --------------------
# Florida uses surveillance cluster_id directly. Add depth quartiles, year bin.
fl_meta <- fl_meta %>%
  mutate(
    depth_q  = ntile(depth_med, 4),
    year_bin = case_when(
      Year <= 2018           ~ "<=2018",
      Year %in% 2019:2020    ~ "2019_2020",
      Year >= 2021           ~ ">=2021",
      TRUE                   ~ "unk"
    )
  )

# ---- Ghana: spatial info ----------------------------------------------------
gh_meta <- gh_meta %>%
  mutate(
    depth_q  = ntile(depth_med, 4),
    year_bin = "GH_single_year"  # Ghana cohort cross-sectional 2022-2023
  )

stopifnot(all(c("patientId","avgLatitude","avgLongitude") %in% names(gh_gis)))

gh_gis_sf <- gh_gis %>%
  st_as_sf(coords = c("avgLongitude","avgLatitude"), crs = 4326) %>%
  st_transform(GH_UTM_EPSG)

gh_coords <- bind_cols(
  patientId = gh_gis_sf$patientId,
  as_tibble(st_coordinates(gh_gis_sf))
) %>% rename(easting_m = X, northing_m = Y)

miss_gis <- setdiff(gh_meta$patientId, gh_coords$patientId)
if (length(miss_gis) > 0) {
  message(sprintf("GH samples without GIS coordinates: %d (dropping)",
                  length(miss_gis)))
  gh_meta <- gh_meta %>% filter(!patientId %in% miss_gis)
}

gh_dist_mat <- gh_coords %>%
  select(easting_m, northing_m) %>% as.matrix() %>%
  dist() %>% as.matrix()
rownames(gh_dist_mat) <- gh_coords$patientId
colnames(gh_dist_mat) <- gh_coords$patientId

# ---- Build pair tables -------------------------------------------------------
build_fl_pairs <- function(fl_det) {
  keys <- fl_det %>% select(Sample, cluster_id, depth_q, lineage, year_bin)
  make_pairs(fl_det$Sample) %>%
    left_join(keys %>% rename_with(~ paste0(.x,"_a"), -Sample) %>% rename(sample_a = Sample),
              by = "sample_a") %>%
    left_join(keys %>% rename_with(~ paste0(.x,"_b"), -Sample) %>% rename(sample_b = Sample),
              by = "sample_b") %>%
    mutate(
      same_cluster    = as.integer(cluster_id_a == cluster_id_b),
      depth_combo     = combo_num(depth_q_a, depth_q_b),
      lineage_concord = as.integer(lineage_a == lineage_b),
      year_combo      = combo_chr(year_bin_a, year_bin_b)
    )
}

build_gh_pairs <- function(gh_det, snp_mat, dist_mat) {
  keys <- gh_det %>% select(Sample, patientId, depth_q, lineage, year_bin)
  make_pairs(gh_det$Sample) %>%
    left_join(keys %>% rename_with(~ paste0(.x,"_a"), -Sample) %>% rename(sample_a = Sample),
              by = "sample_a") %>%
    left_join(keys %>% rename_with(~ paste0(.x,"_b"), -Sample) %>% rename(sample_b = Sample),
              by = "sample_b") %>%
    filter(patientId_a != patientId_b) %>%
    mutate(
      snp_dist           = lookup_pair(snp_mat,  sample_a,    sample_b),
      spatial_dist_m     = lookup_pair(dist_mat, patientId_a, patientId_b),
      # Cluster definitions (canonical primary: <=12 SNP-only)
      same_cluster_12_snp_only      = as.integer(snp_dist <= GH_SNP_PRIMARY),
      same_cluster_5_snp_only       = as.integer(snp_dist <= GH_SNP_SENSITIVE),
      same_cluster_12_snp_plus_spat = as.integer(snp_dist <= GH_SNP_PRIMARY   &
                                                  spatial_dist_m <= GH_SPATIAL_M),
      same_cluster_5_snp_plus_spat  = as.integer(snp_dist <= GH_SNP_SENSITIVE &
                                                  spatial_dist_m <= GH_SPATIAL_M),
      depth_combo     = combo_num(depth_q_a, depth_q_b),
      lineage_concord = as.integer(lineage_a == lineage_b),
      year_combo      = combo_chr(year_bin_a, year_bin_b)
    )
}

# ---- Main inference: Florida (regression + standardization) ----------------
fit_fl_standardize <- function(pairs) {
  m <- glm(shared_iSNV ~ same_cluster + depth_combo + lineage_concord + year_combo,
           data = pairs, family = poisson(link = "log"))

  V_a  <- vcovCL(m, cluster = pairs$sample_a, type = "HC1")
  V_b  <- vcovCL(m, cluster = pairs$sample_b, type = "HC1")
  V_ab <- vcovCL(m, cluster = interaction(pairs$sample_a, pairs$sample_b),
                 type = "HC1")
  V_dyadic <- V_a + V_b - V_ab

  within  <- pairs %>% mutate(same_cluster = 1L)
  between <- pairs %>% mutate(same_cluster = 0L)
  pred_within  <- predict(m, newdata = within  %>% filter(pairs$same_cluster == 1),
                          type = "response")
  pred_between <- predict(m, newdata = between %>% filter(pairs$same_cluster == 0),
                          type = "response")

  pr_log <- coef(m)["same_cluster"]
  pr_se  <- sqrt(V_dyadic["same_cluster","same_cluster"])

  # Raw 2x2 (PPV/NPV/enrichment from observed counts, not regression)
  m2x2 <- compute_2x2_metrics(pairs, same_col = "same_cluster")

  list(
    sens       = mean(pred_within),
    spec       = 1 - mean(pred_between),
    pr         = exp(pr_log),
    pr_lo      = exp(pr_log - 1.96 * pr_se),
    pr_hi      = exp(pr_log + 1.96 * pr_se),
    j_w        = mean(pairs$jaccard[pairs$same_cluster == 1]),
    j_b        = mean(pairs$jaccard[pairs$same_cluster == 0]),
    n_w        = sum(pairs$same_cluster == 1),
    n_b        = sum(pairs$same_cluster == 0),
    # 2x2 metrics
    ppv        = m2x2$ppv,
    npv        = m2x2$npv,
    enrichment = m2x2$ppv_enrichment,
    prevalence = m2x2$prevalence,
    sens_raw   = m2x2$sens_raw,
    spec_raw   = m2x2$spec_raw,
    TP = m2x2$TP, FP = m2x2$FP, FN = m2x2$FN, TN = m2x2$TN
  )
}

# ---- Main inference: Ghana (unmatched, modified Poisson) -------------------
fit_gh_unmatched <- function(pairs, same_col) {
  pairs <- pairs %>% mutate(.same = .data[[same_col]])

  n_w <- sum(pairs$.same == 1L)
  n_b <- sum(pairs$.same == 0L)
  if (n_w < 2) {
    return(list(sens=NA, spec=NA, pr=NA, pr_lo=NA, pr_hi=NA,
                j_w=NA, j_b=NA, n_w=n_w, n_b=n_b,
                ppv=NA, npv=NA, enrichment=NA, prevalence=NA,
                sens_raw=NA, spec_raw=NA,
                TP=NA, FP=NA, FN=NA, TN=NA))
  }

  sens <- mean(pairs$shared_iSNV[pairs$.same == 1L])
  spec <- 1 - mean(pairs$shared_iSNV[pairs$.same == 0L])
  j_w  <- mean(pairs$jaccard[pairs$.same == 1L])
  j_b  <- mean(pairs$jaccard[pairs$.same == 0L])

  m <- tryCatch(
    glm(shared_iSNV ~ .same, data = pairs, family = poisson(link = "log")),
    error = function(e) NULL
  )

  if (!is.null(m) && !is.na(coef(m)[".same"])) {
    V_a  <- vcovCL(m, cluster = pairs$sample_a, type = "HC1")
    V_b  <- vcovCL(m, cluster = pairs$sample_b, type = "HC1")
    V_ab <- vcovCL(m, cluster = interaction(pairs$sample_a, pairs$sample_b),
                   type = "HC1")
    V <- V_a + V_b - V_ab
    pr_log <- coef(m)[".same"]
    se     <- sqrt(V[".same",".same"])
    pr     <- exp(pr_log)
    lo     <- exp(pr_log - 1.96*se)
    hi     <- exp(pr_log + 1.96*se)
  } else {
    pr <- NA_real_; lo <- NA_real_; hi <- NA_real_
  }

  # Raw 2x2 (PPV/NPV/enrichment)
  m2x2 <- compute_2x2_metrics(pairs, same_col = ".same")

  list(sens=sens, spec=spec, pr=pr, pr_lo=lo, pr_hi=hi,
       j_w=j_w, j_b=j_b, n_w=n_w, n_b=n_b,
       ppv        = m2x2$ppv,
       npv        = m2x2$npv,
       enrichment = m2x2$ppv_enrichment,
       prevalence = m2x2$prevalence,
       sens_raw   = m2x2$sens_raw,
       spec_raw   = m2x2$spec_raw,
       TP = m2x2$TP, FP = m2x2$FP, FN = m2x2$FN, TN = m2x2$TN)
}

# ---- Bootstrap helper -------------------------------------------------------
# Returns 95% percentile CIs for PR, PPV, and enrichment jointly. Resamples
# sample IDs (Cameron-Gelbach-Miller-style cluster bootstrap), reconstructs
# pairs, and recomputes the fit. Degenerate same-sample pairs are dropped.
bootstrap_metrics <- function(pairs, fit_fn, n_reps = N_BOOTSTRAP, ...) {
  samples <- unique(c(pairs$sample_a, pairs$sample_b))
  pr_b   <- numeric(n_reps)
  ppv_b  <- numeric(n_reps)
  enr_b  <- numeric(n_reps)
  sens_b <- numeric(n_reps)
  spec_b <- numeric(n_reps)

  message(sprintf("    n_reps=%d, n_samples=%d", n_reps, length(samples)))

  for (b in seq_len(n_reps)) {
    samp_b <- sample(samples, length(samples), replace = TRUE)
    boot_pairs <- pairs %>%
      filter(sample_a %in% samp_b, sample_b %in% samp_b) %>%
      filter(sample_a != sample_b)
    if (nrow(boot_pairs) < 10) {
      pr_b[b] <- ppv_b[b] <- enr_b[b] <- sens_b[b] <- spec_b[b] <- NA_real_
      next
    }
    fit <- tryCatch(fit_fn(boot_pairs, ...), error = function(e) NULL)
    if (is.null(fit)) {
      pr_b[b] <- ppv_b[b] <- enr_b[b] <- sens_b[b] <- spec_b[b] <- NA_real_
    } else {
      pr_b[b]   <- fit$pr
      ppv_b[b]  <- fit$ppv
      enr_b[b]  <- fit$enrichment
      sens_b[b] <- fit$sens
      spec_b[b] <- fit$spec
    }
    if (b %% 100 == 0) message(sprintf("    bootstrap rep %d/%d", b, n_reps))
  }

  q <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE)
  q_pr   <- q(pr_b)
  q_ppv  <- q(ppv_b)
  q_enr  <- q(enr_b)
  q_sens <- q(sens_b)
  q_spec <- q(spec_b)

  list(
    pr_boot_lo   = q_pr[1],   pr_boot_hi   = q_pr[2],
    ppv_boot_lo  = q_ppv[1],  ppv_boot_hi  = q_ppv[2],
    enr_boot_lo  = q_enr[1],  enr_boot_hi  = q_enr[2],
    sens_boot_lo = q_sens[1], sens_boot_hi = q_sens[2],
    spec_boot_lo = q_spec[1], spec_boot_hi = q_spec[2],
    n_valid      = sum(!is.na(pr_b))
  )
}

# ---- Sensitivity grid runner ------------------------------------------------
run_one_cell <- function(threshold_tier, cluster_def, zero_drop,
                          fl_det_full, gh_det_full,
                          fl_vsets, gh_vsets) {

  if (zero_drop) {
    fl_det <- fl_det_full %>% filter(n_calibrated_isnv >= 1)
    gh_det <- gh_det_full %>% filter(n_calibrated_isnv >= 1)
  } else {
    fl_det <- fl_det_full
    gh_det <- gh_det_full
  }

  fl_pairs <- build_fl_pairs(fl_det) %>% pair_isnv_metrics(fl_vsets)
  gh_pairs <- build_gh_pairs(gh_det, gh_snp_mat, gh_dist_mat) %>%
                pair_isnv_metrics(gh_vsets)

  # Florida: surveillance cluster (one definition only)
  if (cluster_def == "fl_surveillance") {
    fit <- fit_fl_standardize(fl_pairs)
    cohort <- "Florida"
  } else {
    # Ghana: pick the right cluster column
    same_col <- switch(
      cluster_def,
      "gh_12_snp_only"      = "same_cluster_12_snp_only",
      "gh_5_snp_only"       = "same_cluster_5_snp_only",
      "gh_12_snp_plus_spat" = "same_cluster_12_snp_plus_spat",
      "gh_5_snp_plus_spat"  = "same_cluster_5_snp_plus_spat"
    )
    fit <- fit_gh_unmatched(gh_pairs, same_col)
    cohort <- "Ghana"
  }

  tibble(
    cohort           = cohort,
    threshold_tier   = threshold_tier,
    cluster_def      = cluster_def,
    zero_drop        = zero_drop,
    apply_tag        = TAG$apply_tag,
    sens             = fit$sens,
    spec             = fit$spec,
    pr               = fit$pr,
    pr_lo            = fit$pr_lo,
    pr_hi            = fit$pr_hi,
    jaccard_within   = fit$j_w,
    jaccard_between  = fit$j_b,
    n_within         = fit$n_w,
    n_between        = fit$n_b,
    # Predictive values and 2x2 counts
    prevalence       = fit$prevalence,
    ppv              = fit$ppv,
    npv              = fit$npv,
    enrichment       = fit$enrichment,
    sens_raw         = fit$sens_raw,
    spec_raw         = fit$spec_raw,
    TP               = fit$TP,
    FP               = fit$FP,
    FN               = fit$FN,
    TN               = fit$TN
  )
}

# ---- Build variant sets for each tier --------------------------------------
load_vsets_for_tier <- function(fl_var_path, gh_var_path, fl_samples, gh_samples) {
  if (!file.exists(fl_var_path) || !file.exists(gh_var_path)) {
    return(list(fl = NULL, gh = NULL))
  }
  list(
    fl = build_variant_sets(readRDS(fl_var_path), fl_samples),
    gh = build_variant_sets(readRDS(gh_var_path), gh_samples)
  )
}

# ---- Run the grid ----------------------------------------------------------
cluster_defs <- c("fl_surveillance",
                  "gh_12_snp_only",
                  "gh_5_snp_only",
                  "gh_12_snp_plus_spat",
                  "gh_5_snp_plus_spat")

threshold_tiers <- c("loose", "primary", "tighter")
zero_drop_levels <- c(FALSE, TRUE)

results <- list()
i <- 0L

for (tier in threshold_tiers) {
  vsets <- load_vsets_for_tier(
    fl_var_paths[[tier]], gh_var_paths[[tier]],
    fl_meta$Sample, gh_meta$Sample
  )
  if (is.null(vsets$fl) || is.null(vsets$gh)) {
    message(sprintf("Skipping tier '%s' (variants missing)", tier))
    next
  }

  for (cd in cluster_defs) {
    for (zd in zero_drop_levels) {
      i <- i + 1L
      message(sprintf("[%d] tier=%s cluster=%s zero_drop=%s",
                      i, tier, cd, zd))
      results[[i]] <- run_one_cell(
        threshold_tier = tier,
        cluster_def    = cd,
        zero_drop      = zd,
        fl_det_full    = fl_meta,
        gh_det_full    = gh_meta,
        fl_vsets       = vsets$fl,
        gh_vsets       = vsets$gh
      )
    }
  }
}

grid_results <- bind_rows(results)

readr::write_csv(grid_results,
                 file.path(PATHS$external, "external_validity_grid.csv"))

# ---- Bootstrap canonical primary cells -------------------------------------
# Bootstrap on canonical primary cells (FL surveillance + GH 12 SNP-only)
# under primary threshold + zero_drop = FALSE. Returns CIs for PR, PPV,
# enrichment, sens, and spec jointly.
if (N_BOOTSTRAP > 0) {
  message("\nBootstrapping canonical primary cells...")

  fl_vsets_prim <- build_variant_sets(readRDS(fl_var_paths$primary), fl_meta$Sample)
  gh_vsets_prim <- build_variant_sets(readRDS(gh_var_paths$primary), gh_meta$Sample)

  fl_pairs_prim <- build_fl_pairs(fl_meta) %>% pair_isnv_metrics(fl_vsets_prim)
  gh_pairs_prim <- build_gh_pairs(gh_meta, gh_snp_mat, gh_dist_mat) %>%
                     pair_isnv_metrics(gh_vsets_prim)

  message("  Florida bootstrap...")
  fl_boot <- bootstrap_metrics(fl_pairs_prim, fit_fl_standardize)

  message("  Ghana bootstrap...")
  gh_boot <- bootstrap_metrics(gh_pairs_prim, fit_gh_unmatched,
                               same_col = "same_cluster_12_snp_only")

  message(sprintf("\nFL primary bootstrap: PR (%.3f, %.3f) | PPV (%.3f, %.3f) | enrichment (%.2f, %.2f) [n_valid=%d]",
                  fl_boot$pr_boot_lo, fl_boot$pr_boot_hi,
                  fl_boot$ppv_boot_lo, fl_boot$ppv_boot_hi,
                  fl_boot$enr_boot_lo, fl_boot$enr_boot_hi,
                  fl_boot$n_valid))
  message(sprintf("GH 12-SNP-only bootstrap: PR (%.3f, %.3f) | PPV (%.3f, %.3f) | enrichment (%.2f, %.2f) [n_valid=%d]",
                  gh_boot$pr_boot_lo, gh_boot$pr_boot_hi,
                  gh_boot$ppv_boot_lo, gh_boot$ppv_boot_hi,
                  gh_boot$enr_boot_lo, gh_boot$enr_boot_hi,
                  gh_boot$n_valid))

  # Persist pair tables for 06 diagnostics
  saveRDS(fl_pairs_prim, file.path(PATHS$external, "fl_pairs.rds"))
  saveRDS(gh_pairs_prim, file.path(PATHS$external, "gh_pairs.rds"))

  bootstrap_summary <- tibble(
    cell        = c("FL_surveillance_primary", "GH_12_snp_only_primary"),
    pr_boot_lo  = c(fl_boot$pr_boot_lo,   gh_boot$pr_boot_lo),
    pr_boot_hi  = c(fl_boot$pr_boot_hi,   gh_boot$pr_boot_hi),
    ppv_boot_lo = c(fl_boot$ppv_boot_lo,  gh_boot$ppv_boot_lo),
    ppv_boot_hi = c(fl_boot$ppv_boot_hi,  gh_boot$ppv_boot_hi),
    enr_boot_lo = c(fl_boot$enr_boot_lo,  gh_boot$enr_boot_lo),
    enr_boot_hi = c(fl_boot$enr_boot_hi,  gh_boot$enr_boot_hi),
    sens_boot_lo = c(fl_boot$sens_boot_lo, gh_boot$sens_boot_lo),
    sens_boot_hi = c(fl_boot$sens_boot_hi, gh_boot$sens_boot_hi),
    spec_boot_lo = c(fl_boot$spec_boot_lo, gh_boot$spec_boot_lo),
    spec_boot_hi = c(fl_boot$spec_boot_hi, gh_boot$spec_boot_hi),
    n_valid     = c(fl_boot$n_valid,    gh_boot$n_valid)
  )
  readr::write_csv(bootstrap_summary,
                   file.path(PATHS$external, "bootstrap_summary.csv"))
}

# ---- Build "main" table (canonical row + monotonicity check) ----------------
canonical <- grid_results %>%
  filter(zero_drop == FALSE,
         (cohort == "Florida" & cluster_def == "fl_surveillance") |
         (cohort == "Ghana"   & cluster_def == "gh_12_snp_only"))

readr::write_csv(canonical, file.path(PATHS$external, "external_validity_main.csv"))

cat("\nCanonical table (zero_drop=FALSE, primary cluster definitions):\n")
print(canonical %>% select(cohort, threshold_tier, sens, spec, ppv, enrichment,
                           pr, pr_lo, pr_hi,
                           jaccard_within, jaccard_between, n_within, n_between))

# Monotonicity check: PR_loose <= PR_primary <= PR_tighter expected in Florida
cat("\nThreshold-ladder PR monotonicity check (Florida, canonical):\n")
fl_ladder <- canonical %>% filter(cohort == "Florida") %>%
  arrange(match(threshold_tier, c("loose","primary","tighter")))
print(fl_ladder %>% select(threshold_tier, pr, pr_lo, pr_hi, ppv, enrichment))

if (nrow(fl_ladder) == 3) {
  monotonic <- fl_ladder$pr[1] <= fl_ladder$pr[2] &
               fl_ladder$pr[2] <= fl_ladder$pr[3]
  cat(sprintf("Monotonicity holds: %s\n", monotonic))
}

cat(sprintf("\n[04_external_validity] Done. Outputs written to: %s/\n",
            PATHS$external))
