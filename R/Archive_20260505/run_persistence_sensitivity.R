# R/run_persistence_sensitivity.R
# P-B sensitivity reruns: persistence application at two regime cells.
# Run AFTER the canonical pipeline (anchor) completes.

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")
suppressPackageStartupMessages({ library(dplyr); library(readr) })

sens_cells <- list(
  list(label = "sens_70_6_005",
       DP_min = 70, AD1_min = 6, MAF_min = 0.05, MAF_max = MAF_MAX),
  list(label = "sens_70_6_015",
       DP_min = 70, AD1_min = 6, MAF_min = 0.15, MAF_max = MAF_MAX)
)

gh_variants_raw <- readRDS(file.path(PATHS$variants, "gh_variants.rds")) %>%
  mutate(Sample = as.character(Sample))

PATHS_canonical <- PATHS  # backup for restoration

for (cell in sens_cells) {
  message(sprintf("\n=== Sensitivity rerun: %s (DP=%g, AD1=%g, MAF=%g) ===",
                  cell$label, cell$DP_min, cell$AD1_min, cell$MAF_min))
  
  # 1. Apply rule
  variants_thr <- apply_thresholds(
    gh_variants_raw,
    DP_min = cell$DP_min, AD1_min = cell$AD1_min,
    MAF_min = cell$MAF_min, MAF_max = cell$MAF_max,
    snv_only = SNV_ONLY, drop_ppe = DROP_PPE_APPLY
  )
  
  # 2. Sub-directories under canonical per-tag dirs
  sens_thresh_dir <- file.path(PATHS_canonical$thresholded, cell$label)
  sens_app_dir    <- file.path(PATHS_canonical$application, cell$label)
  dir.create(sens_thresh_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(sens_app_dir,    recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(variants_thr,
          file.path(sens_thresh_dir, "gh_variants_thr_primary.rds"))
  
  # 3. Redirect PATHS for the sourced application; the script reads
  #    PATHS$thresholded and PATHS$application from .GlobalEnv.
  PATHS$thresholded <<- sens_thresh_dir
  PATHS$application <<- sens_app_dir
  
  source("R/application_persistence_analysis.R")
}

# 4. Restore PATHS
PATHS <<- PATHS_canonical
message("\n[run_persistence_sensitivity] PATHS restored to canonical.")