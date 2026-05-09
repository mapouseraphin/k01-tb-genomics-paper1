# =============================================================================
# 00_pipeline_config.R
#
# Single source of truth for canonical specification flags and tag construction.
# Sourced at the top of every pipeline script (00 -> 08) so a single edit here
# propagates through the entire run.
#
# Canonical specification for Paper 1:
#   - iSNV definition: SNVs only (indels excluded)
#   - PE/PPE: retained in primary; excluded as sensitivity
#   - Calibration set: M0 within-visit replicate pairs only
#   - Guards: original (MIN_EVAL_PAIRS=10, MAX_PROP_BOTH0=0.70)
#
# To run a sensitivity variant: flip a flag below, re-run 00->07. Output
# directories AND tables/figures carry the active spec tag automatically;
# previous runs are not overwritten.
# =============================================================================

# ---- Canonical flags --------------------------------------------------------
SNV_ONLY        <- TRUE    # canonical: SNVs only. Sensitivity: FALSE (indels in)
DROP_PPE_CAL    <- FALSE   # canonical: PE/PPE retained in calibration step
DROP_PPE_APPLY  <- FALSE   # canonical: PE/PPE retained in application step
M0_ONLY         <- TRUE    # canonical: restrict cal_pairs to M0 within-visit
ZERO_ISNV_DROP  <- FALSE   # canonical: include all samples in 04
                           # Sensitivity: TRUE (restrict to >=1 calibrated iSNV)

# ---- Calibration guards (script 02) -----------------------------------------
MIN_EVAL_PAIRS <- 10L      # original stability floor; relax only with explicit
                           # documentation
MAX_PROP_BOTH0 <- 0.70     # original ceiling on fraction of zero-zero pairs

# ---- Variant-class restriction grid (script 01) -----------------------------
DP_GRID  <- seq(40, 80, by = 10)
AD1_GRID <- 1:6
MAF_GRID <- c(0.02, 0.03, 0.05, 0.10, 0.15, 0.20)
MAF_MAX  <- 0.45

# ---- Cohort QC ---------------------------------------------------------------
MIN_COV_INCLUDE <- 50L     # medianCov floor for analytic frame inclusion

# ---- Ghana cluster definitions (script 04) ----------------------------------
GH_SNP_PRIMARY    <- 12L   # primary cluster definition: <=12 SNP-only (Asare 2020)
GH_SNP_SENSITIVE  <- 5L    # sensitivity: <=5 SNP-only (stricter)
GH_SPATIAL_M      <- 1500L # used only in spatial sensitivity analyses
GH_UTM_EPSG       <- 32630 # Ghana UTM zone (Korle-Bu)
GH_MATCH_RATIO    <- 3L    # 1:3 cell-level matching for sensitivity panels

# ---- Bootstrap settings (script 04) -----------------------------------------
SEED        <- 20260428L
N_BOOTSTRAP <- 1000L       # production: 1000. Set lower (e.g., 200) for dev runs.

# ---- Tag construction --------------------------------------------------------
# Output directories carry a 4-part tag describing variant-class and PPE handling.
# Calibration vs application tags are kept separate because PE/PPE handling can
# differ across the two steps (e.g., calibration with PPE retained, application
# sensitivity with PPE excluded).
build_tag <- function(snv_only = SNV_ONLY,
                      drop_ppe_cal = DROP_PPE_CAL,
                      drop_ppe_apply = DROP_PPE_APPLY) {
  variant_part   <- if (snv_only) "snv_only" else "with_indels"
  ppe_cal_part   <- if (drop_ppe_cal)   "ppe_excluded" else "ppe_retained"
  ppe_apply_part <- if (drop_ppe_apply) "ppe_excluded" else "ppe_retained"
  list(
    variant   = variant_part,
    ppe_cal   = ppe_cal_part,
    ppe_apply = ppe_apply_part,
    cal_tag   = paste(variant_part, ppe_cal_part,   sep = "_"),
    apply_tag = paste(variant_part, ppe_apply_part, sep = "_")
  )
}

TAG <- build_tag()

# ---- Output paths -----------------------------------------------------------
# All per-spec outputs (data_derived AND outputs/tables, outputs/figures) carry
# the active spec's tag. Cross-spec contrast outputs (script 08) sit at top-
# level outputs/contrast/.
PATHS <- list(
  data_raw        = "data_raw",
  meta            = "data_derived/00_metadata",
  variants        = "data_derived/00_variants",
  calibration     = file.path("data_derived",
                              paste0("01_calibration_",  TAG$cal_tag)),
  thresholded     = file.path("data_derived",
                              paste0("03_thresholded_", TAG$apply_tag)),
  external        = file.path("data_derived",
                              paste0("04_external_",    TAG$apply_tag)),
  diagnostics     = file.path("data_derived",
                              paste0("06_diagnostics_", TAG$apply_tag)),
  tables          = file.path("outputs", "tables",  TAG$apply_tag),
  figures         = file.path("outputs", "figures", TAG$apply_tag),
  contrast_tables  = file.path("outputs", "contrast", "tables"),
  contrast_figures = file.path("outputs", "contrast", "figures")
)

# Create output directories ---------------------------------------------------
for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- Cross-spec path resolver (used by script 08) ---------------------------
# Returns the per-spec output paths for an arbitrary spec without touching the
# active TAG. Script 08 calls this once per spec to enumerate available outputs.
build_paths_for_spec <- function(snv_only, drop_ppe_cal, drop_ppe_apply) {
  tag <- build_tag(snv_only      = snv_only,
                   drop_ppe_cal  = drop_ppe_cal,
                   drop_ppe_apply = drop_ppe_apply)
  list(
    cal_tag      = tag$cal_tag,
    apply_tag    = tag$apply_tag,
    calibration  = file.path("data_derived",
                             paste0("01_calibration_",  tag$cal_tag)),
    thresholded  = file.path("data_derived",
                             paste0("03_thresholded_", tag$apply_tag)),
    external     = file.path("data_derived",
                             paste0("04_external_",    tag$apply_tag)),
    diagnostics  = file.path("data_derived",
                             paste0("06_diagnostics_", tag$apply_tag)),
    tables       = file.path("outputs", "tables",  tag$apply_tag),
    figures      = file.path("outputs", "figures", tag$apply_tag)
  )
}

# Diagnostic banner -----------------------------------------------------------
banner <- function() {
  cat(strrep("=", 78), "\n", sep = "")
  cat(sprintf("Pipeline canonical specification\n"))
  cat(sprintf("  SNV_ONLY         = %s\n", SNV_ONLY))
  cat(sprintf("  DROP_PPE_CAL     = %s  (calibration step)\n", DROP_PPE_CAL))
  cat(sprintf("  DROP_PPE_APPLY   = %s  (application step)\n", DROP_PPE_APPLY))
  cat(sprintf("  M0_ONLY          = %s\n", M0_ONLY))
  cat(sprintf("  ZERO_ISNV_DROP   = %s\n", ZERO_ISNV_DROP))
  cat(sprintf("  MIN_EVAL_PAIRS   = %d\n", MIN_EVAL_PAIRS))
  cat(sprintf("  MAX_PROP_BOTH0   = %.2f\n", MAX_PROP_BOTH0))
  cat(sprintf("  N_BOOTSTRAP      = %d\n", N_BOOTSTRAP))
  cat(sprintf("  Calibration tag  = %s\n", TAG$cal_tag))
  cat(sprintf("  Application tag  = %s\n", TAG$apply_tag))
  cat(strrep("=", 78), "\n", sep = "")
}
