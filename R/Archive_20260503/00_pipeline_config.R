# =============================================================================
# 00_pipeline_config.R
#
# Single source of truth for canonical specification flags and tag construction.
# Sourced at the top of every pipeline script (00 -> 06) so a single edit here
# propagates through the entire run.
#
# Canonical specification for Paper 1 (LOCKED 2026-04-29):
#   - iSNV definition: SNVs only (indels excluded)
#   - PE/PPE: retained in primary; excluded as sensitivity (Sens A)
#   - Universe: Definition C -- positions where >=1 PASS GATK call exists
#               in the pair (Definition A available as Sens B)
#   - MAF binning: 4 bins [0.02,0.05), [0.05,0.10), [0.10,0.25), [0.25,0.45]
#                  (3-bin merged-low binning available as Sens C)
#   - Calibration set: M0 within-visit replicate pairs only
#   - Guards: original (MIN_EVAL_PAIRS=10, MAX_PROP_BOTH0=0.70)
#
# To run a sensitivity variant: flip the relevant flag below, re-run the
# pipeline. Output directories carry the active spec tag automatically;
# previous runs are not overwritten. Canonical run keeps the original tag
# `snv_only_ppe_retained`; sensitivities append disambiguating suffixes.
# =============================================================================

# ---- Canonical flags --------------------------------------------------------
# CRITICAL: these defaults reproduce the locked Paper 1 results. Changing any
# default below changes the headline run.
SNV_ONLY        <- TRUE    # canonical: SNVs only. 
DROP_PPE_CAL    <- TRUE   # canonical: PE/PPE retained in calibration step
DROP_PPE_APPLY  <- TRUE   # canonical: PE/PPE retained in application step
M0_ONLY         <- TRUE    # canonical: restrict cal_pairs to M0 within-visit

# ---- LCA universe + binning flags (added 2026-04-30) -----------------------
# Sensitivity LCAs B and C operate on the canonical Stage 1 calibration but
# rebuild the LCA dataset under alternative universe / binning choices.
UNIVERSE        <- "C"        # "C" = >=1 PASS GATK call in pair (canonical)
                              # "A" = >=1 GATK call (PASS or fail) in pair
MAF_BINNING     <- "default"  # "default" = [0.02,0.05) [0.05,0.10) [0.10,0.25) [0.25,0.45]
                              # "merged_low" = [0.02,0.10) [0.10,0.25) [0.25,0.45]

# NOTE: ZERO_ISNV_DROP removed 2026-04-30. It was a sample-level filter for
# the retired pair-level external validity framework. Inactive in the LCA
# pipeline; the LCA likelihood handles all-zero pairs correctly.

# ---- Calibration guards (script 02) -----------------------------------------
MIN_EVAL_PAIRS <- 6L      # original stability floor; relax only with explicit
                           # documentation
MAX_PROP_BOTH0 <- 0.95     # original ceiling on fraction of zero-zero pairs

# ---- Variant-class restriction grid (script 01) -----------------------------
DP_GRID  <- seq(40, 80, by = 10)
AD1_GRID <- 1:6
MAF_GRID <- c(0.02, 0.03, 0.05, 0.10, 0.15, 0.20)
MAF_MAX  <- 0.45

# ---- Cohort QC ---------------------------------------------------------------
MIN_COV_INCLUDE <- 50L     # medianCov floor for analytic frame inclusion

# ---- Bootstrap settings ------------------------------------------------------
SEED                 <- 20260428L
N_BOOTSTRAP          <- 1000L  # external-validity bootstrap (legacy; kept for archived 04)
N_BOOTSTRAP_LEX_SEL  <- 1000L  # Stage 1 lexicographic-selection stability bootstrap (script 02b)

# Legacy: Florida cluster definitions, kept commented out as a marker that
# external validity is no longer a Paper 1 deliverable. Restore only if a
# Paper 1 v2 reincorporates Florida.
# GH_SNP_PRIMARY    <- 12L
# GH_SNP_SENSITIVE  <- 5L
# GH_SPATIAL_M      <- 1500L
# GH_UTM_EPSG       <- 32630
# GH_MATCH_RATIO    <- 3L

# ---- Tag construction --------------------------------------------------------
# Stage 1 (calibration) output depends on SNV_ONLY and DROP_PPE_CAL only.
# Stage 2a (apply_thresholds) output depends on SNV_ONLY and DROP_PPE_APPLY.
# Stage 2b (LCA dataset / fits) additionally depends on UNIVERSE and MAF_BINNING;
#   when both are at default, lca_tag = apply_tag (preserves the original
#   canonical tag exactly). When non-default, suffixes "_univA" / "_mafmerge"
#   are appended so canonical and sensitivity outputs do not collide.
build_tag <- function(snv_only       = SNV_ONLY,
                      drop_ppe_cal   = DROP_PPE_CAL,
                      drop_ppe_apply = DROP_PPE_APPLY,
                      universe       = UNIVERSE,
                      maf_binning    = MAF_BINNING) {

  variant_part   <- if (snv_only)       "snv_only"      else "with_indels"
  ppe_cal_part   <- if (drop_ppe_cal)   "ppe_excluded"  else "ppe_retained"
  ppe_apply_part <- if (drop_ppe_apply) "ppe_excluded"  else "ppe_retained"

  cal_tag   <- paste(variant_part, ppe_cal_part,   sep = "_")
  apply_tag <- paste(variant_part, ppe_apply_part, sep = "_")

  # LCA tag: only append suffixes when UNIVERSE / MAF_BINNING differ from
  # canonical default. Keeps the canonical lca_tag identical to apply_tag,
  # which preserves the locked Paper 1 output directory name exactly.
  lca_suffix <- ""
  if (!identical(universe, "C"))            lca_suffix <- paste0(lca_suffix, "_univ", universe)
  if (!identical(maf_binning, "default"))   lca_suffix <- paste0(lca_suffix, "_mafmerge")
  lca_tag <- paste0(apply_tag, lca_suffix)

  list(
    variant     = variant_part,
    ppe_cal     = ppe_cal_part,
    ppe_apply   = ppe_apply_part,
    universe    = universe,
    maf_binning = maf_binning,
    cal_tag     = cal_tag,
    apply_tag   = apply_tag,
    lca_tag     = lca_tag
  )
}

TAG <- build_tag()

# ---- Output paths -----------------------------------------------------------
# All per-spec outputs (data_derived AND outputs/tables, outputs/figures) carry
# the active spec's tag. Calibration paths use cal_tag; thresholded paths use
# apply_tag; LCA paths use lca_tag.
PATHS <- list(
  data_raw         = "data_raw",
  meta             = "data_derived/00_metadata",
  variants         = "data_derived/00_variants",
  calibration      = file.path("data_derived",
                               paste0("01_calibration_",  TAG$cal_tag)),
  bootstrap        = file.path("data_derived",
                               paste0("01_calibration_",  TAG$cal_tag),
                               "bootstrap"),
  thresholded      = file.path("data_derived",
                               paste0("03_thresholded_", TAG$apply_tag)),
  lca              = file.path("data_derived",
                               paste0("05_lca_",         TAG$lca_tag)),
  lca_fits         = file.path("data_derived",
                               paste0("05_lca_",         TAG$lca_tag), "fits"),
  lca_diag         = file.path("data_derived",
                               paste0("05_lca_",         TAG$lca_tag), "diagnostics"),
  tables           = file.path("outputs", "tables",       TAG$lca_tag),
  figures          = file.path("outputs", "figures",      TAG$lca_tag),
  supplemental     = file.path("outputs", "supplemental", TAG$lca_tag),
  cross_spec       = file.path("outputs", "cross_spec_lca")  # NEW: holds
                                                             # canonical-vs-
                                                             # sensitivity LCA
                                                             # comparisons
  # Legacy paths kept commented to mark the architectural decision:
  # external         = file.path("data_derived", paste0("04_external_",    TAG$apply_tag)),
  # diagnostics_inv  = file.path("data_derived", paste0("06_diagnostics_", TAG$apply_tag)),
  # contrast_tables  = file.path("outputs", "contrast", "tables"),
  # contrast_figures = file.path("outputs", "contrast", "figures")
)

# Create output directories ---------------------------------------------------
for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- Cross-spec path resolver -----------------------------------------------
# Returns the per-spec output paths for an arbitrary spec without touching the
# active TAG. Used by sensitivity-LCA aggregation in script 11.
build_paths_for_spec <- function(snv_only, drop_ppe_cal, drop_ppe_apply,
                                 universe = "C", maf_binning = "default") {
  tag <- build_tag(snv_only       = snv_only,
                   drop_ppe_cal   = drop_ppe_cal,
                   drop_ppe_apply = drop_ppe_apply,
                   universe       = universe,
                   maf_binning    = maf_binning)
  list(
    cal_tag      = tag$cal_tag,
    apply_tag    = tag$apply_tag,
    lca_tag      = tag$lca_tag,
    calibration  = file.path("data_derived",
                             paste0("01_calibration_",  tag$cal_tag)),
    thresholded  = file.path("data_derived",
                             paste0("03_thresholded_", tag$apply_tag)),
    lca          = file.path("data_derived",
                             paste0("05_lca_",         tag$lca_tag)),
    lca_fits     = file.path("data_derived",
                             paste0("05_lca_",         tag$lca_tag), "fits"),
    tables       = file.path("outputs", "tables",      tag$lca_tag),
    figures      = file.path("outputs", "figures",     tag$lca_tag),
    supplemental = file.path("outputs", "supplemental", tag$lca_tag)
  )
}

# Diagnostic banner -----------------------------------------------------------
banner <- function() {
  cat(strrep("=", 78), "\n", sep = "")
  cat(sprintf("Pipeline specification\n"))
  cat(sprintf("  SNV_ONLY         = %s\n", SNV_ONLY))
  cat(sprintf("  DROP_PPE_CAL     = %s  (calibration step)\n", DROP_PPE_CAL))
  cat(sprintf("  DROP_PPE_APPLY   = %s  (application step)\n", DROP_PPE_APPLY))
  cat(sprintf("  M0_ONLY          = %s\n", M0_ONLY))
  cat(sprintf("  UNIVERSE         = %s\n", UNIVERSE))
  cat(sprintf("  MAF_BINNING      = %s\n", MAF_BINNING))
  cat(sprintf("  MIN_EVAL_PAIRS   = %d\n", MIN_EVAL_PAIRS))
  cat(sprintf("  MAX_PROP_BOTH0   = %.2f\n", MAX_PROP_BOTH0))
  cat(sprintf("  N_BOOTSTRAP_LEX  = %d  (Stage 1 selection-stability bootstrap)\n",
              N_BOOTSTRAP_LEX_SEL))
  cat(sprintf("  Calibration tag  = %s\n", TAG$cal_tag))
  cat(sprintf("  Application tag  = %s\n", TAG$apply_tag))
  cat(sprintf("  LCA tag          = %s\n", TAG$lca_tag))
  is_canonical <- SNV_ONLY && !DROP_PPE_CAL && !DROP_PPE_APPLY &&
                  identical(UNIVERSE, "C") && identical(MAF_BINNING, "default")
  cat(sprintf("  Spec status      = %s\n",
              if (is_canonical) "CANONICAL (locked Paper 1 spec)" else "SENSITIVITY"))
  cat(strrep("=", 78), "\n", sep = "")
}
