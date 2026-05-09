# =============================================================================
# 00_pipeline_config.R
#
# Single source of truth for canonical specification flags and tag construction.
# Sourced at the top of every pipeline script (00 -> 06) so a single edit here
# propagates through the entire run.
#
# Canonical specification for Paper 1 (LOCKED 2026-05-03):
#   - iSNV definition: SNVs only (indels excluded)
#   - PE/PPE: retained in primary; excluded as Sens A (full pipeline rerun)
#   - FILTER: PASS-only (universal in input; Sens B retired 2026-05-03)
#   - MAF binning: 4 bins [0.02,0.05), [0.05,0.10), [0.10,0.25), [0.25,0.45]
#                  Sens C (merged_low) retired -- bin-drop policy makes it moot
#   - Calibration set: M0 within-visit replicate pairs only
#   - Calibration selection: lex (canonical) and Pareto (parallel sensitivity)
#   - Stage 2: bootstrap stabilization (B = 1000)
#   - Three-tier ladder: looser / primary / tighter (no sentinel)
#   - Guards: MIN_EVAL_PAIRS = 5, MAX_PROP_BOTH0 = 0.90
#
# To run a sensitivity variant: flip the relevant flag below, re-run the
# pipeline. Output directories carry the active spec tag automatically;
# previous runs are not overwritten. Canonical run keeps the original tag
# `snv_only_ppe_retained_lex`; sensitivities append disambiguating suffixes.
# =============================================================================

# ---- Canonical flags --------------------------------------------------------
# CRITICAL: these defaults reproduce the locked Paper 1 results. Changing any
# default below changes the headline run.
#
# IMPORTANT: each flag is set ONLY IF NOT ALREADY DEFINED. This lets
# run_paper1.R override flags before re-sourcing this file (the orchestrator
# pattern is: assign with <<- in .GlobalEnv, then source("R/00_pipeline_config.R")
# to rebuild TAG and PATHS for that spec). Without the `exists()` guard,
# unconditional assignment here would clobber the orchestrator's override.
if (!exists("SNV_ONLY",       envir = .GlobalEnv, inherits = FALSE)) SNV_ONLY       <- TRUE
if (!exists("DROP_PPE_CAL",   envir = .GlobalEnv, inherits = FALSE)) DROP_PPE_CAL   <- TRUE
if (!exists("DROP_PPE_APPLY", envir = .GlobalEnv, inherits = FALSE)) DROP_PPE_APPLY <- TRUE
if (!exists("M0_ONLY",        envir = .GlobalEnv, inherits = FALSE)) M0_ONLY        <- TRUE

# ---- Selection method (added 2026-05-03) -----------------------------------
# Stage 1 calibration selection method. Two-stage workflow:
#   Stage 1 (this method): observed-data primary
#   Stage 2 (02b):         bootstrap stabilization (modal cell)
# "lex"   = lexicographic ranking on six concordance criteria (canonical)
# "pareto"= hybrid Pareto + closest-to-ideal (parallel sensitivity)
if (!exists("SELECTION_METHOD", envir = .GlobalEnv, inherits = FALSE))
  SELECTION_METHOD <- "lex"

# ---- Calibration guards (script 02) -----------------------------------------
MIN_EVAL_PAIRS <- 5L      # canonical stability floor
MAX_PROP_BOTH0 <- 0.95     # canonical ceiling on fraction of zero-zero pairs

# ---- Variant-class restriction grid (script 01) -----------------------------
DP_GRID  <- seq(40, 200, by = 10)
AD1_GRID <- 3:8
MAF_GRID <- c(0.02, 0.03, 0.05, 0.10, 0.15, 0.20)
MAF_MAX  <- 0.50            # symmetric MAF: pmin(AD1, AD2) / DP <= 0.5 by
                            # construction; 0.50 retains essentially all
                            # iSNV-range calls without crossing the 0.5 boundary

# ---- Cohort QC ---------------------------------------------------------------
MIN_COV_INCLUDE <- 50L      # medianCov floor for analytic frame inclusion

# ---- Bootstrap settings ------------------------------------------------------
SEED            <- 20260428L
N_BOOTSTRAP     <- 1000L    # external-validity bootstrap (legacy; kept for
                            # archived 04_*.R scripts; not used in Paper 1)
N_BOOTSTRAP_SEL <- 1000L    # Stage 2 selection-stability bootstrap (script 02b)

# ---- Tag construction --------------------------------------------------------
# Stage 1 (calibration) output depends on SNV_ONLY, DROP_PPE_CAL, and
#   SELECTION_METHOD (so lex and Pareto runs do not collide).
# Stage 2a (apply_thresholds) output depends on SNV_ONLY, DROP_PPE_APPLY,
#   and SELECTION_METHOD (so the calibrated rule that's applied is method-
#   aware; Sens A under lex differs from canonical under lex by PE/PPE).
# Stage 2b (LCA dataset / fits) inherits apply_tag (canonical Universe = C
#   and default MAF binning are the only options after 2026-05-03 retirements).
build_tag <- function(snv_only         = SNV_ONLY,
                      drop_ppe_cal     = DROP_PPE_CAL,
                      drop_ppe_apply   = DROP_PPE_APPLY,
                      selection_method = SELECTION_METHOD) {

  variant_part   <- if (snv_only)       "snv_only"      else "with_indels"
  ppe_cal_part   <- if (drop_ppe_cal)   "ppe_excluded"  else "ppe_retained"
  ppe_apply_part <- if (drop_ppe_apply) "ppe_excluded"  else "ppe_retained"
  method_part    <- selection_method  # "lex" or "pareto"

  cal_tag   <- paste(variant_part, ppe_cal_part,   method_part, sep = "_")
  apply_tag <- paste(variant_part, ppe_apply_part, method_part, sep = "_")
  lca_tag   <- apply_tag

  list(
    variant          = variant_part,
    ppe_cal          = ppe_cal_part,
    ppe_apply        = ppe_apply_part,
    selection_method = selection_method,
    cal_tag          = cal_tag,
    apply_tag        = apply_tag,
    lca_tag          = lca_tag
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
  cross_spec       = file.path("outputs", "cross_spec_lca")  # canonical-vs-
                                                             # SensA LCA
                                                             # comparisons
)

# Create output directories ---------------------------------------------------
for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

# ---- Path overrides (post-build, persistent across re-sources) -------------
# To redirect specific PATHS entries across script re-sources (e.g. point
# PATHS$variants at a sensitivity-specific variant table), set
#   PATHS_OVERRIDES <- list(<key> = <new_path>, ...)
# in .GlobalEnv BEFORE the next source() call. Clear with
#   rm(PATHS_OVERRIDES, envir = .GlobalEnv).
if (exists("PATHS_OVERRIDES", envir = .GlobalEnv, inherits = FALSE) &&
    is.list(PATHS_OVERRIDES) && length(PATHS_OVERRIDES) > 0) {
  for (key in names(PATHS_OVERRIDES)) {
    if (!is.null(PATHS_OVERRIDES[[key]])) {
      PATHS[[key]] <- PATHS_OVERRIDES[[key]]
      message(sprintf("[PATHS override] %s -> %s", key, PATHS_OVERRIDES[[key]]))
      dir.create(PATHS[[key]], recursive = TRUE, showWarnings = FALSE)
    }
  }
}

# ---- Cross-spec path resolver -----------------------------------------------
# Returns the per-spec output paths for an arbitrary spec without touching the
# active TAG. Used by sensitivity-LCA aggregation in cross_spec_sens_a_comparison.R.
build_paths_for_spec <- function(snv_only, drop_ppe_cal, drop_ppe_apply,
                                 selection_method = "lex") {
  tag <- build_tag(snv_only         = snv_only,
                   drop_ppe_cal     = drop_ppe_cal,
                   drop_ppe_apply   = drop_ppe_apply,
                   selection_method = selection_method)
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
  cat(sprintf("  SNV_ONLY          = %s\n", SNV_ONLY))
  cat(sprintf("  DROP_PPE_CAL      = %s  (calibration step)\n", DROP_PPE_CAL))
  cat(sprintf("  DROP_PPE_APPLY    = %s  (application step)\n", DROP_PPE_APPLY))
  cat(sprintf("  M0_ONLY           = %s\n", M0_ONLY))
  cat(sprintf("  SELECTION_METHOD  = %s  (Stage 1 calibration selector)\n",
              SELECTION_METHOD))
  cat(sprintf("  MIN_EVAL_PAIRS    = %d\n", MIN_EVAL_PAIRS))
  cat(sprintf("  MAX_PROP_BOTH0    = %.2f\n", MAX_PROP_BOTH0))
  cat(sprintf("  N_BOOTSTRAP_SEL   = %d  (Stage 2 selection-stability bootstrap)\n",
              N_BOOTSTRAP_SEL))
  cat(sprintf("  Calibration tag   = %s\n", TAG$cal_tag))
  cat(sprintf("  Application tag   = %s\n", TAG$apply_tag))
  cat(sprintf("  LCA tag           = %s\n", TAG$lca_tag))
  is_canonical <- SNV_ONLY && !DROP_PPE_CAL && !DROP_PPE_APPLY &&
                  M0_ONLY && identical(SELECTION_METHOD, "lex")
  cat(sprintf("  Spec status       = %s\n",
              if (is_canonical) "CANONICAL (locked Paper 1 spec)" else "SENSITIVITY"))
  cat(strrep("=", 78), "\n", sep = "")
}
