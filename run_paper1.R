# =============================================================================
# run_paper1.R
#
# Paper 1 (Microbial Genomics) reproducibility entry point.
#
# Reproduces the calibration and application pipeline presented in the
# manuscript. Deterministic given the seed pinned in 00_pipeline_config.R.
#
# Manuscript primary spec: SNV-only, PE/PPE excluded, M0 within-visit pairs,
# lexicographic threshold selection. Tag: snv_only_ppe_excluded_lex.
#
# Pipeline stages:
#   prep         build_sample_metadata_v2 -> 00_prep_metadata
#   calibrate    01_build_cal_pairs -> 02_calibrate_lexicographic ->
#                02b_calibration_bootstrap -> 02c_depth_concordance_diagnostic
#   apply        03_apply_thresholds
#   figures      09_make_calibration_selection_figure_table   (Fig 2, Table 2)
#   application  application_per_patient_prevalence           (Table 4)
#                application_logistic_regression              (Table 5 logistic)
#                application_persistence_analysis             (Table 5 persistence,
#                                                              Tables S4, S5, Fig 5)
#  
# Usage:
#   Rscript run_paper1.R
#   Rscript run_paper1.R --from application --to application
#   Rscript run_paper1.R --scripts-dir R
#   Rscript run_paper1.R --dry-run
#
# Must be invoked from the repository root. PATHS in 00_pipeline_config.R
# are relative (data_raw, data_derived, outputs) and anchor to the CWD.
# Pipeline scripts may live in the repo root or in R/; auto-detected, or
# specified with --scripts-dir.
#
# To capture stderr (sub-script messages, warnings, errors) in the log file:
#   Rscript run_paper1.R 2>&1 | tee logs/run_paper1_$(date +%Y%m%d_%H%M%S).log
# =============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  }
  library(optparse)
})

# ---- Argument parsing -------------------------------------------------------
opt_list <- list(
  make_option("--from", type = "character", default = "prep",
              help = "First stage: prep, calibrate, apply, figures, application. [default: %default]"),
  make_option("--to", type = "character", default = "application",
              help = "Last stage: prep, calibrate, apply, figures, application. [default: %default]"),
  make_option("--scripts-dir", type = "character", default = NULL,
              help = "Path to pipeline scripts directory. Auto-detected if omitted."),
  make_option("--dry-run", action = "store_true", default = FALSE,
              help = "Print plan without executing."),
  make_option("--log-dir", type = "character", default = "logs",
              help = "Directory for per-run log files. [default: %default]")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# ---- Helpers ----------------------------------------------------------------
detect_scripts_dir <- function() {
  candidates <- c(".", "R", "code", "scripts")
  for (d in candidates) {
    if (file.exists(file.path(d, "00_pipeline_config.R"))) return(d)
  }
  stop("Cannot find 00_pipeline_config.R in ., R/, code/, or scripts/. ",
       "Run from the repo root, or pass --scripts-dir. Current cwd: ", getwd())
}

stage_idx <- function(s, order) {
  i <- match(s, order)
  if (is.na(i)) {
    stop(sprintf("Unknown stage: %s. Valid: %s", s, paste(order, collapse = ", ")))
  }
  i
}

log_msg <- function(...) cat(paste0(..., "\n"))

run_stage <- function(stage_name, scripts_dir, stages, dry_run) {
  for (script in stages[[stage_name]]) {
    src_path <- file.path(scripts_dir, script)
    if (!file.exists(src_path)) {
      stop(sprintf("Missing script: %s (cwd=%s)", src_path, getwd()))
    }
    log_msg(sprintf("[%s] sourcing %s", stage_name, src_path))
    if (dry_run) next
    t0 <- Sys.time()
    source(src_path, echo = FALSE)
    log_msg(sprintf("  done in %s", format(round(Sys.time() - t0, 2))))
  }
}

# ---- Configuration ----------------------------------------------------------
SPEC <- list(
  SNV_ONLY         = TRUE,
  DROP_PPE_CAL     = TRUE,
  DROP_PPE_APPLY   = TRUE,
  M0_ONLY          = TRUE,
  SELECTION_METHOD = "lex"
)

STAGES <- list(
  prep        = c("build_sample_metadata_v2.R",
                  "00_prep_metadata.R"),
  calibrate   = c("01_build_cal_pairs.R",
                  "02_calibrate_lexicographic.R",
                  "02b_calibration_bootstrap.R",
                  "02c_depth_concordance_diagnostic.R"),
  apply       = c("03_apply_thresholds.R"),
  figures     = c("09_make_calibration_selection_figure_table.R"),
  application = c("application_per_patient_prevalence.R",
                  "application_logistic_regression.R",
                  "application_persistence_analysis.R")
)

STAGE_ORDER <- c("prep", "calibrate", "apply", "figures", "application")

# ---- Resolve scripts directory ---------------------------------------------
SCRIPTS_DIR <- if (is.null(opt$`scripts-dir`)) {
  detect_scripts_dir()
} else {
  opt$`scripts-dir`
}
if (!file.exists(file.path(SCRIPTS_DIR, "00_pipeline_config.R"))) {
  stop(sprintf("Scripts directory '%s' does not contain 00_pipeline_config.R",
               SCRIPTS_DIR))
}

# ---- Resolve stage range ----------------------------------------------------
from_i <- stage_idx(opt$from, STAGE_ORDER)
to_i   <- stage_idx(opt$to,   STAGE_ORDER)
if (from_i > to_i) stop("--from stage must precede --to stage")
stages_to_run <- STAGE_ORDER[from_i:to_i]

# ---- Logging ----------------------------------------------------------------
dir.create(opt$`log-dir`, recursive = TRUE, showWarnings = FALSE)
run_id   <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(opt$`log-dir`, sprintf("run_paper1_%s.log", run_id))

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE, type = "output")

on.exit({
  log_msg("---- sessionInfo() ----")
  print(sessionInfo())
  log_msg(sprintf("Run finished: %s", Sys.time()))
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

# ---- Banner -----------------------------------------------------------------
log_msg(sprintf("Run started:  %s", Sys.time()))
log_msg(sprintf("Run ID:       %s", run_id))
log_msg(sprintf("Log file:     %s", log_file))
log_msg(sprintf("Working dir:  %s", getwd()))
log_msg(sprintf("Scripts dir:  %s", SCRIPTS_DIR))
log_msg(sprintf("Stages:       %s", paste(stages_to_run, collapse = " -> ")))
log_msg(sprintf("Dry run:      %s", opt$`dry-run`))
log_msg("")

# ---- Apply spec overrides and load config -----------------------------------
log_msg("[spec] Applying manuscript primary spec overrides:")
for (k in names(SPEC)) {
  assign(k, SPEC[[k]], envir = .GlobalEnv)
  log_msg(sprintf("  %s = %s", k, format(SPEC[[k]])))
}
source(file.path(SCRIPTS_DIR, "00_pipeline_config.R"), local = FALSE)
if (exists("banner", envir = .GlobalEnv)) banner()

# ---- Execute ----------------------------------------------------------------
for (stg in stages_to_run) {
  log_msg(sprintf("\n==== STAGE: %s ====", stg))
  run_stage(stg, SCRIPTS_DIR, STAGES, opt$`dry-run`)
}

log_msg("\n==== Paper 1 pipeline complete. ====")
