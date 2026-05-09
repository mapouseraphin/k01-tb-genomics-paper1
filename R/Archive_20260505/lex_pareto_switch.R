# Clean R session at project root.
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
setwd("~/Library/CloudStorage/OneDrive-UniversityofFlorida/Research/manuscripts/Submitted/within_between_host_isnv_calibration")

# 1. Override flags BEFORE sourcing config. The exists() guards in
#    00_pipeline_config.R will preserve these.
SELECTION_METHOD <- "pareto"
DROP_PPE_CAL     <- FALSE
DROP_PPE_APPLY   <- FALSE
SNV_ONLY         <- TRUE
M0_ONLY          <- TRUE

# 2. Source config; TAG and PATHS rebuild for the Pareto tag.
source("R/00_pipeline_config.R")

# 3. Sanity check — confirm the tag is the Pareto one and PATHS is fresh.
print(TAG$cal_tag)            # expect: "snv_only_ppe_retained_pareto"
print(PATHS$calibration)      # expect: data_derived/01_calibration_snv_only_ppe_retained_pareto

# 4. Pareto needs its own cal_pairs.rds at this tag. cal_pairs depends on
#    DROP_PPE_CAL + SNV_ONLY only (not on SELECTION_METHOD), but the per-tag
#    PATHS isolates output, so 01 must rerun to populate the Pareto directory.
source("R/01_build_cal_pairs.R")

# 5. Pareto Stage 1 — closest-to-ideal on the (J, C) frontier.
source("R/02_calibrate_hybrid_pareto.R")

# 6. Stage 2 bootstrap — the script auto-dispatches on SELECTION_METHOD.
source("R/02b_calibration_bootstrap.R")
