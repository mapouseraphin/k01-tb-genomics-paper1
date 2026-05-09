# =============================================================================
# application_per_patient_prevalence.R   (NEW 2026-05-04)
#
# Per-patient iSNV detection prevalence under the calibrated rule, stratified
# by configurable epidemiologic covariates. Reports prevalence at all three
# tiers (Looser, Primary, Tighter) with Wilson 95% binomial confidence
# intervals.
#
#
# CONFIG
# ------
# COVARIATES (vector of column names in gh_meta): which patient-level
#   covariates to stratify by. The script will produce a row per (covariate,
#   level, tier) combination plus an overall row per tier. Default covariates
#   are HIV, sex, age (banded), and baseline smear grade if present; the
#   list can be expanded by passing a different vector before sourcing.
#   Covariates must be discrete (factor or character). Continuous covariates
#   (e.g. age, year of enrollment) are auto-banded; banding is configurable
#   via the AGE_BANDS / CONT_BANDS arguments to band_covariate().
#
# UNIT OF ANALYSIS
# ----------------
# Per patient: a patient is "detected at tier t" if at least one of their
# samples (across all timepoints M0/M1/M2 in the analytic cohort) has at
# least one variant call passing tier t. Patients with no samples are
# excluded (none expected; gh_meta is the analytic cohort).
#
# INPUTS (per spec)
# -----------------
#   PATHS$meta/gh_meta.rds                          (sample-level metadata)
#   PATHS$thresholded/gh_variants_thr_looser.rds
#   PATHS$thresholded/gh_variants_thr_primary.rds
#   PATHS$thresholded/gh_variants_thr_tighter.rds
#
# OUTPUTS (per spec, written to PATHS$tables)
# -------------------------------------------
#   table4_per_patient_prevalence_<spec_tag>.csv    (long-format)
#   table4_per_patient_prevalence_<spec_tag>.tex    (paper-ready)
#   per_patient_detection_<spec_tag>.rds            (patient-level dataset
#                                                     for downstream use,
#                                                     e.g. shared-iSNV
#                                                     analysis)
#
# DEPENDENCIES
# ------------
#   00_pipeline_config.R must be sourced first (defines PATHS, TAG).
#   gh_calls_*.rds must exist (run 03_apply_thresholds.R first).
#
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

stopifnot(exists("PATHS"), exists("TAG"))

cat(sprintf(
  "[application_per_patient_prevalence] Spec tag: %s\nOutput dir: %s\n\n",
  TAG$apply_tag, PATHS$tables))

# ---- Configurable covariate list (defaults; override by setting before source) ----
if (!exists("COVARIATES", envir = .GlobalEnv, inherits = FALSE)) {
  COVARIATES <- c("HIV", "Sex", "age_band", "smear_grade")
}
if (!exists("AGE_BANDS", envir = .GlobalEnv, inherits = FALSE)) {
  AGE_BANDS <- c(0, 25, 35, 45, 55, Inf)  # 5 bands; final = "55+"
}
if (!exists("WILSON_CONF_LEVEL", envir = .GlobalEnv, inherits = FALSE)) {
  WILSON_CONF_LEVEL <- 0.95
}

# ---- Helpers ---------------------------------------------------------------

#' Wilson score interval for a binomial proportion.
#' Returns a tibble with columns prev, lo, hi.
wilson_ci <- function(x, n, conf_level = WILSON_CONF_LEVEL) {
  if (n == 0) return(tibble(prev = NA_real_, lo = NA_real_, hi = NA_real_))
  z <- qnorm(1 - (1 - conf_level) / 2)
  p <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p * (1 - p) + z^2 / (4 * n)) / n) / denom
  tibble(prev = p, lo = max(0, centre - margin), hi = min(1, centre + margin))
}

#' Auto-band a continuous covariate.
band_continuous <- function(x, breaks, prefix = "") {
  if (!is.numeric(x)) return(as.character(x))
  bands <- cut(x, breaks = breaks, right = FALSE, include.lowest = TRUE)
  # Re-label to be reader-friendly
  labels <- levels(bands)
  pretty <- sapply(labels, function(lbl) {
    bounds <- as.numeric(strsplit(gsub("[\\[\\]\\(\\)]", "", lbl), ",")[[1]])
    if (is.infinite(bounds[2])) sprintf("%s%g+", prefix, bounds[1])
    else                        sprintf("%s%g-%g", prefix, bounds[1], bounds[2] - 1)
  })
  pretty[as.integer(bands)]
}

#' Coerce a covariate column for stratification. Returns character vector.
prepare_covariate <- function(meta_df, covar_name) {
  if (!covar_name %in% names(meta_df)) {
    warning(sprintf("Covariate '%s' not found in gh_meta — skipping.",
                    covar_name))
    return(NULL)
  }
  x <- meta_df[[covar_name]]
  # Auto-band continuous age
  if (covar_name == "age_band" && !"age_band" %in% names(meta_df) &&
      "Age" %in% names(meta_df)) {
    x <- band_continuous(meta_df$Age, AGE_BANDS, prefix = "")
  } else if (covar_name == "age_band" && !is.null(x) && is.numeric(x)) {
    x <- band_continuous(x, AGE_BANDS, prefix = "")
  } else if (is.numeric(x)) {
    # Catch-all for any other continuous covariate
    cat(sprintf(
      "  NOTE: covariate '%s' is numeric; auto-banding into quartiles.\n",
      covar_name))
    qbreaks <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    x <- band_continuous(x, qbreaks)
  } else if (is.factor(x)) {
    x <- as.character(x)
  }
  x[is.na(x)] <- "Unknown"
  x
}

# ---- Step 1: Load metadata + thresholded calls -----------------------------

gh_meta <- readRDS(file.path(PATHS$meta, "gh_meta.rds"))

# Pre-aggregate to one row per patient with all covariates
patient_tbl <- gh_meta |>
  distinct(patientId, .keep_all = TRUE)  # one row per patient (covariates
                                         # must be patient-invariant; if not,
                                         # this takes the first row)

# Add age_band if Age present
if ("Age" %in% names(patient_tbl) && !"age_band" %in% names(patient_tbl)) {
  patient_tbl$age_band <- band_continuous(patient_tbl$Age, AGE_BANDS)
}

cat(sprintf("  Cohort: %d patients across %d samples\n\n",
            n_distinct(patient_tbl$patientId), nrow(gh_meta)))

# Read all three tier call tables
tier_files <- list(
  Looser  = file.path(PATHS$thresholded, "gh_variants_thr_looser.rds"),
  Primary = file.path(PATHS$thresholded, "gh_variants_thr_primary.rds"),
  Tighter = file.path(PATHS$thresholded, "gh_variants_thr_tighter.rds")
)

# Defensive: catch empty paths AND missing files (avoid vacuous-all bug)
stopifnot(
  length(unlist(tier_files)) == 3L,
  all(nzchar(unlist(tier_files))),
  all(file.exists(unlist(tier_files)))
)

tier_calls <- lapply(tier_files, readRDS)

cat(sprintf("  Calls per tier (sample-level rows):\n"))
for (t in names(tier_calls)) {
  cat(sprintf("    %s: %d calls across %d samples\n",
              t, nrow(tier_calls[[t]]),
              n_distinct(tier_calls[[t]]$Sample)))
}
cat("\n")

# ---- Step 2: Compute per-patient detection at each tier --------------------

# A sample is "detected at tier t" if it has >= 1 call in the tier's call table.
# A patient is "detected at tier t" if any of their samples are detected.

per_patient_detection <- patient_tbl |>
  select(patientId, all_of(intersect(names(patient_tbl), COVARIATES))) |>
  arrange(patientId)

for (t in names(tier_calls)) {
  detected_samples <- unique(tier_calls[[t]]$Sample)
  detected_patients <- gh_meta |>
    filter(Sample %in% detected_samples) |>
    pull(patientId) |>
    unique()
  per_patient_detection[[paste0("detected_", t)]] <-
    per_patient_detection$patientId %in% detected_patients
}

# Save per-patient dataset for downstream (shared-iSNV analysis, etc.)
out_rds <- file.path(PATHS$tables,
                     sprintf("per_patient_detection_%s.rds", TAG$apply_tag))
saveRDS(per_patient_detection, out_rds)
cat(sprintf("  Per-patient dataset saved to: %s\n\n", out_rds))

# ---- Step 3: Compute prevalence by stratum ---------------------------------

compute_stratum_row <- function(df, tier_name, stratum_label, stratum_value) {
  detected_col <- paste0("detected_", tier_name)
  n_total <- nrow(df)
  n_det   <- sum(df[[detected_col]], na.rm = TRUE)
  ci      <- wilson_ci(n_det, n_total)
  tibble(
    tier            = tier_name,
    covariate       = stratum_label,
    level           = stratum_value,
    n_patients      = n_total,
    n_detected      = n_det,
    prevalence      = ci$prev,
    ci_lo           = ci$lo,
    ci_hi           = ci$hi
  )
}

results <- list()

# Overall row per tier
for (t in names(tier_calls)) {
  results[[length(results) + 1]] <- compute_stratum_row(
    per_patient_detection, t, "Overall", "All patients")
}

# Stratified rows per (covariate, level, tier)
for (covar in COVARIATES) {
  if (!covar %in% names(per_patient_detection)) {
    # Allow age_band to be derived from Age if not directly present
    if (covar == "age_band" && "Age" %in% names(patient_tbl)) {
      per_patient_detection$age_band <- band_continuous(
        patient_tbl$Age[match(per_patient_detection$patientId,
                              patient_tbl$patientId)],
        AGE_BANDS)
    } else {
      next
    }
  }
  per_patient_detection[[covar]] <- prepare_covariate(per_patient_detection,
                                                      covar)
  levs <- sort(unique(per_patient_detection[[covar]]))
  for (lev in levs) {
    sub <- per_patient_detection[per_patient_detection[[covar]] == lev, ,
                                 drop = FALSE]
    for (t in names(tier_calls)) {
      results[[length(results) + 1]] <- compute_stratum_row(
        sub, t, covar, lev)
    }
  }
}

table4 <- bind_rows(results)

# ---- Step 4: Format and write outputs --------------------------------------

# Long-format CSV
out_csv <- file.path(PATHS$tables,
                     sprintf("table4_per_patient_prevalence_%s.csv",
                             TAG$apply_tag))
write_csv(table4, out_csv)

# Paper-ready wide format: one row per (covariate, level), columns per tier
table4_wide <- table4 |>
  mutate(prevalence_pct = sprintf("%.1f", 100 * prevalence),
         ci_str         = sprintf("(%.1f, %.1f)", 100 * ci_lo, 100 * ci_hi)) |>
  unite("prev_with_ci", prevalence_pct, ci_str, sep = " ", remove = FALSE) |>
  select(covariate, level, tier, n_patients, n_detected, prev_with_ci) |>
  pivot_wider(names_from = tier,
              values_from = c(n_detected, prev_with_ci)) |>
  arrange(covariate, level)

out_csv_wide <- file.path(PATHS$tables,
                          sprintf("table4_wide_%s.csv", TAG$apply_tag))
write_csv(table4_wide, out_csv_wide)

cat(sprintf("[application_per_patient_prevalence] Outputs:\n"))
cat(sprintf("  Long-format: %s\n", out_csv))
cat(sprintf("  Wide-format: %s\n", out_csv_wide))

# Brief printout
cat("\nOverall per-patient detection prevalence:\n")
overall <- table4 |> filter(covariate == "Overall")
print(overall |> select(tier, n_detected, n_patients, prevalence,
                        ci_lo, ci_hi) |> as.data.frame())

cat("\nBy HIV status (Primary tier):\n")
hiv_primary <- table4 |>
  filter(covariate == "HIV", tier == "Primary")
print(hiv_primary |> select(level, n_detected, n_patients, prevalence,
                            ci_lo, ci_hi) |> as.data.frame())

invisible(table4)
