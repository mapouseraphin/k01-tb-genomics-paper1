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
#   PATHS$meta/gh_meta.rds                                (sample-level metadata)
#   PATHS$thresholded/gh_variants_thr_looser.rds          (thresholded calls;
#   PATHS$thresholded/gh_variants_thr_primary.rds          written by
#   PATHS$thresholded/gh_variants_thr_tighter.rds          03_apply_thresholds.R)
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
#   gh_variants_thr_*.rds must exist (run 03_apply_thresholds.R first).
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
  COVARIATES <- c("HIV", "Sex", "age_band", "Fever", "Hemoptysis", "lineage")
}
if (!exists("AGE_BANDS", envir = .GlobalEnv, inherits = FALSE)) {
  # Binary banding (locked 2026-05-05): <25 vs >=25. The Ghana cohort skews
  # young; a 5-band scheme produced unstable strata. Override by setting
  # AGE_BANDS in .GlobalEnv before sourcing the script.
  AGE_BANDS <- c(0, 25, Inf)
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
#'
#' Builds reader-friendly band labels directly from the breaks vector, rather
#' than parsing cut() level strings (the previous implementation had a regex
#' bug in TRE: the character class `[\\[\\]\\(\\)]` did not reliably strip
#' brackets, leaving stray characters that propagated to "NA-NA" labels).
#'
#' Label conventions:
#'   - First band, when the lower bound is 0 or -Inf:   "<hi"
#'   - Last band, when the upper bound is Inf:           ">=lo"
#'   - Otherwise (closed-open interval [lo, hi)):        "lo-(hi-1)"
#' "(hi-1)" is appropriate for integer-valued covariates (e.g., Age in years);
#' for non-integer covariates the same labels still parse unambiguously.
band_continuous <- function(x, breaks, prefix = "") {
  if (!is.numeric(x)) return(as.character(x))

  n_bands <- length(breaks) - 1L
  if (n_bands < 1L) stop("band_continuous: need at least 2 breaks.")

  band_labels <- character(n_bands)
  for (i in seq_len(n_bands)) {
    lo <- breaks[i]
    hi <- breaks[i + 1L]
    is_first <- (i == 1L)
    is_last  <- (i == n_bands)

    band_labels[i] <-
      if (is_last && is.infinite(hi)) {
        sprintf("%s>=%g", prefix, lo)
      } else if (is_first && (is.infinite(lo) || lo == 0)) {
        sprintf("%s<%g", prefix, hi)
      } else {
        sprintf("%s%g-%g", prefix, lo, hi - 1)
      }
  }

  bands <- cut(x, breaks = breaks, right = FALSE, include.lowest = TRUE,
               labels = band_labels)
  as.character(bands)
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

# Read all three tier call tables. Source: 03_apply_thresholds.R writes
# gh_variants_thr_<tier>.rds to PATHS$thresholded for each tier.
tier_files <- list(
  Looser  = file.path(PATHS$thresholded, "gh_variants_thr_looser.rds"),
  Primary = file.path(PATHS$thresholded, "gh_variants_thr_primary.rds"),
  Tighter = file.path(PATHS$thresholded, "gh_variants_thr_tighter.rds")
)

# A given tier may be NULL if 02b found no axis-wise dominant candidate
# (see R/02b_calibration_bootstrap.R). In that case the corresponding
# gh_variants_thr_<tier>.rds is absent; we skip the tier silently and
# proceed with the rest. Primary is mandatory.
tier_present <- vapply(tier_files, file.exists, logical(1))
if (!tier_present["Primary"]) {
  stop("application_per_patient_prevalence: required Primary tier file missing: ",
       tier_files$Primary)
}
tier_files <- tier_files[tier_present]
tier_calls <- lapply(tier_files, readRDS)

if (any(!tier_present)) {
  cat(sprintf("  NOTE: tier(s) absent (no axis-wise dominant candidate): %s\n",
              paste(names(tier_present)[!tier_present], collapse = ", ")))
}

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
