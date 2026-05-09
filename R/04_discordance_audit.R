# =============================================================================
# 04_discordance_audit.R
#
# PRE-INFERENCE GATING ANALYSIS for Bayesian misclassification model.
#
# Goal: For each calibrated tier (looser, primary, tighter), compute the
# per-patient within-patient pattern of the detection indicator Y_ij and
# determine whether observed discordance is sufficient to identify the
# detection-correction parameters in Stage 2.
#
# Per-sample indicator:
#   Y_ij = 1 if sample i of patient j has >= 1 variant call passing tier t
#   Y_ij = 0 otherwise
#
# Per-patient pattern (multi-sample patients only):
#   all_zero : sum_i Y_ij == 0
#   all_one  : sum_i Y_ij == n_j
#   mixed    : 0 < sum_i Y_ij < n_j   <-- informative for theta
#
# Singleton patients (n_j == 1) are tabulated separately and excluded from
# the discordance fraction (uninformative for replicate-based identification).
#
# Gating rule on discordant_fraction = n_mixed / n_eligible:
#   < 0.05            -> uninformative; report naive prevalence only
#   0.05 - 0.20       -> weakly informative; proceed with caution
#   >= 0.20           -> well-powered; proceed with default specification
#
# Inputs:
#   PATHS$meta/gh_meta.rds
#   PATHS$thresholded/gh_variants_thr_{looser,primary,tighter}.rds
#
# Outputs:
#   PATHS$thresholded/discordance_audit_summary_<tag>.csv  (tier-level)
#   PATHS$thresholded/discordance_audit_<tag>.rds          (full structure)
#   Console summary of pattern counts and gating decision per tier.
# =============================================================================

source("R/00_pipeline_config.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(purrr)
})

cat(sprintf("[04_discordance_audit] Spec tag: %s\n\n", TAG$apply_tag))

# ---- Inputs ----------------------------------------------------------------
gh_meta <- readRDS(file.path(PATHS$meta, "gh_meta.rds"))

tier_files <- list(
  Looser  = file.path(PATHS$thresholded, "gh_variants_thr_looser.rds"),
  Primary = file.path(PATHS$thresholded, "gh_variants_thr_primary.rds"),
  Tighter = file.path(PATHS$thresholded, "gh_variants_thr_tighter.rds")
)
tier_present <- vapply(tier_files, file.exists, logical(1))
if (!tier_present["Primary"])
  stop("Primary tier file missing: ", tier_files$Primary)
tier_files <- tier_files[tier_present]
tier_calls <- lapply(tier_files, readRDS)

cat("Tiers present:", paste(names(tier_files), collapse = ", "), "\n")

# ---- Resolve depth field defensively ---------------------------------------
depth_field <- intersect(c("depth_med", "coverage_median"),
                         names(gh_meta))[1]
if (is.na(depth_field))
  stop("No depth field in gh_meta (expected depth_med or coverage_median).")
cat(sprintf("Using depth field: %s\n", depth_field))

# ---- Restrict to baseline (M0) samples -------------------------------------
# Bayesian model is on the baseline cross-section. M1/M2 samples are out
# of scope for prevalence-of-baseline-diversity inference.
if ("visit" %in% names(gh_meta)) {
  m0 <- gh_meta %>% filter(visit == "M0")
  cat(sprintf("Restricted to baseline (M0): %d samples / %d patients\n\n",
              nrow(m0), n_distinct(m0$patientId)))
} else {
  m0 <- gh_meta
  warning("No visit field; using all samples in gh_meta.")
  cat(sprintf("Using all samples: %d samples / %d patients\n\n",
              nrow(m0), n_distinct(m0$patientId)))
}

# ---- Per-tier audit --------------------------------------------------------
audit_one_tier <- function(tier_name, calls, meta) {
  detected_samples <- unique(calls$Sample)

  per_sample <- meta %>%
    transmute(
      Sample, patientId,
      depth     = suppressWarnings(as.numeric(.data[[depth_field]])),
      log_depth = log(suppressWarnings(as.numeric(.data[[depth_field]]))),
      Y         = as.integer(Sample %in% detected_samples)
    )

  per_patient <- per_sample %>%
    group_by(patientId) %>%
    summarise(
      n_samples      = n(),
      n_positive     = sum(Y),
      mean_log_depth = mean(log_depth, na.rm = TRUE),
      min_depth      = min(depth, na.rm = TRUE),
      max_depth      = max(depth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      pattern = case_when(
        n_samples  == 1               ~ "singleton",
        n_positive == 0               ~ "all_zero",
        n_positive == n_samples       ~ "all_one",
        TRUE                          ~ "mixed"
      )
    )

  # Eligibility: only multi-sample patients identify theta via discordance
  eligible <- per_patient %>% filter(n_samples >= 2)

  pattern_counts <- eligible %>%
    count(pattern, name = "n_patients") %>%
    complete(pattern = c("all_zero", "mixed", "all_one"),
             fill = list(n_patients = 0))

  n_eligible  <- nrow(eligible)
  n_mixed     <- pattern_counts$n_patients[pattern_counts$pattern == "mixed"]
  n_all_zero  <- pattern_counts$n_patients[pattern_counts$pattern == "all_zero"]
  n_all_one   <- pattern_counts$n_patients[pattern_counts$pattern == "all_one"]
  n_singleton <- sum(per_patient$n_samples == 1)

  discordant_frac <- if (n_eligible == 0) NA_real_ else n_mixed / n_eligible

  gate <- case_when(
    is.na(discordant_frac)        ~ "no_eligible",
    discordant_frac <  0.05       ~ "uninformative",
    discordant_frac <  0.20       ~ "weakly_informative",
    TRUE                          ~ "well_powered"
  )

  # Diagnostic: depth distribution within mixed-pattern patients
  # (low-depth driving discordance is the expected, model-friendly pattern)
  mixed_pids <- eligible$patientId[eligible$pattern == "mixed"]
  depth_within_mixed <- per_sample %>%
    filter(patientId %in% mixed_pids) %>%
    group_by(Y) %>%
    summarise(
      n            = n(),
      median_depth = median(depth, na.rm = TRUE),
      q25_depth    = quantile(depth, 0.25, na.rm = TRUE),
      q75_depth    = quantile(depth, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(tier = tier_name)

  list(
    summary = tibble(
      tier             = tier_name,
      n_patients_total = nrow(per_patient),
      n_singleton      = n_singleton,
      n_eligible       = n_eligible,
      n_all_zero       = n_all_zero,
      n_mixed          = n_mixed,
      n_all_one        = n_all_one,
      naive_prev       = (n_all_one + n_mixed) / max(n_eligible, 1L),
      discordant_frac  = discordant_frac,
      gate             = gate
    ),
    per_sample         = per_sample %>% mutate(tier = tier_name),
    per_patient        = per_patient %>% mutate(tier = tier_name),
    depth_within_mixed = depth_within_mixed
  )
}

audits <- imap(tier_calls, ~ audit_one_tier(.y, .x, m0))

# ---- Combine and write -----------------------------------------------------
summary_tbl         <- bind_rows(map(audits, "summary"))
per_patient_tbl     <- bind_rows(map(audits, "per_patient"))
per_sample_tbl      <- bind_rows(map(audits, "per_sample"))
depth_within_mixed  <- bind_rows(map(audits, "depth_within_mixed"))

out_summary <- file.path(PATHS$thresholded,
                         sprintf("discordance_audit_summary_%s.csv",
                                 TAG$apply_tag))
out_full    <- file.path(PATHS$thresholded,
                         sprintf("discordance_audit_%s.rds",
                                 TAG$apply_tag))
write_csv(summary_tbl, out_summary)
saveRDS(list(summary            = summary_tbl,
             per_patient        = per_patient_tbl,
             per_sample         = per_sample_tbl,
             depth_within_mixed = depth_within_mixed),
        out_full)

# ---- Console report --------------------------------------------------------
cat("\n", strrep("=", 78), "\n", sep = "")
cat("DISCORDANCE AUDIT SUMMARY (per tier)\n")
cat(strrep("=", 78), "\n", sep = "")
print(as.data.frame(summary_tbl), row.names = FALSE, digits = 3)

cat("\nDEPTH WITHIN MIXED-PATTERN PATIENTS (Y=0 vs Y=1 samples within same patient):\n")
print(as.data.frame(depth_within_mixed), row.names = FALSE, digits = 3)

cat("\n", strrep("-", 78), "\n", sep = "")
cat("GATING DECISION\n")
cat(strrep("-", 78), "\n", sep = "")
for (i in seq_len(nrow(summary_tbl))) {
  s <- summary_tbl[i, ]
  cat(sprintf("  %-7s: discordant fraction = %.3f (%d/%d eligible) -> %s\n",
              s$tier, s$discordant_frac, s$n_mixed, s$n_eligible, s$gate))
}

cat("\nGate codes:\n")
cat("  uninformative      : <5% mixed   -> report naive prevalence; misclassification\n")
cat("                                       model adds nothing\n")
cat("  weakly_informative : 5-20% mixed -> proceed; expect wide CrIs on theta\n")
cat("  well_powered       : >=20% mixed -> proceed with default specification\n")

cat(sprintf("\nOutputs:\n  %s\n  %s\n", out_summary, out_full))
