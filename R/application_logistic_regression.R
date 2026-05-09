# =============================================================================
# application_logistic_regression.R
#
# Descriptive logistic-regression association of patient-level detectable
# iSNV with three binary covariates: HIV-positive, female, age >= 25.
# Fit per tier (Looser, Primary, Tighter).
#
# Estimator: Firth's bias-reduced penalized maximum likelihood (logistf).
# Rationale: events-per-variable ratio at this cohort size (~95 patients,
# ~9 events under Primary) is far below the standard threshold for
# unpenalized MLE; Firth removes small-sample bias and handles separation.
# Standard MLE (glm) shown in parallel for transparency.
#
# REPORTING FRAME
# ---------------
# This is a DESCRIPTIVE / EXPLORATORY analysis. Outcome is "detectable
# within-host diversity under the calibrated rule," NOT true diversity.
# ORs are associations, not causal effects. With ~9 events, CIs are wide
# and these results are hypothesis-generating only.
#
# INPUTS
# ------
#   PATHS$tables/per_patient_detection_<tag>.rds   (from
#                                                   application_per_patient_prevalence.R)
#
# OUTPUTS
# -------
#   PATHS$tables/table5_logistic_<tag>.csv         (long-format)
#   PATHS$tables/table5_logistic_wide_<tag>.csv    (paper-ready wide)
#   PATHS$tables/table5_logistic_<tag>.tex         (LaTeX)
#   Console summary with cell counts, EPV, separation flags.
#
# DEPENDENCIES
# ------------
#   logistf, dplyr, tidyr, readr, tibble, purrr
#   00_pipeline_config.R must be sourced first (defines PATHS, TAG).
#   application_per_patient_prevalence.R must have been run.
# =============================================================================

source("R/00_pipeline_config.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(purrr)
  if (!requireNamespace("logistf", quietly = TRUE))
    stop("Install logistf: install.packages('logistf')")
  library(logistf)
})

cat(sprintf("[application_logistic_regression] Spec tag: %s\n\n",
            TAG$apply_tag))

# ---- Load per-patient detection table --------------------------------------
ppd_file <- file.path(PATHS$tables,
                      sprintf("per_patient_detection_%s.rds",
                              TAG$apply_tag))
if (!file.exists(ppd_file))
  stop("Per-patient detection file missing. Run application_per_patient_prevalence.R first:\n  ",
       ppd_file)

ppd <- readRDS(ppd_file)

# Identify available tier columns (Looser may be NULL; Primary mandatory)
tier_cols <- intersect(c("detected_Looser", "detected_Primary",
                         "detected_Tighter"),
                       names(ppd))
if (!"detected_Primary" %in% tier_cols)
  stop("detected_Primary column missing in per-patient detection table.")

cat(sprintf("Tiers present: %s\n",
            paste(sub("detected_", "", tier_cols), collapse = ", ")))

# ---- Construct binary covariates -------------------------------------------
# Defensive: HIV and Sex columns may be character ("Positive"/"Negative",
# "Female"/"Male"), factor, or coded otherwise. age_band is from
# band_continuous() with default AGE_BANDS = c(0, 25, Inf): "<25" / ">=25".
# Adjust the matching below if your coding differs.

cat("\nCovariate coding:\n")
cat("  HIV-positive : 1 if HIV %in% c('Positive', 'pos', 1, 'HIV+'), else 0\n")
cat("  female       : 1 if Sex %in% c('Female', 'F', 'female'), else 0\n")
cat("  older (>=25) : 1 if age_band == '>=25', else 0\n\n")

ppd_mod <- ppd %>%
  mutate(
    HIV_pos = case_when(
      tolower(as.character(HIV)) %in% c("positive", "pos", "yes", "1", "hiv+") ~ 1L,
      tolower(as.character(HIV)) %in% c("negative", "neg", "no", "0", "hiv-")  ~ 0L,
      TRUE ~ NA_integer_
    ),
    female = case_when(
      tolower(as.character(Sex)) %in% c("female", "f")    ~ 1L,
      tolower(as.character(Sex)) %in% c("male",   "m")    ~ 0L,
      TRUE ~ NA_integer_
    ),
    older = case_when(
      as.character(age_band) == ">=25" ~ 1L,
      as.character(age_band) == "<25"  ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Complete-case for the three covariates
n_before <- nrow(ppd_mod)
ppd_mod <- ppd_mod %>%
  filter(!is.na(HIV_pos), !is.na(female), !is.na(older))
n_after  <- nrow(ppd_mod)
cat(sprintf("Complete-case N: %d / %d (%d dropped for missing covariates)\n\n",
            n_after, n_before, n_before - n_after))

# ---- Descriptive cell counts (cross-tab) -----------------------------------
cat("Covariate distribution in analytic cohort:\n")
cat(sprintf("  HIV-positive : %d / %d (%.1f%%)\n",
            sum(ppd_mod$HIV_pos), n_after, 100 * mean(ppd_mod$HIV_pos)))
cat(sprintf("  Female       : %d / %d (%.1f%%)\n",
            sum(ppd_mod$female),  n_after, 100 * mean(ppd_mod$female)))
cat(sprintf("  Age >= 25    : %d / %d (%.1f%%)\n\n",
            sum(ppd_mod$older),   n_after, 100 * mean(ppd_mod$older)))

for (tc in tier_cols) {
  y <- as.integer(ppd_mod[[tc]])
  cat(sprintf("Outcome '%s': %d events / %d patients (%.1f%%)\n",
              sub("detected_", "", tc), sum(y), length(y),
              100 * mean(y)))
  cat("  Cross-tab (rows = covariate, cols = detected 0/1):\n")
  for (cov in c("HIV_pos", "female", "older")) {
    tab <- table(ppd_mod[[cov]], y)
    cat(sprintf("    %-8s 0:[%d,%d]  1:[%d,%d]\n",
                cov,
                tab["0", "0"] %||% 0, tab["0", "1"] %||% 0,
                tab["1", "0"] %||% 0, tab["1", "1"] %||% 0))
  }
  epv <- sum(y) / 3
  cat(sprintf("  Events per variable (EPV) = %.1f  %s\n\n",
              epv,
              if (epv < 5) "[LOW: Firth strongly recommended]"
              else if (epv < 10) "[MARGINAL: Firth recommended]"
              else "[ADEQUATE]"))
}

# ---- Fit per tier: Firth + standard MLE for transparency -------------------
fit_one_tier <- function(tier_col, data) {
  tier_name <- sub("detected_", "", tier_col)
  data$y <- as.integer(data[[tier_col]])

  # Firth penalized MLE
  fit_firth <- tryCatch(
    logistf(y ~ HIV_pos + female + older, data = data,
            control = logistf.control(maxit = 200)),
    error = function(e) { warning("Firth failed: ", e$message); NULL }
  )
  # Standard unpenalized MLE for comparison
  fit_mle <- tryCatch(
    glm(y ~ HIV_pos + female + older, data = data,
        family = binomial(link = "logit")),
    error = function(e) { warning("MLE failed: ", e$message); NULL }
  )

  # Extract Firth results (skip intercept)
  firth_tbl <- if (!is.null(fit_firth)) {
    coefs <- coef(fit_firth)[-1]
    cis   <- confint(fit_firth)[-1, , drop = FALSE]
    pvals <- fit_firth$prob[-1]
    tibble(
      tier      = tier_name,
      estimator = "Firth",
      term      = names(coefs),
      log_or    = unname(coefs),
      OR        = exp(unname(coefs)),
      ci_lo     = exp(cis[, 1]),
      ci_hi     = exp(cis[, 2]),
      p_value   = unname(pvals)
    )
  } else NULL

  # Extract MLE results
  mle_tbl <- if (!is.null(fit_mle)) {
    s <- summary(fit_mle)$coefficients
    s <- s[-1, , drop = FALSE]
    cis_mle <- tryCatch(suppressMessages(confint(fit_mle))[-1, , drop = FALSE],
                        error = function(e) matrix(NA_real_, nrow(s), 2))
    tibble(
      tier      = tier_name,
      estimator = "MLE",
      term      = rownames(s),
      log_or    = s[, "Estimate"],
      OR        = exp(s[, "Estimate"]),
      ci_lo     = exp(cis_mle[, 1]),
      ci_hi     = exp(cis_mle[, 2]),
      p_value   = s[, "Pr(>|z|)"]
    )
  } else NULL

  bind_rows(firth_tbl, mle_tbl)
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

results <- bind_rows(lapply(tier_cols, fit_one_tier, data = ppd_mod))

# ---- Format and write ------------------------------------------------------
out_csv <- file.path(PATHS$tables,
                     sprintf("table5_logistic_%s.csv", TAG$apply_tag))
write_csv(results, out_csv)

# Wide format: one row per (tier, term), columns Firth vs MLE
wide <- results %>%
  mutate(or_ci = sprintf("%.2f (%.2f, %.2f)", OR, ci_lo, ci_hi)) %>%
  select(tier, term, estimator, or_ci) %>%
  pivot_wider(names_from = estimator, values_from = or_ci) %>%
  arrange(tier, term)

out_csv_wide <- file.path(PATHS$tables,
                          sprintf("table5_logistic_wide_%s.csv",
                                  TAG$apply_tag))
write_csv(wide, out_csv_wide)

# LaTeX (Firth only, as primary reportable estimator)
firth_only <- results %>%
  filter(estimator == "Firth") %>%
  mutate(
    or_ci = sprintf("%.2f (%.2f, %.2f)", OR, ci_lo, ci_hi),
    p_str = sprintf("%.3f", p_value)
  ) %>%
  select(tier, term, or_ci, p_str)

build_latex_table <- function(df) {
  hdr <- c("Tier", "Term", "OR (95\\% CI)", "$p$")
  body <- apply(df, 1, function(r) paste(r, collapse = " & "))
  c("\\begin{tabular}{llcr}",
    "\\hline",
    paste(paste(hdr, collapse = " & "), "\\\\"),
    "\\hline",
    paste(body, "\\\\"),
    "\\hline",
    "\\end{tabular}",
    "\\\\\\footnotesize\\textit{Firth's penalized maximum likelihood;",
    "OR = odds of detectable within-host diversity under the calibrated rule;",
    "exploratory descriptive analysis at small N.}"
  )
}
writeLines(paste(build_latex_table(firth_only), collapse = "\n"),
           file.path(PATHS$tables,
                     sprintf("table5_logistic_%s.tex", TAG$apply_tag)))

# ---- Console summary -------------------------------------------------------
cat(strrep("=", 78), "\n", sep = "")
cat("LOGISTIC REGRESSION RESULTS (Firth)\n")
cat(strrep("=", 78), "\n", sep = "")
print(as.data.frame(firth_only), row.names = FALSE)

cat("\nMLE comparison (transparency check; large gaps from Firth ",
    "indicate sparse-data bias in MLE):\n", sep = "")
print(as.data.frame(
  results %>% filter(estimator == "MLE") %>%
    mutate(or_ci = sprintf("%.2f (%.2f, %.2f)", OR, ci_lo, ci_hi)) %>%
    select(tier, term, or_ci, p_value)
), row.names = FALSE)

cat(sprintf("\nOutputs:\n  %s\n  %s\n  %s\n",
            out_csv, out_csv_wide,
            file.path(PATHS$tables,
                      sprintf("table5_logistic_%s.tex", TAG$apply_tag))))

cat("\nReporting frame: descriptive / exploratory association under the\n")
cat("calibrated detection rule. Outcome is DETECTABLE diversity, not true\n")
cat("diversity. Wide CIs reflect cohort size; results are hypothesis-\n")
cat("generating only.\n")
