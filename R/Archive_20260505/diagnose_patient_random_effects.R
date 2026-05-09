# =============================================================================
# diagnose_patient_random_effects.R
#
# Sanity check on M2 primary's patient random effect u_p posterior.
#
# Concern: M2 primary fits sigma_u ~ 1.4, meaning patient-level logit detection
# spans roughly +/- 3 (95% interval). A reasonable normal distribution; but
# is it actually shaped like a normal, or are a few patients at extreme values
# encoding patient-level T-clustering rather than per-call latent T?
#
# Tests:
#   1. Per-patient u_p posterior median + 95% CrI (sorted)
#   2. Distribution shape: mean, SD, range, skewness, count beyond +/- 2*sigma_u
#   3. Patient-level W=1 rate vs u_p median (should be monotone if u captures
#      "patient detectability" rather than something stranger)
#
# Outputs:
#   <LCA_PATH>/diagnostics/m2_patient_random_effects.csv
#   <LCA_PATH>/diagnostics/m2_patient_random_effects.md
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr); library(readr); library(tibble); library(tidyr)
})

LCA_PATH  <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
FITS_PATH <- file.path(LCA_PATH, "fits")
DIAG_PATH <- file.path(LCA_PATH, "diagnostics")
dir.create(DIAG_PATH, recursive = TRUE, showWarnings = FALSE)

print_via_message <- function(df) {
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}

# ---- Load ------------------------------------------------------------------
fit_path <- file.path(FITS_PATH, "fit_hierarchical_primary.rds")
stopifnot(file.exists(fit_path))
fit <- readRDS(fit_path)

lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))
lca_fit <- lca %>% filter(maf_bin != paste0(">", 0.45))

# Reconstruct patient_int -> patientId mapping (must match fit script's
# as.integer(factor(patientId)))
patient_map <- lca_fit %>%
  mutate(patient_int = as.integer(factor(patientId))) %>%
  distinct(patient_int, patientId) %>%
  arrange(patient_int)
P <- nrow(patient_map)
message(sprintf("Loaded M2 primary fit; %d patients", P))

# ---- Extract u_p posteriors -------------------------------------------------
u_draws <- fit$draws("u", format = "draws_df")
u_cols  <- grep("^u\\[", colnames(u_draws), value = TRUE)
stopifnot(length(u_cols) == P)

u_summ <- tibble(
  patient_int = seq_len(P),
  patientId   = patient_map$patientId,
  u_median    = vapply(u_cols, function(c) median(u_draws[[c]]), numeric(1)),
  u_q025      = vapply(u_cols, function(c) quantile(u_draws[[c]], 0.025), numeric(1)),
  u_q975      = vapply(u_cols, function(c) quantile(u_draws[[c]], 0.975), numeric(1))
)

# ---- Patient-level marginal W=1 rate (descriptive) -------------------------
patient_W <- bind_rows(
  lca_fit %>% transmute(patientId, slot = slot_a, W = W1),
  lca_fit %>% transmute(patientId, slot = slot_b, W = W2)
) %>%
  group_by(patientId) %>%
  summarise(n_obs = n(),
            W_rate = round(mean(W), 3),
            .groups = "drop")

u_summ <- u_summ %>% left_join(patient_W, by = "patientId") %>%
  mutate(across(c(u_median, u_q025, u_q975), ~round(.x, 2)))

# ---- Distribution stats ----------------------------------------------------
sigma_u_est <- median(fit$draws("sigma_u", format = "draws_matrix"))

shape_stats <- tibble(
  n_patients   = P,
  sigma_u_est  = round(sigma_u_est, 3),
  u_mean       = round(mean(u_summ$u_median), 3),
  u_sd         = round(sd(u_summ$u_median), 3),
  u_min        = round(min(u_summ$u_median), 3),
  u_max        = round(max(u_summ$u_median), 3),
  range_in_SD  = round((max(u_summ$u_median) - min(u_summ$u_median)) / sigma_u_est, 2),
  skewness     = round(mean((u_summ$u_median - mean(u_summ$u_median))^3) /
                         sd(u_summ$u_median)^3, 3),
  n_outside_2SD = sum(abs(u_summ$u_median) > 2 * sigma_u_est),
  pct_outside_2SD = round(100 * sum(abs(u_summ$u_median) > 2 * sigma_u_est) / P, 1)
)

# ---- Console summary -------------------------------------------------------
message("\n--- Per-patient u_p posterior (sorted by median) ---")
print_via_message(u_summ %>% arrange(u_median))

message("\n--- Distribution shape stats ---")
print_via_message(shape_stats)

# Interpretation
message("\n--- Interpretation ---")
range_sd_val <- shape_stats$range_in_SD
skew_val     <- shape_stats$skewness
pct_out      <- shape_stats$pct_outside_2SD

if (range_sd_val > 6) {
  message(sprintf("- Range in SD = %.2f: WIDE (>6 SD). May indicate a few outlier patients.",
                  range_sd_val))
} else {
  message(sprintf("- Range in SD = %.2f: reasonable for %d patients under Normal(0, sigma_u).",
                  range_sd_val, P))
}
if (abs(skew_val) > 1) {
  message(sprintf("- Skewness = %.2f: SUBSTANTIAL (|skew|>1). Distribution is asymmetric.",
                  skew_val))
} else {
  message(sprintf("- Skewness = %.2f: roughly symmetric.", skew_val))
}
if (pct_out > 10) {
  message(sprintf("- %.1f%% of patients beyond +/- 2 SD: HIGH (expected ~5%% under Normal).",
                  pct_out))
} else {
  message(sprintf("- %.1f%% of patients beyond +/- 2 SD: consistent with Normal (expected ~5%%).",
                  pct_out))
}

# Correlation of u_median with marginal W rate
cor_uW <- cor(u_summ$u_median, u_summ$W_rate)
message(sprintf("- Correlation(u_median, marginal W_rate): %.2f", cor_uW))
message("  Expectation: positive correlation. u_p should track patient-level")
message("  detection rate. If near zero, u is not capturing what it should.")

# ---- Save -------------------------------------------------------------------
write_csv(u_summ,       file.path(DIAG_PATH, "m2_patient_random_effects.csv"))
write_csv(shape_stats,  file.path(DIAG_PATH, "m2_patient_random_effects_shape.csv"))

md <- c(
  sprintf("# M2 patient random effect diagnostic (%s)",
          format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  sprintf("M2 primary fit: sigma_u median = %.3f", sigma_u_est),
  "",
  "## Distribution shape",
  "```",
  capture.output(print(as.data.frame(shape_stats), row.names = FALSE)),
  "```",
  "",
  sprintf("Correlation(u_median, marginal patient W_rate) = %.2f", cor_uW),
  "",
  "## Per-patient u_p (sorted by posterior median)",
  "",
  "```",
  capture.output(print(as.data.frame(u_summ %>% arrange(u_median)),
                       row.names = FALSE)),
  "```"
)
writeLines(md, file.path(DIAG_PATH, "m2_patient_random_effects.md"))

message(sprintf("\n[diagnose_patient_random_effects] Done. %s/", DIAG_PATH))
