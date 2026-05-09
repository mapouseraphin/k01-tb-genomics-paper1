# =============================================================================
# diagnose_within_patient_icc.R
#
# Quantifies within-patient correlation of W, conditional on slot. Motivates
# (or de-motivates) the patient random effect in M2.
#
# Method:
#   For each (patient, slot, position) -- after deduplicating positions that
#   appear in multiple pairs of the same patient -- compute the empirical
#   W rate per patient per slot. Then for each slot:
#     - report mean W rate, between-patient SD of rates, and a
#       method-of-moments ICC for binary data.
#
# Method-of-moments ICC for binary outcomes:
#   V_obs        = var(p_hat_p)             observed between-patient variance
#   V_within_avg = mean(p_hat_p * (1 - p_hat_p) / n_p)   sampling component
#   V_between    = max(0, V_obs - V_within_avg)          true between-patient
#   p_bar        = mean(W)                  overall rate
#   ICC          = V_between / (V_between + p_bar * (1 - p_bar))
#
# Interpretation:
#   ICC near 0  -> within-patient W observations behave independently;
#                  conditional independence may hold; M2 sigma_u likely small.
#   ICC > 0.1   -> substantial within-patient correlation; M2 is needed and
#                  conditional independence in M1 is violated.
#
# Caveat:
#   This diagnostic measures *marginal* W correlation within patient x slot.
#   It does not directly test conditional independence W1 ⊥ W2 | T because
#   T is latent and may itself cluster by patient. Either source -- shared
#   noise (technical) or shared T (biological) -- motivates the patient
#   random effect.
#
# 2026-04-30 PATCH: NA-safe. Drops rows with NA slot before computing ICC
# (with a warning). Guards interpretation loop against NA ICC and against
# slots with n_patients < 2 (where var() is undefined).
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

LCA_PATH  <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
DIAG_PATH <- file.path(LCA_PATH, "diagnostics")
dir.create(DIAG_PATH, recursive = TRUE, showWarnings = FALSE)

lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))

# Drop sentinel rows for consistency with M1/M2 fitting set
lca <- lca %>% filter(maf_bin != paste0(">", 0.45))

print_via_message <- function(df) {
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}

# ---- Build long form: one row per (patient, slot, position) ---------------
# Each variant call has W1 (slot_a side) and W2 (slot_b side). Reshape to
# long, then deduplicate by (patient, slot, CHROM, POS, REF, ALT) so that a
# sample appearing in multiple pair_types is counted once per position.
lca_long <- bind_rows(
  lca %>% transmute(patientId, slot = slot_a,
                    CHROM, POS, REF, ALT, W = W1),
  lca %>% transmute(patientId, slot = slot_b,
                    CHROM, POS, REF, ALT, W = W2)
) %>%
  distinct(patientId, slot, CHROM, POS, REF, ALT, .keep_all = TRUE)

# ---- PATCH: drop NA slots upfront, with a diagnostic summary ----------------
n_total <- nrow(lca_long)
n_na_slot <- sum(is.na(lca_long$slot))
if (n_na_slot > 0) {
  message(sprintf("WARNING: %d / %d (%.2f%%) long-form rows have NA slot.",
                  n_na_slot, n_total, 100 * n_na_slot / n_total))
  affected <- lca_long %>%
    filter(is.na(slot)) %>%
    distinct(patientId) %>%
    pull(patientId)
  message(sprintf("  Affected patient IDs (%d): %s", length(affected),
                  paste(head(affected, 10), collapse = ", "),
                  if (length(affected) > 10) " ..." else ""))
  message("  Dropping these rows from ICC computation; investigate parse_slot() ",
          "upstream if this is unexpected.")
  lca_long <- lca_long %>% filter(!is.na(slot))
}
n_na_pid <- sum(is.na(lca_long$patientId))
if (n_na_pid > 0) {
  message(sprintf("WARNING: %d rows with NA patientId; dropping.", n_na_pid))
  lca_long <- lca_long %>% filter(!is.na(patientId))
}

message(sprintf("Long-form (patient x slot x position) observations: %s",
                format(nrow(lca_long), big.mark = ",")))
message(sprintf("Distinct patients: %d  |  distinct slots: %s",
                n_distinct(lca_long$patientId),
                paste(sort(unique(lca_long$slot)), collapse = ", ")))

# ---- Per-patient W rate by slot --------------------------------------------
patient_rates <- lca_long %>%
  group_by(slot, patientId) %>%
  summarise(n_obs   = n(),
            n_W1    = sum(W),
            W_rate  = mean(W),
            .groups = "drop")

message("\n--- Per-patient W rate distribution by slot ---")
patient_rate_summary <- patient_rates %>%
  group_by(slot) %>%
  summarise(
    n_patients     = n(),
    n_obs_total    = sum(n_obs),
    obs_per_pt_med = round(median(n_obs), 1),
    obs_per_pt_iqr = paste0(quantile(n_obs, 0.25), "-", quantile(n_obs, 0.75)),
    p_bar          = round(mean(W_rate), 3),
    sd_rates       = round(sd(W_rate), 3),
    .groups = "drop"
  )
print_via_message(patient_rate_summary)

# ---- Method-of-moments ICC by slot -----------------------------------------
icc_by_slot <- patient_rates %>%
  group_by(slot) %>%
  summarise(
    n_patients     = n(),
    p_bar          = mean(W_rate),
    V_obs          = if (n() >= 2) var(W_rate) else NA_real_,
    V_within_avg   = mean(W_rate * (1 - W_rate) / pmax(n_obs, 1)),
    V_between      = pmax(0, V_obs - V_within_avg),
    ICC            = V_between / (V_between + p_bar * (1 - p_bar)),
    .groups = "drop"
  ) %>%
  mutate(across(c(p_bar, V_obs, V_within_avg, V_between, ICC), ~round(.x, 4)))

message("\n--- Method-of-moments ICC of W by slot ---")
print_via_message(icc_by_slot)

# ---- Interpretation (NA-safe) ----------------------------------------------
message("\n--- Interpretation ---")
for (i in seq_len(nrow(icc_by_slot))) {
  s       <- icc_by_slot$slot[i]
  icc_val <- icc_by_slot$ICC[i]
  n_pt    <- icc_by_slot$n_patients[i]

  if (is.na(s)) {
    next  # already filtered, but defensive
  }
  if (n_pt < 2) {
    message(sprintf("  %s: ICC = NA -- only %d patient(s); ICC undefined.",
                    s, n_pt))
    next
  }
  if (is.na(icc_val)) {
    message(sprintf("  %s: ICC = NA -- variance components undefined ",
                    "(likely degenerate W_rate distribution).", s))
    next
  }

  msg <- if (icc_val < 0.05) {
    sprintf("  %s: ICC = %.3f -- LOW; within-patient observations behave nearly independently.",
            s, icc_val)
  } else if (icc_val < 0.15) {
    sprintf("  %s: ICC = %.3f -- MODERATE; some within-patient clustering.",
            s, icc_val)
  } else {
    sprintf("  %s: ICC = %.3f -- HIGH; substantial within-patient clustering. M2 needed.",
            s, icc_val)
  }
  message(msg)
}
message("")
message("Overall guidance:")
message("  - Any slot with ICC > 0.05 motivates including a patient random")
message("    effect in M2.")
message("  - All ICCs near 0 -> M1's conditional independence assumption may")
message("    be defensible; the M1 PPC failure may have other sources.")

# ---- Save -------------------------------------------------------------------
write_csv(patient_rate_summary, file.path(DIAG_PATH, "icc_patient_rate_summary.csv"))
write_csv(icc_by_slot,         file.path(DIAG_PATH, "icc_by_slot.csv"))
saveRDS(patient_rates,         file.path(DIAG_PATH, "icc_patient_rates.rds"))

message(sprintf("\n[diagnose_within_patient_icc] Done. Tables written to %s/",
                DIAG_PATH))
