# =============================================================================
# diagnose_w_asymmetry.R
#
# Characterizes the (W1, W2) cell asymmetry observed in the M1 smoke test
# (pp_01 = 0.250 by model, obs_prop = 0.524 in data; pp_10 = 0.250, obs = 0.250).
#
# Hypotheses to discriminate:
#   - H1 (Option B feasible): asymmetry is a labeling artifact -- sampleA is
#     deterministically the same slot in each pair_type, and that slot has
#     systematically lower detection. Random shuffling W1 <-> W2 would restore
#     exchangeability.
#   - H2 (Option C/D required): asymmetry reflects real differences between
#     sample slots (E0/EM vs S0.x). Slot-specific detection parameters needed.
#   - H3 (other): asymmetry varies across pair_types in ways neither H1 nor
#     H2 explains.
#
# Outputs (all to stderr via message()):
#   D1. Slot composition of sampleA and sampleB by pair_type
#   D2. Marginal W=1 rate in sampleA vs sampleB (overall and by pair_type)
#   D3. Joint (W1, W2) cell counts by pair_type, with asymmetry ratio
#   D4. Same as D3 stratified by maf_bin (drives whether asymmetry depends on
#       MAF -- relevant for M3)
#   D5. Conclusion table mapping observed pattern to recommended option (B/C/D)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

LCA_PATH <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))

print_via_message <- function(df) {
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}

# ---- D1: slot composition by pair_type --------------------------------------
message("\n--- D1: Slot composition of sampleA / sampleB by pair_type ---")
D1 <- lca %>%
  distinct(pair_id, pair_type, slot_a, slot_b) %>%
  count(pair_type, slot_a, slot_b, name = "n_pairs") %>%
  arrange(pair_type)
print_via_message(D1)
message("Interpretation: if (slot_a, slot_b) is always the same within each")
message("pair_type, sampleA assignment is deterministic by slot order.")

# ---- D2: marginal W=1 rate in sampleA vs sampleB ----------------------------
message("\n--- D2: Marginal W=1 rate, sampleA vs sampleB ---")
D2_overall <- tibble(
  scope = "overall",
  pair_type = "all",
  n_calls   = nrow(lca),
  W1_eq1    = sum(lca$W1),
  W2_eq1    = sum(lca$W2),
  W1_rate   = round(mean(lca$W1), 3),
  W2_rate   = round(mean(lca$W2), 3),
  ratio_W2_to_W1 = round(mean(lca$W2) / mean(lca$W1), 2)
)

D2_by_type <- lca %>%
  group_by(pair_type) %>%
  summarise(
    scope     = "by_pair_type",
    n_calls   = n(),
    W1_eq1    = sum(W1),
    W2_eq1    = sum(W2),
    W1_rate   = round(mean(W1), 3),
    W2_rate   = round(mean(W2), 3),
    ratio_W2_to_W1 = round(mean(W2) / pmax(mean(W1), 1e-9), 2),
    .groups = "drop"
  ) %>%
  select(scope, pair_type, n_calls, W1_eq1, W2_eq1, W1_rate, W2_rate, ratio_W2_to_W1)

D2 <- bind_rows(D2_overall, D2_by_type)
print_via_message(D2)
message("Interpretation: ratio_W2_to_W1 > 1 means sampleB is MORE likely to be")
message("rule-positive than sampleA. If consistent across pair_types -> sample")
message("position is the driver. If varies by pair_type -> slot is the driver.")

# ---- D3: joint cells by pair_type with asymmetry ratio ----------------------
message("\n--- D3: Joint (W1, W2) cells by pair_type, with asymmetry ratio ---")
D3 <- lca %>%
  count(pair_type, W1, W2, name = "n") %>%
  pivot_wider(names_from = c(W1, W2), values_from = n,
              names_prefix = "cell_", values_fill = 0L) %>%
  rename(cell_00 = cell_0_0, cell_01 = cell_0_1,
         cell_10 = cell_1_0, cell_11 = cell_1_1) %>%
  mutate(
    n_total  = cell_00 + cell_01 + cell_10 + cell_11,
    asym_ratio = round(cell_01 / pmax(cell_10, 1e-9), 2),  # >1 => sampleB-favored
    asym_diff  = cell_01 - cell_10
  ) %>%
  select(pair_type, n_total, cell_00, cell_01, cell_10, cell_11,
         asym_ratio, asym_diff) %>%
  arrange(pair_type)
print_via_message(D3)
message("asym_ratio = (W1=0,W2=1) / (W1=1,W2=0). Symmetric LCA assumes ratio=1.")
message("Ratio > 1 means sampleB is more often the W=1 of the pair.")

# ---- D4: same stratified by MAF bin -----------------------------------------
message("\n--- D4: Asymmetry by pair_type x MAF bin ---")
D4 <- lca %>%
  filter(!is.na(maf_bin)) %>%
  count(pair_type, maf_bin, W1, W2, name = "n") %>%
  pivot_wider(names_from = c(W1, W2), values_from = n,
              names_prefix = "cell_", values_fill = 0L) %>%
  rename_with(~ gsub("cell_", "c", .x), starts_with("cell_")) %>%
  rename(c00 = c0_0, c01 = c0_1, c10 = c1_0, c11 = c1_1) %>%
  mutate(asym_ratio = round(c01 / pmax(c10, 1e-9), 2)) %>%
  select(pair_type, maf_bin, c00, c01, c10, c11, asym_ratio) %>%
  arrange(pair_type, maf_bin)
print_via_message(D4)
message("Interpretation: if asym_ratio is stable across MAF bins within a")
message("pair_type, the asymmetry is sample-position-driven and MAF-invariant.")
message("If asym_ratio varies with MAF, slot-specific MAF-dependent detection")
message("is needed (more complex M3 sensitivity).")

# ---- D5: decision rule mapping ---------------------------------------------
message("\n--- D5: Decision rule ---")
message("")
message("Examine D1, D2, D3, D4 together:")
message("")
message("OPTION B (random W1/W2 shuffle in build):")
message("  - All pair_types have similar asym_ratio (within 30%% of each other)")
message("  - asym_ratio same direction across pair_types (all >1 or all <1)")
message("  - D1 shows sampleA is deterministically a single slot per pair_type")
message("  - => asymmetry is a labeling artifact; shuffling restores exchangeability")
message("")
message("OPTION C/D (slot-specific detection):")
message("  - asym_ratio differs substantially across pair_types")
message("  - OR asym_ratio in S0.1_S0.2 is approximately 1 but EM_S0.x are not")
message("    (this implicates the EM/E0 slot specifically)")
message("  - OR D4 shows MAF-dependence in the asymmetry")
message("  - => slot is a real driver; need slot-specific or randomised model")
message("")
message("Reference: M1 v1 overall (W1=0,W2=1) prop = 0.524, (W1=1,W2=0) = 0.250")
message("           overall asym_ratio = ~2.1")

message("\n[diagnose_w_asymmetry] Done. No artifacts written; this is a console")
message("diagnostic only.")
