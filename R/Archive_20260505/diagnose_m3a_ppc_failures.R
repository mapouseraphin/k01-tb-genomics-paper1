# =============================================================================
# diagnose_m3a_ppc_failures.R
#
# Investigates which PPC cells M3a fails. Output: failure pattern by
# pair_type, MAF bin, and observation count, to determine whether failures
# are concentrated in expected places (low-n cells in [0.05, 0.10)) or
# scattered randomly (model-fit concern).
#
# Inputs:
#   <LCA_PATH>/fits/fit_strata_ppc.csv   (already produced by fit_lca_models.R)
#   <LCA_PATH>/lca_dataset.rds           (for n per cell)
#
# Outputs (console only):
#   D1. PPC pass rate by pair_type
#   D2. PPC pass rate by MAF bin
#   D3. PPC pass rate by cell type (00, 01, 10, 11)
#   D4. Failed cells with observation counts
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

LCA_PATH  <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
FITS_PATH <- file.path(LCA_PATH, "fits")

print_via_message <- function(df) {
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}

ppc <- read_csv(file.path(FITS_PATH, "fit_strata_ppc.csv"),
                show_col_types = FALSE)
lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))
lca_fit <- lca %>% filter(maf_bin != paste0(">", 0.45))

# ---- Compute per-cell n ----------------------------------------------------
cell_n <- lca_fit %>%
  count(pair_type, maf_bin, W1, W2, name = "n_cell") %>%
  mutate(cell = paste0(W1, W2)) %>%
  select(pair_type, maf_bin = maf_bin, cell, n_cell) %>%
  mutate(maf_bin = as.character(maf_bin))

ppc_n <- ppc %>%
  left_join(cell_n, by = c("pair_type", "maf_bin", "cell"))

message("\n=== M3a PPC failure investigation ===")
message(sprintf("Overall: %d / %d cells within 95%% CrI (%.0f%%)",
                sum(ppc_n$within_CrI, na.rm = TRUE), nrow(ppc_n),
                100 * mean(ppc_n$within_CrI, na.rm = TRUE)))

# ---- D1: by pair_type ------------------------------------------------------
D1 <- ppc_n %>% group_by(pair_type) %>%
  summarise(n_cells = n(),
            n_pass  = sum(within_CrI),
            pct     = round(100 * n_pass / n_cells), .groups = "drop")
message("\n--- D1: by pair_type ---")
print_via_message(D1)

# ---- D2: by MAF bin --------------------------------------------------------
D2 <- ppc_n %>% group_by(maf_bin) %>%
  summarise(n_cells = n(),
            n_pass  = sum(within_CrI),
            pct     = round(100 * n_pass / n_cells), .groups = "drop")
message("\n--- D2: by MAF bin ---")
print_via_message(D2)

# ---- D3: by cell type ------------------------------------------------------
D3 <- ppc_n %>% group_by(cell) %>%
  summarise(n_cells = n(),
            n_pass  = sum(within_CrI),
            pct     = round(100 * n_pass / n_cells), .groups = "drop")
message("\n--- D3: by cell type (00 = both reject; 11 = both pass) ---")
print_via_message(D3)

# ---- D4: failed cells with n -----------------------------------------------
D4 <- ppc_n %>% filter(!within_CrI) %>%
  arrange(pair_type, maf_bin, cell)
message("\n--- D4: Failed cells, with observed n in each ---")
print_via_message(D4 %>%
  select(pair_type, maf_bin, cell, n_cell, obs_prop, pp_median, pp_q025, pp_q975))

# ---- Concentration check ---------------------------------------------------
message("\n--- Concentration check ---")
fails_by_size <- ppc_n %>%
  mutate(size_bucket = cut(n_cell, c(0, 5, 10, 20, 50, Inf),
                          labels = c("<=5", "6-10", "11-20", "21-50", ">50"),
                          right = TRUE)) %>%
  group_by(size_bucket) %>%
  summarise(n_cells = n(),
            n_pass  = sum(within_CrI),
            pct     = round(100 * n_pass / n_cells), .groups = "drop")
print_via_message(fails_by_size)
message("If failures concentrate in small-n buckets, the 63%% PPC reflects")
message("expected statistical variability in low-count cells, not model misfit.")
