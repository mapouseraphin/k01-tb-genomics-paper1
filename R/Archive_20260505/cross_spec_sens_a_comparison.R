# =============================================================================
# cross_spec_sens_a_comparison.R   (NEW 2026-05-03)
#
# Replaces the retired M3a-vs-M3b comparison block in fit_lca_models.R.
#
# Compares the M3a (slot x MAF) sens/spec posteriors from the canonical
# specification (PE/PPE retained) against the Sens A specification
# (PE/PPE excluded). Produces:
#
#   PATHS$cross_spec/supp_table_S3_canonical_vs_sensA.csv
#   PATHS$cross_spec/supp_fig_S3_canonical_vs_sensA.pdf
#
# These are the cross-spec sensitivity outputs replacing the old Supp Fig S3
# / Supp Table S3 (which were M3a-vs-M3b under the retired primary-selection
# logic).
#
# Inputs:
#   <canonical>: build_paths_for_spec(snv_only=TRUE,  drop_ppe_cal=FALSE,
#                                     drop_ppe_apply=FALSE,
#                                     selection_method=SELECTION_METHOD)
#                $lca_fits/fit_strata_summary.csv
#   <sens_A>   : build_paths_for_spec(snv_only=TRUE,  drop_ppe_cal=TRUE,
#                                     drop_ppe_apply=TRUE,
#                                     selection_method=SELECTION_METHOD)
#                $lca_fits/fit_strata_summary.csv
#
# Behaviour:
#   - If either fit summary is missing, exit cleanly with a message.
#     (This script can be run after only the canonical run completes; in
#     that case it does nothing -- it requires both runs.)
#
# Outputs:
#   PATHS$cross_spec/supp_table_S3_canonical_vs_sensA.csv
#   PATHS$cross_spec/supp_fig_S3_canonical_vs_sensA.pdf
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(ggplot2)
})

dir.create(PATHS$cross_spec, recursive = TRUE, showWarnings = FALSE)

# ---- Resolve per-spec paths ------------------------------------------------
canon_paths <- build_paths_for_spec(snv_only         = TRUE,
                                    drop_ppe_cal     = FALSE,
                                    drop_ppe_apply   = FALSE,
                                    selection_method = SELECTION_METHOD)
sensA_paths <- build_paths_for_spec(snv_only         = TRUE,
                                    drop_ppe_cal     = TRUE,
                                    drop_ppe_apply   = TRUE,
                                    selection_method = SELECTION_METHOD)

canon_summ_path <- file.path(canon_paths$lca_fits, "fit_strata_summary.csv")
sensA_summ_path <- file.path(sensA_paths$lca_fits, "fit_strata_summary.csv")

# ---- Existence check -------------------------------------------------------
if (!file.exists(canon_summ_path) || !file.exists(sensA_summ_path)) {
  message("[cross_spec_sens_a_comparison] One or both M3a fits not found; ",
          "skipping cross-spec comparison.")
  message("  canonical (expected): ", canon_summ_path,
          "  -> ", if (file.exists(canon_summ_path)) "present" else "MISSING")
  message("  Sens A    (expected): ", sensA_summ_path,
          "  -> ", if (file.exists(sensA_summ_path)) "present" else "MISSING")
  message("  This script requires BOTH the canonical (PE/PPE retained) and ",
          "the Sens A (PE/PPE excluded) runs to have completed under ",
          sprintf("SELECTION_METHOD = %s.", SELECTION_METHOD))
  message("  No outputs written.")
  invisible(return(NULL))
}

message(sprintf("[cross_spec_sens_a_comparison] Loading M3a posteriors:"))
message(sprintf("  canonical : %s", canon_summ_path))
message(sprintf("  Sens A    : %s", sensA_summ_path))

canon_summ <- read_csv(canon_summ_path, show_col_types = FALSE)
sensA_summ <- read_csv(sensA_summ_path, show_col_types = FALSE)

# ---- Build comparison table -----------------------------------------------
mk_block <- function(df, spec_label) {
  df %>%
    filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, maf_bin,
              spec        = spec_label,
              median      = round(median, 3),
              q2.5        = round(q2.5, 3),
              q97.5       = round(q97.5, 3),
              CrI_width   = round(q97.5 - q2.5, 3),
              rhat        = round(rhat, 4),
              ess_bulk    = round(ess_bulk, 0))
}

cross_tbl <- bind_rows(
  mk_block(canon_summ, "canonical (PE/PPE retained)"),
  mk_block(sensA_summ, "Sens A (PE/PPE excluded)")
) %>%
  arrange(param, slot, maf_bin, spec)

write_csv(cross_tbl,
          file.path(PATHS$cross_spec, "supp_table_S3_canonical_vs_sensA.csv"))

# ---- Compute cell-level shifts (Sens A median - canonical median) ---------
shifts <- cross_tbl %>%
  select(param, slot, maf_bin, spec, median) %>%
  pivot_wider(names_from = spec, values_from = median) %>%
  mutate(delta = `Sens A (PE/PPE excluded)` - `canonical (PE/PPE retained)`,
         abs_delta = abs(delta))

write_csv(shifts,
          file.path(PATHS$cross_spec, "supp_table_S3b_canonical_vs_sensA_deltas.csv"))

# ---- Comparison figure -----------------------------------------------------
SLOT_LEVELS <- c("EM", "S0.1", "S0.2")
SLOT_COLORS <- c("EM" = "#0072B2", "S0.1" = "#D55E00", "S0.2" = "#009E73")

# Order MAF bins as they appear in the data, with [0.02,0.05) absent (dropped
# from M3a fitting) per locked decision.
maf_levels_observed <- cross_tbl %>%
  pull(maf_bin) %>%
  unique() %>%
  na.omit() %>%
  sort()

cross_plot <- cross_tbl %>%
  filter(!is.na(slot), !is.na(maf_bin)) %>%
  mutate(slot    = factor(slot,    levels = SLOT_LEVELS),
         maf_bin = factor(maf_bin, levels = maf_levels_observed),
         param_label = factor(ifelse(param == "sens",
                                     "Sensitivity", "Specificity"),
                              levels = c("Sensitivity", "Specificity")))

supp_fig_S3 <- ggplot(cross_plot,
                     aes(x = maf_bin, y = median, color = slot, shape = spec)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5),
                position = position_dodge(width = 0.7), width = 0.2,
                linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.7), size = 2.4,
             fill = "white", stroke = 0.9) +
  facet_grid(param_label ~ slot, scales = "free_x") +
  scale_color_manual(values = SLOT_COLORS, guide = "none") +
  scale_shape_manual(
    values = c("canonical (PE/PPE retained)" = 21,
               "Sens A (PE/PPE excluded)"    = 24),
    name   = "Specification") +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(x = "MAF bin", y = "Posterior median (95% CrI)",
       caption = paste(
         "Cross-spec M3a comparison: canonical (circles) vs Sens A (triangles).",
         "Both are M3a (slot x MAF detection + patient random effect) under",
         sprintf("SELECTION_METHOD = %s.", SELECTION_METHOD),
         "Posterior shifts indicate the contribution of PE/PPE-region calls to",
         "the LCA's measurement-property estimates.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top",
        strip.background = element_rect(fill = "grey95", color = NA),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave(file.path(PATHS$cross_spec, "supp_fig_S3_canonical_vs_sensA.pdf"),
       supp_fig_S3, width = 7, height = 5)

# ---- Console summary -------------------------------------------------------
message("\n[cross_spec_sens_a_comparison] Outputs written to ",
        PATHS$cross_spec, "/")
message("  supp_table_S3_canonical_vs_sensA.csv")
message("  supp_table_S3b_canonical_vs_sensA_deltas.csv")
message("  supp_fig_S3_canonical_vs_sensA.pdf")

# Largest shifts (top-5 by absolute delta)
top_shifts <- shifts %>%
  arrange(desc(abs_delta)) %>%
  head(5)
message("\nTop-5 cell-level shifts (Sens A - canonical):")
message(paste(capture.output(print(as.data.frame(top_shifts), row.names = FALSE)),
              collapse = "\n"))
