# =============================================================================
# cross_spec_sens_a_comparison.R
#
# Builds Supplemental Table S3 / Figure S3: side-by-side comparison of M3a
# per-cell (sens, spec) posteriors fitted under the canonical specification
# (PE/PPE retained) vs the Sens A specification (PE/PPE excluded).
#
# Replaces the prior M3a-vs-M3b comparison, which was retired 2026-05-03 after
# M3b's catastrophic convergence failure.
#
# Inputs (resolved via build_paths_for_spec):
#   Canonical: data_derived/05_lca_snv_only_ppe_retained_<sel>/fits/
#                fit_strata.rds, fit_strata_summary.csv
#   Sens A   : data_derived/05_lca_snv_only_ppe_excluded_<sel>/fits/
#                fit_strata.rds, fit_strata_summary.csv
#
# Where <sel> is the active SELECTION_METHOD (lex or pareto).
#
# Outputs (PATHS$cross_spec):
#   supp_table_S3_canonical_vs_sensA.csv
#   supp_fig_S3_canonical_vs_sensA.pdf
#
# Run-prereq: this script requires both canonical and Sens A pipeline runs
# to have completed. If either fit file is missing, the script exits cleanly
# with a message; it does not fail the orchestrator.
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

dir.create(PATHS$cross_spec, recursive = TRUE, showWarnings = FALSE)

paths_canonical <- build_paths_for_spec(
  snv_only         = TRUE,
  drop_ppe_cal     = FALSE,    # canonical: PE/PPE retained
  drop_ppe_apply   = FALSE,
  selection_method = SELECTION_METHOD
)
paths_sens_a    <- build_paths_for_spec(
  snv_only         = TRUE,
  drop_ppe_cal     = TRUE,     # Sens A: PE/PPE excluded
  drop_ppe_apply   = TRUE,
  selection_method = SELECTION_METHOD
)

summ_canonical <- file.path(paths_canonical$lca_fits, "fit_strata_summary.csv")
summ_sens_a    <- file.path(paths_sens_a$lca_fits,    "fit_strata_summary.csv")

if (!file.exists(summ_canonical) || !file.exists(summ_sens_a)) {
  message("[cross_spec_sens_a_comparison] One or both fits not found:")
  message(sprintf("  canonical: %s (exists=%s)", summ_canonical, file.exists(summ_canonical)))
  message(sprintf("  sens_a   : %s (exists=%s)", summ_sens_a,    file.exists(summ_sens_a)))
  message("Run both canonical and Sens A pipelines before invoking this script.")
  return(invisible(NULL))
}

read_summ <- function(path, label) {
  read_csv(path, show_col_types = FALSE) %>%
    filter(param %in% c("sens", "spec")) %>%
    transmute(
      param,
      slot,
      maf_bin,
      spec_label = label,
      median     = round(median, 3),
      q2.5       = round(q2.5,  3),
      q97.5      = round(q97.5, 3),
      CrI_width  = round(q97.5 - q2.5, 3),
      rhat,
      ess_bulk
    )
}

s_canon <- read_summ(summ_canonical, "Canonical (PE/PPE retained)")
s_sensA <- read_summ(summ_sens_a,    "Sens A (PE/PPE excluded)")

cmp <- bind_rows(s_canon, s_sensA) %>%
  arrange(param, slot, maf_bin, spec_label)

write_csv(cmp, file.path(PATHS$cross_spec,
                         "supp_table_S3_canonical_vs_sensA.csv"))
message(sprintf("[cross_spec] Wrote: %s/supp_table_S3_canonical_vs_sensA.csv",
                PATHS$cross_spec))

# Per-cell delta summary (Sens A median - canonical median)
deltas <- s_canon %>%
  select(param, slot, maf_bin, median_canonical = median) %>%
  inner_join(s_sensA %>% select(param, slot, maf_bin, median_sensA = median),
             by = c("param", "slot", "maf_bin")) %>%
  mutate(delta = round(median_sensA - median_canonical, 3)) %>%
  arrange(param, slot, maf_bin)

write_csv(deltas, file.path(PATHS$cross_spec,
                            "supp_table_S3_deltas.csv"))
message(sprintf("[cross_spec] Wrote: %s/supp_table_S3_deltas.csv",
                PATHS$cross_spec))

# ---- Figure S3: forest plot, canonical vs Sens A by (param, slot, MAF) -----
plot_df <- bind_rows(s_canon, s_sensA) %>%
  mutate(
    param      = factor(param, levels = c("sens", "spec"),
                        labels = c("Sensitivity", "Specificity")),
    slot       = factor(slot,  levels = c("EM", "S0.1", "S0.2")),
    cell_label = sprintf("%s | %s", slot, maf_bin),
    spec_label = factor(spec_label,
                        levels = c("Canonical (PE/PPE retained)",
                                   "Sens A (PE/PPE excluded)"))
  )

p <- ggplot(plot_df,
            aes(x = median, xmin = q2.5, xmax = q97.5,
                y = cell_label, colour = spec_label,
                shape = spec_label)) +
  geom_pointrange(position = position_dodge(width = 0.5),
                  size = 0.3) +
  facet_wrap(~ param, ncol = 2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Posterior median (95% CrI)", y = NULL,
       colour = "Specification", shape = "Specification",
       title = "M3a per-cell posteriors: canonical vs Sens A (PE/PPE excluded)") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank())

ggsave(file.path(PATHS$cross_spec, "supp_fig_S3_canonical_vs_sensA.pdf"),
       p, width = 9, height = 6)
message(sprintf("[cross_spec] Wrote: %s/supp_fig_S3_canonical_vs_sensA.pdf",
                PATHS$cross_spec))

message("\n[cross_spec_sens_a_comparison] Done.")
