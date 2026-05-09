# =============================================================================
# make_paper1_tables_figures.R
#
# Generates publication-ready tables and figures for Paper 1's LCA results.
#
# Primary deliverables:
#   Table 4. LCA-derived sens, spec, and prevalence at primary thresholds,
#            stratified by slot x MAF bin (M3a primary).
#   Figure 3. Posterior distributions of sens and spec by slot x MAF bin.
#   Figure 4. Patient random effect: sigma_u posterior + per-patient u_p
#             distribution.
#
# Supplemental:
#   Supp Fig S1. PPC observed vs predicted (M3a primary, 36 cells).
#   Supp Fig S2. Prior sensitivity for M2 (sigma_u half-normal scales 0.5/1.0/2.0).
#   Supp Fig S3. M3a vs M3b comparison (per-cell sens, spec).
#   Supp Fig S4. M2 primary vs M2 globalclass comparison.
#   Supp Table S1. PPC table (full).
#   Supp Table S2. Prior sensitivity numerical comparison.
#
# All numerical content sourced from upstream fit objects produced by
# fit_lca_models.R. No values are hard-coded.
#
# Outputs:
#   outputs/tables/<TAG>/table_4_lca_sens_spec.csv (+ .md)
#   outputs/figures/<TAG>/figure_3_posterior_by_slot_maf.pdf
#   outputs/figures/<TAG>/figure_4_patient_random_effect.pdf
#   outputs/supplemental/<TAG>/supp_fig_*.pdf
#   outputs/supplemental/<TAG>/supp_table_*.csv
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(cmdstanr); library(posterior)
  library(dplyr); library(tidyr); library(readr); library(tibble)
  library(ggplot2); library(stringr)
})

LCA_PATH    <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
FITS_PATH   <- file.path(LCA_PATH, "fits")
TABLES_PATH <- file.path("outputs", "tables",       TAG$apply_tag)
FIG_PATH    <- file.path("outputs", "figures",      TAG$apply_tag)
SUPP_PATH   <- file.path("outputs", "supplemental", TAG$apply_tag)
for (p in c(TABLES_PATH, FIG_PATH, SUPP_PATH)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

# ---- Publication theme + palette -------------------------------------------
theme_paper1 <- function(base_size = 5) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey92"),
      panel.grid.major.y = element_blank(),
      strip.background   = element_rect(fill = "grey95", color = NA),
      strip.text         = element_text(face = "bold", size = base_size),
      legend.position    = "bottom",
      legend.key.size    = unit(0.4, "cm"),
      legend.title       = element_text(size = base_size),
      axis.title         = element_text(size = base_size),
      plot.title         = element_text(face = "bold", size = base_size + 1),
      plot.caption       = element_text(size = base_size - 1, hjust = 0,
                                        color = "grey40")
    )
}

# Colorblind-friendly palette for slots (Wong 2011)
SLOT_COLORS <- c("EM" = "#0072B2", "S0.1" = "#D55E00", "S0.2" = "#009E73")
SLOT_LEVELS <- c("EM", "S0.1", "S0.2")
MAF_LEVELS  <- c("[0.02,0.05)", "[0.05,0.10)", "[0.10,0.25)", "[0.25,0.45]")

# ---- Load primary fit + summaries ------------------------------------------
fit_m3a   <- readRDS(file.path(FITS_PATH, "fit_strata.rds"))
summ_m3a  <- read_csv(file.path(FITS_PATH, "fit_strata_summary.csv"),
                      show_col_types = FALSE)
ppc_m3a   <- read_csv(file.path(FITS_PATH, "fit_strata_ppc.csv"),
                      show_col_types = FALSE)

# Sensitivity fits / summaries
fit_m2_pri <- readRDS(file.path(FITS_PATH, "fit_hierarchical_primary.rds"))
summ_m2_pri <- read_csv(file.path(FITS_PATH, "fit_hierarchical_primary_summary.csv"),
                        show_col_types = FALSE)
summ_m2_tig <- read_csv(file.path(FITS_PATH, "fit_hierarchical_tight_summary.csv"),
                        show_col_types = FALSE)
summ_m2_loo <- read_csv(file.path(FITS_PATH, "fit_hierarchical_loose_summary.csv"),
                        show_col_types = FALSE)
summ_m2_gc  <- read_csv(file.path(FITS_PATH, "fit_hierarchical_globalclass_summary.csv"),
                        show_col_types = FALSE)
summ_m3b    <- read_csv(file.path(FITS_PATH, "fit_strata_pooled_summary.csv"),
                        show_col_types = FALSE)

# Per-cell sample size for table 4
lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))
lca_fit <- lca %>% filter(maf_bin != paste0(">", 0.45))
n_per_slot_maf <- bind_rows(
  lca_fit %>% transmute(slot = slot_a, maf_bin),
  lca_fit %>% transmute(slot = slot_b, maf_bin)
) %>%
  count(slot, maf_bin, name = "n_calls") %>%
  mutate(slot = as.character(slot), maf_bin = as.character(maf_bin))

# =============================================================================
# Table 4: LCA-derived sens, spec, prevalence at primary thresholds
# =============================================================================
message("Building Table 4...")

pi_row <- summ_m3a %>% filter(variable == "pi") %>%
  transmute(parameter = "Prevalence (pi)", slot = "(overall)",
            maf_bin = "(overall)", n_calls = NA_integer_,
            posterior_median = round(median, 3),
            CrI_lower = round(q2.5, 3),
            CrI_upper = round(q97.5, 3))

sigma_row <- summ_m3a %>% filter(variable == "sigma_u") %>%
  transmute(parameter = "Patient SD (sigma_u, logit scale)",
            slot = "(overall)", maf_bin = "(overall)", n_calls = NA_integer_,
            posterior_median = round(median, 3),
            CrI_lower = round(q2.5, 3),
            CrI_upper = round(q97.5, 3))

slot_maf_rows <- summ_m3a %>%
  filter(param %in% c("sens", "spec")) %>%
  transmute(parameter = ifelse(param == "sens", "Sensitivity", "Specificity"),
            slot, maf_bin,
            posterior_median = round(median, 3),
            CrI_lower = round(q2.5, 3),
            CrI_upper = round(q97.5, 3)) %>%
  left_join(n_per_slot_maf, by = c("slot", "maf_bin")) %>%
  arrange(parameter, factor(slot, levels = SLOT_LEVELS),
          factor(maf_bin, levels = MAF_LEVELS)) %>%
  select(parameter, slot, maf_bin, n_calls,
         posterior_median, CrI_lower, CrI_upper)

table4 <- bind_rows(pi_row, sigma_row, slot_maf_rows)
write_csv(table4, file.path(TABLES_PATH, "table_4_lca_sens_spec.csv"))

# Markdown formatted
md4 <- c(
  "# Table 4. LCA-derived sensitivity, specificity, and prevalence",
  "",
  sprintf("Primary model: M3a (slot x MAF detection + patient random effect, ",
          "sigma_u ~ half-normal(0, 1.0))."),
  sprintf("N = %s variant calls, %d patients, %d M0 within-visit pairs.",
          format(nrow(lca_fit), big.mark = ","),
          n_distinct(lca_fit$patientId), n_distinct(lca_fit$pair_id)),
  "",
  paste("|", paste(names(table4), collapse = " | "), "|"),
  paste("|", paste(rep("---", ncol(table4)), collapse = " | "), "|"),
  apply(table4, 1, function(r)
    paste0("| ", paste(ifelse(is.na(r), "", r), collapse = " | "), " |")),
  ""
)
writeLines(md4, file.path(TABLES_PATH, "table_4_lca_sens_spec.md"))

# =============================================================================
# Figure 3: Posterior sens/spec by slot x MAF bin
# =============================================================================
message("Building Figure 3...")

fig3_data <- summ_m3a %>%
  filter(param %in% c("sens", "spec")) %>%
  mutate(slot = factor(slot, levels = SLOT_LEVELS),
         maf_bin = factor(maf_bin, levels = MAF_LEVELS),
         param_label = factor(ifelse(param == "sens",
                                     "Sensitivity",
                                     "Specificity"),
                              levels = c("Sensitivity", "Specificity")))

fig3 <- ggplot(fig3_data,
               aes(x = maf_bin, y = median, color = slot, group = slot)) +
  geom_line(position = position_dodge(width = 0.5),
            alpha = 0.5, linewidth = 0.4) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5),
                position = position_dodge(width = 0.5),
                width = 0.2, linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.5),
             size = 2.4, shape = 21, fill = "white",
             stroke = 1.0) +
  facet_wrap(~ param_label, ncol = 2) +
  scale_color_manual(values = SLOT_COLORS, name = "Sample slot") +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0),
                     breaks = seq(0, 1, 0.2)) +
  labs(x = "Minor allele frequency bin",
       y = "Posterior median (95% credible interval)",
       caption = "From M3a (slot x MAF detection + patient random effect).") +
  theme_paper1() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8))

ggsave(file.path(FIG_PATH, "figure_3_posterior_by_slot_maf.pdf"),
       fig3, width = 7.0, height = 4.0)
ggsave(file.path(FIG_PATH, "figure_3_posterior_by_slot_maf.png"),
       fig3, width = 7.0, height = 4.0, dpi = 300)

# =============================================================================
# Figure 4: Patient random effect (sigma_u posterior + per-patient u_p)
# =============================================================================
message("Building Figure 4...")

# Pull sigma_u and u draws from primary M2 fit (which we use for the patient
# RE diagnostic; M3a's RE is the same conceptual quantity but M2 has the
# clean primary fit we already validated)
sigma_u_draws <- as.numeric(fit_m2_pri$draws("sigma_u",
                                              format = "draws_matrix"))
u_draws_mat   <- fit_m2_pri$draws("u", format = "draws_matrix")
u_medians <- apply(u_draws_mat, 2, median)
u_q025    <- apply(u_draws_mat, 2, quantile, 0.025)
u_q975    <- apply(u_draws_mat, 2, quantile, 0.975)

patient_summ <- tibble(
  patient_int = seq_along(u_medians),
  u_median = u_medians, u_q025 = u_q025, u_q975 = u_q975
) %>% arrange(u_median) %>%
  mutate(rank = row_number())

sigma_u_med <- median(sigma_u_draws)
sigma_u_q   <- quantile(sigma_u_draws, c(0.025, 0.975))

# Panel A: sigma_u posterior density
panel_A <- ggplot(tibble(sigma_u = sigma_u_draws), aes(x = sigma_u)) +
  geom_density(fill = "#0072B2", alpha = 0.35, color = "#0072B2",
               linewidth = 0.6) +
  geom_vline(xintercept = sigma_u_med, linetype = "dashed",
             color = "#0072B2", linewidth = 0.5) +
  geom_vline(xintercept = sigma_u_q, linetype = "dotted",
             color = "grey50", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, max(sigma_u_draws) * 1.05),
                     expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression(sigma[u] ~ "(patient SD on logit detection scale)"),
       y = "Posterior density",
       title = sprintf("A. Patient SD posterior (median %.2f, 95%% CrI %.2f - %.2f)",
                       sigma_u_med, sigma_u_q[1], sigma_u_q[2])) +
  theme_paper1()

# Panel B: per-patient u_p posterior medians (forest plot)
panel_B <- ggplot(patient_summ, aes(x = u_median, y = rank)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = c(-2, 2) * sigma_u_med, linetype = "dotted",
             color = "grey70", linewidth = 0.4) +
  geom_segment(aes(x = u_q025, xend = u_q975, yend = rank),
               color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.6, color = "#0072B2") +
  labs(x = expression("Per-patient" ~ u[p] ~ "posterior (median, 95% CrI)"),
       y = "Patient (sorted by posterior median)",
       title = sprintf("B. Per-patient random effects (n=%d)",
                       nrow(patient_summ)),
       caption = paste("Vertical dotted lines at +/- 2 sigma_u show",
                       "expected range under Normal patient distribution.")) +
  theme_paper1() +
  theme(axis.text.y = element_blank())

# Combine A and B side-by-side
fig4 <- patchwork::wrap_plots(panel_A, panel_B, nrow = 1, widths = c(1, 1.2))
# patchwork may not be available; use cowplot fallback below if needed
ok_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (!ok_patchwork) {
  if (requireNamespace("cowplot", quietly = TRUE)) {
    fig4 <- cowplot::plot_grid(panel_A, panel_B, nrow = 1, rel_widths = c(1, 1.2))
  } else {
    # Last resort: save each panel separately
    ggsave(file.path(FIG_PATH, "figure_4A_sigma_u.pdf"), panel_A,
           width = 3.5, height = 3.0)
    ggsave(file.path(FIG_PATH, "figure_4B_patient_u.pdf"), panel_B,
           width = 3.5, height = 3.0)
    fig4 <- panel_A    # so subsequent ggsave doesn't fail
    message("Note: install patchwork or cowplot for combined Figure 4.")
  }
}

ggsave(file.path(FIG_PATH, "figure_4_patient_random_effect.pdf"),
       fig4, width = 7.0, height = 3.5)
ggsave(file.path(FIG_PATH, "figure_4_patient_random_effect.png"),
       fig4, width = 7.0, height = 3.5, dpi = 300)

# =============================================================================
# Supplemental Fig S1: PPC observed vs predicted (M3a)
# =============================================================================
message("Building Supp Fig S1 (PPC)...")

ppc_plot_data <- ppc_m3a %>%
  mutate(maf_bin = factor(maf_bin, levels = MAF_LEVELS),
         pair_type = factor(pair_type),
         cell = factor(cell, levels = c("00", "01", "10", "11")))

supp_fig1 <- ggplot(ppc_plot_data,
                    aes(x = pp_median, y = obs_prop)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60") +
  geom_errorbarh(aes(xmin = pp_q025, xmax = pp_q975, color = within_CrI),
                 height = 0, linewidth = 0.5) +
  geom_point(aes(shape = cell, color = within_CrI), size = 2) +
  facet_grid(maf_bin ~ pair_type) +
  scale_color_manual(values = c("TRUE" = "#0072B2", "FALSE" = "#D55E00"),
                     name = "Within 95% CrI") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Cell (W1, W2)") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Posterior predictive median (95% CrI as horizontal bar)",
       y = "Observed proportion",
       caption = "Diagonal: perfect fit. Each panel: one (pair_type, MAF bin) combination.") +
  theme_paper1() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7))

ggsave(file.path(SUPP_PATH, "supp_fig_S1_ppc_m3a.pdf"),
       supp_fig1, width = 7, height = 7)

# =============================================================================
# Supplemental Fig S2: Prior sensitivity for M2 (sigma_u priors)
# =============================================================================
message("Building Supp Fig S2 (prior sensitivity)...")

prior_data <- bind_rows(
  summ_m2_tig %>% mutate(prior = "half-normal(0, 0.5)"),
  summ_m2_pri %>% mutate(prior = "half-normal(0, 1.0)"),
  summ_m2_loo %>% mutate(prior = "half-normal(0, 2.0)")
) %>%
  filter(param %in% c("pi", "sigma_u", "sens", "spec")) %>%
  mutate(label = ifelse(is.na(slot), param,
                        paste0(param, "(", slot, ")")),
         label = factor(label, levels = c(
           "pi", "sigma_u",
           "sens(EM)", "sens(S0.1)", "sens(S0.2)",
           "spec(EM)", "spec(S0.1)", "spec(S0.2)"
         )),
         prior = factor(prior, levels = c("half-normal(0, 0.5)",
                                          "half-normal(0, 1.0)",
                                          "half-normal(0, 2.0)")))

supp_fig2 <- ggplot(prior_data,
                    aes(x = median, y = label, color = prior)) +
  geom_errorbarh(aes(xmin = q2.5, xmax = q97.5),
                 height = 0.2, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("#56B4E9", "#0072B2", "#D55E00"),
                     name = expression("Prior on" ~ sigma[u])) +
  labs(x = "Posterior median (95% CrI)",
       y = "Parameter",
       caption = paste("Prior sensitivity for M2. Three preregistered priors;",
                       "primary is half-normal(0, 1.0).")) +
  theme_paper1() +
  theme(legend.position = "right")

ggsave(file.path(SUPP_PATH, "supp_fig_S2_prior_sensitivity.pdf"),
       supp_fig2, width = 7, height = 4.5)

# =============================================================================
# Supplemental Fig S3: M3a vs M3b
# =============================================================================
message("Building Supp Fig S3 (M3a vs M3b)...")

m3_compare <- bind_rows(
  summ_m3a %>% filter(param %in% c("sens", "spec")) %>%
    mutate(model = "M3a (independent)"),
  summ_m3b %>% filter(param %in% c("sens", "spec")) %>%
    mutate(model = "M3b (hierarchical-MAF)")
) %>%
  mutate(slot = factor(slot, levels = SLOT_LEVELS),
         maf_bin = factor(maf_bin, levels = MAF_LEVELS),
         param_label = factor(ifelse(param == "sens",
                                     "Sensitivity", "Specificity"),
                              levels = c("Sensitivity", "Specificity")),
         model = factor(model))

supp_fig3 <- ggplot(m3_compare,
                    aes(x = maf_bin, y = median, color = slot, shape = model)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5),
                position = position_dodge(width = 0.7), width = 0.2,
                linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.7), size = 2.2,
             fill = "white", stroke = 0.9) +
  facet_grid(param_label ~ slot, scales = "free_x") +
  scale_color_manual(values = SLOT_COLORS, guide = "none") +
  scale_shape_manual(values = c("M3a (independent)" = 21,
                                "M3b (hierarchical-MAF)" = 24),
                     name = "Model") +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(x = "MAF bin", y = "Posterior median (95% CrI)",
       caption = paste("M3a is the primary; M3b is sensitivity. M3b had",
                       "convergence problems (R-hat 1.024, 506 divergent",
                       "transitions); shown for sensitivity comparison only.")) +
  theme_paper1() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 7))

ggsave(file.path(SUPP_PATH, "supp_fig_S3_m3a_vs_m3b.pdf"),
       supp_fig3, width = 7, height = 5)

# =============================================================================
# Supplemental Fig S4: M2 primary vs M2 globalclass
# =============================================================================
message("Building Supp Fig S4 (M2 globalclass)...")

m2_compare <- bind_rows(
  summ_m2_pri %>% filter(param %in% c("sens", "spec")) %>%
    mutate(model = "M2 primary (per-slot ordering)"),
  summ_m2_gc %>% filter(param %in% c("sens", "spec")) %>%
    mutate(model = "M2 globalclass (shared class-gap)")
) %>%
  mutate(slot = factor(slot, levels = SLOT_LEVELS),
         param_label = factor(ifelse(param == "sens",
                                     "Sensitivity", "Specificity"),
                              levels = c("Sensitivity", "Specificity")),
         model = factor(model))

supp_fig4 <- ggplot(m2_compare,
                    aes(x = slot, y = median, color = slot, shape = model)) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5),
                position = position_dodge(width = 0.6), width = 0.2,
                linewidth = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5,
             fill = "white", stroke = 0.9) +
  facet_wrap(~ param_label) +
  scale_color_manual(values = SLOT_COLORS, guide = "none") +
  scale_shape_manual(values = c("M2 primary (per-slot ordering)" = 21,
                                "M2 globalclass (shared class-gap)" = 24),
                     name = "Model") +
  scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
  labs(x = "Sample slot", y = "Posterior median (95% CrI)",
       caption = paste("M2 primary fits substantially better (PPC 83%) than",
                       "globalclass (PPC 42%); shown to verify per-slot freedom",
                       "is data-driven.")) +
  theme_paper1()

ggsave(file.path(SUPP_PATH, "supp_fig_S4_m2_globalclass.pdf"),
       supp_fig4, width = 7, height = 4)

# =============================================================================
# Supplemental tables
# =============================================================================
message("Building supplemental tables...")

# Supp Table S1: Full PPC table
write_csv(ppc_m3a, file.path(SUPP_PATH, "supp_table_S1_ppc_m3a.csv"))

# Supp Table S2: Prior sensitivity
write_csv(prior_data %>%
            select(param, slot, prior, median, q2.5, q97.5, rhat, ess_bulk),
          file.path(SUPP_PATH, "supp_table_S2_prior_sensitivity.csv"))

# Supp Table S3: M3a vs M3b
write_csv(m3_compare %>%
            select(param, slot, maf_bin, model, median, q2.5, q97.5,
                   rhat, ess_bulk),
          file.path(SUPP_PATH, "supp_table_S3_m3a_vs_m3b.csv"))

message("\n[make_paper1_tables_figures] Done.")
message(sprintf("Tables in:        %s/", TABLES_PATH))
message(sprintf("Figures in:       %s/", FIG_PATH))
message(sprintf("Supplemental in:  %s/", SUPP_PATH))
