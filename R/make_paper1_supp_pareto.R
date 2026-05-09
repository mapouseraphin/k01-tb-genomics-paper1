# =============================================================================
# make_paper1_supp_pareto.R
#
# Build supplemental Pareto stringency-sensitivity outputs for Paper 1.
# Reads from BOTH per-tag LCA trees (lex canonical + Pareto sensitivity) and
# writes into the canonical lex supplemental directory so the supplement
# travels with the canonical manuscript.
#
# Run AFTER both LCA pipelines have completed (lex + Pareto), e.g. as the
# final source() call inside run_paper1_dual.R. Active TAG must be lex
# (canonical) when this script runs; it resolves the Pareto tree
# explicitly via build_paths_for_spec().
#
# Outputs (under PATHS$supplemental, i.e. lex tag):
#   supp_table_S6_pareto_lca_posterior.csv
#   supp_table_S7_lex_vs_pareto_lca.csv
#   supp_fig_S5_pareto_lca_sens_spec.pdf
#   supp_fig_S5b_lex_vs_pareto_sens_spec.pdf
#   supp_pareto_log.txt
#
# Robust to:
#   - Pareto LCA tree missing or empty (script exits cleanly with a logged
#     reason; lex outputs unaffected).
#   - fit_strata_summary.csv schema differences (uses select(any_of(...))).
#   - data-thin cells in the Pareto fit (footnotes in S6/S7 flag them).
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
  library(ggplot2)
})

# ---- Resolve both per-method paths -----------------------------------------
# Active spec must be lex for SUPP_OUT to point at the canonical supplement.
paths_lex <- build_paths_for_spec(
  snv_only         = SNV_ONLY,
  drop_ppe_cal     = DROP_PPE_CAL,
  drop_ppe_apply   = DROP_PPE_APPLY,
  selection_method = "lex"
)
paths_par <- build_paths_for_spec(
  snv_only         = SNV_ONLY,
  drop_ppe_cal     = DROP_PPE_CAL,
  drop_ppe_apply   = DROP_PPE_APPLY,
  selection_method = "pareto"
)

SUPP_OUT <- paths_lex$supplemental
dir.create(SUPP_OUT, recursive = TRUE, showWarnings = FALSE)

# ---- Logger ----------------------------------------------------------------
LOG <- character()
log_msg <- function(...) {
  args <- list(...)
  s <- if (length(args) == 1) as.character(args[[1]]) else do.call(sprintf, args)
  message(s); LOG <<- c(LOG, s); invisible(s)
}
write_log <- function() {
  writeLines(LOG, file.path(SUPP_OUT, "supp_pareto_log.txt"))
}

log_msg("[make_paper1_supp_pareto] Active spec (canonical): %s", paths_lex$lca_tag)
log_msg("                          Sensitivity spec       : %s", paths_par$lca_tag)
log_msg("                          Supplement output      : %s", SUPP_OUT)

# ---- Locate Pareto LCA fit summary -----------------------------------------
par_summ_path <- file.path(paths_par$lca_fits, "fit_strata_summary.csv")
lex_summ_path <- file.path(paths_lex$lca_fits, "fit_strata_summary.csv")

if (!file.exists(par_summ_path)) {
  log_msg("Pareto LCA summary not found at: %s", par_summ_path)
  log_msg("Pareto sensitivity outputs will NOT be built. Lex outputs unaffected.")
  log_msg("To produce these outputs, run the Pareto pipeline through fit_lca_models.R first.")
  write_log()
  invisible(return(NULL))
}
if (!file.exists(lex_summ_path)) {
  log_msg("Lex LCA summary not found at: %s", lex_summ_path)
  log_msg("Cannot build comparison Table S7 without lex summary. Aborting.")
  write_log()
  invisible(return(NULL))
}

par_summ <- read_csv(par_summ_path, show_col_types = FALSE)
lex_summ <- read_csv(lex_summ_path, show_col_types = FALSE)

log_msg("Loaded fit_strata_summary.csv from both trees:")
log_msg("  lex   : %d rows", nrow(lex_summ))
log_msg("  pareto: %d rows", nrow(par_summ))

# ---- Sanity: required columns ----------------------------------------------
required_cols <- c("param", "slot", "maf_bin",
                   "median", "q2.5", "q97.5", "rhat", "ess_bulk")
missing_lex <- setdiff(required_cols, names(lex_summ))
missing_par <- setdiff(required_cols, names(par_summ))
if (length(missing_lex) || length(missing_par)) {
  log_msg("Schema mismatch in fit_strata_summary.csv:")
  if (length(missing_lex)) log_msg("  lex    missing: %s", paste(missing_lex, collapse = ", "))
  if (length(missing_par)) log_msg("  pareto missing: %s", paste(missing_par, collapse = ", "))
  log_msg("Aborting; verify fit_lca_models.R wrote expected columns.")
  write_log()
  invisible(return(NULL))
}

# ---- Tag data-thin cells from the M3a identification check ------------------
# These are slot x MAF cells with < 10 replicate disagreements in the Pareto
# fitting set; their posteriors are wide and prior-influenced. We do not
# auto-detect from the summary (no n_disagree column); we hard-code the cells
# observed in the Pareto run log. Adjust here if a future Pareto rerun changes
# the data-thin cell list.
PARETO_DATA_THIN_CELLS <- tribble(
  ~slot,   ~maf_bin,
  "EM",    "[0.05,0.10)",
  "S0.1",  "[0.05,0.10)",
  "S0.2",  "[0.05,0.10)"
  # Add more rows if future runs flag additional <10-disagreement cells.
)

mark_thin <- function(df) {
  df %>%
    left_join(PARETO_DATA_THIN_CELLS %>% mutate(data_thin = TRUE),
              by = c("slot", "maf_bin")) %>%
    mutate(data_thin = ifelse(is.na(data_thin), FALSE, TRUE))
}

# ---- TABLE S6: Pareto M3a posterior summary --------------------------------
fmt_par <- function(df) {
  df %>%
    mutate(
      posterior_median = round(median, 3),
      CrI_lower        = round(q2.5,   3),
      CrI_upper        = round(q97.5,  3),
      rhat             = round(rhat,   4),
      ess_bulk         = round(ess_bulk, 0)
    )
}

# Sens/spec rows (per slot x MAF bin)
ssens_par <- par_summ %>%
  filter(param %in% c("sens", "spec")) %>%
  fmt_par() %>%
  mark_thin() %>%
  transmute(
    parameter = ifelse(param == "sens", "Sensitivity", "Specificity"),
    slot, maf_bin,
    posterior_median, CrI_lower, CrI_upper, rhat, ess_bulk,
    note = ifelse(data_thin,
                  "data-thin (n_disagree < 10); posterior prior-influenced",
                  "")
  ) %>%
  arrange(parameter, slot, maf_bin)

# Prevalence + sigma_u rows
extras_par <- par_summ %>%
  filter(param %in% c("pi", "sigma_u")) %>%
  fmt_par() %>%
  transmute(
    parameter = case_when(
      param == "pi"      ~ "Prevalence (pi)",
      param == "sigma_u" ~ "Patient SD (sigma_u, logit scale)"
    ),
    slot     = "(overall)",
    maf_bin  = "(overall)",
    posterior_median, CrI_lower, CrI_upper, rhat, ess_bulk,
    note = ifelse(parameter == "Prevalence (pi)" & CrI_upper - CrI_lower > 0.5,
                  "uninformative (CrI width > 0.5; not transportable as prior)",
                  "")
  )

table_S6 <- bind_rows(extras_par, ssens_par)
write_csv(table_S6, file.path(SUPP_OUT, "supp_table_S6_pareto_lca_posterior.csv"))
log_msg("Wrote: supp_table_S6_pareto_lca_posterior.csv  (%d rows)", nrow(table_S6))

# ---- TABLE S7: lex vs Pareto side-by-side ----------------------------------
mk_side <- function(df, lbl) {
  df %>%
    filter(param %in% c("sens", "spec")) %>%
    transmute(
      parameter = ifelse(param == "sens", "Sensitivity", "Specificity"),
      slot, maf_bin,
      !!paste0("median_", lbl)   := round(median, 3),
      !!paste0("CrI_low_", lbl)  := round(q2.5,   3),
      !!paste0("CrI_high_", lbl) := round(q97.5,  3)
    )
}

# Lex and Pareto bin labels may differ in their renderings; normalize before join
normalize_bin <- function(x) {
  x <- gsub("\\s+", "", x)
  x <- gsub("0\\.5\\]",  "0.45]", x)  # tolerate 0.45 vs 0.50 upper-bin labels
  x
}

lex_side <- mk_side(lex_summ, "lex") %>%
  mutate(maf_bin = normalize_bin(maf_bin))
par_side <- mk_side(par_summ, "par") %>%
  mutate(maf_bin = normalize_bin(maf_bin))

table_S7 <- lex_side %>%
  full_join(par_side, by = c("parameter", "slot", "maf_bin")) %>%
  left_join(PARETO_DATA_THIN_CELLS %>%
              mutate(maf_bin = normalize_bin(maf_bin),
                     pareto_data_thin = TRUE),
            by = c("slot", "maf_bin")) %>%
  mutate(
    pareto_data_thin = ifelse(is.na(pareto_data_thin), FALSE, TRUE),
    delta_median     = round(median_lex - median_par, 3),
    note = case_when(
      pareto_data_thin ~ "Pareto cell data-thin; delta dominated by Pareto prior",
      TRUE             ~ ""
    )
  ) %>%
  arrange(parameter, slot, maf_bin)

write_csv(table_S7, file.path(SUPP_OUT, "supp_table_S7_lex_vs_pareto_lca.csv"))
log_msg("Wrote: supp_table_S7_lex_vs_pareto_lca.csv     (%d rows)", nrow(table_S7))

# ---- FIGURE S5: Pareto sens/spec by slot x MAF -----------------------------
fig_data <- par_summ %>%
  filter(param %in% c("sens", "spec")) %>%
  mark_thin() %>%
  mutate(
    parameter = factor(ifelse(param == "sens", "Sensitivity", "Specificity"),
                       levels = c("Sensitivity", "Specificity")),
    slot      = factor(slot, levels = c("EM", "S0.1", "S0.2"))
  )

# Explicit ordering for x-axis (alphabetical sort would put [0.10,...) first)
maf_levels_present <- unique(fig_data$maf_bin)
maf_order <- c("[0.05,0.10)", "[0.05, 0.10)",
               "[0.10,0.25)", "[0.10, 0.25)",
               "[0.25,0.5]",  "[0.25, 0.5]",
               "[0.25,0.45]", "[0.25, 0.45]")
maf_order <- intersect(maf_order, maf_levels_present)
if (length(maf_order) < length(maf_levels_present)) {
  # Fallback: keep whatever ordering is in the data
  maf_order <- maf_levels_present
}
fig_data$maf_bin <- factor(fig_data$maf_bin, levels = maf_order)

p_S5 <- ggplot(fig_data,
               aes(x = maf_bin, y = median, ymin = q2.5, ymax = q97.5,
                   color = slot, group = slot)) +
  geom_pointrange(position = position_dodge(0.45), size = 0.5) +
  geom_line(position = position_dodge(0.45), alpha = 0.4) +
  # Mark data-thin cells
  geom_point(data = fig_data %>% filter(data_thin),
             aes(x = maf_bin, y = median),
             position = position_dodge(0.45),
             shape = 8, size = 2.6, color = "black",
             show.legend = FALSE) +
  facet_wrap(~ parameter) +
  scale_color_manual(values = c(EM = "#0072B2", S0.1 = "#D55E00", S0.2 = "#009E73")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "MAF bin (M0)", y = "Posterior median (95% CrI)",
       color = "Slot",
       title = "Pareto stringency-sensitivity LCA (M3a)",
       subtitle = sprintf(
         "Strict primary: DP >= 80, AD1 >= 6, MAF in [0.20, %.2f]; n_evaluable = 5 pairs",
         MAF_MAX),
       caption = paste(
         "Black asterisks mark data-thin cells (replicate disagreements < 10);",
         "their posteriors are wide and prior-influenced.",
         "Stringency sensitivity to lex canonical (DP >= 40, AD1 >= 6, MAF >= 0.02).",
         sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 8, color = "gray30"))

ggsave(file.path(SUPP_OUT, "supp_fig_S5_pareto_lca_sens_spec.pdf"),
       p_S5, width = 8, height = 4.8)
log_msg("Wrote: supp_fig_S5_pareto_lca_sens_spec.pdf")

# ---- FIGURE S5b: lex vs Pareto direct overlay ------------------------------
side_data <- bind_rows(
  lex_summ %>% filter(param %in% c("sens", "spec")) %>% mutate(method = "lex (canonical)"),
  par_summ %>% filter(param %in% c("sens", "spec")) %>% mutate(method = "pareto (sensitivity)")
) %>%
  mutate(
    parameter = factor(ifelse(param == "sens", "Sensitivity", "Specificity"),
                       levels = c("Sensitivity", "Specificity")),
    slot      = factor(slot, levels = c("EM", "S0.1", "S0.2")),
    maf_bin   = normalize_bin(maf_bin)
  )

# Common MAF ordering
maf_order_b <- c("[0.05,0.10)", "[0.10,0.25)", "[0.25,0.45]", "[0.25,0.5]")
maf_order_b <- intersect(maf_order_b, unique(side_data$maf_bin))
side_data$maf_bin <- factor(side_data$maf_bin, levels = maf_order_b)

p_S5b <- ggplot(side_data,
                aes(x = maf_bin, y = median, ymin = q2.5, ymax = q97.5,
                    color = method, group = method)) +
  geom_pointrange(position = position_dodge(0.55), size = 0.45) +
  geom_line(position = position_dodge(0.55), alpha = 0.45) +
  facet_grid(parameter ~ slot) +
  scale_color_manual(values = c("lex (canonical)" = "#0072B2",
                                "pareto (sensitivity)" = "#D55E00")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "MAF bin (M0)", y = "Posterior median (95% CrI)",
       color = "Method",
       title = "Lex canonical vs Pareto stringency-sensitivity (M3a)",
       caption = paste(
         "Lex canonical: DP >= 40, AD1 >= 6, MAF >= 0.02.",
         "Pareto sensitivity: DP >= 80, AD1 >= 6, MAF >= 0.20.",
         "Pareto cells in [0.05,0.10) are data-thin (omit/footnote in main text).",
         sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 8, color = "gray30"))

ggsave(file.path(SUPP_OUT, "supp_fig_S5b_lex_vs_pareto_sens_spec.pdf"),
       p_S5b, width = 9, height = 6)
log_msg("Wrote: supp_fig_S5b_lex_vs_pareto_sens_spec.pdf")

# ---- Console summary -------------------------------------------------------
log_msg("")
log_msg("[make_paper1_supp_pareto] Done. Supplemental Pareto outputs:")
log_msg("  %s/", SUPP_OUT)
log_msg("    supp_table_S6_pareto_lca_posterior.csv")
log_msg("    supp_table_S7_lex_vs_pareto_lca.csv")
log_msg("    supp_fig_S5_pareto_lca_sens_spec.pdf")
log_msg("    supp_fig_S5b_lex_vs_pareto_sens_spec.pdf")
log_msg("    supp_pareto_log.txt")

write_log()
