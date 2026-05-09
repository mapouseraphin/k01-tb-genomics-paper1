# =============================================================================
# 07_make_tables_figures.R
#
# Produce all main-text and supplementary tables and figures for Paper 1
# from the pipeline outputs.
#
# IMPORTANT: This script writes to PATHS$tables and PATHS$figures, which are
# now spec-tagged (see 00_pipeline_config.R). Running under a different spec
# writes to a different output subdirectory; previous spec runs are NOT
# overwritten.
#
# Reads:
#   PATHS$meta/{gh,fl}_meta.rds                                (Table 1)
#   PATHS$calibration/cal_lexicographic_thresholds.csv         (Table 2)
#   PATHS$calibration/cal_pairs.rds                            (Table 3)
#   PATHS$calibration/evaluable_pairs_at_primary_by_slot.csv   (Table 2 footnote)
#   PATHS$calibration/cal_lexicographic_heatmap.pdf            (Figure 2)
#   PATHS$external/external_validity_main.csv                  (Table 4, 5)
#   PATHS$external/external_validity_grid.csv                  (Table S1, Figure S1)
#   PATHS$external/bootstrap_summary.csv                       (Table 4)
#   PATHS$external/{fl,gh}_pairs.rds                           (Figure 4)
#   PATHS$diagnostics/q1_*.csv                                 (Table S4, Figure S2)
#   PATHS$diagnostics/q2_lineage_stratified_jaccard.csv        (Table S3)
#   PATHS$diagnostics/q5_gh_cluster_def_sensitivity.csv        (Table S2)
#
# Writes (under PATHS$tables and PATHS$figures, both spec-tagged):
#   table1_cohort_characteristics.{csv,md}
#   table2_calibrated_thresholds.{csv,md}
#   table3_internal_validity.{csv,md}
#   table4_external_validity_main.{csv,md}      [NOW WITH PPV / enrichment]
#   table5_threshold_ladder_pr.{csv,md}
#   tableS1_sensitivity_grid.{csv,md}            [NOW WITH PPV / enrichment]
#   tableS2_gh_cluster_def.{csv,md}
#   tableS3_lineage_stratified.{csv,md}
#   tableS4_position_concentration.{csv,md}
#   figure3_pr_ladder_forest.{pdf,png}
#   figure4_jaccard_distributions.{pdf,png}
#   figureS1_sensitivity_grid_forest.{pdf,png}
#   figureS2_position_concentration.{pdf,png}
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

suppressPackageStartupMessages({
  library(tidyverse)
  library(scales)
})

# Ggplot theme
theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "gray92"),
      panel.grid.major.y = element_line(color = "gray92"),
      strip.text         = element_text(face = "bold"),
      plot.title         = element_text(face = "bold"),
      plot.caption       = element_text(size = base_size - 2,
                                        color = "gray40", hjust = 0)
    )
}

# Helper: write CSV + markdown side-by-side
write_table <- function(df, slug, caption = NULL) {
  csv_path <- file.path(PATHS$tables, paste0(slug, ".csv"))
  md_path  <- file.path(PATHS$tables, paste0(slug, ".md"))

  readr::write_csv(df, csv_path)

  md_lines <- character(0)
  if (!is.null(caption)) md_lines <- c(md_lines, paste0("**", caption, "**"), "")
  md_lines <- c(md_lines, knitr::kable(df, format = "pipe"))
  writeLines(md_lines, md_path)

  message(sprintf("  wrote %s (%d rows)", slug, nrow(df)))
}

# Helper: format estimate with CI
fmt_est_ci <- function(est, lo, hi, digits = 2) {
  ifelse(is.na(est), "—",
         sprintf(sprintf("%%.%df (%%.%df, %%.%df)", digits, digits, digits),
                 est, lo, hi))
}

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "—",
         sprintf(sprintf("%%.%df%%%%", digits), 100 * x))
}

# Caption suffix describing canonical run
spec_caption <- sprintf(
  "Canonical specification: SNV_ONLY=%s, DROP_PPE_CAL=%s, DROP_PPE_APPLY=%s, M0_ONLY=%s. Tag: %s.",
  SNV_ONLY, DROP_PPE_CAL, DROP_PPE_APPLY, M0_ONLY, TAG$apply_tag
)

# =============================================================================
# Table 1 -- cohort characteristics (patient-level for Ghana, sample-level
#            for Florida by convention; Ghana also reports sample-level n).
# =============================================================================
make_table1 <- function() {
  message("\n[Table 1] Cohort characteristics")

  gh_meta <- readRDS(file.path(PATHS$meta, "gh_meta.rds"))
  fl_meta <- readRDS(file.path(PATHS$meta, "fl_meta.rds"))

  # Ghana: collapse to patient-level for demographics
  gh_pt <- gh_meta %>%
    group_by(patientId) %>%
    summarise(
      Age       = first(Age),
      Sex       = first(Sex),
      HIV       = first(HIV),
      lineage   = first(lineage),
      depth_med = median(depth_med, na.rm = TRUE),
      n_samples = n(),
      .groups   = "drop"
    )

  fmt_n_pct <- function(x, level) {
    n <- sum(x == level, na.rm = TRUE)
    pct <- 100 * n / sum(!is.na(x))
    sprintf("%d (%.1f%%)", n, pct)
  }

  fmt_med_iqr <- function(x) {
    sprintf("%.0f (%.0f, %.0f)",
            median(x, na.rm = TRUE),
            quantile(x, 0.25, na.rm = TRUE),
            quantile(x, 0.75, na.rm = TRUE))
  }

  rows <- tibble(
    Characteristic = c(
      "Patients (n)",
      "Samples (n)",
      "Age, years (median, IQR)",
      "Sex: Female",
      "Sex: Male",
      "HIV: Negative",
      "HIV: Positive",
      "Lineage: L4",
      "Lineage: non-L4",
      "Median sequencing depth, x (median, IQR)"
    ),
    Ghana = c(
      sprintf("%d", n_distinct(gh_meta$patientId)),
      sprintf("%d", nrow(gh_meta)),
      fmt_med_iqr(gh_pt$Age),
      fmt_n_pct(gh_pt$Sex, "Female"),
      fmt_n_pct(gh_pt$Sex, "Male"),
      fmt_n_pct(gh_pt$HIV, "HIV-"),
      fmt_n_pct(gh_pt$HIV, "HIV+"),
      fmt_n_pct(gh_pt$lineage, "L4"),
      sprintf("%d (%.1f%%)",
              sum(gh_pt$lineage != "L4" & !is.na(gh_pt$lineage)),
              100 * mean(gh_pt$lineage != "L4", na.rm = TRUE)),
      fmt_med_iqr(gh_meta$depth_med)
    ),
    Florida = c(
      sprintf("%d", nrow(fl_meta)),
      sprintf("%d", nrow(fl_meta)),
      if ("Age" %in% names(fl_meta) && any(!is.na(fl_meta$Age)))
        fmt_med_iqr(fl_meta$Age) else "—",
      if ("Sex" %in% names(fl_meta)) fmt_n_pct(fl_meta$Sex, "Female") else "—",
      if ("Sex" %in% names(fl_meta)) fmt_n_pct(fl_meta$Sex, "Male") else "—",
      if ("HIV" %in% names(fl_meta)) fmt_n_pct(fl_meta$HIV, "HIV-") else "—",
      if ("HIV" %in% names(fl_meta)) fmt_n_pct(fl_meta$HIV, "HIV+") else "—",
      fmt_n_pct(fl_meta$lineage, "L4"),
      fmt_n_pct(fl_meta$lineage, "L_non_L4"),
      fmt_med_iqr(fl_meta$depth_med)
    )
  )

  write_table(rows, "table1_cohort_characteristics",
              caption = paste(
                "Table 1. Cohort characteristics, Ghana and Florida.",
                "Demographics for Ghana are summarized at the patient level;",
                "sequencing depth is summarized over all included samples."))
}

# =============================================================================
# Table 2 -- calibrated thresholds + diagnostic columns
# =============================================================================
make_table2 <- function() {
  message("\n[Table 2] Calibrated thresholds")

  thr <- readr::read_csv(file.path(PATHS$calibration,
                                   "cal_lexicographic_thresholds.csv"),
                         show_col_types = FALSE)

  out <- thr %>%
    transmute(
      Selection            = factor(selection,
                                    levels = c("loose", "primary", "tighter")),
      `DP_min`             = DP_min,
      `AD1_min`            = AD1_min,
      `MAF_min`            = MAF_min,
      `MAF_max`            = MAF_MAX,
      `n_eval_pairs`       = pairs_evaluable,
      `prop_both0`         = round(prop_both0, 3),
      `median_Jaccard`     = round(jacc_med_eval, 3),
      `mean_confirm_rate`  = round(conf_med_eval, 3)
    ) %>%
    arrange(Selection)

  write_table(out, "table2_calibrated_thresholds",
              caption = paste("Table 2. Calibrated decision rules for iSNV detection.",
                              spec_caption,
                              "Loose is a no-filter sentinel; primary and tighter are grid-picked under lexicographic prioritization."))
}

# =============================================================================
# Table 3 -- internal validity: within-visit replicate concordance, by pair type
# =============================================================================
make_table3 <- function() {
  message("\n[Table 3] Internal validity by pair type")

  cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))
  thr <- readRDS(file.path(PATHS$calibration, "thresholds.rds"))
  prim <- thr$primary

  pp <- cal_pairs %>%
    filter(DP_min == prim$DP_min, AD1_min == prim$AD1_min,
           MAF_min == prim$MAF_min) %>%
    mutate(slot_a = parse_slot(sampleA),
           slot_b = parse_slot(sampleB),
           pair_type = case_when(
             (slot_a == "EM"   & slot_b %in% c("S0.1","S0.2")) |
             (slot_b == "EM"   & slot_a %in% c("S0.1","S0.2"))   ~ "EM<->morning_spot",
             (slot_a == "S0.1" & slot_b == "S0.2") |
             (slot_a == "S0.2" & slot_b == "S0.1")               ~ "S0.1<->S0.2",
             TRUE                                                  ~ "other"
           ))

  out <- pp %>%
    group_by(pair_type) %>%
    summarise(
      n_pairs           = n(),
      n_evaluable       = sum(nA > 0 & nB > 0),
      median_jaccard    = round(median(jaccard[nA > 0 & nB > 0], na.rm = TRUE), 3),
      iqr_low_jaccard   = round(quantile(jaccard[nA > 0 & nB > 0], 0.25, na.rm = TRUE), 3),
      iqr_high_jaccard  = round(quantile(jaccard[nA > 0 & nB > 0], 0.75, na.rm = TRUE), 3),
      mean_confirm_both = round(mean(confirm_both[nA > 0 & nB > 0], na.rm = TRUE), 3),
      .groups           = "drop"
    ) %>%
    arrange(pair_type)

  write_table(out, "table3_internal_validity",
              caption = paste("Table 3. Internal validity: within-visit replicate concordance at primary thresholds, by sample-slot pair type.",
                              spec_caption,
                              "Evaluable pairs are those with at least one calibrated iSNV in both replicates."))
}

# =============================================================================
# Table 4 -- main external validity (canonical primary, FL + GH 12 SNP-only)
#
# UPDATED: now includes within-cluster prevalence, PPV, and PPV enrichment
# (PPV / prevalence). Bootstrap 95% CIs reported for PR, PPV, and enrichment.
# =============================================================================
make_table4 <- function() {
  message("\n[Table 4] External validity, canonical primary")

  main <- readr::read_csv(file.path(PATHS$external,
                                    "external_validity_main.csv"),
                          show_col_types = FALSE)

  boot_path <- file.path(PATHS$external, "bootstrap_summary.csv")
  if (file.exists(boot_path)) {
    boot <- readr::read_csv(boot_path, show_col_types = FALSE)
  } else {
    boot <- tibble(cell = character(0))
  }

  # Restrict to primary tier for Table 4
  primary_only <- main %>% filter(threshold_tier == "primary")

  # Attach bootstrap CIs
  primary_only <- primary_only %>%
    mutate(
      cell = case_when(
        cohort == "Florida" & cluster_def == "fl_surveillance" ~ "FL_surveillance_primary",
        cohort == "Ghana"   & cluster_def == "gh_12_snp_only"  ~ "GH_12_snp_only_primary",
        TRUE                                                    ~ NA_character_
      )
    ) %>%
    left_join(boot, by = "cell")

  out <- primary_only %>%
    transmute(
      Cohort               = cohort,
      `Cluster definition` = case_when(
        cluster_def == "fl_surveillance" ~ "FL surveillance clusters",
        cluster_def == "gh_12_snp_only"  ~ "<=12 SNP-only (Ghana)"
      ),
      `n within`           = n_within,
      `n between`          = n_between,
      `Within-cluster prevalence` = round(prevalence, 3),
      `Sensitivity`        = round(sens, 3),
      `Specificity`        = round(spec, 3),
      `PPV`                = round(ppv, 3),
      `Enrichment (PPV/prev)` = round(enrichment, 2),
      `Mean Jaccard within`  = round(jaccard_within, 3),
      `Mean Jaccard between` = round(jaccard_between, 3),
      `PR (analytic 95% CI)` = fmt_est_ci(pr, pr_lo, pr_hi, 3),
      `PR (bootstrap 95% CI)` = ifelse(is.na(pr_boot_lo), "—",
                                       sprintf("(%.3f, %.3f)",
                                               pr_boot_lo, pr_boot_hi)),
      `PPV (bootstrap 95% CI)` = ifelse(is.na(ppv_boot_lo), "—",
                                        sprintf("(%.3f, %.3f)",
                                                ppv_boot_lo, ppv_boot_hi)),
      `Enrichment (bootstrap 95% CI)` = ifelse(is.na(enr_boot_lo), "—",
                                               sprintf("(%.2f, %.2f)",
                                                       enr_boot_lo, enr_boot_hi))
    )

  write_table(out, "table4_external_validity_main",
              caption = paste("Table 4. External validity at primary thresholds.",
                              spec_caption,
                              "Florida uses surveillance cluster membership as reference; Ghana uses pairwise SNP distance <=12 (Asare et al. 2020).",
                              "PR is from modified Poisson regression with dyadic-clustered standard errors.",
                              "PPV = P(same cluster | shared iSNV); Enrichment = PPV / within-cluster prevalence (>1: test enriches above baseline; ~1: negative-control behavior).",
                              "Bootstrap 95% CIs use the percentile method on cluster-resampled sample IDs."))
}

# =============================================================================
# Table 5 -- threshold-ladder PR (loose/primary/tighter)
# =============================================================================
make_table5 <- function() {
  message("\n[Table 5] Threshold-ladder PR")

  main <- readr::read_csv(file.path(PATHS$external,
                                    "external_validity_main.csv"),
                          show_col_types = FALSE)

  out <- main %>%
    mutate(threshold_tier = factor(threshold_tier,
                                   levels = c("loose","primary","tighter"))) %>%
    arrange(cohort, threshold_tier) %>%
    transmute(
      Cohort         = cohort,
      `Cluster definition` = case_when(
        cluster_def == "fl_surveillance" ~ "FL surveillance",
        cluster_def == "gh_12_snp_only"  ~ "<=12 SNP-only (Ghana)"
      ),
      Threshold      = threshold_tier,
      `n within`     = n_within,
      `n between`    = n_between,
      Sensitivity    = round(sens, 3),
      Specificity    = round(spec, 3),
      PPV            = round(ppv, 3),
      Enrichment     = round(enrichment, 2),
      `PR (95% CI)`  = fmt_est_ci(pr, pr_lo, pr_hi, 3)
    )

  write_table(out, "table5_threshold_ladder_pr",
              caption = paste("Table 5. PR across the threshold ladder.",
                              spec_caption,
                              "Under classical measurement-error theory with non-differential misclassification, PR is expected to increase as thresholds tighten (PR_loose <= PR_primary <= PR_tighter) in settings with genuine transmission signal.",
                              "PPV and enrichment columns added to track whether tightening also moves predictive value above baseline."))
}

# =============================================================================
# Supplementary tables
# =============================================================================
make_tableS1 <- function() {
  message("\n[Table S1] Full sensitivity grid")
  grid <- readr::read_csv(file.path(PATHS$external,
                                    "external_validity_grid.csv"),
                          show_col_types = FALSE)

  out <- grid %>%
    transmute(
      Cohort           = cohort,
      Threshold        = threshold_tier,
      `Cluster def`    = cluster_def,
      `Zero-iSNV dropped` = zero_drop,
      `Apply tag`      = apply_tag,
      `n within`       = n_within,
      `n between`      = n_between,
      `Within-cluster prevalence` = round(prevalence, 3),
      Sensitivity      = round(sens, 3),
      Specificity      = round(spec, 3),
      PPV              = round(ppv, 3),
      NPV              = round(npv, 3),
      Enrichment       = round(enrichment, 2),
      `Jaccard within` = round(jaccard_within, 3),
      `Jaccard between` = round(jaccard_between, 3),
      `PR (95% CI)`    = fmt_est_ci(pr, pr_lo, pr_hi, 3),
      TP               = TP,
      FP               = FP,
      FN               = FN,
      TN               = TN
    )

  write_table(out, "tableS1_sensitivity_grid",
              caption = paste("Table S1. Full sensitivity grid: cohort x threshold x cluster definition x zero-iSNV handling.",
                              spec_caption,
                              "PPV = P(same cluster | shared iSNV); NPV = P(different cluster | no shared iSNV); Enrichment = PPV / within-cluster prevalence."))
}

make_tableS2 <- function() {
  message("\n[Table S2] GH cluster-definition sensitivity")
  src <- file.path(PATHS$diagnostics, "q5_gh_cluster_def_sensitivity.csv")
  if (!file.exists(src)) {
    message("  skipped (q5 output not found; run 06 first)"); return(invisible())
  }
  q5 <- readr::read_csv(src, show_col_types = FALSE)

  out <- q5 %>%
    transmute(
      `Cluster definition` = definition,
      `n within`           = n_within,
      `n between`          = n_between,
      Sensitivity          = round(sens, 3),
      Specificity          = round(spec, 3),
      `Jaccard within`     = round(jaccard_within, 3),
      `Jaccard between`    = round(jaccard_between, 3),
      `PR (95% CI)`        = fmt_est_ci(pr, pr_lo, pr_hi, 3)
    )

  write_table(out, "tableS2_gh_cluster_def",
              caption = "Table S2. Ghana cluster-definition sensitivity. PR estimated via modified Poisson regression with dyadic-clustered standard errors.")
}

make_tableS3 <- function() {
  message("\n[Table S3] Lineage-stratified Jaccard")
  src <- file.path(PATHS$diagnostics, "q2_lineage_stratified_jaccard.csv")
  if (!file.exists(src)) {
    message("  skipped (q2 output not found; run 06 first)"); return(invisible())
  }
  q2 <- readr::read_csv(src, show_col_types = FALSE)

  out <- q2 %>%
    transmute(
      `Cluster definition` = definition,
      `Cluster class`      = cluster_class,
      `Lineage pair`       = lineage_pair,
      `n pairs`            = n_pairs,
      `Mean Jaccard`       = round(mean_jaccard, 3),
      `% pairs with shared iSNV` = round(pct_shared, 1)
    )

  write_table(out, "tableS3_lineage_stratified",
              caption = "Table S3. Lineage stratification of pair-level Jaccard concordance, Ghana (negative-control supporting evidence).")
}

make_tableS4 <- function() {
  message("\n[Table S4] Position concentration in between-cluster pool")
  src_files <- c(
    GH12 = file.path(PATHS$diagnostics, "q1_gh12_between_position_freq.csv"),
    GH5  = file.path(PATHS$diagnostics, "q1_gh5_between_position_freq.csv"),
    FL   = file.path(PATHS$diagnostics, "q1_fl_between_position_freq.csv")
  )
  if (any(!file.exists(src_files))) {
    message("  skipped (q1 outputs not all present; run 06 first)")
    return(invisible())
  }

  read_q1 <- function(p, src_label) {
    readr::read_csv(p, show_col_types = FALSE) %>%
      mutate(source = src_label)
  }

  q1_all <- bind_rows(
    read_q1(src_files["GH12"], "GH <=12 SNP between"),
    read_q1(src_files["GH5"],  "GH <=5 SNP between"),
    read_q1(src_files["FL"],   "FL surveillance between")
  )

  # Top-K cumulative concentration table
  out <- q1_all %>%
    group_by(source) %>%
    summarise(
      `Distinct positions in between-cluster pool` = max(rank, na.rm = TRUE),
      `% in top 5`  = round(100 * cum_pct[rank == 5][1], 1),
      `% in top 10` = round(100 * cum_pct[rank == 10][1], 1),
      `% in top 20` = round(100 * cum_pct[rank == 20][1], 1),
      `% in top 50` = round(100 * cum_pct[rank == 50][1], 1),
      .groups = "drop"
    )

  write_table(out, "tableS4_position_concentration",
              caption = "Table S4. Position-level concentration in the between-cluster shared-iSNV pool, by setting and cluster definition. Reading: in Ghana, a small number of recurrent positions account for the majority of between-cluster sharing, consistent with lineage-typical homoplasy under weak transmission signal (negative-control).")
}

# =============================================================================
# Figure 3 -- threshold-ladder PR forest plot
# =============================================================================
make_figure3 <- function() {
  message("\n[Figure 3] Threshold-ladder PR forest plot")

  main <- readr::read_csv(file.path(PATHS$external,
                                    "external_validity_main.csv"),
                          show_col_types = FALSE)

  pdat <- main %>%
    mutate(
      threshold_tier = factor(threshold_tier,
                              levels = c("loose","primary","tighter")),
      cohort_label = case_when(
        cohort == "Florida"                              ~ "Florida\n(surveillance clusters)",
        cohort == "Ghana" & cluster_def == "gh_12_snp_only" ~ "Ghana\n(<=12 SNP-only)"
      )
    )

  p <- ggplot(pdat,
              aes(x = pr, y = threshold_tier, color = cohort_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
    geom_pointrange(aes(xmin = pr_lo, xmax = pr_hi),
                    size = 0.6, fatten = 4) +
    facet_wrap(~ cohort_label, ncol = 1, scales = "free_y") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)) +
    scale_color_manual(values = c("#0072B2", "#D55E00"), guide = "none") +
    labs(
      x = "Prevalence ratio (PR), log2 scale",
      y = NULL,
      title = "PR across the threshold ladder",
      caption = paste(
        "PR for shared calibrated iSNVs: same cluster vs different cluster.",
        "Reference line: PR = 1 (no enrichment).",
        sprintf("Tag: %s.", TAG$apply_tag)
      )
    ) +
    theme_paper()

  ggsave(file.path(PATHS$figures, "figure3_pr_ladder_forest.pdf"),
         p, width = 6, height = 5)
  ggsave(file.path(PATHS$figures, "figure3_pr_ladder_forest.png"),
         p, width = 6, height = 5, dpi = 300)
  message("  wrote figure3_pr_ladder_forest")
}

# =============================================================================
# Figure 4 -- within vs between Jaccard distribution by cohort
# =============================================================================
make_figure4 <- function() {
  message("\n[Figure 4] Jaccard distribution within vs between")

  fl_pairs_path <- file.path(PATHS$external, "fl_pairs.rds")
  gh_pairs_path <- file.path(PATHS$external, "gh_pairs.rds")

  if (!file.exists(fl_pairs_path) || !file.exists(gh_pairs_path)) {
    message("  skipped (fl_pairs / gh_pairs not present; run 04 first)")
    return(invisible())
  }

  fl_pairs <- readRDS(fl_pairs_path)
  gh_pairs <- readRDS(gh_pairs_path)

  fl_long <- fl_pairs %>%
    transmute(
      cohort = "Florida (surveillance clusters)",
      jaccard,
      cluster_class = if_else(same_cluster == 1, "within", "between")
    )

  gh_long <- gh_pairs %>%
    transmute(
      cohort = "Ghana (<=12 SNP-only)",
      jaccard,
      cluster_class = if_else(same_cluster_12_snp_only == 1, "within", "between")
    )

  pdat <- bind_rows(fl_long, gh_long) %>%
    mutate(cluster_class = factor(cluster_class,
                                  levels = c("between", "within")))

  p <- ggplot(pdat, aes(x = jaccard, fill = cluster_class)) +
    geom_density(alpha = 0.4, color = NA) +
    facet_wrap(~ cohort, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(between = "#999999", within = "#D55E00"),
                      name = "Cluster class") +
    labs(
      x = "Pair-level Jaccard concordance",
      y = "Density",
      title = "Pair-level Jaccard distribution: within vs between cluster",
      caption = paste(
        "Florida: within > between concordance reflects recovered transmission signal.",
        "Ghana: similar within and between distributions reflect graceful failure of the calibrated measurement under weak transmission signal (negative-control).",
        sprintf("Tag: %s.", TAG$apply_tag),
        sep = " "
      )
    ) +
    theme_paper()

  ggsave(file.path(PATHS$figures, "figure4_jaccard_distributions.pdf"),
         p, width = 7, height = 5.5)
  ggsave(file.path(PATHS$figures, "figure4_jaccard_distributions.png"),
         p, width = 7, height = 5.5, dpi = 300)
  message("  wrote figure4_jaccard_distributions")
}

# =============================================================================
# Figure S1 -- full sensitivity grid forest plot
# =============================================================================
make_figureS1 <- function() {
  message("\n[Figure S1] Sensitivity grid forest plot")

  grid <- readr::read_csv(file.path(PATHS$external,
                                    "external_validity_grid.csv"),
                          show_col_types = FALSE)

  pdat <- grid %>%
    filter(!is.na(pr)) %>%
    mutate(
      threshold_tier = factor(threshold_tier,
                              levels = c("loose","primary","tighter")),
      zero_lab = if_else(zero_drop, "Excl. zero-iSNV", "All samples"),
      cluster_def_label = case_when(
        cluster_def == "fl_surveillance"      ~ "FL surveillance",
        cluster_def == "gh_12_snp_only"       ~ "GH <=12 SNP-only",
        cluster_def == "gh_5_snp_only"        ~ "GH <=5 SNP-only",
        cluster_def == "gh_12_snp_plus_spat"  ~ "GH <=12 SNP + spatial",
        cluster_def == "gh_5_snp_plus_spat"   ~ "GH <=5 SNP + spatial"
      )
    )

  p <- ggplot(pdat, aes(x = pr, y = threshold_tier,
                        color = zero_lab, shape = zero_lab)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
    geom_pointrange(aes(xmin = pr_lo, xmax = pr_hi),
                    position = position_dodge(width = 0.5),
                    size = 0.4, fatten = 3) +
    facet_wrap(~ cluster_def_label, ncol = 1) +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.25, 0.5, 1.0, 2.0, 4.0)) +
    scale_color_manual(values = c("All samples" = "#0072B2",
                                  "Excl. zero-iSNV" = "#D55E00"),
                       name = NULL) +
    scale_shape_manual(values = c("All samples" = 16,
                                  "Excl. zero-iSNV" = 17),
                       name = NULL) +
    labs(
      x = "Prevalence ratio (PR), log2 scale",
      y = NULL,
      title = "Sensitivity grid: PR across threshold tier x cluster definition x zero-iSNV handling",
      caption = sprintf(
        "Each cell shows analytic 95%% CI from modified Poisson regression with dyadic-clustered SE. Tag: %s.",
        TAG$apply_tag
      )
    ) +
    theme_paper(base_size = 9)

  ggsave(file.path(PATHS$figures, "figureS1_sensitivity_grid_forest.pdf"),
         p, width = 8, height = 11)
  ggsave(file.path(PATHS$figures, "figureS1_sensitivity_grid_forest.png"),
         p, width = 8, height = 11, dpi = 300)
  message("  wrote figureS1_sensitivity_grid_forest")
}

# =============================================================================
# Figure S2 -- position concentration curves (06 Q1)
# =============================================================================
make_figureS2 <- function() {
  message("\n[Figure S2] Position concentration curves")

  src_files <- list(
    "GH <=12 SNP between"     = file.path(PATHS$diagnostics, "q1_gh12_between_position_freq.csv"),
    "GH <=5 SNP between"      = file.path(PATHS$diagnostics, "q1_gh5_between_position_freq.csv"),
    "FL surveillance between" = file.path(PATHS$diagnostics, "q1_fl_between_position_freq.csv")
  )

  if (any(!file.exists(unlist(src_files)))) {
    message("  skipped (q1 outputs not all present; run 06 first)")
    return(invisible())
  }

  read_q1 <- function(p, lbl) {
    readr::read_csv(p, show_col_types = FALSE) %>%
      mutate(source = lbl)
  }

  pdat <- bind_rows(lapply(names(src_files),
                            function(lbl) read_q1(src_files[[lbl]], lbl)))

  p <- ggplot(pdat %>% filter(rank <= 50),
              aes(x = rank, y = 100 * cum_pct,
                  color = source, group = source)) +
    geom_step(linewidth = 0.7) +
    geom_hline(yintercept = c(50, 70, 90),
               linetype = "dotted", color = "gray70") +
    scale_color_manual(
      values = c(
        "GH <=12 SNP between"     = "#D55E00",
        "GH <=5 SNP between"      = "#E69F00",
        "FL surveillance between" = "#0072B2"
      ),
      name = NULL
    ) +
    labs(
      x = "Rank of position (top-K)",
      y = "Cumulative % of all between-cluster shared-iSNV instances",
      title = "Position-level concentration in between-cluster shared-iSNV pool",
      caption = paste(
        "Steeper curves indicate stronger concentration: a few recurrent positions account for most between-cluster sharing.",
        "GH concentration > FL is consistent with lineage-typical homoplasy under weak transmission signal (negative-control).",
        sprintf("Tag: %s.", TAG$apply_tag),
        sep = " "
      )
    ) +
    theme_paper()

  ggsave(file.path(PATHS$figures, "figureS2_position_concentration.pdf"),
         p, width = 7, height = 5)
  ggsave(file.path(PATHS$figures, "figureS2_position_concentration.png"),
         p, width = 7, height = 5, dpi = 300)
  message("  wrote figureS2_position_concentration")
}

# =============================================================================
# Run all
# =============================================================================
required_inputs <- c(
  meta_gh        = file.path(PATHS$meta,        "gh_meta.rds"),
  meta_fl        = file.path(PATHS$meta,        "fl_meta.rds"),
  thresholds_csv = file.path(PATHS$calibration, "cal_lexicographic_thresholds.csv"),
  cal_pairs      = file.path(PATHS$calibration, "cal_pairs.rds"),
  thresholds_rds = file.path(PATHS$calibration, "thresholds.rds"),
  ext_main       = file.path(PATHS$external,    "external_validity_main.csv"),
  ext_grid       = file.path(PATHS$external,    "external_validity_grid.csv")
)
miss <- required_inputs[!file.exists(required_inputs)]
if (length(miss)) {
  cat("Missing required inputs:\n")
  for (n in names(miss)) cat(sprintf("  %-20s %s\n", n, miss[[n]]))
  stop("Run pipeline 00 -> 04 (and optionally 06) before 07.")
}

# Main tables
make_table1()
make_table2()
make_table3()
make_table4()
make_table5()

# Supplementary tables (06 outputs may be optional)
make_tableS1()
make_tableS2()
make_tableS3()
make_tableS4()

# Figures
make_figure3()
make_figure4()
make_figureS1()
make_figureS2()

# Note: Figure 1 (DAG) is not generated here -- typically authored externally
# in a vector graphics tool. Figure 2 (calibration heatmap) is produced by 02.

message(sprintf("\n[07_make_tables_figures] Done."))
message(sprintf("  Tables:  %s/", PATHS$tables))
message(sprintf("  Figures: %s/", PATHS$figures))
