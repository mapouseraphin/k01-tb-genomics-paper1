# =============================================================================
# 08_cross_spec_contrast.R
#
# Cross-spec contrast tables and figures. Reads outputs from each available
# specification subdirectory and produces side-by-side comparisons across
# the three Paper 1 analytic specifications:
#
#   Spec A — CANONICAL                  : SNV-only, PE/PPE retained
#   Spec B — INDEL-EFFECT SENSITIVITY   : with-indels, PE/PPE retained
#   Spec C — FULL NON-CANONICAL         : with-indels, PE/PPE excluded
#
# Spec A vs Spec B isolates the indel-exclusion effect (PE/PPE held constant).
# Spec B vs Spec C isolates the PE/PPE-handling effect (variant class held
# constant). Spec A vs Spec C is the full canonical-vs-non-canonical contrast.
#
# Outputs (PATHS$contrast_tables and PATHS$contrast_figures, NOT spec-tagged
# because they combine results from multiple specs):
#
#   tableC1_cross_spec_thresholds.{csv,md}
#       Calibrated DP/AD1/MAF thresholds at primary, by spec.
#
#   tableC2_cross_spec_external_validity_primary.{csv,md}
#       External validity at primary thresholds, by cohort x spec. Includes
#       PR (analytic + bootstrap), PPV, enrichment, sens, spec, prevalence.
#
#   tableC3_cross_spec_threshold_ladder.{csv,md}
#       Threshold ladder PR by cohort x spec x tier, plus PPV / enrichment.
#
#   figureC1_cross_spec_pr_forest.{pdf,png}
#       Forest plot of canonical-primary PR across the three specs, faceted
#       by cohort.
#
#   figureC2_cross_spec_threshold_ladder.{pdf,png}
#       Threshold-ladder PR forest plot, dodged by spec.
#
#   figureC3_cross_spec_ppv_enrichment.{pdf,png}
#       Forest plot of PPV enrichment (PPV / within-cluster prevalence) at
#       canonical primary, by cohort x spec.
#
# Behavior: for each spec, the script checks whether expected upstream
# artifacts exist; specs without complete outputs are silently skipped. If
# fewer than 2 specs are available, the script writes a warning and exits.
#
# This script is run AFTER all desired specs have been completed end-to-end
# (00 -> 07). It is NOT part of the per-spec orchestrator.
# =============================================================================

source("R/00_pipeline_config.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

cat(strrep("=", 78), "\n", sep = "")
cat("Cross-spec contrast — script 08\n")
cat(strrep("=", 78), "\n", sep = "")

# ---- Specs to attempt --------------------------------------------------------
specs <- list(
  list(label = "Canonical (SNV-only, PE/PPE retained)",
       short = "canonical",
       snv_only = TRUE,  drop_ppe_cal = FALSE, drop_ppe_apply = FALSE),
  list(label = "Sensitivity: with indels, PE/PPE retained",
       short = "indels_ppe_retained",
       snv_only = FALSE, drop_ppe_cal = FALSE, drop_ppe_apply = FALSE),
  list(label = "Sensitivity: with indels, PE/PPE excluded",
       short = "indels_ppe_excluded",
       snv_only = FALSE, drop_ppe_cal = TRUE,  drop_ppe_apply = TRUE)
)

# Add tag and resolved paths to each spec entry
for (i in seq_along(specs)) {
  s <- specs[[i]]
  paths <- build_paths_for_spec(s$snv_only, s$drop_ppe_cal, s$drop_ppe_apply)
  specs[[i]]$paths     <- paths
  specs[[i]]$apply_tag <- paths$apply_tag
  specs[[i]]$cal_tag   <- paths$cal_tag
}

# ---- Determine which specs have complete outputs -----------------------------
has_complete_outputs <- function(spec) {
  required <- c(
    file.path(spec$paths$calibration, "cal_lexicographic_thresholds.csv"),
    file.path(spec$paths$external,    "external_validity_grid.csv"),
    file.path(spec$paths$external,    "external_validity_main.csv")
  )
  all(file.exists(required))
}

cat("\nSpecs detected:\n")
status_rows <- character(0)
for (s in specs) {
  status <- if (has_complete_outputs(s)) "AVAILABLE" else "missing"
  msg <- sprintf("  %-46s [%s]", s$label, status)
  cat(msg, "\n")
  status_rows <- c(status_rows, msg)
}

available_idx <- vapply(specs, has_complete_outputs, logical(1))
specs_avail <- specs[available_idx]

if (length(specs_avail) < 2) {
  cat("\n[08_cross_spec_contrast] Fewer than 2 specs available; skipping contrast.\n")
  cat("To produce contrasts, run the pipeline under each desired spec configuration.\n")
  cat("(Edit R/00_pipeline_config.R flags, then re-run run_pipeline.R.)\n")
  quit(save = "no", status = 0)
}

cat(sprintf("\nProceeding with %d available spec(s).\n", length(specs_avail)))

# ---- Output directories ------------------------------------------------------
dir.create(PATHS$contrast_tables,  recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$contrast_figures, recursive = TRUE, showWarnings = FALSE)

# ---- Helpers -----------------------------------------------------------------
write_table_md <- function(df, slug, caption) {
  csv_path <- file.path(PATHS$contrast_tables, paste0(slug, ".csv"))
  md_path  <- file.path(PATHS$contrast_tables, paste0(slug, ".md"))

  readr::write_csv(df, csv_path)

  md_lines <- c(paste0("**", caption, "**"), "",
                knitr::kable(df, format = "pipe"))
  writeLines(md_lines, md_path)

  message(sprintf("  wrote %s (%d rows)", slug, nrow(df)))
}

fmt_est_ci <- function(est, lo, hi, digits = 3) {
  ifelse(is.na(est), "—",
         sprintf(sprintf("%%.%df (%%.%df, %%.%df)", digits, digits, digits),
                 est, lo, hi))
}

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

# Color palette consistent across all contrast figures (Okabe-Ito subset)
spec_colors <- c(
  "canonical"           = "#0072B2",
  "indels_ppe_retained" = "#E69F00",
  "indels_ppe_excluded" = "#999999"
)

spec_shapes <- c(
  "canonical"           = 16,
  "indels_ppe_retained" = 17,
  "indels_ppe_excluded" = 15
)

# ---- Read all available spec outputs into a single list -------------------
read_spec_outputs <- function(spec) {
  thr_path  <- file.path(spec$paths$calibration, "cal_lexicographic_thresholds.csv")
  grid_path <- file.path(spec$paths$external,    "external_validity_grid.csv")
  main_path <- file.path(spec$paths$external,    "external_validity_main.csv")
  boot_path <- file.path(spec$paths$external,    "bootstrap_summary.csv")

  thr  <- readr::read_csv(thr_path,  show_col_types = FALSE)
  grid <- readr::read_csv(grid_path, show_col_types = FALSE)
  main <- readr::read_csv(main_path, show_col_types = FALSE)
  boot <- if (file.exists(boot_path)) {
    readr::read_csv(boot_path, show_col_types = FALSE)
  } else {
    tibble(cell = character(0))
  }

  list(
    spec_short = spec$short,
    spec_label = spec$label,
    apply_tag  = spec$apply_tag,
    thresholds = thr,
    grid       = grid,
    main       = main,
    bootstrap  = boot
  )
}

all_outputs <- lapply(specs_avail, read_spec_outputs)

# Helper: spec_label factor with levels in canonical -> sensitivities order
spec_label_levels <- vapply(specs, function(s) s$label, character(1))

# =============================================================================
# Table C1 — cross-spec calibrated thresholds
# =============================================================================
make_tableC1 <- function() {
  message("\n[Table C1] Cross-spec calibrated thresholds")

  rows <- map_dfr(all_outputs, function(out) {
    out$thresholds %>%
      filter(selection == "primary") %>%
      transmute(
        Spec        = out$spec_label,
        `Apply tag` = out$apply_tag,
        `DP min`    = DP_min,
        `AD1 min`   = AD1_min,
        `MAF min`   = MAF_min,
        `MAF max`   = MAF_MAX,
        `n eval pairs` = pairs_evaluable,
        `prop both0`   = round(prop_both0, 3),
        `median Jaccard at primary` = round(jacc_med_eval, 3),
        `mean confirm rate`         = round(conf_med_eval, 3)
      )
  }) %>%
    mutate(Spec = factor(Spec, levels = spec_label_levels)) %>%
    arrange(Spec)

  write_table_md(rows, "tableC1_cross_spec_thresholds",
                 paste("Table C1. Calibrated decision rules at the primary tier across analytic specifications.",
                       "Each row is one spec; columns describe the calibrated DP/AD1/MAF thresholds",
                       "and within-pair concordance metrics at the lexicographic primary."))
}

# =============================================================================
# Table C2 — cross-spec external validity at primary
# =============================================================================
make_tableC2 <- function() {
  message("\n[Table C2] Cross-spec external validity at primary")

  rows <- map_dfr(all_outputs, function(out) {
    main_p <- out$main %>% filter(threshold_tier == "primary")

    main_p <- main_p %>%
      mutate(cell = case_when(
        cohort == "Florida" & cluster_def == "fl_surveillance" ~ "FL_surveillance_primary",
        cohort == "Ghana"   & cluster_def == "gh_12_snp_only"  ~ "GH_12_snp_only_primary",
        TRUE                                                    ~ NA_character_
      )) %>%
      left_join(out$bootstrap, by = "cell")

    main_p %>%
      transmute(
        Spec        = out$spec_label,
        Cohort      = cohort,
        `Cluster definition` = case_when(
          cluster_def == "fl_surveillance" ~ "FL surveillance",
          cluster_def == "gh_12_snp_only"  ~ "<=12 SNP-only (Ghana)"
        ),
        `n within`  = n_within,
        `n between` = n_between,
        `Within-cluster prev.` = round(prevalence, 3),
        Sens        = round(sens, 3),
        Spec_       = round(spec, 3),
        PPV         = round(ppv, 3),
        Enrichment  = round(enrichment, 2),
        `PR (analytic 95% CI)`     = fmt_est_ci(pr, pr_lo, pr_hi, 3),
        `PR (bootstrap 95% CI)`    = ifelse(is.na(pr_boot_lo), "—",
                                            sprintf("(%.3f, %.3f)",
                                                    pr_boot_lo, pr_boot_hi)),
        `PPV (bootstrap 95% CI)`   = ifelse(is.na(ppv_boot_lo), "—",
                                            sprintf("(%.3f, %.3f)",
                                                    ppv_boot_lo, ppv_boot_hi)),
        `Enrichment (bootstrap 95% CI)` = ifelse(is.na(enr_boot_lo), "—",
                                                  sprintf("(%.2f, %.2f)",
                                                          enr_boot_lo, enr_boot_hi))
      )
  }) %>%
    rename(Specificity = Spec_) %>%
    mutate(Spec = factor(Spec, levels = spec_label_levels)) %>%
    arrange(Cohort, Spec)

  write_table_md(rows, "tableC2_cross_spec_external_validity_primary",
                 paste("Table C2. External validity at primary thresholds across analytic specifications.",
                       "Each cohort appears in one row per available spec.",
                       "Bootstrap 95% CIs use the percentile method on cluster-resampled sample IDs.",
                       "Enrichment = PPV / within-cluster prevalence; values >1 indicate the calibrated test enriches above baseline for true within-cluster pairs (transmission signal recovery); values ~1 indicate negative-control behavior."))
}

# =============================================================================
# Table C3 — cross-spec threshold ladder
# =============================================================================
make_tableC3 <- function() {
  message("\n[Table C3] Cross-spec threshold ladder")

  rows <- map_dfr(all_outputs, function(out) {
    out$main %>%
      transmute(
        Spec      = out$spec_label,
        Cohort    = cohort,
        `Cluster definition` = case_when(
          cluster_def == "fl_surveillance" ~ "FL surveillance",
          cluster_def == "gh_12_snp_only"  ~ "<=12 SNP-only (Ghana)"
        ),
        Threshold = factor(threshold_tier, levels = c("loose","primary","tighter")),
        `n within`  = n_within,
        `n between` = n_between,
        Sens      = round(sens, 3),
        Spec_     = round(spec, 3),
        PPV       = round(ppv, 3),
        Enrichment = round(enrichment, 2),
        `PR (analytic 95% CI)` = fmt_est_ci(pr, pr_lo, pr_hi, 3)
      )
  }) %>%
    rename(Specificity = Spec_) %>%
    mutate(Spec = factor(Spec, levels = spec_label_levels)) %>%
    arrange(Cohort, Spec, Threshold)

  write_table_md(rows, "tableC3_cross_spec_threshold_ladder",
                 paste("Table C3. Threshold-ladder prevalence ratio (PR) across analytic specifications.",
                       "Under classical measurement-error theory with non-differential misclassification,",
                       "PR is expected to increase from loose to tighter in cohorts with genuine transmission signal."))
}

# =============================================================================
# Figure C1 — cross-spec PR forest at canonical primary
# =============================================================================
make_figureC1 <- function() {
  message("\n[Figure C1] Cross-spec PR forest")

  pdat <- map_dfr(all_outputs, function(out) {
    main_p <- out$main %>% filter(threshold_tier == "primary")
    main_p %>%
      mutate(spec_short = out$spec_short,
             spec_label = out$spec_label)
  }) %>%
    mutate(
      cohort_label = case_when(
        cohort == "Florida"                              ~ "Florida (surveillance)",
        cohort == "Ghana" & cluster_def == "gh_12_snp_only" ~ "Ghana (<=12 SNP-only)"
      ),
      spec_label = factor(spec_label, levels = spec_label_levels)
    )

  p <- ggplot(pdat,
              aes(x = pr, y = spec_label, color = spec_short)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
    geom_pointrange(aes(xmin = pr_lo, xmax = pr_hi),
                    size = 0.5, fatten = 4) +
    facet_wrap(~ cohort_label, ncol = 1, scales = "free_y") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.25, 0.5, 1.0, 2.0, 4.0)) +
    scale_color_manual(values = spec_colors, guide = "none") +
    labs(
      x = "Prevalence ratio (PR), log2 scale",
      y = NULL,
      title = "Cross-spec PR at primary thresholds",
      caption = paste(
        "Analytic 95% CIs from modified Poisson regression with dyadic-clustered SE.",
        "Reference line: PR=1.",
        sep = "\n"
      )
    ) +
    theme_paper()

  ggsave(file.path(PATHS$contrast_figures, "figureC1_cross_spec_pr_forest.pdf"),
         p, width = 7, height = 5)
  ggsave(file.path(PATHS$contrast_figures, "figureC1_cross_spec_pr_forest.png"),
         p, width = 7, height = 5, dpi = 300)
  message("  wrote figureC1_cross_spec_pr_forest")
}

# =============================================================================
# Figure C2 — cross-spec threshold-ladder PR
# =============================================================================
make_figureC2 <- function() {
  message("\n[Figure C2] Cross-spec threshold ladder")

  pdat <- map_dfr(all_outputs, function(out) {
    out$main %>%
      mutate(spec_short = out$spec_short,
             spec_label = out$spec_label)
  }) %>%
    mutate(
      threshold_tier = factor(threshold_tier,
                              levels = c("loose","primary","tighter")),
      cohort_label = case_when(
        cohort == "Florida"                              ~ "Florida (surveillance)",
        cohort == "Ghana" & cluster_def == "gh_12_snp_only" ~ "Ghana (<=12 SNP-only)"
      ),
      spec_label = factor(spec_label, levels = spec_label_levels)
    )

  p <- ggplot(pdat,
              aes(x = pr, y = threshold_tier,
                  color = spec_short, shape = spec_short)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
    geom_pointrange(aes(xmin = pr_lo, xmax = pr_hi),
                    position = position_dodge(width = 0.6),
                    size = 0.4, fatten = 3) +
    facet_wrap(~ cohort_label, ncol = 1, scales = "free_y") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.25, 0.5, 1.0, 2.0, 4.0)) +
    scale_color_manual(values = spec_colors, name = "Specification",
                       labels = function(x) {
                         m <- setNames(spec_label_levels,
                                       vapply(specs, function(s) s$short, character(1)))
                         unname(m[x])
                       }) +
    scale_shape_manual(values = spec_shapes, name = "Specification",
                       labels = function(x) {
                         m <- setNames(spec_label_levels,
                                       vapply(specs, function(s) s$short, character(1)))
                         unname(m[x])
                       }) +
    labs(
      x = "Prevalence ratio (PR), log2 scale",
      y = "Threshold tier",
      title = "Threshold-ladder PR across analytic specifications",
      caption = paste(
        "Framework prediction: monotonic PR increase from loose to tighter in cohorts with genuine transmission signal.",
        "Reference line: PR=1.",
        sep = "\n"
      )
    ) +
    theme_paper(base_size = 9) +
    theme(legend.position = "bottom")

  ggsave(file.path(PATHS$contrast_figures, "figureC2_cross_spec_threshold_ladder.pdf"),
         p, width = 8, height = 6)
  ggsave(file.path(PATHS$contrast_figures, "figureC2_cross_spec_threshold_ladder.png"),
         p, width = 8, height = 6, dpi = 300)
  message("  wrote figureC2_cross_spec_threshold_ladder")
}

# =============================================================================
# Figure C3 — cross-spec PPV enrichment forest at primary
# =============================================================================
make_figureC3 <- function() {
  message("\n[Figure C3] Cross-spec PPV enrichment forest")

  pdat <- map_dfr(all_outputs, function(out) {
    main_p <- out$main %>% filter(threshold_tier == "primary") %>%
      mutate(cell = case_when(
        cohort == "Florida" & cluster_def == "fl_surveillance" ~ "FL_surveillance_primary",
        cohort == "Ghana"   & cluster_def == "gh_12_snp_only"  ~ "GH_12_snp_only_primary",
        TRUE                                                    ~ NA_character_
      )) %>%
      left_join(out$bootstrap, by = "cell")

    main_p %>%
      mutate(spec_short = out$spec_short,
             spec_label = out$spec_label)
  }) %>%
    mutate(
      cohort_label = case_when(
        cohort == "Florida"                              ~ "Florida (surveillance)",
        cohort == "Ghana" & cluster_def == "gh_12_snp_only" ~ "Ghana (<=12 SNP-only)"
      ),
      spec_label = factor(spec_label, levels = spec_label_levels)
    )

  p <- ggplot(pdat,
              aes(x = enrichment, y = spec_label, color = spec_short)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
    geom_pointrange(aes(xmin = enr_boot_lo, xmax = enr_boot_hi),
                    size = 0.5, fatten = 4) +
    facet_wrap(~ cohort_label, ncol = 1, scales = "free_y") +
    scale_x_continuous(trans = "log2",
                       breaks = c(0.5, 0.75, 1.0, 1.5, 2.0, 3.0)) +
    scale_color_manual(values = spec_colors, guide = "none") +
    labs(
      x = "PPV enrichment (PPV / within-cluster prevalence), log2 scale",
      y = NULL,
      title = "Cross-spec PPV enrichment at primary thresholds",
      caption = paste(
        "Enrichment >1: calibrated test enriches above baseline within-cluster prevalence (transmission signal recovery).",
        "Enrichment ~1: test fires no more informatively than chance given local cluster prevalence (negative-control).",
        "Bootstrap 95% CIs (percentile method).",
        sep = "\n"
      )
    ) +
    theme_paper()

  ggsave(file.path(PATHS$contrast_figures, "figureC3_cross_spec_ppv_enrichment.pdf"),
         p, width = 7, height = 5)
  ggsave(file.path(PATHS$contrast_figures, "figureC3_cross_spec_ppv_enrichment.png"),
         p, width = 7, height = 5, dpi = 300)
  message("  wrote figureC3_cross_spec_ppv_enrichment")
}

# =============================================================================
# Run all
# =============================================================================
make_tableC1()
make_tableC2()
make_tableC3()

make_figureC1()
make_figureC2()
make_figureC3()

cat(sprintf("\n[08_cross_spec_contrast] Done.\n"))
cat(sprintf("  Tables:  %s/\n", PATHS$contrast_tables))
cat(sprintf("  Figures: %s/\n", PATHS$contrast_figures))
