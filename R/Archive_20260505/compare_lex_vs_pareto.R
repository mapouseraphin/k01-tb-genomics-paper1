# =============================================================================
# compare_lex_vs_pareto.R
#
# Cross-method calibration comparison: lex vs hybrid Pareto.
#
# Reads from BOTH per-tag calibration trees populated by:
#   - canonical lex pipeline: 01 -> 02_calibrate_lexicographic.R -> 02b
#   - Pareto switch:          01 -> 02_calibrate_hybrid_pareto.R -> 02b
#
# Produces:
#   compareA_selected_cells.csv         - all tiers x both methods x stages
#   compareB_concordance_at_primary.csv - pair-level metrics at each primary
#   compareC_bootstrap_stability.csv    - B, joint + marginal stability
#   compare_report.md                   - narrative + decision checklist
#
# Run AFTER both calibration pipelines have completed (each through 02b).
# Comparison is calibration-layer only; LCA is run on the chosen primary
# downstream of this decision.
#
# Usage:
#   source("R/compare_lex_vs_pareto.R")
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr)
})

# ---- Resolve per-method paths ----------------------------------------------
# Use the active spec flags (SNV_ONLY, DROP_PPE_*) so this works for canonical
# AND any sensitivity spec (e.g., Sens A) provided both methods were run there.
paths_lex <- build_paths_for_spec(
  snv_only         = SNV_ONLY,
  drop_ppe_cal     = DROP_PPE_CAL,
  drop_ppe_apply   = DROP_PPE_APPLY,
  selection_method = "lex"
)
paths_pareto <- build_paths_for_spec(
  snv_only         = SNV_ONLY,
  drop_ppe_cal     = DROP_PPE_CAL,
  drop_ppe_apply   = DROP_PPE_APPLY,
  selection_method = "pareto"
)

OUT <- file.path("outputs", "cross_method_lex_vs_pareto",
                 sprintf("%s_vs_%s", paths_lex$cal_tag, paths_pareto$cal_tag))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

message("[compare_lex_vs_pareto] Comparing:")
message(sprintf("  lex tree    : %s", paths_lex$calibration))
message(sprintf("  pareto tree : %s", paths_pareto$calibration))
message(sprintf("  Output      : %s", OUT))

# ---- Helpers ----------------------------------------------------------------
safe_rds <- function(path) if (file.exists(path)) readRDS(path) else NULL
safe_csv <- function(path) if (file.exists(path))
  read_csv(path, show_col_types = FALSE) else NULL

load_method <- function(p, label) {
  out <- list(label = label, paths = p)
  out$thr_obs       <- safe_rds(file.path(p$calibration, "thr_primary_observed.rds"))
  out$thr_boot      <- safe_rds(file.path(p$calibration, "thr_primary_bootstrap.rds"))
  out$thresholds    <- safe_rds(file.path(p$calibration, "thresholds.rds"))
  grid_name <- if (label == "lex") "cal_lexicographic_grid_summary.csv"
               else                "cal_hybrid_pareto_grid_summary.csv"
  out$grid_summary  <- safe_csv(file.path(p$calibration, grid_name))
  out$threshold_csv <- safe_csv(file.path(p$calibration, "cal_lexicographic_thresholds.csv"))
  out$boot_stab     <- safe_rds(file.path(p$calibration, "bootstrap",
                                          "bootstrap_stability.rds"))
  out$boot_summary  <- safe_csv(file.path(p$calibration, "bootstrap",
                                          "bootstrap_summary_table.csv"))
  if (label == "pareto") {
    out$frontier <- safe_csv(file.path(p$calibration, "cal_hybrid_pareto_frontier.csv"))
  }

  missing <- character()
  for (k in c("thr_obs", "thr_boot", "thresholds", "grid_summary", "boot_stab")) {
    if (is.null(out[[k]])) missing <- c(missing, k)
  }
  if (length(missing)) {
    message(sprintf("  %s: missing artifacts -> %s",
                    label, paste(missing, collapse = ", ")))
  }
  out
}

L <- load_method(paths_lex,    "lex")
P <- load_method(paths_pareto, "pareto")

# Hard fail if either method incomplete on the core artifacts
incomplete <- c(
  if (is.null(L$thr_obs))     "lex thr_primary_observed.rds",
  if (is.null(L$thr_boot))    "lex thr_primary_bootstrap.rds",
  if (is.null(L$grid_summary))"lex grid summary",
  if (is.null(L$boot_stab))   "lex bootstrap_stability.rds",
  if (is.null(P$thr_obs))     "pareto thr_primary_observed.rds",
  if (is.null(P$thr_boot))    "pareto thr_primary_bootstrap.rds",
  if (is.null(P$grid_summary))"pareto grid summary",
  if (is.null(P$boot_stab))   "pareto bootstrap_stability.rds"
)
if (length(incomplete)) {
  stop("Comparison cannot proceed; missing artifacts:\n  - ",
       paste(incomplete, collapse = "\n  - "),
       "\nRun both calibration pipelines (canonical lex + lex_pareto_switch.R) ",
       "end-to-end through 02b before comparing.")
}

# ---- Cell tibble helpers ----------------------------------------------------
cell_row <- function(cell, source_lbl, method_lbl, tier) {
  if (is.null(cell)) return(NULL)
  tibble(method = method_lbl, tier = tier, source = source_lbl,
         DP_min  = cell$DP_min,
         AD1_min = cell$AD1_min,
         MAF_min = cell$MAF_min,
         MAF_max = cell$MAF_max %||% NA_real_)
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

cells_agree <- function(a, b) {
  !is.null(a) && !is.null(b) &&
    isTRUE(a$DP_min  == b$DP_min) &&
    isTRUE(a$AD1_min == b$AD1_min) &&
    isTRUE(a$MAF_min == b$MAF_min)
}

fmt_cell <- function(cell) {
  if (is.null(cell)) return("(missing)")
  sprintf("DP=%g, AD1=%g, MAF=%g", cell$DP_min, cell$AD1_min, cell$MAF_min)
}

# ---- TABLE A: selected cells (all tiers, both methods) ---------------------
selected_cells <- bind_rows(
  cell_row(L$thresholds$looser,  "Stage 1 (top-K min stringency)",   "lex",    "looser"),
  cell_row(L$thr_obs,            "Stage 1 (lex primary)",             "lex",    "primary"),
  cell_row(L$thr_boot,           "Stage 2 (bootstrap modal)",         "lex",    "primary"),
  cell_row(L$thresholds$tighter, "Stage 1 (top-K max stringency)",    "lex",    "tighter"),
  cell_row(P$thresholds$looser,  "Stage 1 (frontier min stringency)", "pareto", "looser"),
  cell_row(P$thr_obs,            "Stage 1 (closest-to-ideal)",        "pareto", "primary"),
  cell_row(P$thr_boot,           "Stage 2 (bootstrap modal)",         "pareto", "primary"),
  cell_row(P$thresholds$tighter, "Stage 1 (frontier max stringency)", "pareto", "tighter")
)

# Mark Stage 1 vs Stage 2 within primary explicitly
selected_cells <- selected_cells %>%
  mutate(stage = case_when(
    grepl("Stage 1", source) ~ "Stage 1",
    grepl("Stage 2", source) ~ "Stage 2",
    TRUE ~ NA_character_
  )) %>%
  select(method, tier, stage, source, DP_min, AD1_min, MAF_min, MAF_max)

write_csv(selected_cells, file.path(OUT, "compareA_selected_cells.csv"))

# ---- TABLE B: concordance metrics at each primary cell ---------------------
lookup_metrics <- function(grid, cell, method_lbl, tier_lbl) {
  if (is.null(grid) || is.null(cell)) return(NULL)
  grid %>%
    filter(DP_min  == cell$DP_min,
           AD1_min == cell$AD1_min,
           MAF_min == cell$MAF_min) %>%
    select(any_of(c("DP_min", "AD1_min", "MAF_min",
                    "pairs_evaluable", "prop_both0", "prop_overlap",
                    "jacc_med_eval", "overlap_med_eval",
                    "conf_med_eval", "r2_med_eval"))) %>%
    mutate(method = method_lbl, tier = tier_lbl) %>%
    select(method, tier, everything())
}

concordance_tbl <- bind_rows(
  lookup_metrics(L$grid_summary, L$thr_obs,  "lex",    "primary (Stage 1)"),
  lookup_metrics(L$grid_summary, L$thr_boot, "lex",    "primary (Stage 2)"),
  lookup_metrics(P$grid_summary, P$thr_obs,  "pareto", "primary (Stage 1)"),
  lookup_metrics(P$grid_summary, P$thr_boot, "pareto", "primary (Stage 2)")
)

# Round for human readability while preserving CSV float
write_csv(concordance_tbl, file.path(OUT, "compareB_concordance_at_primary.csv"))

# ---- TABLE C: bootstrap stability side by side -----------------------------
stab_row <- function(stab) {
  if (is.null(stab)) return(rep(NA_character_, 12))
  c(
    stab$selection_method,
    sprintf("%d", stab$n_reps),
    sprintf("%d (%.1f%%)", stab$n_failed, 100 * stab$n_failed / stab$n_reps),
    sprintf("%d", stab$n_ok),
    sprintf("(DP=%g, AD1=%g, MAF=%g)",
            stab$observed_primary$DP_min,
            stab$observed_primary$AD1_min,
            stab$observed_primary$MAF_min),
    sprintf("%.1f%%", 100 * stab$prop_observed_primary),
    sprintf("(DP=%g, AD1=%g, MAF=%g)",
            stab$modal_cell$DP_min,
            stab$modal_cell$AD1_min,
            stab$modal_cell$MAF_min),
    sprintf("%.1f%%", 100 * stab$prop_modal_joint),
    if (stab$modal_eq_observed) "yes" else "NO",
    sprintf("%.1f%%", 100 * stab$marginal_stability$prop_at_modal[
      stab$marginal_stability$axis == "DP_min"]),
    sprintf("%.1f%%", 100 * stab$marginal_stability$prop_at_modal[
      stab$marginal_stability$axis == "AD1_min"]),
    sprintf("%.1f%%", 100 * stab$marginal_stability$prop_at_modal[
      stab$marginal_stability$axis == "MAF_min"])
  )
}

stab_tbl <- tibble(
  metric = c(
    "Selection method",
    "Bootstrap reps (B)",
    "Failed reps",
    "Reps yielding primary (n_ok)",
    "Stage 1 observed primary",
    "P(rep == observed primary)",
    "Stage 2 modal cell",
    "P(modal cell)",
    "Modal == observed?",
    "Marginal P(DP_min = modal)",
    "Marginal P(AD1_min = modal)",
    "Marginal P(MAF_min = modal)"
  ),
  lex    = stab_row(L$boot_stab),
  pareto = stab_row(P$boot_stab)
)
write_csv(stab_tbl, file.path(OUT, "compareC_bootstrap_stability.csv"))

# ---- Cross-method agreement diagnostics ------------------------------------
stage1_agree <- cells_agree(L$thr_obs,  P$thr_obs)
stage2_agree <- cells_agree(L$thr_boot, P$thr_boot)

# Concordance deltas at the active (Stage 2) primary
delta_at_primary <- function() {
  lex_row <- concordance_tbl %>% filter(method == "lex",    tier == "primary (Stage 2)")
  par_row <- concordance_tbl %>% filter(method == "pareto", tier == "primary (Stage 2)")
  if (nrow(lex_row) == 0 || nrow(par_row) == 0) return(NULL)
  tibble(
    metric = c("pairs_evaluable", "prop_both0", "jacc_med_eval",
               "conf_med_eval", "r2_med_eval"),
    lex    = c(lex_row$pairs_evaluable, lex_row$prop_both0,
               lex_row$jacc_med_eval,   lex_row$conf_med_eval,
               lex_row$r2_med_eval),
    pareto = c(par_row$pairs_evaluable, par_row$prop_both0,
               par_row$jacc_med_eval,   par_row$conf_med_eval,
               par_row$r2_med_eval)
  ) %>% mutate(delta_lex_minus_pareto = lex - pareto)
}
delta_tbl <- delta_at_primary()

# ---- Markdown decision report ---------------------------------------------
fmt_pct <- function(x) if (is.na(x)) "NA" else sprintf("%.1f%%", 100 * x)

md <- c(
  "# Lex vs Pareto calibration comparison",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Spec: lex tag = `%s`, pareto tag = `%s`",
          paths_lex$cal_tag, paths_pareto$cal_tag),
  sprintf("Output dir: `%s`", OUT),
  "",
  "## 1. Headline: do the methods agree?",
  "",
  sprintf("- **Stage 1 cells agree:** %s", if (stage1_agree) "**YES**" else "**NO**"),
  sprintf("- **Stage 2 cells agree:** %s", if (stage2_agree) "**YES**" else "**NO**"),
  "",
  if (stage2_agree) {
    "When Stage 2 cells agree, the selection method is operationally moot for downstream LCA. Pick lex as canonical for naming consistency; report Pareto in the supplement as method-equivalence evidence."
  } else {
    "Stage 2 cells differ. Decision criteria below resolve which becomes canonical."
  },
  "",
  "## 2. Selected cells",
  "",
  sprintf("- **lex Stage 1:**     %s", fmt_cell(L$thr_obs)),
  sprintf("- **lex Stage 2:**     %s", fmt_cell(L$thr_boot)),
  sprintf("- **Pareto Stage 1:**  %s", fmt_cell(P$thr_obs)),
  sprintf("- **Pareto Stage 2:**  %s", fmt_cell(P$thr_boot)),
  "",
  "## 3. Concordance at each primary cell",
  "",
  "Pair-level metrics computed on the calibration sample at each cell.",
  "Higher Jaccard / Confirm and lower prop_both0 are preferred.",
  "",
  "```",
  paste(capture.output(print(as.data.frame(concordance_tbl), row.names = FALSE)),
        collapse = "\n"),
  "```",
  "",
  if (!is.null(delta_tbl)) {
    c("### Stage 2 head-to-head delta (lex - pareto)",
      "",
      "Positive = lex better; negative = Pareto better.",
      "",
      "```",
      paste(capture.output(print(as.data.frame(delta_tbl), row.names = FALSE)),
            collapse = "\n"),
      "```",
      "")
  } else "",
  "## 4. Bootstrap stability (Stage 2)",
  "",
  "```",
  paste(capture.output(print(as.data.frame(stab_tbl), row.names = FALSE)),
        collapse = "\n"),
  "```",
  "",
  "## 5. Decision criteria (in priority order)",
  "",
  "1. **Joint stability at the active primary** (P(modal cell)).",
  "   This is the bootstrap-stabilized Stage 2 selection probability.",
  "   Higher = the primary cell is reproducibly selected under resampling.",
  "   Use as primary tiebreaker.",
  "",
  "2. **Concordance at the active primary** (jacc_med_eval, conf_med_eval).",
  "   Higher = the rule recovers consistent calls across replicates.",
  "   If methods are within ~0.02 on each, treat as a tie on this axis.",
  "",
  "3. **Marginal axis stability**. A method that stabilizes 2 of 3 axes",
  "   above ~70% is more defensible than one that stabilizes only 1.",
  "   Look for which axis (DP, AD1, MAF) is unstable under each method.",
  "",
  "4. **Operational properties**. Lower DP / AD1 cells are more inclusive",
  "   (more calls retained, more LCA evaluable cells); higher cells are",
  "   stricter (fewer calls, narrower posteriors but possibly thinner data).",
  "   This should match Paper 2's needs: an inclusive rule survives",
  "   transport to lower-coverage Florida data; a strict rule produces",
  "   tighter LCA but may fail on cohort transferability.",
  "",
  "5. **Cell agreement override**. If Stage 2 cells agree across methods,",
  "   skip 1-4 and pick lex (canonical).",
  "",
  "## 6. Suggested Methods text (insert after decision)",
  "",
  "> We compared two threshold-selection methods at Stage 1 (lexicographic",
  "> ranking on six concordance criteria vs. hybrid Pareto + closest-to-ideal",
  "> on the (Jaccard, confirmation rate) plane) and ran each through identical",
  sprintf("> Stage 2 bootstrap stabilization (B = %d). The two methods selected",
          L$boot_stab$n_reps),
  sprintf("> %s Stage 2 cell%s. We adopted [CHOSEN_METHOD] as canonical for",
          if (stage2_agree) "the same"           else "different",
          if (stage2_agree) ""                   else "s"),
  "> downstream LCA characterization based on [STABILITY/CONCORDANCE/AGREEMENT]",
  "> (Supp Table SX). The non-canonical method is reported as a parallel",
  "> sensitivity in Supp Table SY.",
  "",
  "## 7. Files in this comparison",
  "",
  "- `compareA_selected_cells.csv` - all tiers, both methods, both stages",
  "- `compareB_concordance_at_primary.csv` - pair-level metrics at each primary",
  "- `compareC_bootstrap_stability.csv` - B, joint and marginal stability",
  ""
)

writeLines(md, file.path(OUT, "compare_report.md"))

# ---- Console summary --------------------------------------------------------
message("")
message("============================================================")
message("[compare_lex_vs_pareto] Comparison complete.")
message("============================================================")
message(sprintf("  Stage 1 agree : %s", stage1_agree))
message(sprintf("  Stage 2 agree : %s", stage2_agree))
message("")
message(sprintf("  lex    Stage 1: %s", fmt_cell(L$thr_obs)))
message(sprintf("  lex    Stage 2: %s", fmt_cell(L$thr_boot)))
message(sprintf("  pareto Stage 1: %s", fmt_cell(P$thr_obs)))
message(sprintf("  pareto Stage 2: %s", fmt_cell(P$thr_boot)))
message("")
message(sprintf("  P(modal cell) lex    : %s",
                fmt_pct(L$boot_stab$prop_modal_joint)))
message(sprintf("  P(modal cell) pareto : %s",
                fmt_pct(P$boot_stab$prop_modal_joint)))
message("")
message(sprintf("  Outputs in: %s/", OUT))
message("    compareA_selected_cells.csv")
message("    compareB_concordance_at_primary.csv")
message("    compareC_bootstrap_stability.csv")
message("    compare_report.md")
message("")
