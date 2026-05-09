# =============================================================================
# 02_calibrate_hybrid_pareto.R
#
# Hybrid Pareto calibration -- parallel sensitivity Stage 1 selector (post
# 2026-05-03 lex revert; Pareto runs alongside lex, both lock independently).
#
# Drop-in alternative to 02_calibrate_lexicographic.R; writes the same output
# filenames (thr_looser.rds, thr_primary.rds, thr_tighter.rds, thresholds.rds,
# cal_lexicographic_thresholds.csv) so all downstream consumers work without
# modification.
#
# Three-tier evidentiary ladder (no sentinel; all real-grid cells):
#   looser  : Pareto-frontier point with MIN stringency among non-primary
#             frontier points. Less stringent than primary on >= 1 axis.
#             (Fallback: top-K min stringency from guard-passing cells if
#             frontier has only 1 point.)
#   primary : Pareto frontier point closest to (1, 1) on (Jaccard, Confirm)
#             tiebreak: highest r2_med_eval, then highest stringency
#   tighter : Pareto-frontier point with MAX stringency among non-primary
#             frontier points. (Fallback: top-K max stringency from
#             guard-passing cells if frontier has only 1 point.)
#
# The "loose" no-filter sentinel was retired 2026-05-03. Looser is now a
# real-grid frontier-resident tier with real concordance metrics.
#
# PARETO METHOD:
#   Objectives are the QUALITY METRICS (jacc_med_eval, conf_med_eval), where
#   conf_med_eval is the median symmetric confirmation rate (revised
#   2026-05-05 -- previously sourced from confirm_both = Jaccard, which
#   collapsed the (J, C) plane to a 1-D line and made closest-to-(1,1)
#   degenerate). conf_med_eval is now sourced from confirm_sym =
#   0.5 * (confirm_a_to_b + confirm_b_to_a), which is genuinely distinct
#   from Jaccard.
#
#   The (DP, AD1, MAF) grid is the search space; each grid cell is one
#   candidate point on the (J, C) plane. Cell A dominates cell B iff
#   A_J >= B_J and A_C >= B_C with at least one strict; the Pareto frontier
#   is the set of non-dominated cells.
#
# Pair-level filters when computing each median (same as lex):
#   - jacc_med_eval, conf_med_eval, overlap_med_eval: nA > 0 AND nB > 0
#   - r2_med_eval: nA > 0 AND nB > 0 AND n_inter >= 2
#
# Grid-level guards (canonical):
#   prop_both0      <= MAX_PROP_BOTH0  (canonical: 0.90)
#   pairs_evaluable >= MIN_EVAL_PAIRS  (canonical: 5)
#
# NOTE: 02b_calibration_bootstrap.R OVERWRITES the looser/tighter cells
#       written here, using axis-wise dominance against the stabilized
#       Primary. The frontier-resident looser/tighter values written here
#       are provisional.
#
# Inputs:  PATHS$calibration/cal_pairs.rds
# Outputs (PATHS$calibration):
#   cal_hybrid_pareto_grid_summary.csv
#   cal_lexicographic_thresholds.csv      (kept name for downstream compat)
#   cal_hybrid_pareto_frontier.csv        (frontier-specific output)
#   cal_lexicographic_heatmap.pdf         (kept name for downstream compat)
#   cal_hybrid_pareto_frontier.pdf        (frontier-specific output)
#   thr_looser.rds, thr_primary.rds, thr_tighter.rds
#   thr_primary_observed.rds  (Stage 1 backup; preserved through 02b)
#   thresholds.rds  (consolidated list with looser/primary/tighter)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(tidyr)
})

TOP_K          <- 10L
TOP_K_FALLBACK <- 30L

# ---- Load -------------------------------------------------------------------
cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))

req_cols <- c("DP_min", "AD1_min", "MAF_min",
              "sampleA", "sampleB",
              "nA", "nB", "n_union", "n_inter",
              "jaccard", "overlap", "confirm_both", "confirm_sym", "maf_cor")
miss <- setdiff(req_cols, names(cal_pairs))
if (length(miss)) stop("cal_pairs missing required columns: ",
                       paste(miss, collapse = ", "))

# ---- Aggregate per grid point ----------------------------------------------
cal_grid <- cal_pairs %>%
  group_by(DP_min, AD1_min, MAF_min) %>%
  summarise(
    pairs_total      = n(),
    pairs_both0      = sum(nA == 0 & nB == 0),
    pairs_any        = sum(n_union > 0),
    pairs_evaluable  = sum(nA > 0 & nB > 0),
    pairs_r2_elig    = sum(n_inter >= 2 & nA > 0 & nB > 0),
    prop_both0       = pairs_both0 / pairs_total,
    prop_overlap     = mean(jaccard > 0, na.rm = TRUE),
    jacc_med_eval    = median(jaccard     [n_union > 0 & nA > 0 & nB > 0], na.rm = TRUE),
    overlap_med_eval = median(overlap     [n_union > 0 & nA > 0 & nB > 0], na.rm = TRUE),
    conf_med_eval    = median(confirm_sym [n_union > 0 & nA > 0 & nB > 0], na.rm = TRUE),
    r2_med_eval      = median((maf_cor^2) [n_inter >= 2 & nA > 0 & nB > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    guard_prop_both0 = prop_both0 <= MAX_PROP_BOTH0,
    guard_n_eval     = pairs_evaluable >= MIN_EVAL_PAIRS,
    pass_guards      = guard_prop_both0 & guard_n_eval
  )

n_passing <- sum(cal_grid$pass_guards)
message(sprintf(
  "Grid-level guards: %d of %d grid points pass (prop_both0 <= %.2f AND n_evaluable >= %d).",
  n_passing, nrow(cal_grid), MAX_PROP_BOTH0, MIN_EVAL_PAIRS
))

if (n_passing == 0) {
  message("\n*** No grid points pass the canonical guards. ***")
  message("Diagnostic: top 10 grid points by pairs_evaluable:")
  print(cal_grid %>% arrange(desc(pairs_evaluable)) %>%
        select(DP_min, AD1_min, MAF_min,
               pairs_evaluable, prop_both0, jacc_med_eval) %>%
        head(10))
  stop("Calibration cannot proceed under canonical guards. Inspect diagnostic above.")
}

# ---- Pareto frontier on (jacc_med_eval, conf_med_eval) --------------------
is_pareto <- function(J, C) {
  ok <- !is.na(J) & !is.na(C)
  keep <- rep(FALSE, length(J))
  for (i in which(ok)) {
    dom <- (J >= J[i] & C >= C[i]) & (J > J[i] | C > C[i]) & ok
    if (!any(dom)) keep[i] <- TRUE
  }
  keep
}

cal_eligible <- cal_grid %>% filter(pass_guards)
cal_eligible$is_frontier <- is_pareto(cal_eligible$jacc_med_eval,
                                      cal_eligible$conf_med_eval)

frontier <- cal_eligible %>% filter(is_frontier)

message(sprintf("\nPareto frontier on (Jaccard, Confirm): %d of %d guard-passing grid points.",
                nrow(frontier), nrow(cal_eligible)))

# ---- Primary = closest-to-ideal on the frontier ----------------------------
add_stringency <- function(df) {
  df %>% mutate(S = min_rank(DP_min) + min_rank(AD1_min) + min_rank(MAF_min))
}

frontier <- frontier %>%
  mutate(dist_to_ideal = sqrt((1 - jacc_med_eval)^2 + (1 - conf_med_eval)^2)) %>%
  add_stringency()

cal_primary_hp <- frontier %>%
  arrange(dist_to_ideal,
          desc(r2_med_eval),
          desc(S),
          desc(DP_min), desc(AD1_min), desc(MAF_min)) %>%
  slice(1)

message("\nPrimary threshold (hybrid Pareto, closest-to-(1,1) on (J, C)):")
print(cal_primary_hp %>% select(DP_min, AD1_min, MAF_min,
                                 dist_to_ideal,
                                 jacc_med_eval, conf_med_eval, r2_med_eval,
                                 prop_both0, pairs_evaluable))

# ---- Looser and tighter via stringency on the frontier ---------------------
not_primary <- function(df) {
  df %>% filter(!(DP_min  == cal_primary_hp$DP_min &
                  AD1_min == cal_primary_hp$AD1_min &
                  MAF_min == cal_primary_hp$MAF_min))
}

# Candidate pool: prefer the frontier (>=2 points). Fallback to top-K
# guard-passing cells under lex-style ranking if the frontier has only one
# point (then looser/tighter aren't available from the frontier).
if (nrow(frontier) >= 2) {
  cand_pool  <- frontier
  pool_label <- sprintf("Pareto frontier (%d points)", nrow(frontier))
} else {
  cal_ranked <- cal_eligible %>%
    arrange(prop_both0,
            desc(prop_overlap),
            desc(jacc_med_eval),
            desc(overlap_med_eval),
            desc(conf_med_eval),
            desc(r2_med_eval))
  cand_pool  <- cal_ranked %>%
    slice_head(n = min(TOP_K_FALLBACK, nrow(cal_ranked))) %>%
    add_stringency()
  pool_label <- sprintf("frontier degenerate (1 point); fallback top-K=%d guard-passers",
                        TOP_K_FALLBACK)
}

pick_tighter <- function(cands) {
  cands %>% not_primary() %>%
    arrange(desc(S), desc(DP_min), desc(AD1_min), desc(MAF_min)) %>%
    slice_head(n = 1)
}
pick_looser <- function(cands) {
  cands %>% not_primary() %>%
    arrange(S, DP_min, AD1_min, MAF_min) %>%
    slice_head(n = 1)
}

cal_tighter_hp <- pick_tighter(cand_pool)
if (nrow(cal_tighter_hp) == 0) {
  cal_tighter_hp <- NULL
  message("\nNOTE: no tighter candidate available; pool exhausted (", pool_label, ").")
}

cal_looser_hp <- pick_looser(cand_pool)
if (nrow(cal_looser_hp) == 0) {
  cal_looser_hp <- NULL
  message("\nNOTE: no looser candidate available; pool exhausted (", pool_label, ").")
}

message(sprintf("\nLooser/tighter selected from candidate pool: %s", pool_label))

# ---- Assemble selected thresholds ------------------------------------------
common_cols <- c("selection", "DP_min", "AD1_min", "MAF_min",
                 "prop_both0", "prop_overlap",
                 "jacc_med_eval", "conf_med_eval", "r2_med_eval",
                 "pairs_evaluable", "pass_guards")

selected <- bind_rows(
  if (!is.null(cal_looser_hp))
    cal_looser_hp %>% select(-any_of(c("is_frontier", "dist_to_ideal", "S"))) %>%
      mutate(selection = "looser") %>% select(any_of(common_cols)),
  cal_primary_hp %>% select(-any_of(c("is_frontier", "dist_to_ideal", "S"))) %>%
    mutate(selection = "primary") %>% select(any_of(common_cols)),
  if (!is.null(cal_tighter_hp))
    cal_tighter_hp %>% select(-any_of(c("is_frontier", "dist_to_ideal", "S"))) %>%
      mutate(selection = "tighter") %>% select(any_of(common_cols))
)

message("\nSelected thresholds:")
print(selected)

readr::write_csv(cal_grid,
                 file.path(PATHS$calibration, "cal_hybrid_pareto_grid_summary.csv"))
# Note: filename is cal_lexicographic_thresholds.csv (not cal_hybrid_pareto_*)
# so downstream code that reads selected thresholds works unchanged. The
# selection method is recorded by the per-tag directory name (TAG$cal_tag
# ends in "_pareto") and by the cal_hybrid_pareto_frontier.csv companion file.
readr::write_csv(selected,
                 file.path(PATHS$calibration, "cal_lexicographic_thresholds.csv"))
readr::write_csv(frontier %>% select(DP_min, AD1_min, MAF_min,
                                      jacc_med_eval, conf_med_eval, r2_med_eval,
                                      dist_to_ideal),
                 file.path(PATHS$calibration, "cal_hybrid_pareto_frontier.csv"))

# ---- Visit-stratified evaluable pair diagnostic ----------------------------
prim <- cal_primary_hp
evaluable_at_primary <- cal_pairs %>%
  filter(DP_min  == prim$DP_min,
         AD1_min == prim$AD1_min,
         MAF_min == prim$MAF_min,
         nA > 0 & nB > 0) %>%
  mutate(
    visit_a = parse_visit(sampleA),
    visit_b = parse_visit(sampleB),
    slot_a  = parse_slot(sampleA),
    slot_b  = parse_slot(sampleB)
  )

message("\nEvaluable pairs at primary, by visit composition:")
print(evaluable_at_primary %>% count(visit_a, visit_b))

message("\nEvaluable pairs at primary, by slot composition:")
print(evaluable_at_primary %>% count(slot_a, slot_b))

readr::write_csv(evaluable_at_primary %>% count(slot_a, slot_b),
                 file.path(PATHS$calibration, "evaluable_pairs_at_primary_by_slot.csv"))

# ---- Heatmap (kept name cal_lexicographic_heatmap.pdf for downstream compat)
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {

  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

  ad1_primary <- cal_primary_hp$AD1_min
  slice_df <- cal_grid %>% filter(AD1_min == ad1_primary)

  mark_df <- selected %>%
    filter(AD1_min == ad1_primary) %>%
    mutate(selection = factor(selection,
                              levels = c("looser", "primary", "tighter")))

  off_slice_msg <- selected %>%
    filter(AD1_min != ad1_primary) %>%
    mutate(msg = sprintf("%s at AD1_min=%g", selection, AD1_min)) %>%
    pull(msg) %>% paste(collapse = "; ")

  make_panel <- function(fill_var, fill_label, fill_limits = c(0, 1)) {
    ggplot(slice_df, aes(x = DP_min, y = MAF_min)) +
      geom_tile(aes(fill = .data[[fill_var]]), color = "gray90", linewidth = 0.1) +
      scale_fill_viridis_c(option = "D", name = fill_label, limits = fill_limits) +
      geom_point(data = slice_df %>% filter(!pass_guards),
                 shape = 4, size = 1.8, color = "white", show.legend = FALSE) +
      geom_point(data = mark_df,
                 aes(shape = selection, size = selection),
                 color = "red", stroke = 1.1, fill = NA) +
      scale_shape_manual(values = c(looser = 24, primary = 21, tighter = 25)) +
      scale_size_manual (values = c(looser = 3.2, primary = 4.2, tighter = 3.2)) +
      labs(x = "DP_min", y = "MAF_min",
           title = fill_label,
           subtitle = sprintf("AD1_min = %g (primary slice)", ad1_primary)) +
      theme_minimal(base_size = 10) +
      theme(panel.grid = element_blank())
  }

  p_jac <- make_panel("jacc_med_eval", "Median Jaccard (evaluable pairs)")
  p_con <- make_panel("conf_med_eval", "Median symmetric confirmation rate")
  p_bot <- make_panel("prop_both0",    "Fraction Cat 1 pairs")

  caption_txt <- paste(
    "Three-tier ladder (all real-grid cells; no-filter sentinel retired 2026-05-03).",
    if (nchar(off_slice_msg)) sprintf("Off-slice: %s.", off_slice_msg) else "",
    sprintf("Tag: %s. Selection: hybrid Pareto.", TAG$cal_tag),
    sep = " "
  )

  p <- p_jac + p_con + p_bot +
    plot_annotation(
      title    = "Hybrid Pareto calibration: MAF x DP at primary AD1",
      subtitle = "Looser (up-triangle), Primary (circle), Tighter (down-triangle). White X = failed grid guard.",
      caption  = caption_txt
    )

  ggsave(file.path(PATHS$calibration, "cal_lexicographic_heatmap.pdf"), p,
         width = 15, height = 4.8)

  # ---- Pareto frontier scatter plot ------------------------------------
  pf <- ggplot(cal_eligible,
               aes(x = jacc_med_eval, y = conf_med_eval)) +
    geom_point(aes(alpha = is_frontier), color = "gray40", size = 1.4) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25), guide = "none") +
    geom_point(data = frontier, color = "black", size = 2) +
    geom_point(data = cal_primary_hp,
               color = "red", shape = 21, size = 4.5, stroke = 1.2, fill = NA) +
    annotate("point", x = 1, y = 1, shape = 8, size = 3, color = "blue") +
    annotate("text",  x = 1, y = 1, label = "  ideal (1,1)",
             hjust = 0, vjust = 0.5, size = 3, color = "blue") +
    coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05)) +
    labs(x = "Median Jaccard (evaluable pairs)",
         y = "Median symmetric confirmation rate (evaluable pairs)",
         title = "Pareto frontier on (Jaccard, Confirm)",
         subtitle = sprintf("Black = frontier; red ring = primary; dim gray = dominated. Tag: %s",
                            TAG$cal_tag)) +
    theme_minimal(base_size = 10)

  ggsave(file.path(PATHS$calibration, "cal_hybrid_pareto_frontier.pdf"), pf,
         width = 6, height = 5.5)
}

# ---- Save thresholds for downstream consumption ----------------------------
to_thr_list <- function(row, MAF_max = MAF_MAX) {
  if (is.null(row) || nrow(row) == 0) return(NULL)
  list(
    DP_min  = row$DP_min,
    AD1_min = row$AD1_min,
    MAF_min = row$MAF_min,
    MAF_max = MAF_max
  )
}

thr_looser  <- to_thr_list(selected %>% filter(selection == "looser"))
thr_primary <- to_thr_list(selected %>% filter(selection == "primary"))
thr_tighter <- to_thr_list(selected %>% filter(selection == "tighter"))

if (!is.null(thr_looser))
  saveRDS(thr_looser,  file.path(PATHS$calibration, "thr_looser.rds"))
saveRDS(thr_primary, file.path(PATHS$calibration, "thr_primary.rds"))
# Stage 1 backup; 02b will overwrite thr_primary.rds with the bootstrap-
# stabilized cell, but thr_primary_observed.rds is preserved as the Stage 1
# reference for downstream reporting.
saveRDS(thr_primary, file.path(PATHS$calibration, "thr_primary_observed.rds"))
if (!is.null(thr_tighter))
  saveRDS(thr_tighter, file.path(PATHS$calibration, "thr_tighter.rds"))

thresholds <- list(
  looser  = thr_looser,
  primary = thr_primary,
  tighter = thr_tighter
)
saveRDS(thresholds, file.path(PATHS$calibration, "thresholds.rds"))

message("\n[02_calibrate_hybrid_pareto] Done. Thresholds written to ",
        PATHS$calibration, "/")
if (!is.null(thr_looser))
  message(sprintf("  Looser : DP=%g, AD1=%g, MAF=%g (min stringency on frontier)",
                  thr_looser$DP_min, thr_looser$AD1_min, thr_looser$MAF_min))
message(sprintf("  Primary: DP=%g, AD1=%g, MAF=%g (closest-to-ideal; Stage 1; 02b will stabilize)",
                thr_primary$DP_min, thr_primary$AD1_min, thr_primary$MAF_min))
if (!is.null(thr_tighter))
  message(sprintf("  Tighter: DP=%g, AD1=%g, MAF=%g (max stringency on frontier)",
                  thr_tighter$DP_min, thr_tighter$AD1_min, thr_tighter$MAF_min))
