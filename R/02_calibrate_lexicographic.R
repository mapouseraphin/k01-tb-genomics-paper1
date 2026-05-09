# =============================================================================
# 02_calibrate_lexicographic.R
#
# Lexicographic calibration -- canonical Stage 1 selector (post 2026-05-03 lex
# revert; lex is canonical primary, Pareto is parallel sensitivity).
#
# Three-tier evidentiary ladder (no sentinel; all real-grid cells):
#   looser  : grid-picked by minimum stringency among top-K candidates,
#             excluding primary. Less stringent than primary on >= 1 axis.
#   primary : grid-picked by six-level lexicographic rule.
#   tighter : grid-picked by maximum stringency score among top-K candidates,
#             excluding primary. More stringent than primary on >= 1 axis.
#
# The "loose" no-filter sentinel was retired 2026-05-03. Looser is now a
# real-grid tier that satisfies the same guards as primary and tighter, so
# all three rows in cal_lexicographic_thresholds.csv are real cells with
# real concordance metrics.
#
# Pair-level filters when computing each median:
#   jacc_med_eval, conf_med_eval, overlap_med_eval:
#       both replicates non-empty (nA > 0 AND nB > 0)
#   r2_med_eval:
#       both non-empty AND intersection has >= 2 positions (immune to
#       zero-imputation in MAF correlation)
#
# Lex sort criteria (in order):
#   1. prop_both0        asc   (fewer zero-zero pairs)
#   2. prop_overlap      desc  (more pairs with any overlap)
#   3. jacc_med_eval     desc  (median Jaccard, evaluable pairs)
#   4. overlap_med_eval  desc  (median overlap coefficient)
#   5. conf_med_eval     desc  (median symmetric confirmation rate; REVISED
#                              2026-05-05 -- was median Jaccard previously
#                              under the name confirm_both, which equals
#                              Jaccard on the evaluable filter and added no
#                              discriminating information)
#   6. r2_med_eval       desc  (median R^2 of MAF; intersection >= 2 sites)
#
# Grid-level guards (canonical):
#   prop_both0      <= MAX_PROP_BOTH0  (canonical: 0.90; was 0.70 pre-2026-05-03)
#   pairs_evaluable >= MIN_EVAL_PAIRS  (canonical: 5;    was 10  pre-2026-05-03)
#
# Looser (top-K min stringency):
#   - Top K = 10 guard-passing candidates from the lexicographic sort.
#   - Looser = candidate with MIN stringency S = rank(DP)+rank(AD1)+rank(MAF).
#   - Fallback to K2 = 30 if top-K cannot supply a non-primary cell.
#
# Tighter (top-K max stringency):
#   - Same candidate pool as looser.
#   - Tighter = candidate with MAX stringency.
#   - Fallback to K2 = 30 if top-K cannot supply.
#
# NOTE: 02b_calibration_bootstrap.R OVERWRITES the looser/tighter cells
#       written here, using axis-wise dominance against the stabilized
#       Primary (rather than rank-sum stringency against the Stage 1
#       Primary). The values written here are provisional.
#
# Inputs:  PATHS$calibration/cal_pairs.rds
# Outputs (PATHS$calibration):
#   cal_lexicographic_grid_summary.csv
#   cal_lexicographic_thresholds.csv
#   cal_lexicographic_heatmap.pdf
#   thr_looser.rds, thr_primary.rds, thr_tighter.rds
#   thr_primary_observed.rds         (Stage 1 backup; preserved through 02b)
#   thresholds.rds                   (consolidated list with looser/primary/tighter)
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

# ---- Lexicographic sort -----------------------------------------------------
cal_ranked <- cal_grid %>%
  filter(pass_guards) %>%
  arrange(
    prop_both0,
    desc(prop_overlap),
    desc(jacc_med_eval),
    desc(overlap_med_eval),
    desc(conf_med_eval),
    desc(r2_med_eval)
  )

cal_primary_lex <- cal_ranked[1, ]
message("\nPrimary threshold (lexicographic):")
print(cal_primary_lex %>% select(DP_min, AD1_min, MAF_min,
                                  prop_both0, prop_overlap,
                                  jacc_med_eval, conf_med_eval, r2_med_eval,
                                  pairs_evaluable))

# ---- Looser and tighter via top-K stringency -------------------------------
add_stringency <- function(df) {
  df %>% mutate(S = min_rank(DP_min) + min_rank(AD1_min) + min_rank(MAF_min))
}
not_primary <- function(df) {
  df %>% filter(!(DP_min  == cal_primary_lex$DP_min &
                  AD1_min == cal_primary_lex$AD1_min &
                  MAF_min == cal_primary_lex$MAF_min))
}

cand_k  <- cal_ranked %>% slice_head(n = min(TOP_K,          nrow(cal_ranked))) %>% add_stringency()
cand_k2 <- cal_ranked %>% slice_head(n = min(TOP_K_FALLBACK, nrow(cal_ranked))) %>% add_stringency()

pick_tighter <- function(cands) {
  cands %>%
    not_primary() %>%
    arrange(desc(S), desc(DP_min), desc(AD1_min), desc(MAF_min)) %>%
    slice_head(n = 1)
}

# Looser: minimum stringency among top-K guard-passers, excluding primary.
# Tie-break: lowest DP, then lowest AD1, then lowest MAF (most permissive
# on each axis sequentially).
pick_looser <- function(cands) {
  cands %>%
    not_primary() %>%
    arrange(S, DP_min, AD1_min, MAF_min) %>%
    slice_head(n = 1)
}

cal_tighter_lex <- pick_tighter(cand_k)
if (nrow(cal_tighter_lex) == 0) cal_tighter_lex <- pick_tighter(cand_k2)
if (nrow(cal_tighter_lex) == 0) {
  cal_tighter_lex <- NULL
  message("\nNOTE: no tighter candidate available; top-K and fallback exhausted.")
}

cal_looser_lex <- pick_looser(cand_k)
if (nrow(cal_looser_lex) == 0) cal_looser_lex <- pick_looser(cand_k2)
if (nrow(cal_looser_lex) == 0) {
  cal_looser_lex <- NULL
  message("\nNOTE: no looser candidate available; top-K and fallback exhausted.")
}

# ---- Assemble selected thresholds ------------------------------------------
common_cols <- c("selection", "DP_min", "AD1_min", "MAF_min",
                 "prop_both0", "prop_overlap",
                 "jacc_med_eval", "conf_med_eval", "r2_med_eval",
                 "pairs_evaluable", "pass_guards")

selected <- bind_rows(
  if (!is.null(cal_looser_lex))
    cal_looser_lex %>% select(-any_of("S")) %>%
      mutate(selection = "looser") %>% select(any_of(common_cols)),
  cal_primary_lex %>% select(-any_of("S")) %>%
    mutate(selection = "primary") %>% select(any_of(common_cols)),
  if (!is.null(cal_tighter_lex))
    cal_tighter_lex %>% select(-any_of("S")) %>%
      mutate(selection = "tighter") %>% select(any_of(common_cols))
)

message("\nSelected thresholds:")
print(selected)

readr::write_csv(cal_grid,
                 file.path(PATHS$calibration, "cal_lexicographic_grid_summary.csv"))
readr::write_csv(selected,
                 file.path(PATHS$calibration, "cal_lexicographic_thresholds.csv"))

# ---- Visit-stratified evaluable pair diagnostic ----------------------------
prim <- cal_primary_lex
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

# ---- Heatmap ----------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE) &&
    requireNamespace("patchwork", quietly = TRUE)) {

  suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

  ad1_primary <- cal_primary_lex$AD1_min
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
      scale_shape_manual(values = c(looser  = 24, primary = 21, tighter = 25)) +
      scale_size_manual (values = c(looser  = 3.2, primary = 4.2, tighter = 3.2)) +
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
    sprintf("Tag: %s. Selection: lexicographic.", TAG$cal_tag),
    sep = " "
  )

  p <- p_jac + p_con + p_bot +
    plot_annotation(
      title    = "Lexicographic calibration: MAF x DP at primary AD1",
      subtitle = "Looser (up-triangle), Primary (circle), Tighter (down-triangle). White X = failed grid guard.",
      caption  = caption_txt
    )

  ggsave(file.path(PATHS$calibration, "cal_lexicographic_heatmap.pdf"), p,
         width = 15, height = 4.8)
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
# Stage 1 observed-data backup; 02b will overwrite thr_primary.rds with the
# bootstrap-stabilized cell, but thr_primary_observed.rds is preserved as the
# Stage 1 reference for downstream reporting.
saveRDS(thr_primary, file.path(PATHS$calibration, "thr_primary_observed.rds"))
if (!is.null(thr_tighter))
  saveRDS(thr_tighter, file.path(PATHS$calibration, "thr_tighter.rds"))

thresholds <- list(
  looser  = thr_looser,
  primary = thr_primary,
  tighter = thr_tighter
)
saveRDS(thresholds, file.path(PATHS$calibration, "thresholds.rds"))

message("\n[02_calibrate_lexicographic] Done. Thresholds written to ",
        PATHS$calibration, "/")
if (!is.null(thr_looser))
  message(sprintf("  Looser : DP=%g, AD1=%g, MAF=%g",
                  thr_looser$DP_min, thr_looser$AD1_min, thr_looser$MAF_min))
message(sprintf("  Primary: DP=%g, AD1=%g, MAF=%g (Stage 1; 02b will stabilize)",
                thr_primary$DP_min, thr_primary$AD1_min, thr_primary$MAF_min))
if (!is.null(thr_tighter))
  message(sprintf("  Tighter: DP=%g, AD1=%g, MAF=%g",
                  thr_tighter$DP_min, thr_tighter$AD1_min, thr_tighter$MAF_min))
