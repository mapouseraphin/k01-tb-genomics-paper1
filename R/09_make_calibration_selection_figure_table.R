# =============================================================================
# 09_make_calibration_selection_figure_table.R
#
# Produces Figure 2 (calibration selection) and Table 2 (calibrated thresholds
# with stability metadata) for the main Methods/Results section.
#
# Method-aware (post 2026-05-03):
#   SELECTION_METHOD = "lex"    -> Panel A highlights top-K lexicographic
#                                   candidates; reads cal_lexicographic_*.csv
#   SELECTION_METHOD = "pareto" -> Panel A highlights Pareto frontier;
#                                   reads cal_hybrid_pareto_*.csv
#
# Three-tier ladder (post 2026-05-03; all real-grid cells; no sentinel):
#   Looser  : real-grid less-stringent reference; written by 02_calibrate_*.R
#             as thr_looser.rds, then OVERWRITTEN by 02b_calibration_bootstrap.R
#             with axis-wise dominance against the stabilized Primary.
#   Primary : bootstrap-stabilized modal cell from 02b; written as
#             thr_primary.rds (and thr_primary_bootstrap.rds backup).
#             Stage 1 observed-data primary preserved as thr_primary_observed.rds.
#   Tighter : real-grid more-stringent reference; same dual-write pattern as
#             Looser. May be NULL if no axis-wise dominant cell exists in the
#             top-K=10 / top-K=30 candidate pool (e.g., Pareto runs where the
#             Primary sits at the corner of the threshold grid).
#
# Render-mode logic (post 2026-05-08; corrected):
#   The script renders two distinct figure modes, selected by whether the
#   bootstrap-stabilized modal cell equals the Stage 1 observed-data primary:
#
#     modal_eq_observed = TRUE  -> AGREE mode.
#       Panel A shows a single Primary marker (red circle).
#       Panel B colors that cell red and labels the legend category "Primary".
#       Subtitle and caption state that Stage 1 and Stage 2 agree.
#       Table 2 has no override footnote.
#
#     modal_eq_observed = FALSE -> OVERRIDE mode.
#       Panel A shows two Primary markers: bootstrap (red) and observed (blue).
#       Panel B colors the two cells red and blue with separate legend entries.
#       Subtitle and caption note the override.
#       Table 2 carries an override footnote.
#
#   Prior versions used dual-marker phrasing in the AGREE-mode Panel A subtitle
#   and a separate purple "modal = observed" Panel B category, which read as an
#   implied override. Fixed by branching marker-construction, palette, and
#   subtitle on modal_eq_observed.
#
# Figure 2 -- two panels:
#   Panel A: (Jaccard, Confirm) plane.
#            All guard-passing cells; the eligible cells (Pareto frontier
#            for "pareto"; top-K lex candidates for "lex") in solid black;
#            non-eligible cells in light gray; ideal point (1,1) marked with
#            a star; the named cells annotated by hollow markers (see render
#            mode above).
#   Panel B: Top-10 cells by bootstrap selection frequency. Horizontal bar
#            chart with bars colored by named-cell membership.
#
# Table 2 -- calibrated thresholds with bootstrap selection frequency:
#   Three rows (Looser, Primary, Tighter). Columns: threshold tuple, n
#   evaluable pairs, prop_both0, median Jaccard, median Confirm, and bootstrap
#   selection frequency. A footnote records the observed-data Stage 1 primary
#   only if it differs from the stabilized primary.
#
# Inputs (method-aware):
#   PATHS$calibration/cal_lexicographic_grid_summary.csv   (if lex)
#                  OR cal_hybrid_pareto_grid_summary.csv   (if pareto)
#   PATHS$calibration/cal_hybrid_pareto_frontier.csv       (only if pareto)
#   PATHS$calibration/cal_lexicographic_thresholds.csv     (both methods)
#   PATHS$calibration/thr_primary_observed.rds             (Stage 1 primary)
#   PATHS$calibration/thr_primary_bootstrap.rds            (Stage 2 stabilized primary)
#   PATHS$calibration/thr_looser.rds                       (real-grid looser tier)
#   PATHS$calibration/thr_tighter.rds                      (real-grid tighter tier)
#   PATHS$bootstrap/bootstrap_stability.rds
#
# Outputs:
#   PATHS$figures/figure2_calibration_selection.{pdf,png}
#   PATHS$tables/table2_calibrated_thresholds.{csv,tex}
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr); library(ggplot2)
  library(patchwork); library(scales)
})

dir.create(PATHS$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$tables,  recursive = TRUE, showWarnings = FALSE)

stopifnot(SELECTION_METHOD %in% c("lex", "pareto"))
TOP_K <- 10L  # for "lex": eligible cells = top-K by lexicographic ranking

# ---- Load (method-aware) ---------------------------------------------------
grid_filename <- if (SELECTION_METHOD == "pareto")
  "cal_hybrid_pareto_grid_summary.csv" else
  "cal_lexicographic_grid_summary.csv"

cal_grid <- read_csv(file.path(PATHS$calibration, grid_filename),
                     show_col_types = FALSE)

# Pareto frontier file is only written by the Pareto Stage 1 script
frontier_in <- if (SELECTION_METHOD == "pareto") {
  read_csv(file.path(PATHS$calibration, "cal_hybrid_pareto_frontier.csv"),
           show_col_types = FALSE)
} else NULL

thr_observed  <- readRDS(file.path(PATHS$calibration, "thr_primary_observed.rds"))
thr_bootstrap <- readRDS(file.path(PATHS$calibration, "thr_primary_bootstrap.rds"))

# Real-grid Looser and Tighter (post 2026-05-03; both written by 02_calibrate_*.R
# and overwritten by 02b_calibration_bootstrap.R against the stabilized Primary,
# or removed if no axis-wise dominant candidate exists). Use file.exists() to
# avoid issuing a gzfile warning when 02b deletes a NULL tier file.
looser_path  <- file.path(PATHS$calibration, "thr_looser.rds")
tighter_path <- file.path(PATHS$calibration, "thr_tighter.rds")
thr_looser   <- if (file.exists(looser_path))  readRDS(looser_path)  else NULL
thr_tighter  <- if (file.exists(tighter_path)) readRDS(tighter_path) else NULL

stability <- readRDS(file.path(PATHS$bootstrap, "bootstrap_stability.rds"))

# ---- Helpers ---------------------------------------------------------------
lookup_cell <- function(grid, DP, AD1, MAF) {
  grid %>% filter(DP_min == DP, AD1_min == AD1, MAF_min == MAF) %>% slice(1)
}

# Bootstrap selection frequency for an arbitrary cell. Returns 0 if the cell
# was never selected; returns NA only if cell_summary is missing.
lookup_boot_prop <- function(stab, DP, AD1, MAF) {
  if (is.null(stab$cell_summary)) return(NA_real_)
  hit <- stab$cell_summary %>%
    filter(DP_min == DP, AD1_min == AD1, MAF_min == MAF)
  if (nrow(hit) == 0) 0 else hit$prop[1]
}

fmt_pct <- function(p) {
  if (is.na(p)) "—" else sprintf("%.1f%%", 100 * p)
}

# ---- Annotate "eligible" membership in cal_grid ----------------------------
# Under "pareto", eligible = on the Pareto frontier on (J, C).
# Under "lex",    eligible = top-K by lexicographic ranking (the looser/tighter
#                            candidate pool).
guard_pass <- cal_grid %>% filter(pass_guards)

if (SELECTION_METHOD == "pareto") {
  front_keys <- frontier_in %>%
    transmute(key = paste(DP_min, AD1_min, MAF_min, sep = "|"))
  guard_pass <- guard_pass %>%
    mutate(key = paste(DP_min, AD1_min, MAF_min, sep = "|"),
           is_eligible = key %in% front_keys$key)
  eligible_label <- sprintf("Pareto frontier (%d cells)",
                            sum(guard_pass$is_eligible))
} else {
  # Lexicographic ranking on the same six criteria used in 02_calibrate_lexicographic.R
  cal_ranked <- guard_pass %>%
    arrange(prop_both0,
            desc(prop_overlap),
            desc(jacc_med_eval),
            desc(overlap_med_eval),
            desc(conf_med_eval),
            desc(r2_med_eval)) %>%
    mutate(lex_rank = row_number())
  guard_pass <- guard_pass %>%
    left_join(cal_ranked %>% select(DP_min, AD1_min, MAF_min, lex_rank),
              by = c("DP_min", "AD1_min", "MAF_min")) %>%
    mutate(is_eligible = !is.na(lex_rank) & lex_rank <= TOP_K)
  eligible_label <- sprintf("Top-%d lex candidates", TOP_K)
}

# ---- Render mode -----------------------------------------------------------
# AGREE: bootstrap modal cell == Stage 1 observed-data primary.
# OVERRIDE: bootstrap modal cell differs from Stage 1 observed-data primary.
# Marker construction, palette, and subtitle all branch on this.

modal_eq_observed <- (thr_observed$DP_min  == thr_bootstrap$DP_min) &
                     (thr_observed$AD1_min == thr_bootstrap$AD1_min) &
                     (thr_observed$MAF_min == thr_bootstrap$MAF_min)

# ---- Marker frames for Panel A ---------------------------------------------
mk_looser    <- if (!is.null(thr_looser))
                  lookup_cell(cal_grid, thr_looser$DP_min,
                              thr_looser$AD1_min, thr_looser$MAF_min) else NULL
mk_observed  <- lookup_cell(cal_grid, thr_observed$DP_min,
                            thr_observed$AD1_min, thr_observed$MAF_min)
mk_bootstrap <- lookup_cell(cal_grid, thr_bootstrap$DP_min,
                            thr_bootstrap$AD1_min, thr_bootstrap$MAF_min)
mk_tighter   <- if (!is.null(thr_tighter))
                  lookup_cell(cal_grid, thr_tighter$DP_min,
                              thr_tighter$AD1_min, thr_tighter$MAF_min) else NULL

# ---- Panel A: (Jaccard, Confirm) plane -------------------------------------
# Subtitle branches on render mode.
panel_A_subtitle <- if (modal_eq_observed) {
  sprintf(paste0(
    "Eligible: %s (black).
    Markers: down-triangle = Looser, ",
    "red circle = Primary, up-triangle = Tighter. ",
    "Stage 1 and Stage 2 select the same cell."),
    eligible_label)
} else {
  sprintf(paste0(
    "Eligible: %s (black). Markers: down-triangle = Looser, ",
    "circle = Primary (red = Stage 2 stabilized, blue = Stage 1 observed), ",
    "up-triangle = Tighter."),
    eligible_label)
}

panel_A <- ggplot(guard_pass,
                  aes(x = jacc_med_eval, y = conf_med_eval)) +
  # Background: ineligible cells in light gray
  geom_point(data = guard_pass %>% filter(!is_eligible),
             color = "gray70", size = 1.6, alpha = 0.6) +
  # Foreground: eligible cells in solid black
  geom_point(data = guard_pass %>% filter(is_eligible),
             color = "black", size = 2.4) +
  # Looser: hollow inverted triangle, dark green
  { if (!is.null(mk_looser) && nrow(mk_looser) > 0)
      geom_point(data = mk_looser,
                 shape = 25, size = 4.5, stroke = 1.2,
                 color = "darkgreen", fill = NA)
    else NULL } +
  # Tighter: hollow upright triangle, purple
  { if (!is.null(mk_tighter) && nrow(mk_tighter) > 0)
      geom_point(data = mk_tighter,
                 shape = 24, size = 4.5, stroke = 1.2,
                 color = "purple4", fill = NA)
    else NULL } +
  # Observed primary (Stage 1): hollow circle, blue (only in OVERRIDE mode)
  { if (!modal_eq_observed && nrow(mk_observed) > 0)
      geom_point(data = mk_observed,
                 shape = 21, size = 5.0, stroke = 1.3,
                 color = "steelblue3", fill = NA)
    else NULL } +
  # Primary marker (red): in AGREE mode this is THE Primary cell; in OVERRIDE
  # mode this is the Stage 2 stabilized cell.
  geom_point(data = mk_bootstrap,
             shape = 21, size = 6.0, stroke = 1.6,
             color = "red3", fill = NA) +
  # Ideal point (1, 1)
  annotate("point", x = 1, y = 1, shape = 8, size = 3.5, color = "blue") +
  annotate("text",  x = 1, y = 1, label = "ideal (1,1)",
           hjust = 1.15, vjust = -0.5, size = 3, color = "blue") +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(0, 1.05)) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Median Jaccard (evaluable pairs)",
       y = "Median confirmation rate (evaluable pairs)",
       title = sprintf("A. Calibration selection on (Jaccard, Confirm) plane (%s)",
                       SELECTION_METHOD),
       subtitle = panel_A_subtitle) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# ---- Panel B: Top-10 bootstrap selection frequency -------------------------
top_n_cells <- 10L

# Marker assignment branches on render mode. In AGREE mode the modal cell is
# tagged "Primary" (single category, red). In OVERRIDE mode the modal cell and
# the Stage 1 observed cell are tagged separately (red and blue).
top_cells <- stability$cell_summary %>%
  head(top_n_cells) %>%
  mutate(label = sprintf("(%g, %g, %g)", DP_min, AD1_min, MAF_min),
         is_looser = if (!is.null(thr_looser))
                       (DP_min == thr_looser$DP_min &
                        AD1_min == thr_looser$AD1_min &
                        MAF_min == thr_looser$MAF_min)
                     else FALSE,
         is_observed = (DP_min == thr_observed$DP_min &
                        AD1_min == thr_observed$AD1_min &
                        MAF_min == thr_observed$MAF_min),
         is_modal    = (DP_min == thr_bootstrap$DP_min &
                        AD1_min == thr_bootstrap$AD1_min &
                        MAF_min == thr_bootstrap$MAF_min),
         is_tighter  = if (!is.null(thr_tighter))
                         (DP_min == thr_tighter$DP_min &
                          AD1_min == thr_tighter$AD1_min &
                          MAF_min == thr_tighter$MAF_min)
                       else FALSE)

if (modal_eq_observed) {
  # AGREE mode: single "Primary" category
  top_cells <- top_cells %>%
    mutate(marker = case_when(
      is_modal    ~ "Primary",
      is_looser   ~ "looser tier",
      is_tighter  ~ "tighter tier",
      TRUE        ~ "other"
    ),
    cum_prop = cumsum(prop)) %>%
    arrange(prop) %>%
    mutate(label = factor(label, levels = label))

  bar_palette <- c(
    "Primary"      = "red3",
    "looser tier"  = "darkgreen",
    "tighter tier" = "purple4",
    "other"        = "gray60"
  )
  legend_order <- c("Primary", "looser tier", "tighter tier", "other")

} else {
  # OVERRIDE mode: separate Stage 1 / Stage 2 categories
  top_cells <- top_cells %>%
    mutate(marker = case_when(
      is_modal     ~ "modal (Stage 2)",
      is_observed  ~ "observed (Stage 1)",
      is_looser    ~ "looser tier",
      is_tighter   ~ "tighter tier",
      TRUE         ~ "other"
    ),
    cum_prop = cumsum(prop)) %>%
    arrange(prop) %>%
    mutate(label = factor(label, levels = label))

  bar_palette <- c(
    "modal (Stage 2)"    = "red3",
    "observed (Stage 1)" = "steelblue3",
    "looser tier"        = "darkgreen",
    "tighter tier"       = "purple4",
    "other"              = "gray60"
  )
  legend_order <- c("modal (Stage 2)", "observed (Stage 1)",
                    "looser tier", "tighter tier", "other")
}

present_categories <- unique(top_cells$marker)
legend_breaks <- intersect(legend_order, present_categories)

panel_B <- ggplot(top_cells, aes(x = label, y = prop, fill = marker)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = bar_palette, name = NULL,
                    breaks = legend_breaks, drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     limits = c(0, max(top_cells$prop) * 1.18)) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * prop)),
            hjust = -0.18, size = 3, color = "gray20") +
  coord_flip() +
  labs(x = NULL, y = "Proportion of bootstrap reps selecting cell",
       title = sprintf("B. Top-%d cells by bootstrap selection frequency",
                       nrow(top_cells)),
       subtitle = sprintf(
         "B = %d resamples; %d (%.1f%%) failed; n_ok = %d. Method: %s.",
         stability$n_reps, stability$n_failed,
         100 * stability$n_failed / stability$n_reps,
         stability$n_ok, SELECTION_METHOD)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.justification = c(0, 0),
        plot.title = element_text(face = "bold"))

# ---- Combine panels --------------------------------------------------------
caption_txt <- paste(
  "Three-tier calibration ladder (all real-grid cells; no-filter sentinel retired):",
  if (!is.null(thr_looser))
    sprintf("Looser (DP=%g, AD1=%g, MAF=%g);",
            thr_looser$DP_min, thr_looser$AD1_min, thr_looser$MAF_min)
  else "Looser n/a;",
  sprintf("Primary (DP=%g, AD1=%g, MAF=%g);",
          thr_bootstrap$DP_min, thr_bootstrap$AD1_min, thr_bootstrap$MAF_min),
  if (!is.null(thr_tighter))
    sprintf("Tighter (DP=%g, AD1=%g, MAF=%g).",
            thr_tighter$DP_min, thr_tighter$AD1_min, thr_tighter$MAF_min)
  else "Tighter n/a.",
  sprintf("\nStage 1: %s on observed cal_pairs selects (DP=%g, AD1=%g, MAF=%g) as primary.",
          SELECTION_METHOD,
          thr_observed$DP_min, thr_observed$AD1_min, thr_observed$MAF_min),
  sprintf("\nStage 2: B = %d nonparametric pair-level bootstrap reruns the same procedure;",
          stability$n_reps),
  "the modal cell across reps is taken as the stabilized Primary.",
  if (modal_eq_observed) "Stage 1 and Stage 2 selections agree."
  else sprintf("Stabilized Primary: (DP=%g, AD1=%g, MAF=%g) [override of Stage 1].",
               thr_bootstrap$DP_min, thr_bootstrap$AD1_min, thr_bootstrap$MAF_min),
  sprintf("Joint stability of Primary: %.1f%%.",
          100 * stability$prop_modal_joint),
  sep = " "
)

fig2 <- (panel_A | panel_B) +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    caption = caption_txt,
    theme = theme(plot.caption = element_text(hjust = 0, size = 9,
                                               color = "gray30",
                                               margin = margin(t = 6)))
  )

ggsave(file.path(PATHS$figures, "figure2_calibration_selection.pdf"),
       fig2, width = 12.5, height = 5.4, useDingbats = FALSE)
ggsave(file.path(PATHS$figures, "figure2_calibration_selection.png"),
       fig2, width = 12.5, height = 5.4, dpi = 300)

message("[09] Figure 2 written to ", PATHS$figures,
        "/figure2_calibration_selection.{pdf,png}")

# ---- Table 2: three-tier calibrated thresholds with bootstrap selection % --
looser_row  <- if (!is.null(thr_looser))
                 lookup_cell(cal_grid, thr_looser$DP_min,
                             thr_looser$AD1_min, thr_looser$MAF_min) else NULL
primary_row <- lookup_cell(cal_grid, thr_bootstrap$DP_min,
                           thr_bootstrap$AD1_min, thr_bootstrap$MAF_min)
tighter_row <- if (!is.null(thr_tighter))
                 lookup_cell(cal_grid, thr_tighter$DP_min,
                             thr_tighter$AD1_min, thr_tighter$MAF_min) else NULL

build_row <- function(selection_label, DP_min, AD1_min, MAF_min, MAF_max,
                      grid_row, boot_prop) {
  tibble(
    Selection         = selection_label,
    DP_min            = DP_min,
    AD1_min           = AD1_min,
    MAF_min           = MAF_min,
    MAF_max           = MAF_max,
    n_eval_pairs      = if (!is.null(grid_row) && nrow(grid_row) > 0)
                          grid_row$pairs_evaluable[1] else NA_integer_,
    prop_both0        = if (!is.null(grid_row) && nrow(grid_row) > 0)
                          round(grid_row$prop_both0[1], 3) else NA_real_,
    median_Jaccard    = if (!is.null(grid_row) && nrow(grid_row) > 0)
                          round(grid_row$jacc_med_eval[1], 3) else NA_real_,
    median_Confirm    = if (!is.null(grid_row) && nrow(grid_row) > 0)
                          round(grid_row$conf_med_eval[1], 3) else NA_real_,
    P_bootstrap       = fmt_pct(boot_prop)
  )
}

# Bootstrap selection probabilities
prop_looser  <- if (!is.null(thr_looser)){
                  lookup_boot_prop(stability, thr_looser$DP_min,
                                   thr_looser$AD1_min, thr_looser$MAF_min)
}                else NA_real_
prop_primary <- stability$prop_modal_joint
prop_tighter <- if (!is.null(thr_tighter)){
                  lookup_boot_prop(stability, thr_tighter$DP_min,
                                   thr_tighter$AD1_min, thr_tighter$MAF_min)
}                else NA_real_

# In AGREE mode the Primary row label is plain "Primary"; in OVERRIDE mode it
# is "Primary (bootstrap-stabilized)" to flag that Stage 2 differs from Stage 1.
primary_row_label <- if (modal_eq_observed) "Primary" else "Primary (bootstrap-stabilized)"

table2 <- bind_rows(
  if (!is.null(thr_looser))
    build_row("Looser",
              thr_looser$DP_min, thr_looser$AD1_min,
              thr_looser$MAF_min, thr_looser$MAF_max,
              grid_row = looser_row, boot_prop = prop_looser),
  build_row(primary_row_label,
            thr_bootstrap$DP_min, thr_bootstrap$AD1_min,
            thr_bootstrap$MAF_min, thr_bootstrap$MAF_max,
            grid_row = primary_row, boot_prop = prop_primary),
  if (!is.null(thr_tighter))
    build_row("Tighter",
              thr_tighter$DP_min, thr_tighter$AD1_min,
              thr_tighter$MAF_min, thr_tighter$MAF_max,
              grid_row = tighter_row, boot_prop = prop_tighter)
)

# Add a footnote line at the bottom about Stage 1 observed primary if differs
if (!modal_eq_observed) {
  table2_with_footnote <- table2 %>%
    mutate(across(everything(), as.character))
  footnote_row <- tibble(
    Selection = sprintf(
      paste0("Note: Stage 1 observed-data primary (%s) was (DP=%g, AD1=%g, MAF=%g), ",
             "selected in %.1f%% of bootstrap reps; ",
             "overridden by stabilized modal cell shown."),
      SELECTION_METHOD,
      thr_observed$DP_min, thr_observed$AD1_min, thr_observed$MAF_min,
      100 * stability$prop_observed_primary
    ),
    DP_min = "", AD1_min = "", MAF_min = "", MAF_max = "",
    n_eval_pairs = "", prop_both0 = "", median_Jaccard = "",
    median_Confirm = "", P_bootstrap = ""
  )
  table2_with_footnote <- bind_rows(table2_with_footnote, footnote_row)
} else {
  table2_with_footnote <- table2
}

write_csv(table2_with_footnote,
          file.path(PATHS$tables, "table2_calibrated_thresholds.csv"))

# ---- LaTeX rendering -------------------------------------------------------
build_latex_table <- function(t2, modal_eq_obs, thr_obs, prop_obs, method) {
  hdr <- c("Selection", "DP$_{\\min}$", "AD1$_{\\min}$",
           "MAF$_{\\min}$", "MAF$_{\\max}$",
           "$n_{\\text{eval}}$", "$p_{\\text{both0}}$",
           "Median Jaccard", "Median Confirm",
           "$P(\\text{bootstrap})$")
  body <- apply(t2, 1, function(r) paste(r, collapse = " & "))
  out <- c(
    "\\begin{tabular}{lccccrrrrr}",
    "\\hline",
    paste(paste(hdr, collapse = " & "), "\\\\"),
    "\\hline",
    paste(body, "\\\\"),
    "\\hline",
    "\\end{tabular}"
  )
  if (!modal_eq_obs) {
    out <- c(out,
             sprintf(
               paste0("\\\\\\footnotesize\\textit{Stage 1 observed-data %s primary:} ",
                      "(DP=%g, AD1=%g, MAF=%g), bootstrap selection rate %.1f\\%%; ",
                      "overridden by modal cell."),
               method,
               thr_obs$DP_min, thr_obs$AD1_min, thr_obs$MAF_min,
               100 * prop_obs))
  }
  paste(out, collapse = "\n")
}

writeLines(build_latex_table(table2, modal_eq_observed,
                              thr_observed, stability$prop_observed_primary,
                              SELECTION_METHOD),
           file.path(PATHS$tables, "table2_calibrated_thresholds.tex"))

# ---- Console summary -------------------------------------------------------
message(sprintf("\n[09] Three-tier ladder (method = %s; render mode = %s):",
                SELECTION_METHOD,
                if (modal_eq_observed) "AGREE" else "OVERRIDE"))
if (!is.null(thr_looser))
  message(sprintf("  Looser : DP=%g, AD1=%g, MAF=%g  (%.1f%% bootstrap)",
                  thr_looser$DP_min, thr_looser$AD1_min, thr_looser$MAF_min,
                  100 * prop_looser))
message(sprintf(  "  Primary: DP=%g, AD1=%g, MAF=%g  (%.1f%% bootstrap, modal)",
                  thr_bootstrap$DP_min, thr_bootstrap$AD1_min, thr_bootstrap$MAF_min,
                  100 * prop_primary))
if (!is.null(thr_tighter))
  message(sprintf("  Tighter: DP=%g, AD1=%g, MAF=%g  (%.1f%% bootstrap)",
                  thr_tighter$DP_min, thr_tighter$AD1_min, thr_tighter$MAF_min,
                  100 * prop_tighter))

message("\n[09] Table 2 written:")
message("  CSV: ", PATHS$tables, "/table2_calibrated_thresholds.csv")
message("  TeX: ", PATHS$tables, "/table2_calibrated_thresholds.tex")
message("\nTable 2 contents:")
print(as.data.frame(table2), row.names = FALSE)
if (!modal_eq_observed) {
  message(sprintf(
    paste0("\nFootnote: Stage 1 primary (%s; DP=%g, AD1=%g, MAF=%g) ",
           "selected in %.1f%% of reps; overridden by stabilized modal."),
    SELECTION_METHOD,
    thr_observed$DP_min, thr_observed$AD1_min, thr_observed$MAF_min,
    100 * stability$prop_observed_primary))
}
