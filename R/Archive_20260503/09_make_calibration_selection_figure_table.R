# =============================================================================
# 09_make_calibration_selection_figure_table.R
#
# Produces Figure 2 (calibration selection) and Table 2 (calibrated thresholds
# with stability metadata) for the main Methods/Results section.
#
# THREE-TIER CALIBRATION LADDER (relabel 2026-04-30):
#   - Looser  = (40, 2, 0.05)  -- Pareto-frontier cell with AD1=2, less stringent
#               than Primary on AD1; serves as the less-stringent ladder reference.
#               (On-disk file: thr_tighter.rds, written by 02_calibrate_hybrid_pareto.R
#               under its prior label; the rename to "Looser" is local to this
#               script and reflects honest direction. The Pareto pipeline output
#               remains unchanged.)
#   - Primary = (70, 6, 0.05)  -- bootstrap-stabilized modal cell; headline rule.
#   - Tighter = (80, 6, 0.05)  -- selected from bootstrap top-N as a strictly
#               more stringent reference (DP +10, all other axes equal to Primary;
#               14.7% bootstrap selection frequency).
#
# IMPORTANT: LCA is fit only at the Primary rule. Looser and Tighter are
# calibration-ladder references that show concordance behavior at less / more
# stringent rules; they are not used as inputs to downstream LCA fitting.
#
# Figure 2 -- two panels:
#   Panel A: Pareto frontier on the (Jaccard, Confirm) plane.
#            All guard-passing cells as points; frontier cells solid black;
#            dominated cells light gray; the ideal point (1,1) marked with a
#            star; the four named cells annotated by hollow markers
#            (looser inverted-triangle green; observed-Stage-1 ring blue,
#            shown only if it differs from bootstrap; bootstrap-stabilized
#            ring red, larger; tighter triangle purple).
#   Panel B: Top-10 cells by bootstrap selection frequency. Horizontal bar
#            chart with bars colored by named-cell membership.
#
# Table 2 -- calibrated thresholds with bootstrap selection frequency:
#   Three rows (Looser, Primary, Tighter). Columns: threshold tuple, n
#   evaluable pairs, prop_both0, median Jaccard, median Confirm, and bootstrap
#   selection frequency. A footnote records the observed-data Stage 1 primary
#   if it differs from the stabilized primary.
#
# Inputs:
#   PATHS$calibration/cal_hybrid_pareto_grid_summary.csv
#   PATHS$calibration/cal_hybrid_pareto_frontier.csv
#   PATHS$calibration/thr_primary_observed.rds          (Stage 1 primary)
#   PATHS$calibration/thr_primary_bootstrap.rds         (Stage 2 stabilized primary)
#   PATHS$calibration/thr_tighter.rds                   (read AS the Looser tier)
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

# =============================================================================
# CONFIG: Tighter tier (more stringent than Primary on >=1 axis)
# =============================================================================
# The Tighter tier is selected from the bootstrap top-N as a cell that is
# strictly more stringent than Primary on at least one axis and not less
# stringent on any other. (80, 6, 0.05) was chosen because:
#   (a) it sits at the DP grid maximum,
#   (b) it holds AD1 and MAF equal to Primary, so the ladder reads as a
#       single-axis stringency move (DP only),
#   (c) it was selected in 14.7% of bootstrap reps (top-3 cell by frequency).
# Edit this block to relocate the Tighter tier; no upstream rerun required.
TIGHTER_CELL <- list(
  DP_min  = 80,
  AD1_min = 6,
  MAF_min = 0.05,
  MAF_max = 0.45
)

# ---- Load ------------------------------------------------------------------
cal_grid    <- read_csv(file.path(PATHS$calibration,
                                  "cal_hybrid_pareto_grid_summary.csv"),
                        show_col_types = FALSE)
frontier_in <- read_csv(file.path(PATHS$calibration,
                                  "cal_hybrid_pareto_frontier.csv"),
                        show_col_types = FALSE)

thr_observed  <- readRDS(file.path(PATHS$calibration, "thr_primary_observed.rds"))
thr_bootstrap <- readRDS(file.path(PATHS$calibration, "thr_primary_bootstrap.rds"))

# NOTE: The on-disk file thr_tighter.rds was written by 02_calibrate_hybrid_pareto.R
# as the "max-stringency on Pareto frontier" cell. After the 2026-04-30 relabel,
# this cell is semantically the LOOSER tier (less stringent than Primary on AD1).
# Read into a variable named to reflect the new semantics.
thr_looser <- tryCatch(
  readRDS(file.path(PATHS$calibration, "thr_tighter.rds")),
  error = function(e) NULL
)

thr_tighter <- TIGHTER_CELL  # new: more-stringent reference

stability <- readRDS(file.path(PATHS$bootstrap, "bootstrap_stability.rds"))

# ---- Helpers ---------------------------------------------------------------
lookup_cell <- function(grid, DP, AD1, MAF) {
  grid %>% filter(DP_min == DP, AD1_min == AD1, MAF_min == MAF) %>% slice(1)
}

# Bootstrap selection frequency for an arbitrary cell. Returns 0 if the cell
# was never selected by any rep; returns NA only if cell_summary is missing.
lookup_boot_prop <- function(stab, DP, AD1, MAF) {
  if (is.null(stab$cell_summary)) return(NA_real_)
  hit <- stab$cell_summary %>%
    filter(DP_min == DP, AD1_min == AD1, MAF_min == MAF)
  if (nrow(hit) == 0) 0 else hit$prop[1]
}

# Format bootstrap frequency for a table cell (handles NA gracefully)
fmt_pct <- function(p) {
  if (is.na(p)) "—" else sprintf("%.1f%%", 100 * p)
}

# ---- Annotate frontier membership in cal_grid ------------------------------
guard_pass <- cal_grid %>% filter(pass_guards)
front_keys <- frontier_in %>%
  transmute(key = paste(DP_min, AD1_min, MAF_min, sep = "|"))
guard_pass <- guard_pass %>%
  mutate(key = paste(DP_min, AD1_min, MAF_min, sep = "|"),
         is_frontier = key %in% front_keys$key)

# ---- Marker frames for Panel A ---------------------------------------------
# Use cal_grid (full grid) for marker lookup so cells that fail observed-data
# guards (e.g., Primary at n_eval=9 < MIN_EVAL_PAIRS=10) still receive a marker.
mk_looser <- NULL
if (!is.null(thr_looser)) {
  mk_looser <- lookup_cell(cal_grid, thr_looser$DP_min,
                           thr_looser$AD1_min, thr_looser$MAF_min)
}
mk_observed  <- lookup_cell(cal_grid, thr_observed$DP_min,
                            thr_observed$AD1_min, thr_observed$MAF_min)
mk_bootstrap <- lookup_cell(cal_grid, thr_bootstrap$DP_min,
                            thr_bootstrap$AD1_min, thr_bootstrap$MAF_min)
mk_tighter   <- lookup_cell(cal_grid, thr_tighter$DP_min,
                            thr_tighter$AD1_min, thr_tighter$MAF_min)

modal_eq_observed <- (thr_observed$DP_min  == thr_bootstrap$DP_min) &
                     (thr_observed$AD1_min == thr_bootstrap$AD1_min) &
                     (thr_observed$MAF_min == thr_bootstrap$MAF_min)

# ---- Panel A: Pareto frontier on (Jaccard, Confirm) plane ------------------
panel_A <- ggplot(guard_pass,
                  aes(x = jacc_med_eval, y = conf_med_eval)) +
  # Background: dominated cells in light gray
  geom_point(data = guard_pass %>% filter(!is_frontier),
             color = "gray70", size = 1.6, alpha = 0.6) +
  # Foreground: frontier cells in solid black
  geom_point(data = guard_pass %>% filter(is_frontier),
             color = "black", size = 2.4) +
  # Looser: hollow inverted triangle, dark green
  { if (!is.null(mk_looser) && nrow(mk_looser) > 0)
      geom_point(data = mk_looser,
                 shape = 25, size = 4.5, stroke = 1.2,
                 color = "darkgreen", fill = NA)
    else NULL } +
  # Tighter: hollow upright triangle, purple
  { if (nrow(mk_tighter) > 0)
      geom_point(data = mk_tighter,
                 shape = 24, size = 4.5, stroke = 1.2,
                 color = "purple4", fill = NA)
    else NULL } +
  # Observed primary (Stage 1): hollow circle, blue (only if differs from bootstrap)
  { if (!modal_eq_observed && nrow(mk_observed) > 0)
      geom_point(data = mk_observed,
                 shape = 21, size = 5.0, stroke = 1.3,
                 color = "steelblue3", fill = NA)
    else NULL } +
  # Bootstrap-stabilized primary (Stage 2): hollow circle, red, headline marker
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
       title = "A. Pareto frontier on (Jaccard, Confirm)",
       subtitle = sprintf(
         "Frontier (black) of %d guard-passing cells; gray = dominated. Markers: \u25BD looser, \u25CB primary (red = bootstrap, blue = Stage 1), \u25B3 tighter.",
         nrow(guard_pass))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# ---- Panel B: Top-10 bootstrap selection frequency -------------------------
top_n_cells <- 10L
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
         is_tighter  = (DP_min == thr_tighter$DP_min &
                        AD1_min == thr_tighter$AD1_min &
                        MAF_min == thr_tighter$MAF_min),
         marker = case_when(
           is_modal & is_observed ~ "modal = observed",
           is_modal               ~ "modal (Stage 2)",
           is_observed            ~ "observed (Stage 1)",
           is_looser              ~ "looser tier",
           is_tighter             ~ "tighter tier",
           TRUE                   ~ "other"
         ),
         cum_prop = cumsum(prop)) %>%
  arrange(prop) %>%
  mutate(label = factor(label, levels = label))

bar_palette <- c(
  "modal = observed"   = "purple3",
  "modal (Stage 2)"    = "red3",
  "observed (Stage 1)" = "steelblue3",
  "looser tier"        = "darkgreen",
  "tighter tier"       = "purple4",
  "other"              = "gray60"
)

# Filter legend to categories that actually appear in the top-N
present_categories <- unique(top_cells$marker)
legend_breaks <- intersect(
  c("modal (Stage 2)", "observed (Stage 1)", "modal = observed",
    "looser tier", "tighter tier", "other"),
  present_categories
)

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
         "B = %d resamples; %d (%.1f%%) failed (no Pareto frontier); n_ok = %d",
         stability$n_reps, stability$n_failed,
         100 * stability$n_failed / stability$n_reps,
         stability$n_ok)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "top",
        legend.justification = c(0, 0),
        plot.title = element_text(face = "bold"))

# ---- Combine panels --------------------------------------------------------
caption_txt <- paste(
  "Three-tier calibration ladder:",
  if (!is.null(thr_looser))
    sprintf("Looser (DP=%g, AD1=%g, MAF=%g);",
            thr_looser$DP_min, thr_looser$AD1_min, thr_looser$MAF_min)
  else "Looser n/a;",
  sprintf("Primary (DP=%g, AD1=%g, MAF=%g);",
          thr_bootstrap$DP_min, thr_bootstrap$AD1_min, thr_bootstrap$MAF_min),
  sprintf("Tighter (DP=%g, AD1=%g, MAF=%g).",
          thr_tighter$DP_min, thr_tighter$AD1_min, thr_tighter$MAF_min),
  "\nStage 1: hybrid Pareto on observed cal_pairs selects",
  sprintf("(DP=%g, AD1=%g, MAF=%g) as closest-to-ideal on the frontier.",
          thr_observed$DP_min, thr_observed$AD1_min, thr_observed$MAF_min),
  "\nStage 2: B = 1000 nonparametric pair-level bootstrap reruns the same procedure;",
  "the modal cell across reps is taken as the stabilized Primary.",
  if (modal_eq_observed) "Stage 1 and Stage 2 selections agree."
  else sprintf("Stabilized Primary: (DP=%g, AD1=%g, MAF=%g) [override of Stage 1].",
               thr_bootstrap$DP_min, thr_bootstrap$AD1_min, thr_bootstrap$MAF_min),
  sprintf("Joint stability of Primary: %.1f%%.",
          100 * stability$prop_modal_joint),
  "\nLCA is fit only at the Primary rule; Looser and Tighter are calibration-ladder references.",
  sprintf("Spec tag: %s.", TAG$cal_tag),
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
looser_row <- NULL
if (!is.null(thr_looser)) {
  looser_row <- lookup_cell(cal_grid, thr_looser$DP_min,
                            thr_looser$AD1_min, thr_looser$MAF_min)
}
primary_row <- lookup_cell(cal_grid, thr_bootstrap$DP_min,
                           thr_bootstrap$AD1_min, thr_bootstrap$MAF_min)
tighter_row <- lookup_cell(cal_grid, thr_tighter$DP_min,
                           thr_tighter$AD1_min, thr_tighter$MAF_min)

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
prop_looser <- NA_real_
if (!is.null(thr_looser)) {
  prop_looser <- lookup_boot_prop(stability, thr_looser$DP_min,
                                  thr_looser$AD1_min, thr_looser$MAF_min)
}
prop_primary <- stability$prop_modal_joint
prop_tighter <- lookup_boot_prop(stability, thr_tighter$DP_min,
                                 thr_tighter$AD1_min, thr_tighter$MAF_min)

table2 <- bind_rows(
  if (!is.null(thr_looser))
    build_row("Looser",
              thr_looser$DP_min, thr_looser$AD1_min,
              thr_looser$MAF_min, thr_looser$MAF_max,
              grid_row = looser_row, boot_prop = prop_looser),
  build_row("Primary (bootstrap-stabilized)",
            thr_bootstrap$DP_min, thr_bootstrap$AD1_min,
            thr_bootstrap$MAF_min, thr_bootstrap$MAF_max,
            grid_row = primary_row, boot_prop = prop_primary),
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
      "Note: Stage 1 observed-data primary was (DP=%g, AD1=%g, MAF=%g), selected in %.1f%% of bootstrap reps; overridden by stabilized modal cell shown.",
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
# Builds a clean tabular environment without depending on xtable / kableExtra.
build_latex_table <- function(t2, modal_eq_obs, thr_obs, prop_obs) {
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
               "\\\\\\footnotesize\\textit{Stage 1 observed-data Pareto primary:} (DP=%g, AD1=%g, MAF=%g), bootstrap selection rate %.1f\\%%; overridden by modal cell.",
               thr_obs$DP_min, thr_obs$AD1_min, thr_obs$MAF_min,
               100 * prop_obs))
  }
  paste(out, collapse = "\n")
}

writeLines(build_latex_table(table2, modal_eq_observed,
                              thr_observed, stability$prop_observed_primary),
           file.path(PATHS$tables, "table2_calibrated_thresholds.tex"))

# ---- Console summary -------------------------------------------------------
message("\n[09] Three-tier ladder:")
if (!is.null(thr_looser))
  message(sprintf("  Looser : DP=%g, AD1=%g, MAF=%g  (%.1f%% bootstrap)",
                  thr_looser$DP_min, thr_looser$AD1_min, thr_looser$MAF_min,
                  100 * prop_looser))
message(sprintf(  "  Primary: DP=%g, AD1=%g, MAF=%g  (%.1f%% bootstrap, modal)",
                  thr_bootstrap$DP_min, thr_bootstrap$AD1_min, thr_bootstrap$MAF_min,
                  100 * prop_primary))
message(sprintf(  "  Tighter: DP=%g, AD1=%g, MAF=%g  (%.1f%% bootstrap)",
                  thr_tighter$DP_min, thr_tighter$AD1_min, thr_tighter$MAF_min,
                  100 * prop_tighter))

message("\n[09] Table 2 written:")
message("  CSV: ", PATHS$tables, "/table2_calibrated_thresholds.csv")
message("  TeX: ", PATHS$tables, "/table2_calibrated_thresholds.tex")
message("\nTable 2 contents:")
print(as.data.frame(table2), row.names = FALSE)
if (!modal_eq_observed) {
  message(sprintf(
    "\nFootnote: Stage 1 primary (DP=%g, AD1=%g, MAF=%g) selected in %.1f%% of reps; overridden by stabilized modal.",
    thr_observed$DP_min, thr_observed$AD1_min, thr_observed$MAF_min,
    100 * stability$prop_observed_primary))
}
