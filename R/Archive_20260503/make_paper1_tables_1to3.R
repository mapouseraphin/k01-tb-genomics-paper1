# =============================================================================
# make_paper1_tables_1to3.R
#
# Ghana-only Tables 1, 2, 3 for Paper 1. Adapted from 07_make_tables_figures.R
# with all Florida content removed, since Paper 1's post-pivot scope is
# Ghana-only (Florida moved to Paper 2 for HIV bias correction).
#
# Reads:
#   PATHS$meta/gh_meta.rds                                     (Table 1)
#   PATHS$calibration/cal_lexicographic_thresholds.csv         (Table 2)
#   PATHS$calibration/cal_pairs.rds                            (Table 3)
#   PATHS$calibration/thresholds.rds                           (Table 3)
#
# Writes (under PATHS$tables, spec-tagged):
#   table_1_cohort_characteristics.{csv,md}
#   table_2_calibrated_thresholds.{csv,md}
#   table_3_internal_validity.{csv,md}
#
# For Table 4 (LCA-derived) and Figures 3-4, see make_paper1_tables_figures.R.
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

# ---- helpers ---------------------------------------------------------------
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
write_table <- function(df, name, caption) {
  csv_path <- file.path(PATHS$tables, paste0(name, ".csv"))
  md_path  <- file.path(PATHS$tables, paste0(name, ".md"))
  write_csv(df, csv_path)
  md <- c(
    paste("###", caption), "",
    paste("|", paste(names(df), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(df)), collapse = " | "), "|"),
    apply(df, 1, function(r)
      paste0("| ", paste(ifelse(is.na(r), "—", r), collapse = " | "), " |"))
  )
  writeLines(md, md_path)
  message(sprintf("  -> %s", csv_path))
}

spec_caption <- sprintf(
  "Canonical specification: SNV_ONLY=%s, DROP_PPE_CAL=%s, DROP_PPE_APPLY=%s, M0_ONLY=%s. Tag: %s.",
  SNV_ONLY, DROP_PPE_CAL, DROP_PPE_APPLY, M0_ONLY, TAG$apply_tag
)

# =============================================================================
# Table 1 -- Ghana cohort characteristics
# =============================================================================
message("\n[Table 1] Ghana cohort characteristics")

gh_meta <- readRDS(file.path(PATHS$meta, "gh_meta.rds"))

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

table1 <- tibble(
  Characteristic = c(
    "Patients (n)",
    "Samples (n)",
    "Age, years (median, IQR)",
    "Sex: Female",
    "Sex: Male",
    "HIV: Negative",
    "HIV: Positive",
    "HIV: Unknown",
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
    sprintf("%d (%.1f%%)",
            sum(is.na(gh_pt$HIV)),
            100 * mean(is.na(gh_pt$HIV))),
    fmt_n_pct(gh_pt$lineage, "L4"),
    sprintf("%d (%.1f%%)",
            sum(gh_pt$lineage != "L4" & !is.na(gh_pt$lineage)),
            100 * mean(gh_pt$lineage != "L4", na.rm = TRUE)),
    fmt_med_iqr(gh_meta$depth_med)
  )
)
write_table(table1, "table_1_cohort_characteristics",
            paste("Table 1. Ghana cohort characteristics.",
                  "Demographics summarised at the patient level;",
                  "sequencing depth at the sample level."))

# =============================================================================
# Table 2 -- calibrated thresholds (lexicographic) + diagnostic columns
# =============================================================================
message("\n[Table 2] Calibrated thresholds (lexicographic)")

thr <- read_csv(file.path(PATHS$calibration,
                          "cal_lexicographic_thresholds.csv"),
                show_col_types = FALSE)

table2 <- thr %>%
  transmute(
    Selection           = factor(selection,
                                 levels = c("loose", "primary", "tighter")),
    DP_min              = DP_min,
    AD1_min             = AD1_min,
    MAF_min             = MAF_min,
    MAF_max             = MAF_MAX,
    n_eval_pairs        = pairs_evaluable,
    prop_both0          = round(prop_both0, 3),
    median_Jaccard      = round(jacc_med_eval, 3),
    mean_confirm_rate   = round(conf_med_eval, 3)
  ) %>%
  arrange(Selection)

write_table(table2, "table_2_calibrated_thresholds",
            paste("Table 2. Calibrated decision rules for iSNV detection.",
                  spec_caption,
                  "Loose is a no-filter sentinel; primary and tighter selected by lexicographic ranking on six concordance criteria."))

# =============================================================================
# Table 3 -- internal validity by sample-slot pair type at primary thresholds
# =============================================================================
message("\n[Table 3] Internal validity by pair type")

cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))
thr_obj   <- readRDS(file.path(PATHS$calibration, "thresholds.rds"))
prim      <- thr_obj$primary

pp <- cal_pairs %>%
  filter(DP_min == prim$DP_min, AD1_min == prim$AD1_min,
         MAF_min == prim$MAF_min) %>%
  mutate(slot_a = parse_slot(sampleA),
         slot_b = parse_slot(sampleB),
         pair_type = case_when(
           (slot_a == "EM"   & slot_b %in% c("S0.1","S0.2")) |
           (slot_b == "EM"   & slot_a %in% c("S0.1","S0.2"))   ~ "EM <-> morning_spot",
           (slot_a == "S0.1" & slot_b == "S0.2") |
           (slot_a == "S0.2" & slot_b == "S0.1")               ~ "S0.1 <-> S0.2",
           TRUE                                                 ~ "other"
         ))

table3 <- pp %>%
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

write_table(table3, "table_3_internal_validity",
            paste("Table 3. Within-visit replicate concordance at primary thresholds, by sample-slot pair type.",
                  spec_caption,
                  "Evaluable pairs have at least one calibrated iSNV in both replicates."))

message(sprintf("\n[make_paper1_tables_1to3] Done. Tables in %s/", PATHS$tables))
