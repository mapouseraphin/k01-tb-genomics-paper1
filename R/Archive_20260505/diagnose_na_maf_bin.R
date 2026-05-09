# =============================================================================
# diagnose_na_maf_bin.R
#
# One-page diagnostic on rows in the LCA dataset where maf_bin is NA.
# Purpose: characterise these rows so the Methods can footnote what is being
# excluded, and surface any data-cleaning issues upstream.
#
# Inputs:
#   <LCA_PATH>/lca_dataset.rds  (built by R/build_lca_dataset.R; pre-exclusion)
#
# Outputs:
#   <LCA_PATH>/feasibility/na_maf_bin_diagnostic.md
#   <LCA_PATH>/feasibility/na_maf_bin_rows.csv  (full 33-row dump)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

LCA_PATH  <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
FEAS_PATH <- file.path(LCA_PATH, "feasibility")
dir.create(FEAS_PATH, recursive = TRUE, showWarnings = FALSE)

lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))
thr <- readRDS(file.path(PATHS$calibration, "thr_primary.rds"))

na_rows <- lca %>% filter(is.na(maf_bin))
n_na <- nrow(na_rows)

message(sprintf("NA-maf_bin rows: %d of %s total (%.2f%%)",
                n_na, format(nrow(lca), big.mark = ","),
                100 * n_na / nrow(lca)))

if (n_na == 0) {
  md_empty <- c(
    sprintf("# NA-maf_bin diagnostic (%s)", format(Sys.time(), "%Y-%m-%d %H:%M")),
    "",
    sprintf("Spec tag: `%s`", TAG$apply_tag),
    sprintf("Threshold: DP>=%g, AD1>=%g, MAF in [%g, %g]",
            thr$DP_min, thr$AD1_min, thr$MAF_min, thr$MAF_max),
    "",
    sprintf("**0 of %s rows have NA maf_bin.**",
            format(nrow(lca), big.mark = ",")),
    "",
    "No diagnostic needed. Likely cause: the `>MAF_max` sentinel bin",
    "(added in build_lca_dataset.R) absorbs alt-allele-majority calls",
    "that previously fell to NA, and Exclusion 1 removes calls with",
    "MAF_max_obs below MAF_min."
  )
  writeLines(md_empty, file.path(FEAS_PATH, "na_maf_bin_diagnostic.md"))
  message("[diagnose_na_maf_bin] No NA-maf_bin rows. Empty diagnostic written.")
} else {

# ---- Q1: Is MAF actually NA, or non-NA but outside the expected range? -----
maf_status <- na_rows %>%
  mutate(
    MAF1_state = case_when(
      is.na(MAF1)            ~ "NA",
      MAF1 < thr$MAF_min     ~ sprintf("<%g", thr$MAF_min),
      MAF1 > thr$MAF_max     ~ sprintf(">%g", thr$MAF_max),
      TRUE                    ~ "in_range"
    ),
    MAF2_state = case_when(
      is.na(MAF2)            ~ "NA",
      MAF2 < thr$MAF_min     ~ sprintf("<%g", thr$MAF_min),
      MAF2 > thr$MAF_max     ~ sprintf(">%g", thr$MAF_max),
      TRUE                    ~ "in_range"
    )
  ) %>%
  count(MAF1_state, MAF2_state, name = "n")

# ---- Q2: How many are disagreement cells? ----------------------------------
W_cells <- na_rows %>% count(W1, W2, name = "n")

# ---- Q3: PE/PPE distribution ------------------------------------------------
region_dist <- na_rows %>% count(region, name = "n")

# ---- Q4: Are AD1/AD2/DP present? -------------------------------------------
field_presence <- tibble(
  field   = c("MAF1","MAF2","DP1","DP2","AD1_1","AD1_2"),
  n_NA    = c(sum(is.na(na_rows$MAF1)),  sum(is.na(na_rows$MAF2)),
              sum(is.na(na_rows$DP1)),   sum(is.na(na_rows$DP2)),
              sum(is.na(na_rows$AD1_1)), sum(is.na(na_rows$AD1_2))),
  n_total = n_na
)

# ---- Q5: Which pairs do they come from? ------------------------------------
pair_dist <- na_rows %>%
  count(pair_id, sampleA, sampleB, name = "n_NA_calls") %>%
  arrange(desc(n_NA_calls))

# ---- Q6: HIV distribution ---------------------------------------------------
hiv_dist <- na_rows %>%
  mutate(HIV_grp = ifelse(is.na(HIV), "HIV unknown", as.character(HIV))) %>%
  count(HIV_grp, name = "n")

# ---- Save full row dump for inspection -------------------------------------
write_csv(na_rows, file.path(FEAS_PATH, "na_maf_bin_rows.csv"))

# ---- Markdown summary -------------------------------------------------------
fmt_tbl <- function(df) paste(capture.output(print(as.data.frame(df))),
                              collapse = "\n")

md <- c(
  sprintf("# NA-maf_bin diagnostic (%s)", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  sprintf("Spec tag: `%s`", TAG$apply_tag),
  sprintf("Threshold: DP>=%g, AD1>=%g, MAF in [%g, %g]",
          thr$DP_min, thr$AD1_min, thr$MAF_min, thr$MAF_max),
  "",
  sprintf("**%d rows of %s have NA maf_bin (%.2f%% of dataset).**",
          n_na, format(nrow(lca), big.mark = ","), 100 * n_na / nrow(lca)),
  "",
  "Rationale: maf_bin is computed as `cut(MAF_max_obs, ...)` over",
  sprintf("[%g, %g, %g, %g, %g]. Rows fail this binning when",
          0.02, 0.05, 0.10, 0.25, thr$MAF_max),
  "MAF_max_obs is NA (both replicates lack a derivable MAF) or falls outside",
  "the bin range. The bins start at 0.02 = MAF_min, so rows where both",
  "replicates have MAF below MAF_min should already be in the [0,0.02) bin",
  "(if that bin was added) -- otherwise they fall to NA.",
  "",
  "## Q1. MAF state in each replicate",
  "",
  "Cross-tab of MAF1 x MAF2 status. `NA` means the replicate did not call",
  "the variant; in_range means MAF in [MAF_min, MAF_max].",
  "",
  "```",
  fmt_tbl(maf_status),
  "```",
  "",
  "## Q2. Joint (W1, W2) cells among NA rows",
  "",
  "If any cell other than (0,0) appears, the row is informative for sens/spec",
  "and the exclusion is non-trivial.",
  "",
  "```",
  fmt_tbl(W_cells),
  "```",
  "",
  "## Q3. Region distribution",
  "",
  "```",
  fmt_tbl(region_dist),
  "```",
  "",
  "## Q4. Field presence (NA counts)",
  "",
  "If DP/AD1 are present but MAF is NA, the issue is upstream MAF derivation",
  "(likely a missing AD2). If everything is NA, the replicate genuinely did",
  "not call.",
  "",
  "```",
  fmt_tbl(field_presence),
  "```",
  "",
  "## Q5. Pair distribution",
  "",
  "Concentration in a small number of pairs suggests a sample-level issue;",
  "spread across many pairs suggests a position-level issue.",
  "",
  "```",
  fmt_tbl(pair_dist),
  "```",
  "",
  "## Q6. HIV distribution",
  "",
  "```",
  fmt_tbl(hiv_dist),
  "```",
  "",
  "## Suggested footnote text",
  "",
  sprintf("> %d variant calls (%.2f%% of the LCA universe) had non-derivable",
          n_na, 100 * n_na / nrow(lca)),
  "> minor allele frequency in both replicates and were excluded from the",
  sprintf("> primary analysis. Of these, %d contributed to disagreement cells.",
          sum(W_cells$n[W_cells$W1 != W_cells$W2])),
  "> Sensitivity analysis retaining these rows under [imputation rule] is",
  "> reported in supplementary materials."
)

writeLines(md, file.path(FEAS_PATH, "na_maf_bin_diagnostic.md"))

# All console output via message() so stdout-eating environments preserve it.
print_via_message <- function(df) {
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}

message("\n--- NA maf_bin diagnostic ---")
message("\nQ1 MAF state cross-tab:")
print_via_message(maf_status)
message("\nQ2 (W1,W2) cells:")
print_via_message(W_cells)
message("\nQ3 Region:")
print_via_message(region_dist)
message("\nQ4 Field presence:")
print_via_message(field_presence)
message("\nQ5 Pair distribution (top rows):")
print_via_message(head(pair_dist, 10))
message("\nQ6 HIV distribution:")
print_via_message(hiv_dist)

message("\n[diagnose_na_maf_bin] Done. Summary written to ",
        file.path(FEAS_PATH, "na_maf_bin_diagnostic.md"))

}  # end else (n_na > 0)
