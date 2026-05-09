# =============================================================================
# feasibility_assessment.R
#
# Produces the tabulations that decide whether stratified LCA models are
# feasible. NO modeling decisions are made here. The output is a set of
# tables for human inspection.
#
# Tables produced:
#   T1  Patients & pairs by HIV status
#   T2  Pair counts by pair-type x HIV
#   T3  LCA-eligible variant calls by HIV x MAF bin
#   T4  Disagreement-cell (W1 != W2) counts by HIV x MAF bin
#   T5  LCA-eligible variant calls by region (PE/PPE) x MAF bin
#   T6  Disagreement-cell counts by region x MAF bin
#   T7  Three-way: variant-call counts by HIV x MAF x pair-type
#   T8  Three-way: disagreement-cell counts by HIV x MAF x pair-type
#   T9  Joint (W1,W2) cell counts overall and by stratum
#
# Inputs:
#   <LCA_PATH>/lca_dataset.rds  (built by R/build_lca_dataset.R)
#
# Outputs:
#   <LCA_PATH>/feasibility/T1_patients_by_HIV.csv
#   <LCA_PATH>/feasibility/T2_pairs_by_pairtype_HIV.csv
#   <LCA_PATH>/feasibility/T3_calls_by_HIV_MAF.csv
#   <LCA_PATH>/feasibility/T4_disagree_by_HIV_MAF.csv
#   <LCA_PATH>/feasibility/T5_calls_by_region_MAF.csv
#   <LCA_PATH>/feasibility/T6_disagree_by_region_MAF.csv
#   <LCA_PATH>/feasibility/T7_calls_by_HIV_MAF_pairtype.csv
#   <LCA_PATH>/feasibility/T8_disagree_by_HIV_MAF_pairtype.csv
#   <LCA_PATH>/feasibility/T9_joint_cells.csv
#   <LCA_PATH>/feasibility/feasibility_tables.rds   (all tables in one list)
#   <LCA_PATH>/feasibility/feasibility_report.md    (printable summary)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble)
})

LCA_PATH <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
FEAS_PATH <- file.path(LCA_PATH, "feasibility")
dir.create(FEAS_PATH, recursive = TRUE, showWarnings = FALSE)

lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))

message(sprintf("Feasibility assessment: %s rows, %d pairs, %d patients",
                format(nrow(lca), big.mark = ","),
                n_distinct(lca$pair_id),
                n_distinct(lca$patientId)))

# Helper: factorise HIV with explicit "HIV unknown" level for tabulation
hiv_tab <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- "HIV unknown"
  factor(out, levels = c("HIV-", "HIV+", "HIV unknown"))
}

# Helper: write + print
# All console output goes through message() (stderr). In some environments
# (Rscript with redirection, HPC jobs, knit/render pipelines), stdout from
# print() is swallowed or truncated; stderr from message() is preserved.
write_table <- function(tbl, name, descr) {
  fn <- file.path(FEAS_PATH, paste0(name, ".csv"))
  write_csv(tbl, fn)
  message("\n--- ", name, ": ", descr, " ---")
  message(sprintf("(rows=%d, cols=%d, cols=[%s])",
                  nrow(tbl), ncol(tbl), paste(names(tbl), collapse = ", ")))
  if (nrow(tbl) > 0) {
    rendered <- capture.output(print(as.data.frame(tbl), row.names = FALSE))
    for (ln in rendered) message(ln)
  } else {
    message("(empty table)")
  }
  invisible(tbl)
}

# ---- T1: patients and pairs by HIV ------------------------------------------
# Constructed with base-R table() to avoid factor-level / pivot edge cases.
hiv_recode <- function(x) {
  out <- as.character(x)
  out[is.na(out)] <- "HIV unknown"
  out
}
all_hiv_levels <- c("HIV-", "HIV+", "HIV unknown")

patient_hiv <- unique(as.data.frame(lca)[, c("patientId", "HIV")])
patient_hiv$HIV_grp <- hiv_recode(patient_hiv$HIV)

pair_hiv <- unique(as.data.frame(lca)[, c("pair_id", "HIV")])
pair_hiv$HIV_grp <- hiv_recode(pair_hiv$HIV)

n_pat  <- as.integer(table(factor(patient_hiv$HIV_grp, levels = all_hiv_levels)))
n_prs  <- as.integer(table(factor(pair_hiv$HIV_grp,    levels = all_hiv_levels)))

T1 <- data.frame(
  HIV_grp    = all_hiv_levels,
  n_patients = n_pat,
  n_pairs    = n_prs,
  stringsAsFactors = FALSE
)
# Drop empty rows (e.g., HIV unknown after exclusion)
T1 <- T1[T1$n_patients > 0 | T1$n_pairs > 0, , drop = FALSE]
rownames(T1) <- NULL
write_table(T1, "T1_patients_by_HIV", "Patients and M0 pairs by HIV status")

# ---- T2: pairs by pair-type x HIV -------------------------------------------
pair_info <- unique(as.data.frame(lca)[, c("pair_id", "pair_type", "HIV")])
pair_info$HIV_grp <- factor(hiv_recode(pair_info$HIV), levels = all_hiv_levels)
pair_info$pair_type <- factor(pair_info$pair_type)

T2_mat <- table(pair_info$pair_type, pair_info$HIV_grp)
T2 <- as.data.frame.matrix(T2_mat)
# Drop HIV columns that are entirely zero (e.g., HIV unknown)
T2 <- T2[, colSums(T2) > 0, drop = FALSE]
T2 <- cbind(pair_type = rownames(T2), T2)
rownames(T2) <- NULL
write_table(T2, "T2_pairs_by_pairtype_HIV", "Pair counts by pair-type and HIV")

# ---- T3: variant calls by HIV x MAF bin -------------------------------------
T3 <- lca %>%
  mutate(HIV_grp = hiv_tab(HIV)) %>%
  count(maf_bin, HIV_grp) %>%
  pivot_wider(names_from = HIV_grp, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(maf_bin)
write_table(T3, "T3_calls_by_HIV_MAF", "LCA-eligible variant calls by HIV x MAF bin")

# ---- T4: disagreement cells by HIV x MAF bin (the IDENTIFYING quantity) ----
T4 <- lca %>%
  filter(W1 != W2) %>%
  mutate(HIV_grp = hiv_tab(HIV)) %>%
  count(maf_bin, HIV_grp) %>%
  pivot_wider(names_from = HIV_grp, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(maf_bin)
write_table(T4, "T4_disagree_by_HIV_MAF",
            "Disagreement-cell (W1 != W2) counts by HIV x MAF bin")

# ---- T5: variant calls by region x MAF bin ----------------------------------
T5 <- lca %>%
  count(maf_bin, region) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(maf_bin)
write_table(T5, "T5_calls_by_region_MAF",
            "LCA-eligible variant calls by region (PE/PPE) x MAF bin")

# ---- T6: disagreement cells by region x MAF bin -----------------------------
T6 <- lca %>%
  filter(W1 != W2) %>%
  count(maf_bin, region) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  arrange(maf_bin)
write_table(T6, "T6_disagree_by_region_MAF",
            "Disagreement-cell counts by region x MAF bin")

# ---- T7: three-way variant calls by HIV x MAF x pair-type -------------------
T7 <- lca %>%
  mutate(HIV_grp = hiv_tab(HIV)) %>%
  count(pair_type, maf_bin, HIV_grp, name = "n_calls") %>%
  arrange(pair_type, maf_bin, HIV_grp)
write_table(T7, "T7_calls_by_HIV_MAF_pairtype",
            "Three-way: variant calls by pair-type x MAF x HIV (long form)")

# ---- T8: three-way disagreement cells by HIV x MAF x pair-type --------------
T8 <- lca %>%
  filter(W1 != W2) %>%
  mutate(HIV_grp = hiv_tab(HIV)) %>%
  count(pair_type, maf_bin, HIV_grp, name = "n_disagree") %>%
  arrange(pair_type, maf_bin, HIV_grp)
write_table(T8, "T8_disagree_by_HIV_MAF_pairtype",
            "Three-way: disagreement cells by pair-type x MAF x HIV (long form)")

# ---- T9: joint (W1,W2) cells overall and by HIV / region --------------------
T9_overall <- lca %>%
  count(W1, W2) %>%
  mutate(prop = round(n / sum(n), 3),
         stratum = "overall", level = "all", .before = 1)

T9_hiv <- lca %>%
  mutate(HIV_grp = hiv_tab(HIV)) %>%
  group_by(HIV_grp) %>%
  count(W1, W2) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup() %>%
  rename(level = HIV_grp) %>%
  mutate(stratum = "HIV", .before = 1) %>%
  mutate(level = as.character(level))

T9_region <- lca %>%
  group_by(region) %>%
  count(W1, W2) %>%
  mutate(prop = round(n / sum(n), 3)) %>%
  ungroup() %>%
  rename(level = region) %>%
  mutate(stratum = "region", .before = 1) %>%
  mutate(level = as.character(level))

T9 <- bind_rows(T9_overall, T9_hiv, T9_region) %>%
  select(stratum, level, W1, W2, n, prop)
write_table(T9, "T9_joint_cells",
            "Joint (W1, W2) cell counts: overall, by HIV, by region")

# ---- Bundle all tables ------------------------------------------------------
feasibility <- list(
  T1 = T1, T2 = T2, T3 = T3, T4 = T4, T5 = T5,
  T6 = T6, T7 = T7, T8 = T8, T9 = T9,
  meta = list(
    n_calls    = nrow(lca),
    n_pairs    = n_distinct(lca$pair_id),
    n_patients = n_distinct(lca$patientId),
    maf_bins   = levels(lca$maf_bin),
    timestamp  = Sys.time(),
    apply_tag  = TAG$apply_tag
  )
)
saveRDS(feasibility, file.path(FEAS_PATH, "feasibility_tables.rds"))

# ---- Markdown summary report (no decision rules applied; PI decides) --------
md <- c(
  sprintf("# LCA Feasibility Assessment (%s)", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  sprintf("Spec tag: `%s`", TAG$apply_tag),
  sprintf("- LCA-eligible variant calls: **%s**", format(nrow(lca), big.mark = ",")),
  sprintf("- M0 within-visit pairs:      **%d**", n_distinct(lca$pair_id)),
  sprintf("- Unique patients:            **%d**", n_distinct(lca$patientId)),
  "",
  "## MAF bin levels",
  "",
  paste("Bins:", paste(levels(lca$maf_bin), collapse = ", ")),
  "",
  "The `>MAF_max` level is a sentinel for alt-allele-majority calls",
  "(MAF in either replicate above the rule's upper bound). W1 and W2 are",
  "0 for these calls by rule construction. They are visible in T3/T5/T7 but",
  "absent from T4/T6/T8 (no disagreement cells), and are dropped from",
  "MAF-stratified models.",
  "",
  "## What to look at first",
  "",
  "1. **T4** -- HIV x MAF disagreement cells. This is the identifying",
  "   quantity for HIV-stratified sens. If HIV+ disagreement counts are",
  "   sparse across MAF bins, M5 (HIV-stratified LCA) is not feasible.",
  "2. **T6** -- region x MAF disagreement cells. Drives M4 feasibility.",
  "3. **T8** -- if any HIV+ x MAF x pair-type cell is empty, the three-way",
  "   stratification is foreclosed; reflect that in any sensitivity analysis.",
  "4. **T9** -- joint cells. If the W1=W2=0 cell dominates by orders of",
  "   magnitude, sens identification is data-poor and posteriors will be wide",
  "   regardless of model.",
  "",
  "## Files",
  paste0("- ", list.files(FEAS_PATH, pattern = "\\.csv$"), collapse = "\n")
)
writeLines(md, file.path(FEAS_PATH, "feasibility_report.md"))

message("\n[feasibility_assessment] Done. Tables written to ", FEAS_PATH, "/")
message("Send T1-T8 (or the feasibility_tables.rds) to discuss model selection.")
