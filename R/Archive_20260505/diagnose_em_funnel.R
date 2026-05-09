# =============================================================================
# diagnose_em_funnel.R
#
# Trace EM (early-morning) samples through the Ghana cohort filter funnel
# to identify where they are being lost. Reads sample_tbl.rds (pre-filter),
# D2_longitudinal.rds (post-filter), and the variant table (independent
# inventory). Run this after a Paper 1 pipeline run that produced too few
# EM samples in gh_meta.
#
# Required input files in PATHS$data_raw:
#   sample_tbl.rds                            (Ghana cohort, pre-filter inventory)
#   D2_longitudinal.rds                       (Ghana cohort, post-filter)
#   ghana_gatkVariantTablesFormatted.rds      (variant table snapshot)
# Optional:
#   qualifying_samples_filter.rds             (Ghana cohort audit table)
#
# Run: Rscript R/diagnose_em_funnel.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

source("R/00_pipeline_config.R")

# ---- Load -------------------------------------------------------------------
sample_tbl <- readRDS(file.path(PATHS$data_raw, "sample_tbl.rds"))
D2_long    <- readRDS(file.path(PATHS$data_raw, "D2_longitudinal.rds"))
gh_v_raw   <- readRDS(file.path(PATHS$data_raw,
                                "ghana_gatkVariantTablesFormatted.rds"))

cat(sprintf("sample_tbl:        %d rows, %d patients\n",
            nrow(sample_tbl), n_distinct(sample_tbl$patientId)))
cat(sprintf("D2_longitudinal:   %d rows, %d patients\n",
            nrow(D2_long), n_distinct(D2_long$patientId)))
cat(sprintf("variant table:     %d rows, %d distinct Sample IDs\n",
            nrow(gh_v_raw), n_distinct(gh_v_raw$Sample)))

# ---- 1) Pre-filter EM inventory --------------------------------------------
em_inv <- sample_tbl %>% filter(slot_code == "E0.1")
cat(sprintf("\n[1] sample_tbl EM (slot_code='E0.1') rows: %d\n", nrow(em_inv)))
cat("    (expected: ~1 per patient = ~150)\n")

# ---- 2) EM funnel through each D2 filter -----------------------------------
# Replicates the codebook's qualifying_samples filter, applied to EM only.
em_funnel <- em_inv %>%
  mutate(
    f01_collected     = !is.na(sample_id),
    f02_valid_id      = f01_collected & id_status != "invalid",
    f03_variant_pass  = f02_valid_id  & variant_pass %in% TRUE,
    f04_mtbc          = f03_variant_pass & tbp_is_mtbc %in% TRUE,
    f05_coverage_pass = f04_mtbc      & coverage_pass_50x %in% TRUE
  ) %>%
  summarise(
    n_total               = n(),
    after_01_collected    = sum(f01_collected),
    after_02_valid_id     = sum(f02_valid_id),
    after_03_variant_pass = sum(f03_variant_pass),
    after_04_mtbc         = sum(f04_mtbc),
    after_05_coverage     = sum(f05_coverage_pass)
  )
cat("\n[2] EM filter funnel (cumulative; each row = N surviving up to that step):\n")
print(t(as.matrix(em_funnel)))

# ---- 3) Failure-reason breakdown -------------------------------------------
em_fail <- em_inv %>%
  mutate(
    fail_reason = case_when(
      is.na(sample_id)                      ~ "01_not_collected",
      id_status == "invalid"                ~ "02_invalid_id",
      !(variant_pass %in% TRUE)             ~ "03_no_PASS_variants",
      !(tbp_is_mtbc %in% TRUE)              ~ "04_not_mtbc",
      !(coverage_pass_50x %in% TRUE)        ~ "05_low_coverage",
      TRUE                                   ~ "06_passes_all"
    )
  )
cat("\n[3] EM rows by terminal failure reason (or pass):\n")
print(as.data.frame(em_fail %>% count(fail_reason)))

# ---- 4) Variant-table EM inventory (independent of cohort pipeline) ---------
# Match any sample ending in "-E0", "-E0.1", "-E0.2", etc.
em_in_v <- unique(gh_v_raw$Sample[grepl("-E0(\\.[0-9]+)?$", gh_v_raw$Sample)])
em_in_v_pids <- unique(sub("-E0(\\.[0-9]+)?$", "", em_in_v))
cat(sprintf("\n[4] Variant table: %d distinct EM Sample IDs (across %d patients)\n",
            length(em_in_v), length(em_in_v_pids)))
cat("    Sample-ID forms observed:\n")
print(table(sub(".*-E", "-E", em_in_v)))

# ---- 5) Three-way cross-check ----------------------------------------------
em_in_st_pids <- em_inv %>% filter(!is.na(sample_id)) %>% pull(patientId) %>% unique()
em_in_d2_pids <- D2_long %>% filter(slot_code == "E0.1") %>% pull(patientId) %>% unique()

cat("\n[5] Three-way EM-patient set comparison:\n")
cat(sprintf("    sample_tbl EM (collected):  %d patients\n", length(em_in_st_pids)))
cat(sprintf("    variant-table EM:           %d patients\n", length(em_in_v_pids)))
cat(sprintf("    D2 EM:                       %d patients\n", length(em_in_d2_pids)))

cat(sprintf("\n    In sample_tbl but NOT in variant table:  %d patients\n",
            length(setdiff(em_in_st_pids, em_in_v_pids))))
cat(sprintf("    In variant table but NOT in sample_tbl:  %d patients\n",
            length(setdiff(em_in_v_pids, em_in_st_pids))))
cat(sprintf("    In variant table but NOT in D2:          %d patients\n",
            length(setdiff(em_in_v_pids, em_in_d2_pids))))
cat(sprintf("    In D2 but NOT in variant table:          %d patients\n",
            length(setdiff(em_in_d2_pids, em_in_v_pids))))

# ---- 6) Drill into 5 example "lost" patients --------------------------------
lost <- setdiff(em_in_v_pids, em_in_d2_pids) %>% head(5)
if (length(lost) > 0) {
  cat(sprintf("\n[6] Drilldown on 5 patients with EM in variant table but NOT in D2:\n"))
  for (pid in lost) {
    cat(sprintf("\n  --- %s ---\n", pid))
    cat("    sample_tbl row(s) for slot_code='E0.1':\n")
    pt_st <- sample_tbl %>%
      filter(patientId == pid, slot_code == "E0.1") %>%
      select(any_of(c("patientId", "slot_code", "sample_id", "id_status",
                      "variant_pass", "tbp_lineage", "tbp_is_mtbc",
                      "coverage_median", "coverage_pass_50x")))
    print(as.data.frame(pt_st), row.names = FALSE)

    cat("    variant-table EM Sample IDs:\n")
    pt_v <- gh_v_raw %>%
      filter(grepl(paste0("^", pid, "-E"), Sample)) %>%
      distinct(Sample) %>%
      arrange(Sample)
    print(as.data.frame(pt_v), row.names = FALSE)
  }
}

# ---- 7) Optional: read the cohort pipeline's audit table -------------------
audit_path <- file.path(PATHS$data_raw, "qualifying_samples_filter.rds")
if (file.exists(audit_path)) {
  cat("\n[7] Cohort pipeline filter audit (qualifying_samples_filter.rds):\n")
  print(as.data.frame(readRDS(audit_path)), row.names = FALSE)
} else {
  cat("\n[7] qualifying_samples_filter.rds not present in data_raw/.\n")
  cat("    Copy from <ghana_project>/derived_data/v_YYYYMMDD/ for cohort-level\n")
  cat("    funnel counts (not slot-stratified).\n")
}
