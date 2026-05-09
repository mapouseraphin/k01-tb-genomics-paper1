# =============================================================================
# 00_prep_metadata.R  (Paper 1: Ghana-only, D2-driven)
#
# Build the Ghana analytic metadata frame for Paper 1 from the Ghana cohort
# pipeline's canonical outputs (D2_longitudinal + patient_tbl). The previous
# build relied on ghana_epiData.rds, which had gaps (NA medianCov on 70 E-slot
# samples; missing patientIds; gender encoded inconsistently). Those gaps
# forced a chain of rescues + propagation patches and still produced
# downstream LCA datasets with only ~33 patients. The D2 + patient_tbl
# pathway eliminates all four root causes by construction.
#
# 2026-05-01 SCOPE: Florida code paths removed (Paper 1 is Ghana-only).
#
# Sources (canonical, from the Ghana cohort pipeline):
#   D2_longitudinal.rds   -- sample-unit, n_expected = 97 (pass band 97-115).
#                            Already filtered upstream:
#                              variant_pass == TRUE  (has PASS variants)
#                              tbp_is_mtbc == TRUE   (MTBc lineage call)
#                              coverage_median >= 50 (depth floor)
#                            Has: sample_id, patientId, visit, visit_label,
#                            slot_code, coverage_median, tbp_lineage,
#                            mgit_smear (sample_tbl-inherited).
#
#   patient_tbl.rds       -- patient-unit, n_expected = 150. Has age,
#                            gender (chr; "0"=Male, "1"=Female per cohort
#                            codebook -- raw REDCap codes preserved),
#                            hiv_positive (lgl, derived: hiv_status_curr
#                            wins over hiv_status).
#
#   ghana_gatkVariantTablesFormatted.rds -- variant-call-unit. The Sample
#                            column normalizes E-slots: E0.1->E0, E1.1->E1,
#                            E2.1->E2 (S-slots unchanged). D2 retains the
#                            canonical form (E0.1 etc.). We add a derived
#                            `Sample` column to gh_meta that matches the
#                            variant-table format, and validate the join.
#
# 2026-05-01 IMPORTANT FIX: gender mapping is "0"=Male, "1"=Female per the
# Ghana cohort codebook. The previous Paper 1 script used numeric 1=Male,
# 2=Female, which silently mapped every value to NA (Sex rescue logged
# "165 -> 165 (0 rescued)").
#
# Outputs (PATHS$meta, PATHS$variants):
#   gh_meta.rds                  -- sample-level analytic frame
#                                   (one row per qualifying sample)
#   gh_pairs.rds                 -- M0-only within-visit replicate pairs
#                                   (sampleA, sampleB in variant-table format,
#                                   plus slot_a, slot_b in canonical form)
#   gh_variants.rds              -- variant table filtered to gh_meta$Sample,
#                                   annotated with build_sample_metadata_v2
#                                   columns for backward compatibility
#   gh_metadata_provenance.rds   -- audit log of join counts and NA counts
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

# ---- Load raw inputs --------------------------------------------------------
gh_v_raw <- readRDS(file.path(PATHS$data_raw,
                              "ghana_gatkVariantTablesFormatted.rds")) #%>%
  #dplyr::filter(Sample != "KBTH056-E0")  # excluded sample (D-series provenance)

# Canonical Ghana sample-level analytic frame.
d2_path <- file.path(PATHS$data_raw, "D2_longitudinal.rds")
if (!file.exists(d2_path)) {
  stop("00_prep_metadata: D2_longitudinal.rds not found at ", d2_path, "\n",
       "  Copy from the Ghana cohort project: \n",
       "    cp <ghana_project>/derived_data/v_YYYYMMDD/D2_longitudinal.rds ",
       PATHS$data_raw, "/")
}
D2_long <- readRDS(d2_path) %>% droplevels()
message(sprintf("  D2_longitudinal: %d rows, %d distinct patientId",
                nrow(D2_long), dplyr::n_distinct(D2_long$patientId)))

# Canonical Ghana patient-level table.
patient_tbl_path <- file.path(PATHS$data_raw, "patient_tbl.rds")
if (!file.exists(patient_tbl_path)) {
  stop("00_prep_metadata: patient_tbl.rds not found at ", patient_tbl_path,
       "\n  Copy from the Ghana cohort project.")
}
patient_tbl <- readRDS(patient_tbl_path)
message(sprintf("  patient_tbl: %d rows, %d distinct patientId",
                nrow(patient_tbl), dplyr::n_distinct(patient_tbl$patientId)))

# ---- Validate required columns ---------------------------------------------
required_d2_cols <- c("sample_id", "patientId", "visit", "visit_label",
                      "slot_code", "coverage_median", "tbp_lineage")
missing_d2 <- setdiff(required_d2_cols, names(D2_long))
if (length(missing_d2) > 0) {
  stop("D2_longitudinal missing required columns: ",
       paste(missing_d2, collapse = ", "),
       "\n  Available: ", paste(names(D2_long), collapse = ", "))
}

required_pt_cols <- c("patientId", "age", "gender", "hiv_positive")
missing_pt <- setdiff(required_pt_cols, names(patient_tbl))
if (length(missing_pt) > 0) {
  stop("patient_tbl missing required columns: ",
       paste(missing_pt, collapse = ", "),
       "\n  Available: ", paste(names(patient_tbl), collapse = ", "))
}

# ---- Map D2 sample_id (canonical) -> variant-table Sample format -----------
# build_sample_metadata_v2.R applies this normalization to the variant
# table's Sample column upstream:
#     E0.1 -> E0, E1.1 -> E1, E2.1 -> E2  (S-slots unchanged)
# D2 retains the canonical form. We invert the mapping to produce a Sample
# column that joins cleanly to the variant table.
map_d2_to_variant_sample <- function(sample_id) {
  s <- as.character(sample_id)
  stringr::str_replace(s, "-E([0-9]+)\\.\\d+$", "-E\\1")
}

# Drop DX (diagnostic) slots: not present in the variant table, not used in
# calibration. Should be 0-1 per patient; logged for transparency.
n_dx <- sum(D2_long$slot_code == "DX", na.rm = TRUE)
if (n_dx > 0) {
  message(sprintf("  Dropping %d DX (diagnostic) samples; not used in calibration.",
                  n_dx))
}
D2_use <- D2_long %>%
  dplyr::filter(slot_code != "DX" | is.na(slot_code)) %>%
  dplyr::mutate(Sample = map_d2_to_variant_sample(sample_id))

# Validate: every D2 sample's mapped Sample should exist in the variant
# table (D2 enforces variant_pass=TRUE upstream).
variant_samples <- unique(gh_v_raw$Sample)
d2_not_in_v <- setdiff(D2_use$Sample, variant_samples)
if (length(d2_not_in_v) > 0) {
  message(sprintf(
    "\nWARNING: %d D2 samples have no rows in the variant table after Sample-ID mapping:",
    length(d2_not_in_v)))
  print(as.data.frame(D2_use %>%
                      dplyr::filter(Sample %in% d2_not_in_v) %>%
                      dplyr::select(patientId, sample_id, slot_code, Sample) %>%
                      dplyr::arrange(Sample)),
        row.names = FALSE)
  message("  Possible causes: (a) variant_pass logic upstream differs from\n",
          "  the variant-table Sample inventory; (b) Sample-ID mapping rule\n",
          "  missed an edge case. These samples are DROPPED from gh_meta.")
}
v_not_in_d2 <- setdiff(variant_samples, D2_use$Sample)
if (length(v_not_in_d2) > 0) {
  message(sprintf(
    "  Note: %d variant-table samples are NOT in D2 (intentionally excluded by D2 filter); will be excluded from gh_variants.",
    length(v_not_in_d2)))
}

# ---- Recode demographics from patient_tbl ----------------------------------
# Per Ghana cohort codebook (CODEBOOK.md, patient_tbl section):
#   gender:       chr; "0" = Male, "1" = Female (raw REDCap codes preserved)
#   hiv_positive: lgl; derived (hiv_status_curr wins over hiv_status)
#   age:          int; age in years at registration
recode_gender_ghana <- function(x) {
  # Defensive: handle character or numeric storage
  xc <- as.character(x)
  out <- dplyr::case_when(
    xc == "1" ~ "Male",
    xc == "2" ~ "Female",
    TRUE      ~ NA_character_
  )
  factor(out, levels = c("Female", "Male"))
}
recode_hiv_logical <- function(x) {
  out <- dplyr::case_when(
    x == TRUE  ~ "HIV+",
    x == FALSE ~ "HIV-",
    TRUE       ~ NA_character_
  )
  factor(out, levels = c("HIV-", "HIV+"))
}

pt_demo <- patient_tbl %>%
  dplyr::transmute(
    patientId = as.character(patientId),
    Sex       = recode_gender_ghana(gender),
    HIV       = recode_hiv_logical(hiv_positive),
    Age       = suppressWarnings(as.integer(age))
  ) %>%
  dplyr::distinct(patientId, .keep_all = TRUE)

# Smear: from sample_tbl's mgit_smear (inherited by D2). Optional column.
has_mgit_smear <- "mgit_smear" %in% names(D2_use)
if (!has_mgit_smear) {
  message("  Note: D2 has no mgit_smear column; smear will be NA in gh_meta.")
}

# ---- Build gh_meta ----------------------------------------------------------
gh_meta <- D2_use %>%
  dplyr::filter(Sample %in% variant_samples) %>%
  dplyr::mutate(
    patientId = as.character(patientId),
    depth_med = suppressWarnings(as.numeric(coverage_median)),
    sublineage = factor(tbp_lineage),                                  # full granularity
    lineage    = factor(sub("^(lineage[0-9]+).*", "\\1", tbp_lineage)), # major lineage only
    smear     = if (has_mgit_smear) factor(.data$mgit_smear)
                else factor(rep(NA_character_, dplyr::n())),
    cohort    = "Ghana",
    run_id    = "GH_single_run"
  ) %>%
  dplyr::left_join(pt_demo, by = "patientId") %>%
  dplyr::select(
    Sample, sample_id, patientId, visit, visit_label, slot_code,
    HIV, Sex, Age, smear, lineage, depth_med,
    cohort, run_id, dplyr::everything()
  )

# Mark NA lineage levels explicitly (informative for Table 1).
gh_meta <- gh_meta %>%
  dplyr::mutate(lineage = forcats::fct_explicit_na(lineage,
                                                    na_level = "unknown"))

# Defensive depth re-check (D2 already enforces >=50; this catches drift).
n_pre_floor <- nrow(gh_meta)
gh_meta <- gh_meta %>%
  dplyr::filter(!is.na(depth_med), depth_med >= MIN_COV_INCLUDE)
n_post_floor <- nrow(gh_meta)
if (n_pre_floor != n_post_floor) {
  message(sprintf(
    "\nWARNING: %d D2 samples failed redundant depth_med >= %dx check.",
    n_pre_floor - n_post_floor, MIN_COV_INCLUDE),
    "\n  This should not happen if D2 was built with QC$coverage_median_min == ",
    MIN_COV_INCLUDE, ".")
}

n_samples_gh <- nrow(gh_meta)
n_pts_gh     <- dplyr::n_distinct(gh_meta$patientId)
message(sprintf("\ngh_meta built: %d samples across %d patients",
                n_samples_gh, n_pts_gh))

# Surface patients with residual NA in core demographics.
orphan_patients <- gh_meta %>%
  dplyr::filter(is.na(HIV) | is.na(Sex) | is.na(Age)) %>%
  dplyr::distinct(patientId) %>%
  dplyr::pull(patientId)
if (length(orphan_patients) > 0) {
  message(sprintf(
    "\nNote: %d patients have NA in HIV/Sex/Age after patient_tbl join:",
    length(orphan_patients)))
  print(as.data.frame(gh_meta %>%
                      dplyr::filter(patientId %in% orphan_patients) %>%
                      dplyr::distinct(patientId, HIV, Sex, Age) %>%
                      dplyr::arrange(patientId)),
        row.names = FALSE)
  message("  These patients are NOT excluded; will appear with NA in Table 1.")
} else {
  message("  All patients have complete HIV/Sex/Age after patient_tbl join.")
}

# Visit and slot composition diagnostics
message("\n[gh_meta] Visit composition:")
print(as.data.frame(gh_meta %>% dplyr::count(visit_label)), row.names = FALSE)
message("\n[gh_meta] Slot composition (canonical slot_code):")
print(as.data.frame(gh_meta %>% dplyr::count(slot_code)), row.names = FALSE)

# ---- Build gh_pairs (M0-only within-visit replicate pairs) -----------------
# Pair construction is now driven by D2-derived metadata, not by parsing the
# variant table. Pairs are enumerated within (patientId, visit_label) for
# patients with >= 2 qualifying samples. Slot pair types are read directly
# from canonical slot_code; expect EM(E0.1)<->S0.1, EM<->S0.2, S0.1<->S0.2.

build_pairs_within_visit <- function(meta_df) {
  meta_df %>%
    dplyr::group_by(patientId, visit_label) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::group_modify(~ {
      cs <- utils::combn(seq_len(nrow(.x)), 2)
      tibble::tibble(
        sampleA = .x$Sample[cs[1, ]],
        sampleB = .x$Sample[cs[2, ]],
        slot_a  = .x$slot_code[cs[1, ]],
        slot_b  = .x$slot_code[cs[2, ]]
      )
    }) %>%
    dplyr::ungroup()
}

if (M0_ONLY) {
  m0_meta <- gh_meta %>% dplyr::filter(visit_label == "M0")
  message(sprintf("\nM0 samples: %d across %d patients",
                  nrow(m0_meta), dplyr::n_distinct(m0_meta$patientId)))

  gh_pairs <- build_pairs_within_visit(m0_meta) %>%
    dplyr::mutate(pair_class = "within_visit_replicate") %>%
    dplyr::select(sampleA, sampleB, pair_class, patientId, visit_label,
                  slot_a, slot_b)
  message(sprintf("Calibration pairs (M0 within-visit replicate): %d across %d patients",
                  nrow(gh_pairs), dplyr::n_distinct(gh_pairs$patientId)))
} else {
  gh_pairs <- build_pairs_within_visit(gh_meta) %>%
    dplyr::mutate(pair_class = "within_visit_replicate") %>%
    dplyr::select(sampleA, sampleB, pair_class, patientId, visit_label,
                  slot_a, slot_b)
  message(sprintf("\nM0_ONLY = FALSE; retaining all within-visit pairs: %d (sensitivity).",
                  nrow(gh_pairs)))
}

# Slot pair-type composition diagnostic
message("\n[gh_pairs] Slot pair-type composition:")
print(as.data.frame(gh_pairs %>% dplyr::count(slot_a, slot_b)),
      row.names = FALSE)

# ---- Build gh_variants ------------------------------------------------------
# Re-use build_sample_metadata_v2.R's parser to attach legacy annotation
# columns (person_id, month_num, sample_type, etc.) for backward
# compatibility with downstream scripts. Filter to samples in gh_meta.
source("R/build_sample_metadata_v2.R")

gh_variants_full <- build_sample_metadata(
  variants        = gh_v_raw,
  sample_col      = "Sample",
  enforce_pattern = FALSE,
  verbose         = FALSE
)
gh_variants_annot <- gh_variants_full$variants_annot %>%
  dplyr::filter(Sample %in% gh_meta$Sample)

message(sprintf("\ngh_variants: %d rows -> %d after filtering to gh_meta samples",
                nrow(gh_v_raw), nrow(gh_variants_annot)))

# ---- Build provenance audit -------------------------------------------------
provenance <- list(
  spec = list(
    M0_ONLY         = M0_ONLY,
    MIN_COV_INCLUDE = MIN_COV_INCLUDE,
    sources         = c("D2_longitudinal.rds", "patient_tbl.rds",
                        "ghana_gatkVariantTablesFormatted.rds"),
    timestamp       = Sys.time()
  ),
  d2_n_samples_input            = nrow(D2_long),
  d2_n_dx_dropped               = n_dx,
  d2_n_unmappable_to_variant    = length(d2_not_in_v),
  d2_unmappable_samples         = d2_not_in_v,
  variant_n_samples_input       = length(variant_samples),
  variant_n_excluded_by_d2      = length(v_not_in_d2),
  patient_tbl_n                 = nrow(patient_tbl),
  orphan_patients_post_join     = orphan_patients,
  gh_meta_n_samples             = n_samples_gh,
  gh_meta_n_patients            = n_pts_gh,
  gh_pairs_n                    = nrow(gh_pairs),
  gh_pairs_n_patients           = dplyr::n_distinct(gh_pairs$patientId),
  na_counts_in_gh_meta          = list(
    HIV     = sum(is.na(gh_meta$HIV)),
    Sex     = sum(is.na(gh_meta$Sex)),
    Age     = sum(is.na(gh_meta$Age)),
    smear   = sum(is.na(gh_meta$smear)),
    lineage = sum(as.character(gh_meta$lineage) == "unknown")
  )
)

# ---- Write artifacts --------------------------------------------------------
saveRDS(gh_meta,           file.path(PATHS$meta,     "gh_meta.rds"))
saveRDS(gh_pairs,          file.path(PATHS$meta,     "gh_pairs.rds"))
saveRDS(gh_variants_annot, file.path(PATHS$variants, "gh_variants.rds"))
saveRDS(provenance,        file.path(PATHS$meta,     "gh_metadata_provenance.rds"))

message("\n[00_prep_metadata] Done. Ghana-only artifacts written.")
message(sprintf("  gh_meta:     %d samples, %d patients",
                provenance$gh_meta_n_samples, provenance$gh_meta_n_patients))
message(sprintf("  gh_pairs:    %d pairs across %d patients",
                provenance$gh_pairs_n, provenance$gh_pairs_n_patients))
message(sprintf("  gh_variants: %d rows", nrow(gh_variants_annot)))
message("  Provenance audit: ", PATHS$meta, "/gh_metadata_provenance.rds")
