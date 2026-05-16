# =============================================================================
# 00_prep_metadata.R  (Paper 1: Ghana-only, D2-driven, 3-candidate Sample)
#
# Build the Ghana analytic metadata frame for Paper 1 from the Ghana cohort
# pipeline's canonical outputs (D2_longitudinal + patient_tbl).
#
# Sample-ID handling (the source of historical pain):
#
#   The Ghana cohort lab convention is one EM (early-morning) sample per
#   visit, so the canonical variant-table Sample suffix is "E0" (no ".1").
#   The cohort spec encodes the slot ontology as "E0.1" in slot_code,
#   but normalize_sample_id() in the cohort pipeline does NOT have a
#   typo-correction rule for E-slots (only S-slots), so sample_id is
#   preserved verbatim from REDCap.
#
#   Three Sample-ID forms can occur for a single EM row depending on where
#   the typo was injected:
#
#     (A) No typo anywhere       sample_id = "KBTH001-E0"
#                                variant_table.Sample = "KBTH001-E0"
#                                -> resolved via Sample_canonical
#
#     (B) Typo at REDCap entry,  sample_id = "KBTH001-E0.1"
#         propagated downstream  variant_table.Sample = "KBTH001-E0.1"
#                                -> resolved via Sample_from_id (sample_id
#                                   matches variant table verbatim)
#
#     (C) Typo introduced after  sample_id = "KBTH001-E0"
#         REDCap (sequencing     variant_table.Sample = "KBTH001-E0.1"
#         pipeline used          -> resolved via Sample_slot_form
#         slot_code as the ID)      (paste(patientId, slot_code) reproduces
#                                   the typo'd form)
#
#   This file evaluates all three candidates per row and takes the first
#   that exists in variant_samples (canonical preferred). Failure to resolve
#   is surfaced explicitly in the per-slot-code crosstab.
#
#   When the variant table contains BOTH canonical and typo'd entries for
#   the same patient (1:many join warning territory), canonical wins and
#   the typo'd duplicate is excluded from gh_variants.
#
#
# Sources:
#   D2_longitudinal.rds   -- sample-unit, n_expected = 97 (pass band 97-115).
#                            Pre-filtered upstream: variant_pass=TRUE,
#                            tbp_is_mtbc=TRUE, coverage_median >= 50.
#   patient_tbl.rds       -- patient-unit, n_expected = 150. Provides
#                            gender (chr; "0"=Male, "1"=Female per cohort
#                            codebook), hiv_positive (lgl), age (int).
#   ghana_gatkVariantTablesFormatted.rds -- variant-call-unit. Sample
#                            column in mixed canonical/typo form.
#
# Outputs (PATHS$meta, PATHS$variants):
#   gh_meta.rds                  -- sample-level analytic frame
#   gh_pairs.rds                 -- M0-only within-visit replicate pairs
#   gh_variants.rds              -- variant table filtered to gh_meta$Sample
#   gh_metadata_provenance.rds   -- audit log
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

# ---- Load raw inputs --------------------------------------------------------
gh_v_raw <- readRDS(file.path(PATHS$data_raw,
                              "ghana_gatkVariantTablesFormatted.rds")) %>%
  dplyr::filter(Sample != "KBTH056-E0")  # excluded sample (D-series provenance)

d2_path <- file.path(PATHS$data_raw, "D2_longitudinal.rds")
if (!file.exists(d2_path)) {
  stop("00_prep_metadata: D2_longitudinal.rds not found at ", d2_path)
}
D2_long <- readRDS(d2_path) %>% droplevels()
message(sprintf("  D2_longitudinal: %d rows, %d distinct patientId",
                nrow(D2_long), dplyr::n_distinct(D2_long$patientId)))

patient_tbl_path <- file.path(PATHS$data_raw, "patient_tbl.rds")
if (!file.exists(patient_tbl_path)) {
  stop("00_prep_metadata: patient_tbl.rds not found at ", patient_tbl_path)
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
       paste(missing_d2, collapse = ", "))
}
required_pt_cols <- c("patientId", "age", "gender", "hiv_positive")
missing_pt <- setdiff(required_pt_cols, names(patient_tbl))
if (length(missing_pt) > 0) {
  stop("patient_tbl missing required columns: ",
       paste(missing_pt, collapse = ", "))
}

# ---- Three-candidate Sample construction ------------------------------------
slot_to_variant_suffix <- function(slot_code) {
  s <- as.character(slot_code)
  dplyr::case_when(
    s == "E0.1"                                  ~ "E0",
    s == "E1.1"                                  ~ "E1",
    s == "E2.1"                                  ~ "E2",
    s %in% c("S0.1","S0.2","S1.1","S2.1","S2.2") ~ s,
    s == "DX"                                     ~ NA_character_,
    TRUE                                          ~ NA_character_
  )
}

variant_samples <- unique(gh_v_raw$Sample)

D2_use <- D2_long %>%
  dplyr::mutate(
    variant_suffix    = slot_to_variant_suffix(slot_code),

    # (A) Canonical: slot_code -> variant-suffix mapping.
    Sample_canonical  = ifelse(!is.na(variant_suffix),
                               paste(patientId, variant_suffix, sep = "-"),
                               NA_character_),

    # (B) Sample_id verbatim: typos introduced at REDCap entry.
    Sample_from_id    = as.character(sample_id),

    # (C) Slot-code form: typos introduced after REDCap (e.g. by the
    # sequencing pipeline using slot_code as the BAM/VCF Sample name).
    Sample_slot_form  = ifelse(!is.na(slot_code) & slot_code != "DX",
                               paste(patientId, slot_code, sep = "-"),
                               NA_character_),

    Sample = dplyr::case_when(
      !is.na(Sample_canonical) & Sample_canonical %in% variant_samples ~ Sample_canonical,
      Sample_from_id %in% variant_samples                              ~ Sample_from_id,
      !is.na(Sample_slot_form) & Sample_slot_form %in% variant_samples ~ Sample_slot_form,
      TRUE                                                              ~ NA_character_
    ),
    Sample_source = dplyr::case_when(
      !is.na(Sample_canonical) & Sample_canonical %in% variant_samples ~ "canonical",
      Sample_from_id %in% variant_samples                              ~ "sample_id_fallback",
      !is.na(Sample_slot_form) & Sample_slot_form %in% variant_samples ~ "slot_form_fallback",
      TRUE                                                              ~ "neither"
    )
  )

# Diagnostic: source breakdown overall
message("\n[Sample resolution] Source breakdown:")
print(as.data.frame(D2_use %>% dplyr::count(Sample_source)), row.names = FALSE)

# Diagnostic: source x slot_code crosstab. Single most informative panel
# for verifying that EM rows are being resolved (and via which candidate).
message("\n[Sample resolution] Source x slot_code crosstab:")
print(as.data.frame(
  D2_use %>%
    dplyr::count(slot_code, Sample_source) %>%
    tidyr::pivot_wider(names_from = Sample_source,
                       values_from = n,
                       values_fill = 0)
), row.names = FALSE)

# Surface fallback rows (typo paths)
fallback_rows <- D2_use %>%
  dplyr::filter(Sample_source %in% c("sample_id_fallback", "slot_form_fallback")) %>%
  dplyr::select(patientId, slot_code, sample_id,
                Sample_canonical, Sample_from_id, Sample_slot_form,
                Sample, Sample_source)
if (nrow(fallback_rows) > 0) {
  message(sprintf(
    "\n[Sample resolution] %d D2 samples resolved via a fallback candidate (typo path):",
    nrow(fallback_rows)))
  print(as.data.frame(fallback_rows), row.names = FALSE)
}

# Surface unresolved rows
neither_rows <- D2_use %>%
  dplyr::filter(Sample_source == "neither") %>%
  dplyr::select(patientId, slot_code, sample_id,
                Sample_canonical, Sample_from_id, Sample_slot_form)
if (nrow(neither_rows) > 0) {
  message(sprintf(
    "\n[Sample resolution] %d D2 samples have NO variant-table match across all 3 candidates:",
    nrow(neither_rows)))
  print(as.data.frame(neither_rows), row.names = FALSE)
  message("  These will be dropped from gh_meta. If this includes EM rows,\n",
          "  the variant table uses a fourth Sample-ID convention not yet\n",
          "  enumerated -- inspect a few unresolved patientIds in gh_v_raw.")
}

D2_use <- D2_use %>% dplyr::filter(!is.na(Sample))

# Inverse direction: variant-table samples no D2 row pointed to.
v_unused <- setdiff(variant_samples, D2_use$Sample)
if (length(v_unused) > 0) {
  message(sprintf(
    "\n  Note: %d variant-table samples are not pointed to by any D2 row\n  (intentional D2 filter exclusions OR canonical-vs-typo duplicate exclusions).",
    length(v_unused)))
}

# ---- Recode demographics from patient_tbl ----------------------------------
recode_gender_ghana <- function(x) {
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

recode_sym_logical <- function(x) {
  out <- dplyr::case_when(
    x == TRUE  ~ "Yes",
    x == FALSE ~ "No",
    TRUE       ~ NA_character_
  )
  factor(out, levels = c("Yes", "No"))
}

pt_demo <- patient_tbl %>%
  dplyr::transmute(
    patientId = as.character(patientId),
    Sex       = recode_gender_ghana(gender),
    HIV       = recode_hiv_logical(hiv_positive),
    Age       = suppressWarnings(as.integer(age)),
    Fever     = recode_sym_logical(sym_fever),
    Hemoptysis = recode_sym_logical(sym_hemoptysis),
    Cough     = recode_sym_logical(sym_cough)
  ) %>%
  dplyr::distinct(patientId, .keep_all = TRUE)

has_mgit_smear <- "mgit_smear" %in% names(D2_use)
if (!has_mgit_smear) {
  message("  Note: D2 has no mgit_smear column; smear will be NA in gh_meta.")
}

# ---- Build gh_meta ----------------------------------------------------------
gh_meta <- D2_use %>%
  dplyr::mutate(
    patientId = as.character(patientId),
    depth_med = suppressWarnings(as.numeric(coverage_median)),
    # Extract major lineage as "L[N]" form (e.g., "L4"). Sub-lineage detail
    # is dropped intentionally; Paper 1 reports lineage at the major level only.
    # tbp_lineage values: "lineage4", "lineage4.1.2", "lineage6", "M.bovis",
    # NA, "lineage_BCG". The regex matches only the canonical lineage[N]
    # prefix; non-MTBc / NA values are coerced to NA -> "unknown" downstream.
    lineage   = factor(ifelse(grepl("^lineage[0-9]+", as.character(tbp_lineage)),
                              sub("^lineage([0-9]+).*", "L\\1",
                                  as.character(tbp_lineage)),
                              NA_character_)),
    lineage = binarize_lineage_l4(lineage),
    smear     = if (has_mgit_smear) factor(.data$mgit_smear)
                else factor(rep(NA_character_, dplyr::n())),
    cohort    = "Ghana",
    run_id    = "GH_single_run"
  ) %>%
  dplyr::left_join(pt_demo, by = "patientId") %>%
  dplyr::select(
    Sample, sample_id, Sample_source, patientId, visit, visit_label, slot_code,
    HIV, Sex, Age, smear, Cough, Fever, Hemoptysis, lineage, depth_med,
    cohort, run_id, dplyr::everything()
  ) %>%
  dplyr::mutate(lineage = forcats::fct_explicit_na(lineage,
                                                    na_level = "unknown"))

# Defensive depth re-check
n_pre_floor <- nrow(gh_meta)
gh_meta <- gh_meta %>%
  dplyr::filter(!is.na(depth_med), depth_med >= MIN_COV_INCLUDE)
n_post_floor <- nrow(gh_meta)
if (n_pre_floor != n_post_floor) {
  message(sprintf(
    "\nWARNING: %d D2 samples failed redundant depth_med >= %dx check.",
    n_pre_floor - n_post_floor, MIN_COV_INCLUDE))
}

n_samples_gh <- nrow(gh_meta)
n_pts_gh     <- dplyr::n_distinct(gh_meta$patientId)
message(sprintf("\ngh_meta built: %d samples across %d patients",
                n_samples_gh, n_pts_gh))

message("\n[gh_meta] Visit composition:")
print(as.data.frame(gh_meta %>% dplyr::count(visit_label)), row.names = FALSE)
message("\n[gh_meta] Slot composition (canonical slot_code):")
print(as.data.frame(gh_meta %>% dplyr::count(slot_code)), row.names = FALSE)

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
                      dplyr::distinct(patientId, HIV, Sex, Age, Cough, Fever, Hemoptysis) %>%
                      dplyr::arrange(patientId)),
        row.names = FALSE)
} else {
  message("  All patients have complete HIV/Sex/Age after patient_tbl join.")
}

# ---- Build gh_pairs (M0-only within-visit replicate pairs) -----------------
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

message("\n[gh_pairs] Slot pair-type composition:")
print(as.data.frame(gh_pairs %>% dplyr::count(slot_a, slot_b)),
      row.names = FALSE)

# ---- Build gh_variants ------------------------------------------------------
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

# ---- Provenance audit -------------------------------------------------------
provenance <- list(
  spec = list(
    M0_ONLY         = M0_ONLY,
    MIN_COV_INCLUDE = MIN_COV_INCLUDE,
    sources         = c("D2_longitudinal.rds", "patient_tbl.rds",
                        "ghana_gatkVariantTablesFormatted.rds"),
    mapping_method  = "three_candidate_canonical_id_slotform (v4)",
    timestamp       = Sys.time()
  ),
  d2_n_samples_input            = nrow(D2_long),
  sample_resolution_by_slot     = D2_long %>%
    dplyr::mutate(
      variant_suffix   = slot_to_variant_suffix(slot_code),
      Sample_canonical = ifelse(!is.na(variant_suffix),
                                paste(patientId, variant_suffix, sep = "-"),
                                NA_character_),
      Sample_from_id   = as.character(sample_id),
      Sample_slot_form = ifelse(!is.na(slot_code) & slot_code != "DX",
                                paste(patientId, slot_code, sep = "-"),
                                NA_character_),
      Sample_source = dplyr::case_when(
        !is.na(Sample_canonical) & Sample_canonical %in% variant_samples ~ "canonical",
        Sample_from_id %in% variant_samples                              ~ "sample_id_fallback",
        !is.na(Sample_slot_form) & Sample_slot_form %in% variant_samples ~ "slot_form_fallback",
        TRUE                                                              ~ "neither"
      )
    ) %>%
    dplyr::count(slot_code, Sample_source) %>%
    as.data.frame(),
  variant_n_samples_input       = length(variant_samples),
  variant_n_unused              = length(v_unused),
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
    Cough   = sum(is.na(gh_meta$Cough)),
    Fever   = sum(is.na(gh_meta$Fever)),
    Hemoptysis = sum(is.na(gh_meta$Hemoptysis)),
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

