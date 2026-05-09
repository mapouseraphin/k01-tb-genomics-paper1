# =============================================================================
# 00_prep_metadata.R  (Paper 1: Ghana-only, D2-driven, dual-candidate Sample)
#
# Build the Ghana analytic metadata frame for Paper 1 from the Ghana cohort
# pipeline's canonical outputs (D2_longitudinal + patient_tbl).
#
# Sample-ID handling notes (the source of historical pain):
#
#   The Ghana cohort lab convention is that EM (early-morning) samples have
#   a single instance per visit, so the canonical variant-table Sample
#   suffix is "E0" (no ".1"). The cohort spec encodes the slot ontology as
#   "E0.1" in slot_code, but normalize_sample_id() in the cohort pipeline
#   does NOT have a typo-correction rule for E-slots (only S-slots), so
#   sample_id is preserved verbatim from REDCap.
#
#   Some lab techs typed "KBTH001-E0.1" instead of "KBTH001-E0", and that
#   typo carries through REDCap into both D2$sample_id and the variant
#   table's Sample column. The variant table therefore contains a mix of
#   "KBTH###-E0" (canonical) and "KBTH###-E0.1" (typo'd) Sample IDs.
#
#   A naive slot_code -> "E0" mapping silently drops the typo'd patients
#   (their variant-table row is "E0.1", not "E0"), and a naive sample_id
#   passthrough silently drops the non-typo'd patients (their D2 sample_id
#   is "E0", but my regex would assume "E0.1"). Neither assumption holds
#   uniformly.
#
#   Resolution (this file): for each D2 row, build TWO candidate Sample
#   values and test them against the variant_samples set in priority order:
#     (1) Canonical  = paste(patientId, slot_to_variant_suffix(slot_code))
#     (2) Sample_from_id = as.character(sample_id)
#   Take whichever appears in variant_samples. If both appear (i.e. the
#   variant table has BOTH "KBTH001-E0" and "KBTH001-E0.1" for the same
#   patient -- the source of the 1:many join warning), canonical wins and
#   the typo'd duplicate is automatically excluded from gh_variants.
#   If neither matches, the D2 row is dropped and surfaced in diagnostics.
#
# 2026-05-01 SCOPE: Florida code paths removed (Paper 1 is Ghana-only).
#
# Sources:
#   D2_longitudinal.rds   -- sample-unit, n_expected = 97 (pass band 97-115).
#                            Pre-filtered upstream: variant_pass=TRUE,
#                            tbp_is_mtbc=TRUE, coverage_median >= 50.
#   patient_tbl.rds       -- patient-unit, n_expected = 150. Provides
#                            gender (chr; "0"=Male, "1"=Female per cohort
#                            codebook), hiv_positive (lgl), age (int).
#   ghana_gatkVariantTablesFormatted.rds -- variant-call-unit. Sample
#                            column is in mixed canonical/typo form.
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

# ---- Dual-candidate Sample construction -------------------------------------
# Canonical mapping: slot_code -> variant-table suffix.
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

# Two-candidate resolution: canonical preferred, sample_id fallback for
# typo'd patients whose variant-table entry preserves the typo.
D2_use <- D2_long %>%
  dplyr::mutate(
    variant_suffix    = slot_to_variant_suffix(slot_code),
    Sample_canonical  = ifelse(!is.na(variant_suffix),
                               paste(patientId, variant_suffix, sep = "-"),
                               NA_character_),
    Sample_from_id    = as.character(sample_id),
    Sample = dplyr::case_when(
      !is.na(Sample_canonical) & Sample_canonical %in% variant_samples ~ Sample_canonical,
      Sample_from_id %in% variant_samples                              ~ Sample_from_id,
      TRUE                                                              ~ NA_character_
    ),
    Sample_source = dplyr::case_when(
      !is.na(Sample_canonical) & Sample_canonical %in% variant_samples ~ "canonical",
      Sample_from_id %in% variant_samples                              ~ "sample_id_fallback",
      TRUE                                                              ~ "neither"
    )
  )

# Diagnostic: how each D2 row was resolved.
message("\n[Sample resolution] Source breakdown:")
print(as.data.frame(D2_use %>% dplyr::count(Sample_source)), row.names = FALSE)

fallback_rows <- D2_use %>%
  dplyr::filter(Sample_source == "sample_id_fallback") %>%
  dplyr::select(patientId, slot_code, sample_id, Sample_canonical, Sample)
if (nrow(fallback_rows) > 0) {
  message(sprintf(
    "\n[Sample resolution] %d D2 samples resolved via sample_id fallback (typo preserved into variant table):",
    nrow(fallback_rows)))
  print(as.data.frame(fallback_rows), row.names = FALSE)
}

neither_rows <- D2_use %>%
  dplyr::filter(Sample_source == "neither") %>%
  dplyr::select(patientId, slot_code, sample_id, Sample_canonical)
if (nrow(neither_rows) > 0) {
  message(sprintf(
    "\n[Sample resolution] %d D2 samples have NO variant-table match (canonical and sample_id both miss):",
    nrow(neither_rows)))
  print(as.data.frame(neither_rows), row.names = FALSE)
  message("  These will be dropped from gh_meta. Likely cause: D2 included a\n",
          "  sample whose variant_pass status differs between the cohort\n",
          "  pipeline's view of the variant table and the snapshot loaded here.")
}

# Surface the duplicate-form patients (canonical AND typo both in variant
# table -- this is the source of the 1:many join warning the user observed).
canonical_in_v <- D2_use$Sample_canonical[!is.na(D2_use$Sample_canonical) &
                                          D2_use$Sample_canonical %in% variant_samples]
sampleid_for_those_pats <- D2_use$Sample_from_id[!is.na(D2_use$Sample_canonical) &
                                                 D2_use$Sample_canonical %in% variant_samples]
typo_duplicates <- sampleid_for_those_pats[
  sampleid_for_those_pats != canonical_in_v &
  sampleid_for_those_pats %in% variant_samples
]
if (length(typo_duplicates) > 0) {
  message(sprintf(
    "\n[Sample resolution] %d patients have BOTH canonical and typo'd Sample IDs in the variant table:",
    length(typo_duplicates)))
  message("  Canonical was kept; typo'd duplicates will be excluded from gh_variants.")
  message("  Typo'd Sample IDs being excluded:")
  print(as.data.frame(D2_use %>%
                      dplyr::filter(Sample_from_id %in% typo_duplicates,
                                    Sample_source == "canonical") %>%
                      dplyr::select(patientId, slot_code,
                                    Sample_canonical, Sample_from_id)),
        row.names = FALSE)
}

# Drop unmappable rows.
D2_use <- D2_use %>% dplyr::filter(!is.na(Sample))

# Inverse direction: variant-table samples that no D2 row pointed to (they
# may be intentionally excluded by D2's filter, OR they may be typo'd
# duplicates being correctly excluded by canonical preference).
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
    xc == "0" ~ "Male",
    xc == "1" ~ "Female",
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

has_mgit_smear <- "mgit_smear" %in% names(D2_use)
if (!has_mgit_smear) {
  message("  Note: D2 has no mgit_smear column; smear will be NA in gh_meta.")
}

# ---- Build gh_meta ----------------------------------------------------------
gh_meta <- D2_use %>%
  dplyr::mutate(
    patientId = as.character(patientId),
    depth_med = suppressWarnings(as.numeric(coverage_median)),
    lineage   = factor(sub("^(lineage[0-9]+).*", "\\1", as.character(tbp_lineage))),
    smear     = if (has_mgit_smear) factor(.data$mgit_smear)
                else factor(rep(NA_character_, dplyr::n())),
    cohort    = "Ghana",
    run_id    = "GH_single_run"
  ) %>%
  dplyr::left_join(pt_demo, by = "patientId") %>%
  dplyr::select(
    Sample, sample_id, Sample_source, patientId, visit, visit_label, slot_code,
    HIV, Sex, Age, smear, lineage, depth_med,
    cohort, run_id, dplyr::everything()
  ) %>%
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
                      dplyr::distinct(patientId, HIV, Sex, Age) %>%
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
# Filter directly to gh_meta$Sample. This naturally excludes any typo'd
# duplicate Sample IDs in the variant table -- when canonical "KBTH001-E0"
# was kept in gh_meta, the typo'd "KBTH001-E0.1" rows in the variant table
# do not appear in gh_meta$Sample and are dropped here.
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
    mapping_method  = "dual_candidate_canonical_then_sample_id (v3)",
    timestamp       = Sys.time()
  ),
  d2_n_samples_input            = nrow(D2_long),
  sample_resolution             = D2_long %>%
    dplyr::mutate(
      variant_suffix   = slot_to_variant_suffix(slot_code),
      Sample_canonical = ifelse(!is.na(variant_suffix),
                                paste(patientId, variant_suffix, sep = "-"),
                                NA_character_),
      Sample_from_id   = as.character(sample_id),
      Sample_source = dplyr::case_when(
        !is.na(Sample_canonical) & Sample_canonical %in% variant_samples ~ "canonical",
        Sample_from_id %in% variant_samples                              ~ "sample_id_fallback",
        TRUE                                                              ~ "neither"
      )
    ) %>%
    dplyr::count(Sample_source) %>%
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
