# =============================================================================
# 00_prep_metadata.R
#
# Read raw cohort metadata and variant tables; produce cleaned analytic
# metadata frames for Ghana and Florida. Writes RDS artifacts consumed by
# downstream scripts.
#
# Canonical specification (see 00_pipeline_config.R):
#   - Calibration set: M0 within-visit replicate pairs only
#     Slot pair types retained: EM<->S0.1, EM<->S0.2, S0.1<->S0.2
#     Post-treatment (M1, M2) within-visit pairs excluded because the
#     measurement-error framework requires the latent T to be stable across
#     samples in a calibration pair; treatment-induced selection over weeks
#     between visits violates this.
#
# 2026-04-30 PATCH: depth_med is now derived from a coalesce of:
#   (1) medianCov from ghana_epiData.rds (preferred when populated)
#   (2) call-level median DP from PASS variants in the variant table (fallback
#       when medianCov is NA)
# Motivation: 70 E0-slot samples in the Ghana cohort have NA medianCov in
# ghana_epiData.rds despite carrying PASS variant calls at full sequencing
# depth (median 293x by call-level DP). The previous code dropped these
# samples at line 88's filter(!is.na(depth_med)), inadvertently excluding
# 124 of 178 candidate within-visit pairs from the LCA. The fallback
# rescues these samples without relaxing the >= MIN_COV_INCLUDE criterion;
# samples that fail depth on EITHER source are still dropped.
#
# A new column depth_med_source records provenance for transparency:
#   "metadata"       = depth_med taken from medianCov (preferred path)
#   "calls_fallback" = depth_med computed from call-level median DP because
#                      medianCov was NA
#   "neither"        = both medianCov and call-level DP were unavailable
#                      (sample is dropped by the filter)
#
# Inputs:
#   data_raw/variant_tables_raw/ghana_gatkVariantTablesFormatted.rds
#   data_raw/variant_tables_raw/florida_gatkVariantTablesFormatted.rds
#   data_raw/ghana_epidata.rds
#   data_raw/florida_epiData.rds
#
# Outputs (PATHS$meta, PATHS$variants):
#   gh_meta.rds, fl_meta.rds        -- analytic metadata, depth_med >= 50
#   gh_variants.rds, fl_variants.rds -- raw variant tables (pass-through)
#   gh_pairs.rds                    -- M0-only within-visit replicate pairs
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

# ---- Load raw inputs --------------------------------------------------------
gh_v_raw <- readRDS(file.path(PATHS$data_raw,
                              "ghana_gatkVariantTablesFormatted.rds")) %>%
  dplyr::filter(Sample != "KBTH056-E0")  # excluded sample (D-series provenance)

fl_v_raw <- readRDS(file.path(PATHS$data_raw,
                              "florida_gatkVariantTablesFormatted.rds"))

gh_epi <- readRDS(file.path(PATHS$data_raw, "ghana_epiData.rds")) %>%
  droplevels()
D2_long <- readRDS(file.path(PATHS$data_raw, "D2_longitudinal.rds")) %>%
  droplevels()
fl_epi <- readRDS(file.path(PATHS$data_raw, "florida_epiData.rds")) %>%
  droplevels()

# ---- 2026-04-30 PATCH: load Ghana patient_tbl for demographic rescue --------
# patient_tbl.rds is the canonical Ghana cohort patient-level table produced
# by the Ghana cohort pipeline (derived_data/v_YYYYMMDD/patient_tbl.rds).
# Used as a fallback source for demographic fields (HIV, Sex, Age, smear)
# for samples whose patientId is missing from ghana_epiData.rds.
patient_tbl_path <- file.path(PATHS$data_raw, "patient_tbl.rds")
if (!file.exists(patient_tbl_path)) {
  stop(
    "00_prep_metadata: patient_tbl.rds not found at ", patient_tbl_path, "\n",
    "  Copy it from the Ghana cohort project: \n",
    "    cp <ghana_project>/derived_data/v_20260501/patient_tbl.rds ",
    PATHS$data_raw, "/\n",
    "  This file is required for demographic rescue of the 70 E0 samples ",
    "missing from ghana_epiData.rds."
  )
}
patient_tbl <- readRDS(patient_tbl_path)
message(sprintf("  patient_tbl: %d rows, %d distinct patientId",
                nrow(patient_tbl),
                if ("patientId" %in% names(patient_tbl))
                  dplyr::n_distinct(patient_tbl$patientId) else NA))

# ---- Ghana: parse sample metadata and within-visit pairs --------------------
source("R/build_sample_metadata_v2.R")

gh_variants <- build_sample_metadata(
  variants        = gh_v_raw,
  sample_col      = "Sample",
  enforce_pattern = FALSE,
  verbose         = TRUE
)

gh_pairs_all <- gh_variants$pairs_same_visit %>%
  mutate(pair_class = "within_visit_replicate") %>%
  select(sampleA, sampleB, pair_class)

message(sprintf("Ghana within-visit pairs (all visits): %d", nrow(gh_pairs_all)))

# ---- Restrict pairs to M0 (canonical) ---------------------------------------
if (M0_ONLY) {
  n_pre <- nrow(gh_pairs_all)
  gh_pairs <- gh_pairs_all %>%
    filter(is_m0_sample(sampleA) & is_m0_sample(sampleB))
  n_post <- nrow(gh_pairs)
  message(sprintf(
    "Calibration pairs restricted to M0 (canonical): %d -> %d (dropped %d post-treatment)",
    n_pre, n_post, n_pre - n_post
  ))

  # Slot composition diagnostic
  slot_breakdown <- gh_pairs %>%
    mutate(slot_a = parse_slot(sampleA),
           slot_b = parse_slot(sampleB)) %>%
    count(slot_a, slot_b)
  message("M0 pair composition by slot:")
  print(slot_breakdown)
} else {
  gh_pairs <- gh_pairs_all
  message("M0_ONLY = FALSE; retaining all within-visit pairs (sensitivity).")
}

# ---- Ghana metadata ---------------------------------------------------------
# Build pre-filter metadata first (without the depth filter), then apply the
# depth-source coalesce, then apply the filter. This is the 2026-04-30 patch.
gh_meta_pre <- prep_metadata_ghana_from_epi(gh_epi, gh_variants$meta) %>%
  mutate(lineage = forcats::fct_explicit_na(lineage, na_level = "unknown"))

if (!"person_id" %in% names(gh_meta_pre) && !"patientId" %in% names(gh_meta_pre)) {
  stop("00_prep_metadata: gh_meta missing patientId/person_id. ",
       "Cluster-robust SE on patientId requires this column.")
}
if ("person_id" %in% names(gh_meta_pre) && !"patientId" %in% names(gh_meta_pre)) {
  gh_meta_pre <- gh_meta_pre %>% rename(patientId = person_id)
}

# ---- 2026-04-30 PATCH: call-level depth fallback ---------------------------
# Compute median DP per sample directly from PASS variant calls. This is a
# defensible proxy for sample-level sequencing depth and is available for
# every sample in the variant table (so it covers samples whose medianCov
# field is NA in ghana_epiData.rds).
message("\n[depth_med rescue] Computing call-level median DP from PASS variants...")
depth_from_calls <- gh_v_raw %>%
  { if ("FILTER" %in% names(.)) dplyr::filter(., FILTER == "PASS") else . } %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    n_PASS_calls   = dplyr::n(),
    depth_med_calls = median(DP, na.rm = TRUE),
    .groups = "drop"
  )
message(sprintf("  Computed call-level depth for %d distinct samples.",
                nrow(depth_from_calls)))

# Coalesce depth_med: prefer medianCov-derived value where present, fall back
# to call-level DP where medianCov was NA. Track provenance.
gh_meta_pre <- gh_meta_pre %>%
  dplyr::left_join(depth_from_calls, by = "Sample") %>%
  dplyr::mutate(
    depth_med_metadata = depth_med,    # preserve the original medianCov value
    depth_med = dplyr::coalesce(depth_med, depth_med_calls),
    depth_med_source = dplyr::case_when(
      !is.na(depth_med_metadata)                          ~ "metadata",
      is.na(depth_med_metadata) & !is.na(depth_med_calls) ~ "calls_fallback",
      TRUE                                                 ~ "neither"
    )
  )

source_summary <- gh_meta_pre %>% dplyr::count(depth_med_source)
message("\n[depth_med rescue] Provenance summary:")
print(as.data.frame(source_summary), row.names = FALSE)

n_rescued <- sum(gh_meta_pre$depth_med_source == "calls_fallback")
if (n_rescued > 0) {
  message(sprintf("\n[depth_med rescue] %d samples rescued via call-level DP.",
                  n_rescued))
  rescued <- gh_meta_pre %>%
    dplyr::filter(depth_med_source == "calls_fallback") %>%
    dplyr::select(Sample, patientId, depth_med, n_PASS_calls) %>%
    dplyr::arrange(depth_med)
  message("[depth_med rescue] Rescued samples (sorted by depth_med):")
  print(as.data.frame(rescued), row.names = FALSE)
}

# Now apply the depth filter. Samples with depth_med NA from BOTH sources, OR
# depth_med < MIN_COV_INCLUDE under whichever source is available, are dropped.
n_pre_filter <- nrow(gh_meta_pre)
gh_meta <- gh_meta_pre %>%
  dplyr::filter(!is.na(depth_med), depth_med >= MIN_COV_INCLUDE)
n_post_filter <- nrow(gh_meta)

message(sprintf(
  "\nGhana metadata rows retained (depth_med >= %dx after rescue): %d / %d (%d patients)",
  MIN_COV_INCLUDE, n_post_filter, n_pre_filter, dplyr::n_distinct(gh_meta$patientId)
))

# Diagnostic: how many of the previously-dropped samples are now retained?
n_kept_via_rescue <- gh_meta %>%
  dplyr::filter(depth_med_source == "calls_fallback") %>%
  nrow()
message(sprintf(
  "  Of those, %d entered the analytic frame via the call-level depth rescue.",
  n_kept_via_rescue
))

# Diagnostic: any samples that the rescue could NOT save (depth_med < floor
# under call-level DP)
dropped_post_rescue <- gh_meta_pre %>%
  dplyr::filter(depth_med_source == "calls_fallback",
                !is.na(depth_med),
                depth_med < MIN_COV_INCLUDE)
if (nrow(dropped_post_rescue) > 0) {
  message(sprintf(
    "  %d call-level-rescued samples STILL fall below %dx and are dropped:",
    nrow(dropped_post_rescue), MIN_COV_INCLUDE))
  print(as.data.frame(dropped_post_rescue %>%
                      dplyr::select(Sample, depth_med, n_PASS_calls) %>%
                      dplyr::arrange(depth_med)),
        row.names = FALSE)
}

# ---- 2026-04-30 PATCH: demographic rescue from patient_tbl.rds -------------
# Fields HIV, Sex, Age, smear are clinical records that cannot be
# reconstructed from variant calls. For samples whose patientId is missing
# from ghana_epiData.rds, fall back to patient_tbl.rds (canonical Ghana
# cohort patient-level table) for these fields. Provenance is tracked
# per-field; coalesce is per-field (a patient may have HIV from epiData but
# Sex from patient_tbl, etc., though typically all-or-nothing).
message("\n[demo rescue] Joining patient_tbl for demographic rescue...")

# Helper: pick the first non-null candidate column from a data frame
pick_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(NULL)
  hit[1]
}

# Identify column names in patient_tbl. The Ghana cohort project stores HIV
# in `hiv_status_curr_raw` (NOT hiv_positive, NOT hiv_status); sex is in
# `gender` as numeric 1=Male, 2=Female.
pid_col   <- pick_col(patient_tbl, c("patientId", "person_id", "PatientID"))
hiv_col   <- pick_col(patient_tbl, c("hiv_status_curr_raw"))   # CANONICAL: hiv_status_curr_raw only
sex_col   <- pick_col(patient_tbl, c("gender", "sex", "Sex"))  # CANONICAL: gender (numeric 1/2)
age_col   <- pick_col(patient_tbl, c("age", "Age", "age_at_enrollment", "enrollment_age"))
smear_col <- pick_col(patient_tbl, c("smear_grade", "sputumSmear", "sputum_smear", "smear"))

if (is.null(pid_col)) {
  stop("[demo rescue] patient_tbl has no recognizable patientId column. ",
       "Tried: patientId, person_id, PatientID. Got: ",
       paste(names(patient_tbl), collapse = ", "))
}
if (is.null(hiv_col)) {
  stop("[demo rescue] patient_tbl missing canonical HIV column 'hiv_status_curr_raw'. ",
       "Available columns: ", paste(names(patient_tbl), collapse = ", "))
}

# Numeric-aware sex recoder for the Ghana cohort `gender` encoding (1=Male,
# 2=Female). Falls back to the standard recode_sex() for character inputs
# (so this function is safe to use regardless of what column got picked).
recode_sex_ghana <- function(x) {
  if (is.numeric(x)) {
    out <- dplyr::case_when(
      x == 1L ~ "Male",
      x == 2L ~ "Female",
      TRUE     ~ NA_character_
    )
    return(factor(out, levels = c("Female", "Male")))
  }
  recode_sex(x)
}

message(sprintf("[demo rescue] patient_tbl columns identified: pid=%s, hiv=%s, sex=%s, age=%s, smear=%s",
                pid_col,
                hiv_col   %||% "<missing>",
                sex_col   %||% "<missing>",
                age_col   %||% "<missing>",
                smear_col %||% "<missing>"))

# Build a one-row-per-patient lookup from patient_tbl with standardized names
ptl <- patient_tbl %>%
  dplyr::rename(patientId_lookup = !!pid_col) %>%
  dplyr::mutate(
    HIV_pt   = if (!is.null(hiv_col))   recode_hiv(.data[[hiv_col]])                    else NA,
    Sex_pt   = if (!is.null(sex_col))   recode_sex_ghana(.data[[sex_col]])              else NA,
    Age_pt   = if (!is.null(age_col))   suppressWarnings(as.numeric(.data[[age_col]]))  else NA_real_,
    smear_pt = if (!is.null(smear_col)) as.character(.data[[smear_col]])                else NA_character_
  ) %>%
  dplyr::select(patientId_lookup, HIV_pt, Sex_pt, Age_pt, smear_pt) %>%
  dplyr::distinct(patientId_lookup, .keep_all = TRUE)

# Snapshot pre-rescue NA counts
n_pre_HIV   <- sum(is.na(gh_meta_pre$HIV))
n_pre_Sex   <- sum(is.na(gh_meta_pre$Sex))
n_pre_Age   <- sum(is.na(gh_meta_pre$Age))
n_pre_smear <- sum(is.na(gh_meta_pre$smear))

# Coalesce per field
gh_meta_pre <- gh_meta_pre %>%
  dplyr::left_join(ptl, by = c("patientId" = "patientId_lookup")) %>%
  dplyr::mutate(
    HIV_metadata   = HIV,
    Sex_metadata   = Sex,
    Age_metadata   = Age,
    smear_metadata = smear,
    HIV   = dplyr::coalesce(HIV,   HIV_pt),
    Sex   = dplyr::coalesce(Sex,   Sex_pt),
    Age   = dplyr::coalesce(Age,   Age_pt),
    smear = dplyr::coalesce(as.character(smear), smear_pt),
    smear = factor(smear),
    # Provenance per field: was the value from epiData, patient_tbl, or neither?
    HIV_source = dplyr::case_when(
      !is.na(HIV_metadata) ~ "epiData",
      !is.na(HIV_pt)       ~ "patient_tbl",
      TRUE                  ~ "none"
    ),
    Sex_source = dplyr::case_when(
      !is.na(Sex_metadata) ~ "epiData",
      !is.na(Sex_pt)       ~ "patient_tbl",
      TRUE                  ~ "none"
    ),
    Age_source = dplyr::case_when(
      !is.na(Age_metadata) ~ "epiData",
      !is.na(Age_pt)       ~ "patient_tbl",
      TRUE                  ~ "none"
    ),
    smear_source = dplyr::case_when(
      !is.na(smear_metadata) ~ "epiData",
      !is.na(smear_pt)       ~ "patient_tbl",
      TRUE                    ~ "none"
    )
  )

# Snapshot post-rescue NA counts
n_post_HIV   <- sum(is.na(gh_meta_pre$HIV))
n_post_Sex   <- sum(is.na(gh_meta_pre$Sex))
n_post_Age   <- sum(is.na(gh_meta_pre$Age))
n_post_smear <- sum(is.na(gh_meta_pre$smear))

message("\n[demo rescue] Rescue summary (rows with NA before -> after):")
message(sprintf("  HIV:   %d -> %d  (%d rescued)",
                n_pre_HIV,   n_post_HIV,   n_pre_HIV   - n_post_HIV))
message(sprintf("  Sex:   %d -> %d  (%d rescued)",
                n_pre_Sex,   n_post_Sex,   n_pre_Sex   - n_post_Sex))
message(sprintf("  Age:   %d -> %d  (%d rescued)",
                n_pre_Age,   n_post_Age,   n_pre_Age   - n_post_Age))
message(sprintf("  Smear: %d -> %d  (%d rescued)",
                n_pre_smear, n_post_smear, n_pre_smear - n_post_smear))

# Per-field provenance summary
for (f in c("HIV", "Sex", "Age", "smear")) {
  src_col <- paste0(f, "_source")
  src_summary <- gh_meta_pre %>%
    dplyr::count(.data[[src_col]]) %>%
    dplyr::rename(source = !!src_col)
  message(sprintf("\n[demo rescue] %s provenance:", f))
  print(as.data.frame(src_summary), row.names = FALSE)
}

# Surface any samples STILL missing demographics after rescue
still_missing <- gh_meta_pre %>%
  dplyr::filter(is.na(HIV) | is.na(Sex) | is.na(Age)) %>%
  dplyr::distinct(patientId, HIV_source, Sex_source, Age_source)
if (nrow(still_missing) > 0) {
  message(sprintf(
    "\n[demo rescue] WARNING: %d patients STILL have NA in HIV/Sex/Age after both sources:",
    nrow(still_missing)))
  print(as.data.frame(still_missing %>% dplyr::arrange(patientId)), row.names = FALSE)
  message("[demo rescue] These patients are NOT excluded; they will appear with NA ",
          "in Table 1 cohort summaries.")
}

# ---- 2026-04-30 PATCH (cont.): patient-level metadata propagation ----------
# The same upstream metadata gap that left medianCov NA also left HIV, Sex,
# Age, and lineage NA for some rescued samples. These attributes are
# patient-level constants (one HIV status per patient, one sex, one age at
# enrollment, one consensus lineage), so we can propagate non-NA values
# within patientId. This makes gh_meta internally consistent at the patient
# level and prevents downstream summaries (e.g. Table 1's first(HIV)) from
# spuriously counting the patient as Unknown when one sample's epiData row
# is incomplete.
#
# Propagation rules:
#   HIV, Sex, Age:  fill NA within patient using first non-NA value.
#   lineage:        fill NA or "unknown" within patient using first non-NA
#                   non-"unknown" value.
# For patients where ALL samples have NA / unknown for an attribute, the
# attribute remains NA / unknown (residual truly-missing).
#
# Discordance check: HIV and Sex should be patient-constants. If discordant
# non-NA values are found within a patient, we leave the original values in
# place and emit a warning. This is an upstream data-quality issue, not a
# rescue concern.

fill_within_patient <- function(x, treat_as_na = NULL) {
  if (length(x) == 0) return(x)
  is_missing <- is.na(x)
  if (!is.null(treat_as_na)) {
    is_missing <- is_missing | (as.character(x) %in% as.character(treat_as_na))
  }
  if (all(is_missing)) return(x)
  donor <- x[!is_missing][1]
  out <- x
  out[is_missing] <- donor
  out
}

check_discordance <- function(df, var, var_name) {
  disc <- df %>%
    dplyr::filter(!is.na(.data[[var]])) %>%
    dplyr::group_by(patientId) %>%
    dplyr::summarise(n_distinct_vals = dplyr::n_distinct(.data[[var]]),
                     vals = paste(sort(unique(as.character(.data[[var]]))),
                                  collapse = "|"),
                     .groups = "drop") %>%
    dplyr::filter(n_distinct_vals > 1)
  if (nrow(disc) > 0) {
    message(sprintf("\nWARNING: %d patient(s) have discordant non-NA %s values:",
                    nrow(disc), var_name))
    print(as.data.frame(disc), row.names = FALSE)
    message(sprintf("  Discordant patients are NOT propagated for %s; ",
                    var_name),
            "original sample-level values retained.")
  }
  disc$patientId
}

message("\n[patient-level propagation] Checking discordance...")
discordant_HIV <- check_discordance(gh_meta_pre, "HIV", "HIV")
discordant_Sex <- check_discordance(gh_meta_pre, "Sex", "Sex")

# Preserve originals for diagnostic
gh_meta_pre <- gh_meta_pre %>%
  dplyr::mutate(
    HIV_orig     = HIV,
    Sex_orig     = Sex,
    Age_orig     = Age,
    lineage_orig = lineage
  ) %>%
  dplyr::group_by(patientId) %>%
  dplyr::mutate(
    HIV     = if (patientId[1] %in% discordant_HIV) HIV
              else fill_within_patient(HIV),
    Sex     = if (patientId[1] %in% discordant_Sex) Sex
              else fill_within_patient(Sex),
    Age     = fill_within_patient(Age),
    lineage = fill_within_patient(lineage, treat_as_na = "unknown")
  ) %>%
  dplyr::ungroup()

# Diagnostics
n_filled_HIV <- sum(is.na(gh_meta_pre$HIV_orig)     & !is.na(gh_meta_pre$HIV))
n_filled_Sex <- sum(is.na(gh_meta_pre$Sex_orig)     & !is.na(gh_meta_pre$Sex))
n_filled_Age <- sum(is.na(gh_meta_pre$Age_orig)     & !is.na(gh_meta_pre$Age))
n_filled_lin <- sum((is.na(gh_meta_pre$lineage_orig) |
                     as.character(gh_meta_pre$lineage_orig) == "unknown") &
                    !is.na(gh_meta_pre$lineage) &
                    as.character(gh_meta_pre$lineage) != "unknown")

message("\n[patient-level propagation] Sample-row values filled from same patient:")
message(sprintf("  HIV:     %d samples", n_filled_HIV))
message(sprintf("  Sex:     %d samples", n_filled_Sex))
message(sprintf("  Age:     %d samples", n_filled_Age))
message(sprintf("  lineage: %d samples (was NA or 'unknown')", n_filled_lin))

# Patient-level residual NA counts
patient_summary <- gh_meta_pre %>%
  dplyr::group_by(patientId) %>%
  dplyr::summarise(
    HIV_known     = any(!is.na(HIV)),
    Sex_known     = any(!is.na(Sex)),
    Age_known     = any(!is.na(Age)),
    lineage_known = any(!is.na(lineage) & as.character(lineage) != "unknown"),
    .groups = "drop"
  )
n_pt_total <- nrow(patient_summary)
message(sprintf(
  "\n[patient-level propagation] Patient-level coverage after propagation (n=%d):",
  n_pt_total))
message(sprintf("  HIV known:     %d (%.1f%%)",
                sum(patient_summary$HIV_known),
                100 * mean(patient_summary$HIV_known)))
message(sprintf("  Sex known:     %d (%.1f%%)",
                sum(patient_summary$Sex_known),
                100 * mean(patient_summary$Sex_known)))
message(sprintf("  Age known:     %d (%.1f%%)",
                sum(patient_summary$Age_known),
                100 * mean(patient_summary$Age_known)))
message(sprintf("  Lineage known: %d (%.1f%%)",
                sum(patient_summary$lineage_known),
                100 * mean(patient_summary$lineage_known)))

# Drop the *_orig diagnostic columns
gh_meta_pre <- gh_meta_pre %>%
  dplyr::select(-HIV_orig, -Sex_orig, -Age_orig, -lineage_orig)

# ---- Florida metadata -------------------------------------------------------
fl_meta <- prep_metadata_florida(fl_epi) %>%
  mutate(
    lineage = forcats::fct_explicit_na(lineage, na_level = "unknown"),
    lineage = factor(
      ifelse(as.character(lineage) == "L4", "L4", "L_non_L4"),
      levels = c("L4", "L_non_L4")
    ),
    date_reported = suppressWarnings(as.Date(date_reported)),
    Year          = suppressWarnings(as.integer(format(date_reported, "%Y")))
  ) %>%
  filter(!is.na(depth_med), depth_med >= MIN_COV_INCLUDE)

message(sprintf("Florida metadata rows retained (medianCov >= %dx): %d",
                MIN_COV_INCLUDE, nrow(fl_meta)))

if (!"cluster_id" %in% names(fl_meta)) {
  stop("00_prep_metadata: fl_meta missing cluster_id. ",
       "External validity analysis requires cluster_id.")
}

# ---- Write artifacts --------------------------------------------------------
saveRDS(gh_meta,                    file.path(PATHS$meta,     "gh_meta.rds"))
saveRDS(fl_meta,                    file.path(PATHS$meta,     "fl_meta.rds"))
saveRDS(gh_pairs,                   file.path(PATHS$meta,     "gh_pairs.rds"))
saveRDS(gh_variants$variants_annot, file.path(PATHS$variants, "gh_variants.rds"))
saveRDS(fl_v_raw,                   file.path(PATHS$variants, "fl_variants.rds"))

# Save metadata provenance summary as a separate small file for the supplement
# Includes both depth provenance (depth_med_source) and demographic provenance
# (HIV_source, Sex_source, Age_source, smear_source).
saveRDS(gh_meta_pre %>%
          dplyr::select(Sample, patientId,
                        depth_med, depth_med_metadata, depth_med_calls, depth_med_source,
                        HIV, HIV_metadata, HIV_source,
                        Sex, Sex_metadata, Sex_source,
                        Age, Age_metadata, Age_source,
                        smear, smear_metadata, smear_source),
        file.path(PATHS$meta, "gh_metadata_provenance.rds"))

message("\n[00_prep_metadata] Done. Artifacts written to ",
        PATHS$meta, "/ and ", PATHS$variants, "/")
message("  Metadata provenance saved to ", PATHS$meta, "/gh_metadata_provenance.rds")
