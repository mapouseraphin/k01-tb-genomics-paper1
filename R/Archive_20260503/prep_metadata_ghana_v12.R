# =============================================================================
# prep_metadata_ghana_v12.R
#
# Build Ghana metadata in Florida-compatible schema from v1.2 Ghana pipeline
# outputs (D2_longitudinal, patient_tbl, mixed_infection).
#
# Replaces prep_metadata_ghana_from_epi() (deprecated; consumed v1.0/v1.1
# pipeline's gh_epi + parsed_meta shape).
#
# Decisions encoded:
#   D1.a   M0 24-hour samples only: slot_code in {S0.1, S0.2, E0.1}
#          (DX excluded by A6: Xpert-only, never sequenced; M1/M2 excluded
#          by D1: treatment-induced information bias.)
#   D2.a   Smear: scanty + 1+/2+/3+ -> Positive; negative -> Negative; NA -> Unknown.
#          (Reporting caveat: scanty is grader-inconsistent across labs.)
#   D3     Lineage prep: extract major lineage from Strain. Binarization to
#          L4 vs L_non_L4 happens downstream in orchestration so the same
#          binarize step applies symmetrically to Florida and Ghana.
#   D6     callable = PCT_20X / 100  (callable_frac_20x analog).
#          Defensive: auto-detects whether values are in percent (>1.5) or
#          fraction ([0,1]).
#   D7     cluster_id = NA on Ghana side (no Florida-style transmission
#          clusters in this dataset); patientId is the Ghana cluster unit
#          per D45 (Ghana project).
#
# Dependencies:
#   isnv_helpers.R  (recode_hiv, recode_sex)
#
# Inputs (R objects, loaded by orchestration script from v1.2 derived_data/v_{tag}/):
#   D2_long  -- D2_longitudinal.rds   (sample-unit, n=176 in v1.2)
#   pat_tbl  -- patient_tbl.rds       (patient-unit, n=150)
#   mix      -- mixed_infection.rds   (sample-unit; optional)
#
# Output:
#   tibble in Florida-compatible schema; M0 24-hour samples only.
#   Columns: Sample, patientId, HIV, Sex, Age, smear, lineage (major form),
#            callable, depth_med, run_id, cohort, cluster_id, visit, slot_code,
#            n_mixed_variants, mixed_infection_flag
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
})

# -----------------------------------------------------------------------------
# Helper: recode_smear_pooled (D2 + D2.a)
# Vocab per CODEBOOK $VOCAB_SMEAR_GRADE: negative | scanty | 1+ | 2+ | 3+
# -----------------------------------------------------------------------------
recode_smear_pooled <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    x %in% c("negative", "Negative", "neg", "NEG", "NEGATIVE") ~ "Negative",
    x %in% c("scanty", "Scanty", "SCANTY",
             "1+", "+1", "2+", "+2", "3+", "+3") ~ "Positive",
    TRUE ~ NA_character_
  )
  factor(ifelse(is.na(out), "Unknown", out),
         levels = c("Negative", "Positive", "Unknown"))
}

# -----------------------------------------------------------------------------
# Helper: extract_major_lineage (D3 prep step)
# tbProfiler Strain -> major lineage label.
# Examples:
#   "lineage4.3.4.2" -> "L4"
#   "lineage5.1"     -> "L5"
#   "lineage4;lineage4.3" -> "L4" (mixed; takes first/dominant)
#   "M.bovis", NA, "non_mtb" -> "L_unknown"
# -----------------------------------------------------------------------------
extract_major_lineage <- function(strain) {
  s <- as.character(strain)
  # For mixed (semicolon-separated) calls, take first/dominant
  s_first <- stringr::str_split(s, ";", simplify = TRUE)[, 1]
  s_first <- stringr::str_trim(s_first)
  major   <- stringr::str_match(s_first, "^lineage([0-9]+)")[, 2]
  out     <- ifelse(is.na(major), "L_unknown", paste0("L", major))
  factor(out)
}

# -----------------------------------------------------------------------------
# Helper: binarize_lineage_l4 (D3)
# Major lineage -> L4 vs L_non_L4. Reference = L4.
# NOTE: L_unknown is folded into L_non_L4. This conflates "not L4" with
# "could not classify". Documented as a study limitation per Nancy's directive.
# Apply symmetrically to Florida (which has L1/L2/L4/unknown) and Ghana
# (which has L4-dominant + L5/L6 + L_unknown).
# -----------------------------------------------------------------------------
binarize_lineage_l4 <- function(major_lineage) {
  x <- as.character(major_lineage)
  out <- ifelse(x == "L4", "L4", "L_non_L4")
  factor(out, levels = c("L4", "L_non_L4"))
}

# -----------------------------------------------------------------------------
# Main: prep_metadata_ghana_v12()
# -----------------------------------------------------------------------------
prep_metadata_ghana_v12 <- function(D2_long, pat_tbl, mix = NULL) {

  # ---- 1. D1.a: M0 24-hour filter ------------------------------------------
  if (!"slot_code" %in% names(D2_long)) {
    stop("D2_long missing slot_code column; cannot apply D1.a M0 filter. ",
         "Per CODEBOOK $sample_tbl, slot_code is derived in module 07. ",
         "Verify D2_longitudinal.rds was built with the v1.2 pipeline.")
  }
  M0_SLOTS <- c("S0.1", "S0.2", "E0.1")
  d2_m0 <- D2_long %>% dplyr::filter(slot_code %in% M0_SLOTS)
  message(sprintf(
    "[D1.a] D2 -> M0 24-hour: %d -> %d samples (%d patients; slots kept: %s)",
    nrow(D2_long), nrow(d2_m0),
    dplyr::n_distinct(d2_m0$patientId),
    paste(M0_SLOTS, collapse = ", ")
  ))

  # ---- 2. Sex codebook preprocess (gender 0/1 -> Male/Female) --------------
  # Per CODEBOOK $patient_tbl: gender is chr "0"=Male, "1"=Female (raw codes).
  # recode_sex() in isnv_helpers.R does NOT accept integer-coded values, so
  # preprocess here.
  if (!"gender" %in% names(pat_tbl)) {
    stop("patient_tbl missing gender column.")
  }
  pat_prep <- pat_tbl %>%
    dplyr::mutate(
      sex_chr = dplyr::case_when(
        as.character(gender) %in% c("0", "Male", "male", "M", "m", "1") ~ "Male",
        as.character(gender) %in% c("1", "Female", "female", "F", "f", "2") ~ "Female",
        TRUE ~ NA_character_
      )
    )

  # ---- 3. D6: resolve callable from PCT_20X / callable_frac_20x ------------
  callable_candidates <- c("callable_frac_20x", "callable_frac_10x",
                           "PCT_20X", "pct_20x", "pct_callable_20x")
  callable_col <- intersect(callable_candidates, names(d2_m0))[1]
  if (is.na(callable_col)) {
    stop("No callable column found in D2_long. Looked for: ",
         paste(callable_candidates, collapse = ", "),
         ". If callable is in sample_tbl but not D2_long, join sample_tbl ",
         "to D2 on sample_id BEFORE calling this function.")
  }
  callable_raw <- as.numeric(d2_m0[[callable_col]])
  d2_m0$callable_resolved <- if (max(callable_raw, na.rm = TRUE) > 1.5) {
    callable_raw / 100   # PCT_20X is in percent units
  } else {
    callable_raw         # already in [0,1]
  }
  message(sprintf(
    "[D6] callable resolved from '%s' (range [%.4f, %.4f])",
    callable_col,
    min(d2_m0$callable_resolved, na.rm = TRUE),
    max(d2_m0$callable_resolved, na.rm = TRUE)
  ))

  # ---- 4. Resolve coverage median (depth_med) ------------------------------
  cov_candidates <- c("coverage_median", "MEDIAN_COVERAGE", "medianCov")
  cov_col <- intersect(cov_candidates, names(d2_m0))[1]
  if (is.na(cov_col)) {
    stop("No coverage column found in D2_long. Looked for: ",
         paste(cov_candidates, collapse = ", "))
  }
  d2_m0$depth_med_resolved <- as.numeric(d2_m0[[cov_col]])
  message(sprintf(
    "[depth_med] resolved from '%s' (range [%.0f, %.0f])",
    cov_col,
    min(d2_m0$depth_med_resolved, na.rm = TRUE),
    max(d2_m0$depth_med_resolved, na.rm = TRUE)
  ))

  # ---- 5. Resolve Strain column (for D3 prep step) -------------------------
  strain_candidates <- c("tbp_lineage", "Strain", "strain", "tbp_strain")
  strain_col <- intersect(strain_candidates, names(d2_m0))[1]
  if (is.na(strain_col)) {
    stop("No tbProfiler Strain column found in D2_long. Looked for: ",
         paste(strain_candidates, collapse = ", "))
  }
  message(sprintf("[D3 prep] Strain column resolved: '%s'", strain_col))

  # ---- 6. Resolve smear column ---------------------------------------------
  smear_candidates <- c("smear", "sputumSmear", "sputum_smear")
  smear_col <- intersect(smear_candidates, names(d2_m0))[1]
  if (is.na(smear_col)) {
    stop("No smear column found in D2_long. Looked for: ",
         paste(smear_candidates, collapse = ", "))
  }

  # ---- 7. Resolve sample-id column -----------------------------------------
  sample_candidates <- c("sample_id", "Sample")
  sample_col <- intersect(sample_candidates, names(d2_m0))[1]
  if (is.na(sample_col)) {
    stop("No sample-id column found in D2_long. Looked for: ",
         paste(sample_candidates, collapse = ", "))
  }

  # ---- 8. Join D2 + patient_tbl (patient-level covariates) -----------------
  if (!"patientId" %in% names(d2_m0) || !"patientId" %in% names(pat_prep)) {
    stop("patientId must be present in both D2_long and patient_tbl.")
  }
  joined <- d2_m0 %>%
    dplyr::left_join(
      pat_prep %>% dplyr::select(patientId, hiv_positive, sex_chr, age),
      by = "patientId"
    )

  # ---- 9. Optional join with mixed_infection -------------------------------
  if (!is.null(mix)) {
    mix_join_col <- intersect(c("sample_id", "Sample"), names(mix))[1]
    if (!is.na(mix_join_col)) {
      mix_cols <- intersect(c("n_mixed_variants", "mixed_infection_flag"),
                            names(mix))
      mix_sub <- mix %>%
        dplyr::select(dplyr::all_of(c(mix_join_col, mix_cols)))
      joined <- joined %>%
        dplyr::left_join(
          mix_sub,
          by = setNames(mix_join_col, sample_col)
        )
    }
  }

  # ---- 10. Build Florida-compatible tibble ---------------------------------
  gh_meta <- joined %>%
    dplyr::transmute(
      Sample        = as.character(.data[[sample_col]]),
      patientId     = as.character(patientId),
      HIV           = recode_hiv(hiv_positive),       # from isnv_helpers.R
      Sex           = recode_sex(sex_chr),            # from isnv_helpers.R
      Age           = suppressWarnings(as.numeric(age)),
      smear         = recode_smear_pooled(.data[[smear_col]]),
      lineage       = extract_major_lineage(.data[[strain_col]]),
      callable      = callable_resolved,
      depth_med     = depth_med_resolved,
      run_id        = "GH_single_run",
      cohort        = "Ghana",
      cluster_id    = NA_integer_,                    # D7: Ghana has no cluster
      visit         = if ("visit" %in% names(joined)) as.character(visit) else "M0",
      slot_code     = as.character(slot_code),
      n_mixed_variants = if ("n_mixed_variants" %in% names(joined)) {
                            as.integer(n_mixed_variants)
                          } else NA_integer_,
      mixed_infection_flag = if ("mixed_infection_flag" %in% names(joined)) {
                                as.logical(mixed_infection_flag)
                             } else NA
    )

  # ---- 11. Sanity warnings -------------------------------------------------
  if (all(is.na(gh_meta$HIV))) {
    warning("All HIV values are NA after recoding. ",
            "Check hiv_positive encoding in patient_tbl.")
  }
  if (all(is.na(gh_meta$Sex))) {
    warning("All Sex values are NA after recoding. ",
            "Check gender codebook in patient_tbl (expected 0=Male, 1=Female).")
  }
  if (sum(gh_meta$lineage == "L4", na.rm = TRUE) == 0) {
    warning("No L4 samples after lineage extraction. ",
            "Check tbp_lineage column for sublineage strings.")
  }

  message(sprintf(
    "[prep_metadata_ghana_v12] Output: %d samples across %d patients",
    nrow(gh_meta), dplyr::n_distinct(gh_meta$patientId)
  ))

  gh_meta
}
