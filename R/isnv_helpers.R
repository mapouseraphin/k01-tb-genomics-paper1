# =============================================================================
# R/isnv_helpers.R
#
# Shared helper functions for the iSNV calibration / transportability pipeline.
# Sourced by every pipeline script (00 -> 06).
#
# Conventions:
#   - SNV-only: REF and ALT both length 1 (substitutions; excludes indels)
#   - PE/PPE handling: drop_ppe argument controls whether PE/PPE positions
#     are excluded; FILTER == "PASS" enforced separately via require_pass.
#   - Canonical primary defaults: snv_only = TRUE, drop_ppe = FALSE
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(forcats)
  library(tibble)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------------------------------------------------------------
# Recoders -- shared by both cohorts
# -----------------------------------------------------------------------------
recode_sex <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    x %in% c("M", "Male", "male", "MALE", "m","1") ~ "Male",
    x %in% c("F", "Female", "female", "FEMALE", "f","2") ~ "Female",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Female", "Male"))
}

recode_sex_ghana <- function(x) {
  # Ghana cohort patient_tbl.rds encodes sex as numeric 1 = Male, 2 = Female.
  # Falls through to character variants for cross-cohort safety.
  x_num <- suppressWarnings(as.numeric(as.character(x)))
  x_chr <- as.character(x)
  out <- dplyr::case_when(
    !is.na(x_num) & x_num == 1                       ~ "Male",
    !is.na(x_num) & x_num == 2                       ~ "Female",
    x_chr %in% c("M", "Male", "male", "MALE", "m")   ~ "Male",
    x_chr %in% c("F", "Female", "female", "FEMALE", "f") ~ "Female",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("Female", "Male"))
}

recode_hiv <- function(x) {
  x <- as.character(x)
  out <- dplyr::case_when(
    x %in% c("Positive", "positive", "POS", "Pos", "pos",
             "Yes", "yes", "Y", "1", "TRUE") ~ "HIV+",
    x %in% c("Negative", "negative", "NEG", "Neg", "neg",
             "No",  "no",  "N", "0", "FALSE") ~ "HIV-",
    TRUE ~ NA_character_
  )
  factor(out, levels = c("HIV-", "HIV+"))
}

# -----------------------------------------------------------------------------
# Variant-class filter (D46)
# -----------------------------------------------------------------------------
# SNV-only restriction: keep only single-nucleotide substitutions.
# Excludes insertions, deletions, and complex events.
#
# Motivation:
#   - Distinct error mechanisms (polymerase slippage, alignment ambiguity)
#   - Asymmetry with SNP-distance reference standard (substitutions only)
#   - 06 diagnostics showed indel-driven concentration in GH between-cluster
#     pool (lineage-typical and repeat-region indels).
filter_snv_only <- function(v) {
  stopifnot(all(c("REF", "ALT") %in% names(v)))
  v %>% dplyr::filter(nchar(REF) == 1L, nchar(ALT) == 1L)
}

# -----------------------------------------------------------------------------
# PE/PPE exclusion
# -----------------------------------------------------------------------------
# Excludes variant calls in PE/PPE gene regions. Tries multiple column
# encodings (logical, numeric, character) for robustness.
exclude_ppe <- function(v) {
  if ("ppe_flag" %in% names(v)) {
    pp <- v$ppe_flag
    drop_levels <- c("PPE", "pe/ppe", "PE/PPE", "TRUE", "true", "1")
    if (is.logical(pp)) {
      return(v %>% filter(!ppe_flag | is.na(ppe_flag)))
    } else if (is.numeric(pp)) {
      return(v %>% filter(ppe_flag == 0 | is.na(ppe_flag)))
    } else if (is.character(pp) || is.factor(pp)) {
      return(v %>% filter(!(as.character(ppe_flag) %in% drop_levels) |
                          is.na(ppe_flag)))
    }
  }
  if ("PPE" %in% names(v)) {
    pp <- v$PPE
    if (is.logical(pp)) return(v %>% filter(!PPE | is.na(PPE)))
    if (is.numeric(pp)) return(v %>% filter(PPE == 0 | is.na(PPE)))
    if (is.character(pp))
      return(v %>% filter(!(PPE %in% c("1","TRUE","true","PPE","pe/ppe")) |
                          is.na(PPE)))
  }
  v  # no PPE column found; return unchanged
}

# -----------------------------------------------------------------------------
# MAF derivation (idempotent)
# -----------------------------------------------------------------------------
ensure_maf <- function(v) {
  if ("MAF" %in% names(v)) return(v)
  if (!all(c("AD1", "AD2", "DP") %in% names(v))) {
    stop("ensure_maf: need MAF column or (AD1, AD2, DP) to derive it.")
  }
  v %>% mutate(MAF = pmin(AD1, AD2) / DP)
}

# -----------------------------------------------------------------------------
# apply_thresholds()
# -----------------------------------------------------------------------------
# Apply the calibrated decision rule to a raw variant table. Order of
# operations:
#   1. ensure MAF derived
#   2. require FILTER == "PASS" (gate on quality)
#   3. drop PE/PPE positions (optional; primary = retain)
#   4. SNV-only filter (canonical = TRUE)
#   5. DP/AD1/MAF threshold filtering
#
# Defaults reflect canonical primary specification:
#   - require_pass = TRUE  (always)
#   - drop_ppe     = FALSE (canonical: retain PE/PPE)
#   - snv_only     = TRUE  (canonical: SNVs only)
apply_thresholds <- function(v, DP_min, AD1_min, MAF_min, MAF_max = 0.45,
                             require_pass = TRUE,
                             drop_ppe     = FALSE,
                             snv_only     = TRUE) {
  v2 <- ensure_maf(v)
  if (require_pass && "FILTER" %in% names(v2)) v2 <- v2 %>% filter(FILTER == "PASS")
  if (drop_ppe)  v2 <- exclude_ppe(v2)
  if (snv_only)  v2 <- filter_snv_only(v2)
  v2 %>% filter(DP >= DP_min, AD1 >= AD1_min, MAF >= MAF_min, MAF <= MAF_max)
}

# -----------------------------------------------------------------------------
# Pairwise concordance metrics (for calibration on Ghana replicates)
# -----------------------------------------------------------------------------
site_key <- function(df) paste(df$CHROM, df$POS, df$REF, df$ALT, sep = ":")

jaccard_pair <- function(a_sites, b_sites) {
  inter <- length(intersect(a_sites, b_sites))
  uni   <- length(union(a_sites, b_sites))
  if (uni == 0) 0 else inter / uni
}

overlap_pair <- function(a_sites, b_sites) {
  inter    <- length(intersect(a_sites, b_sites))
  min_size <- min(length(a_sites), length(b_sites))
  if (min_size == 0) 0 else inter / min_size
}

confirm_rate_pair <- function(a_tbl, b_tbl) {
  # Returns four named values per pair:
  #   confirm_both   = |A ∩ B| / |A ∪ B|        (= Jaccard on the evaluable
  #                                                filter; kept for diagnostic
  #                                                transparency only; NOT used
  #                                                in lex sort or Pareto plane)
  #   confirm_a_to_b = |A ∩ B| / |A|            (asymmetric: fraction of A's
  #                                                calls confirmed in B)
  #   confirm_b_to_a = |A ∩ B| / |B|            (asymmetric: fraction of B's
  #                                                calls confirmed in A)
  #   confirm_sym    = 0.5 * (a_to_b + b_to_a)  (SYMMETRIC AVERAGE of the
  #                                                two asymmetric rates; used
  #                                                in lex level 5 and as the
  #                                                Pareto y-axis from 2026-05-05
  #                                                onward; replaces the
  #                                                redundant confirm_both = J)
  #
  # Reference for asymmetric per-replicate confirmation in iSNV calibration:
  #   [CITATION: McCrone JT, Lauring AS. Measurements of intrahost viral
  #    diversity are extremely sensitive to systematic errors in variant
  #    calling. J Virol. 2016;90(15):6884-6895. doi:10.1128/JVI.00667-16.
  # The symmetric average extension (confirm_sym) is our convention for
  # this paper.
  a_sites <- site_key(a_tbl); b_sites <- site_key(b_tbl)
  union_sites <- union(a_sites, b_sites)
  if (length(union_sites) == 0) {
    return(c(confirm_both   = NA_real_,
             confirm_a_to_b = NA_real_,
             confirm_b_to_a = NA_real_,
             confirm_sym    = NA_real_))
  }
  a_has <- union_sites %in% a_sites
  b_has <- union_sites %in% b_sites
  a_to_b <- if (length(a_sites)) mean(a_sites %in% b_sites) else NA_real_
  b_to_a <- if (length(b_sites)) mean(b_sites %in% a_sites) else NA_real_
  c(confirm_both   = mean(a_has & b_has),
    confirm_a_to_b = a_to_b,
    confirm_b_to_a = b_to_a,
    confirm_sym    = if (is.na(a_to_b) || is.na(b_to_a)) NA_real_
                     else 0.5 * (a_to_b + b_to_a))
}

maf_cor_pair <- function(a_tbl, b_tbl) {
  a2 <- a_tbl %>% transmute(k = paste(CHROM, POS, REF, ALT, sep = ":"), MAF = MAF)
  b2 <- b_tbl %>% transmute(k = paste(CHROM, POS, REF, ALT, sep = ":"), MAF = MAF)
  by_site <- full_join(a2, b2, by = "k", suffix = c(".a", ".b")) %>%
    mutate(MAF.a = replace_na(MAF.a, 0),
           MAF.b = replace_na(MAF.b, 0))
  if (nrow(by_site) < 2 || sd(by_site$MAF.a) == 0 || sd(by_site$MAF.b) == 0) {
    return(NA_real_)
  }
  cor(by_site$MAF.a, by_site$MAF.b, use = "complete.obs")
}

# -----------------------------------------------------------------------------
# Visit / slot parsing for Ghana sample IDs
# -----------------------------------------------------------------------------
# Convention: KBTH###-<visit><slot>
#   Visit 0 (M0): E0, S0.1, S0.2
#   Visit 1 (M1): E1, S1.1, S1.2
#   Visit 2 (M2): E2, S2.1, S2.2
parse_visit <- function(Sample) {
  v <- str_extract(as.character(Sample), "(?<=-)[ES][0-9]")
  case_when(
    str_detect(v, "0$") ~ "M0",
    str_detect(v, "1$") ~ "M1",
    str_detect(v, "2$") ~ "M2",
    TRUE                ~ NA_character_
  )
}

parse_slot <- function(Sample) {
  # Conventions:
  #   - Canonical Sample IDs end with E0/E1/E2 (mapped from slot_code E0.1/E1.1/
  #     E2.1 by 00_prep_metadata.R's slot_to_variant_suffix()) or with S0.1/S0.2/
  #     S1.1/S1.2/S2.1/S2.2.
  #   - Typo-path Sample IDs (when the sequencing pipeline used slot_code as the
  #     Sample name) end with E0.1/E1.1/E2.1 verbatim. The regex below matches
  #     both forms by allowing an optional ".1" suffix on E[012].
  s <- as.character(Sample)
  case_when(
    str_detect(s, "E0(\\.1)?$|EM$")                ~ "EM",
    str_detect(s, "S0\\.1$|S01$")                  ~ "S0.1",
    str_detect(s, "S0\\.2$|S02$")                  ~ "S0.2",
    str_detect(s, "E1(\\.1)?$")                    ~ "E1",
    str_detect(s, "S1\\.1$|S11$")                  ~ "S1.1",
    str_detect(s, "S1\\.2$|S12$")                  ~ "S1.2",
    str_detect(s, "E2(\\.1)?$")                    ~ "E2",
    str_detect(s, "S2\\.1$|S21$")                  ~ "S2.1",
    str_detect(s, "S2\\.2$|S22$")                  ~ "S2.2",
    TRUE                                            ~ NA_character_
  )
}

is_m0_sample <- function(Sample) {
  parse_visit(Sample) == "M0" & !is.na(parse_visit(Sample))
}

# -----------------------------------------------------------------------------
# Metadata preparation -- Ghana
# -----------------------------------------------------------------------------
prep_metadata_ghana_from_epi <- function(ghana_epiData, parsed_meta) {
  gh_meta <- parsed_meta %>%
    dplyr::left_join(ghana_epiData,
                     by = c("person_id" = "patientId", "Sample")) %>%
    dplyr::mutate(
      HIV       = recode_hiv(livingWithHIV),
      Sex       = recode_sex(dplyr::coalesce(sex)),
      Age       = suppressWarnings(as.numeric(age)),
      smear     = factor(sputumSmear),
      lineage   = factor(majorLineage),
      callable  = suppressWarnings(as.numeric(callable_frac_10x %||%
                                              callable_frac_20x)),
      depth_med = suppressWarnings(as.numeric(medianCov)),
      cohort    = "Ghana",
      run_id    = "GH_single_run",
      Sample    = as.character(Sample)
    ) %>%
    dplyr::select(
      Sample, HIV, Sex, Age, smear, lineage, callable, depth_med,
      run_id, cohort, dplyr::everything()
    )
  if (all(is.na(gh_meta$HIV))) {
    warning("All HIV values are NA after recoding -- check livingWithHIV encoding.")
  }
  gh_meta
}

# -----------------------------------------------------------------------------
# Metadata preparation -- Florida
# -----------------------------------------------------------------------------
prep_metadata_florida <- function(florida_epiData) {
  fl_meta <- florida_epiData %>%
    dplyr::mutate(
      HIV       = recode_hiv(livingWithHIV %||% HIV),
      Sex       = recode_sex(sex %||% Sex),
      Age       = suppressWarnings(as.numeric(age %||% Age)),
      lineage   = factor(majorLineage %||% lineage),
      callable  = suppressWarnings(as.numeric(callable_frac_10x %||%
                                              callable_frac_20x)),
      depth_med = suppressWarnings(as.numeric(medianCov)),
      cohort    = "Florida",
      Sample    = as.character(Sample)
    ) %>%
    dplyr::select(
      Sample, HIV, Sex, Age, lineage, callable, depth_med,
      cluster_id, cohort, dplyr::everything()
    )
  fl_meta
}

# =============================================================================
# prep_metadata_ghana_v12.R
#
# DEPRECATED 2026-05-05.
# This function is not called by the canonical Paper 1 pipeline (which uses
# 00_prep_metadata.R directly to build gh_meta from D2_longitudinal.rds and
# patient_tbl.rds). It is retained only for archival reference; do not call
# from new code. Known incompatibility: returns gh_meta with a `visit` column
# but no `visit_label`, which would break application_persistence_analysis.R.
.Deprecated_prep_metadata_ghana_v12 <- TRUE  # marker; the function below is
                                             # NOT removed for git-history reasons
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