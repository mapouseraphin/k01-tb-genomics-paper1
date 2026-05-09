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
# Inputs:
#   data_raw/variant_tables_raw/ghana_gatkVariantTablesFormatted.rds
#   data_raw/variant_tables_raw/florida_gatkVariantTablesFormatted.rds
#   data_raw/ghana_epiData.rds
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
fl_epi <- readRDS(file.path(PATHS$data_raw, "florida_epiData.rds")) %>%
  droplevels()

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
gh_meta <- prep_metadata_ghana_from_epi(gh_epi, gh_variants$meta) %>%
  mutate(lineage = forcats::fct_explicit_na(lineage, na_level = "unknown")) %>%
  filter(!is.na(depth_med), depth_med >= MIN_COV_INCLUDE)

if (!"person_id" %in% names(gh_meta) && !"patientId" %in% names(gh_meta)) {
  stop("00_prep_metadata: gh_meta missing patientId/person_id. ",
       "Cluster-robust SE on patientId requires this column.")
}
if ("person_id" %in% names(gh_meta) && !"patientId" %in% names(gh_meta)) {
  gh_meta <- gh_meta %>% rename(patientId = person_id)
}

message(sprintf(
  "Ghana metadata rows retained (medianCov >= %dx): %d rows across %d patients",
  MIN_COV_INCLUDE, nrow(gh_meta), n_distinct(gh_meta$patientId)
))

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

message("\n[00_prep_metadata] Done. Artifacts written to ",
        PATHS$meta, "/ and ", PATHS$variants, "/")
