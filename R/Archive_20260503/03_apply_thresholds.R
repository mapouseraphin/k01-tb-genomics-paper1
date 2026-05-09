# =============================================================================
# 03_apply_thresholds.R
#
# Apply calibrated lexicographic thresholds (from 02_calibrate_lexicographic.R)
# to raw GATK variant tables for both cohorts. Produce per-sample iSNV sets at
# loose, primary, and tighter levels.
#
# Filter chain (controlled by config flags):
#   1. sample in analytic frame
#   2. FILTER == "PASS"
#   3. PE/PPE handling (DROP_PPE_APPLY: canonical FALSE = retain)
#   4. SNV-only filter (canonical TRUE)
#   5. DP, AD1, MAF thresholds
#
# Inputs:
#   PATHS$variants/{gh,fl}_variants.rds
#   PATHS$meta/{gh,fl}_meta.rds
#   PATHS$calibration/thresholds.rds
#
# Outputs (PATHS$thresholded):
#   {gh,fl}_variants_thr_{loose,primary,tighter}.rds
#   {gh,fl}_calibrated_isnv_count.rds  (per-sample counts under primary)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

# ---- Inputs and existence checks --------------------------------------------
inputs <- c(
  fl_variants = file.path(PATHS$variants,    "fl_variants.rds"),
  gh_variants = file.path(PATHS$variants,    "gh_variants.rds"),
  fl_meta     = file.path(PATHS$meta,        "fl_meta.rds"),
  gh_meta     = file.path(PATHS$meta,        "gh_meta.rds"),
  thresholds  = file.path(PATHS$calibration, "thresholds.rds")
)
miss <- inputs[!file.exists(inputs)]
if (length(miss)) {
  stop("[03_apply_thresholds] Missing inputs:\n  ",
       paste(names(miss), miss, sep = ": ", collapse = "\n  "))
}

fl_variants_raw <- readRDS(inputs["fl_variants"]) %>%
  mutate(Sample = as.character(Sample))
gh_variants_raw <- readRDS(inputs["gh_variants"]) %>%
  mutate(Sample = as.character(Sample))
fl_meta    <- readRDS(inputs["fl_meta"])
gh_meta    <- readRDS(inputs["gh_meta"])
thresholds <- readRDS(inputs["thresholds"])

# ---- Threshold structure sanity --------------------------------------------
stopifnot(all(c("primary", "tighter", "loose") %in% names(thresholds)))
required_thr_fields <- c("DP_min", "AD1_min", "MAF_min", "MAF_max")
stopifnot(all(required_thr_fields %in% names(thresholds$primary)))
stopifnot(all(required_thr_fields %in% names(thresholds$loose)))
if (!is.null(thresholds$tighter)) {
  stopifnot(all(required_thr_fields %in% names(thresholds$tighter)))
}

print_thr <- function(label, thr) {
  if (is.null(thr)) {
    cat(sprintf("  %-7s: NULL\n", label))
  } else {
    cat(sprintf("  %-7s: DP=%g  AD1=%g  MAF_min=%.3f  MAF_max=%.3f\n",
                label, thr$DP_min, thr$AD1_min, thr$MAF_min, thr$MAF_max))
  }
}

cat("Calibrated thresholds:\n")
print_thr("loose",   thresholds$loose)
print_thr("primary", thresholds$primary)
print_thr("tighter", thresholds$tighter)
cat(sprintf("\nApplication flags: SNV_ONLY=%s, DROP_PPE_APPLY=%s\n\n",
            SNV_ONLY, DROP_PPE_APPLY))

# ---- Apply thresholds to each cohort x tier --------------------------------
apply_to_cohort <- function(v_raw, meta, thr, label) {
  if (is.null(thr)) return(NULL)
  result <- apply_thresholds(
    v_raw %>% filter(Sample %in% meta$Sample),
    DP_min       = thr$DP_min,
    AD1_min      = thr$AD1_min,
    MAF_min      = thr$MAF_min,
    MAF_max      = thr$MAF_max,
    require_pass = TRUE,
    drop_ppe     = DROP_PPE_APPLY,
    snv_only     = SNV_ONLY
  )
  cat(sprintf("  %-15s: %d calls across %d samples\n",
              label, nrow(result), n_distinct(result$Sample)))
  result
}

cat("Variant calls by cohort x threshold level:\n")
fl_loose   <- apply_to_cohort(fl_variants_raw, fl_meta, thresholds$loose,   "FL loose")
gh_loose   <- apply_to_cohort(gh_variants_raw, gh_meta, thresholds$loose,   "GH loose")
fl_primary <- apply_to_cohort(fl_variants_raw, fl_meta, thresholds$primary, "FL primary")
gh_primary <- apply_to_cohort(gh_variants_raw, gh_meta, thresholds$primary, "GH primary")
fl_tighter <- apply_to_cohort(fl_variants_raw, fl_meta, thresholds$tighter, "FL tighter")
gh_tighter <- apply_to_cohort(gh_variants_raw, gh_meta, thresholds$tighter, "GH tighter")

# ---- Write per-cohort, per-tier variant call tables ------------------------
out_paths <- list(
  fl_loose   = "fl_variants_thr_loose.rds",
  gh_loose   = "gh_variants_thr_loose.rds",
  fl_primary = "fl_variants_thr_primary.rds",
  gh_primary = "gh_variants_thr_primary.rds",
  fl_tighter = "fl_variants_thr_tighter.rds",
  gh_tighter = "gh_variants_thr_tighter.rds"
)

for (key in names(out_paths)) {
  obj <- get(key)
  if (!is.null(obj)) {
    saveRDS(obj, file.path(PATHS$thresholded, out_paths[[key]]))
  }
}

# ---- Per-sample iSNV counts (canonical "detectable" indicator under primary)
build_count <- function(meta, calls) {
  if (is.null(calls)) {
    return(meta %>% select(Sample) %>%
           mutate(n_calibrated_isnv = 0L))
  }
  meta %>%
    select(Sample) %>%
    left_join(calls %>% count(Sample, name = "n_calibrated_isnv"),
              by = "Sample") %>%
    mutate(n_calibrated_isnv = replace_na(n_calibrated_isnv, 0L))
}

fl_isnv_count <- build_count(fl_meta, fl_primary)
gh_isnv_count <- build_count(gh_meta, gh_primary)

saveRDS(fl_isnv_count, file.path(PATHS$thresholded, "fl_calibrated_isnv_count.rds"))
saveRDS(gh_isnv_count, file.path(PATHS$thresholded, "gh_calibrated_isnv_count.rds"))

cat(sprintf("\nFL detectable (>=1 calibrated iSNV under primary): %d / %d (%.1f%%)\n",
            sum(fl_isnv_count$n_calibrated_isnv >= 1), nrow(fl_isnv_count),
            100 * mean(fl_isnv_count$n_calibrated_isnv >= 1)))
cat(sprintf("GH detectable (>=1 calibrated iSNV under primary): %d / %d (%.1f%%)\n",
            sum(gh_isnv_count$n_calibrated_isnv >= 1), nrow(gh_isnv_count),
            100 * mean(gh_isnv_count$n_calibrated_isnv >= 1)))

cat(sprintf("\n[03_apply_thresholds] Done. Outputs written to: %s/\n",
            PATHS$thresholded))
