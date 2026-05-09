# =============================================================================
# 03_apply_thresholds.R
#
# Apply calibrated thresholds (from 02_calibrate_*.R / 02b_calibration_bootstrap.R)
# to the raw GATK variant table. Produce per-sample iSNV sets at looser, primary,
# and tighter levels.
#
# Filter chain (controlled by config flags):
#   1. sample in analytic frame
#   2. FILTER == "PASS"
#   3. PE/PPE handling (DROP_PPE_APPLY: canonical FALSE = retain)
#   4. SNV-only filter (canonical TRUE)
#   5. DP, AD1, MAF thresholds
#
# Three-tier ladder (post 2026-05-03):
#   looser  = real-grid less-stringent reference (selection: looser)
#   primary = bootstrap-stabilized headline rule (selection: primary)
#   tighter = real-grid more-stringent reference (selection: tighter)
#
# Backwards-compat shim: if a stale thresholds.rds with the old `loose` slot
# is encountered, it is treated as `looser` with a deprecation warning.
#
# Inputs:
#   PATHS$variants/gh_variants.rds
#   PATHS$meta/gh_meta.rds
#   PATHS$calibration/thresholds.rds
#
# Outputs (PATHS$thresholded):
#   gh_variants_thr_{looser,primary,tighter}.rds
#   gh_calibrated_isnv_count.rds  (per-sample counts under primary)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

banner()

# ---- Inputs and existence checks --------------------------------------------
inputs <- c(
  gh_variants = file.path(PATHS$variants,    "gh_variants.rds"),
  gh_meta     = file.path(PATHS$meta,        "gh_meta.rds"),
  thresholds  = file.path(PATHS$calibration, "thresholds.rds")
)
miss <- inputs[!file.exists(inputs)]
if (length(miss)) {
  stop("[03_apply_thresholds] Missing inputs:\n  ",
       paste(names(miss), miss, sep = ": ", collapse = "\n  "))
}

gh_variants_raw <- readRDS(inputs["gh_variants"]) %>%
  mutate(Sample = as.character(Sample))
gh_meta    <- readRDS(inputs["gh_meta"])
thresholds <- readRDS(inputs["thresholds"])

# ---- Backwards-compat shim: old `loose` slot -> `looser` -------------------
# If a stale thresholds.rds (from a prior run with the no-filter sentinel)
# is present, normalize the names and warn. Old `loose` = no-filter sentinel,
# which was retired 2026-05-03; treating it as `looser` here allows downstream
# scripts to run, but the upstream calibration script(s) should be re-run
# under the new spec to produce a real-grid `looser` tier.
if ("loose" %in% names(thresholds) && !"looser" %in% names(thresholds)) {
  warning("[03_apply_thresholds] thresholds.rds uses legacy `loose` slot. ",
          "This is the retired no-filter sentinel. Renaming to `looser` for ",
          "compatibility, but you should re-run 02_calibrate_*.R to produce ",
          "a real-grid `looser` tier under the post-2026-05-03 spec.")
  thresholds$looser <- thresholds$loose
  thresholds$loose  <- NULL
}

# ---- Threshold structure sanity --------------------------------------------
required_top <- c("primary")  # primary is mandatory; looser/tighter optional
miss_top <- setdiff(required_top, names(thresholds))
if (length(miss_top))
  stop("[03_apply_thresholds] thresholds.rds missing required slot(s): ",
       paste(miss_top, collapse = ", "))

required_thr_fields <- c("DP_min", "AD1_min", "MAF_min", "MAF_max")
stopifnot(all(required_thr_fields %in% names(thresholds$primary)))
if (!is.null(thresholds$looser))
  stopifnot(all(required_thr_fields %in% names(thresholds$looser)))
if (!is.null(thresholds$tighter))
  stopifnot(all(required_thr_fields %in% names(thresholds$tighter)))

print_thr <- function(label, thr) {
  if (is.null(thr)) {
    cat(sprintf("  %-7s: NULL\n", label))
  } else {
    cat(sprintf("  %-7s: DP=%g  AD1=%g  MAF_min=%.3f  MAF_max=%.3f\n",
                label, thr$DP_min, thr$AD1_min, thr$MAF_min, thr$MAF_max))
  }
}

cat("Calibrated thresholds:\n")
print_thr("looser",  thresholds$looser)
print_thr("primary", thresholds$primary)
print_thr("tighter", thresholds$tighter)
cat(sprintf("\nApplication flags: SNV_ONLY=%s, DROP_PPE_APPLY=%s\n\n",
            SNV_ONLY, DROP_PPE_APPLY))

# ---- Apply thresholds to the cohort, per tier ------------------------------
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

cat("Variant calls by threshold level (Ghana cohort):\n")
gh_looser  <- apply_to_cohort(gh_variants_raw, gh_meta, thresholds$looser,  "GH looser")
gh_primary <- apply_to_cohort(gh_variants_raw, gh_meta, thresholds$primary, "GH primary")
gh_tighter <- apply_to_cohort(gh_variants_raw, gh_meta, thresholds$tighter, "GH tighter")

# ---- Write per-tier variant call tables ------------------------------------
out_paths <- list(
  gh_looser  = "gh_variants_thr_looser.rds",
  gh_primary = "gh_variants_thr_primary.rds",
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

gh_isnv_count <- build_count(gh_meta, gh_primary)

saveRDS(gh_isnv_count, file.path(PATHS$thresholded, "gh_calibrated_isnv_count.rds"))

cat(sprintf("\nGH detectable (>=1 calibrated iSNV under primary): %d / %d (%.1f%%)\n",
            sum(gh_isnv_count$n_calibrated_isnv >= 1), nrow(gh_isnv_count),
            100 * mean(gh_isnv_count$n_calibrated_isnv >= 1)))

cat(sprintf("\n[03_apply_thresholds] Done. Outputs written to: %s/\n",
            PATHS$thresholded))
