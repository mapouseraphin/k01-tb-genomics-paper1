# =============================================================================
# diagnose_patientId_missing.R
#
# Diagnoses why patientId is NA in some rows of lca_dataset.rds.
#
# Hypothesis: gh_meta is coverage-filtered (depth_med >= 50) but gh_pairs is
# not. Pairs containing a low-coverage or coverage-NA sample lose their
# sampleA -> gh_meta join, producing NA patientId in build_lca_dataset.R.
#
# Diagnostic outputs (all to stderr via message()):
#   D1. NA patientId row count and pair count
#   D2. Distinct sampleA values with NA patientId
#   D3. Are those sampleA values present in gh_meta?
#   D4. Coverage (depth_med) of those samples (from raw variant table or
#       Ghana metadata pre-filter, if recoverable)
#   D5. Whether sampleB of the same pair *would* have produced a patientId
#   D6. Whether patientId can be parsed from the Sample ID directly
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble)
})

LCA_PATH <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
lca      <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))
gh_meta  <- readRDS(file.path(PATHS$meta, "gh_meta.rds"))
gh_pairs <- readRDS(file.path(PATHS$meta, "gh_pairs.rds"))

print_via_message <- function(df, n = NULL) {
  if (!is.null(n)) df <- head(df, n)
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}

# ---- D1: scope of the problem ----------------------------------------------
n_total   <- nrow(lca)
n_na_pid  <- sum(is.na(lca$patientId))
pairs_with_na <- lca %>% filter(is.na(patientId)) %>%
  distinct(pair_id, sampleA, sampleB)

message("\n--- D1: Scope ---")
message(sprintf("Total LCA rows:                 %s", format(n_total, big.mark = ",")))
message(sprintf("Rows with NA patientId:         %s (%.1f%%)",
                format(n_na_pid, big.mark = ","), 100 * n_na_pid / n_total))
message(sprintf("Pairs with at least one NA pid: %d", nrow(pairs_with_na)))
message(sprintf("Distinct sampleA in NA rows:    %d",
                n_distinct(pairs_with_na$sampleA)))

if (n_na_pid == 0) {
  message("\nNo NA patientId rows. Diagnostic exits.")
} else {

  # ---- D2: which sampleA values are missing? ------------------------------
  na_samples <- unique(c(pairs_with_na$sampleA))
  message("\n--- D2: Distinct sampleA values producing NA patientId ---")
  message(paste(na_samples, collapse = ", "))

  # ---- D3: are those sampleA values present in gh_meta? -------------------
  in_meta_A <- na_samples %in% gh_meta$Sample
  d3 <- tibble(sampleA = na_samples,
               in_gh_meta = in_meta_A)
  message("\n--- D3: sampleA presence in gh_meta ---")
  print_via_message(d3)
  message(sprintf("Of %d distinct sampleA, %d are in gh_meta, %d are not.",
                  length(na_samples), sum(in_meta_A), sum(!in_meta_A)))

  # ---- D4: coverage of missing sampleA, if recoverable from the Ghana    --
  # epiData (which lives upstream of the depth filter)                     --
  raw_meta <- tryCatch(
    readRDS(file.path(PATHS$data_raw, "ghana_epiData.rds")),
    error = function(e) NULL
  )
  if (!is.null(raw_meta) && "Sample" %in% names(raw_meta)) {
    cov_lookup <- raw_meta %>%
      transmute(Sample = as.character(Sample),
                medianCov_raw = suppressWarnings(as.numeric(medianCov)))
    d4 <- tibble(sampleA = na_samples) %>%
      left_join(cov_lookup, by = c("sampleA" = "Sample"))
    message("\n--- D4: depth_med (medianCov) for missing sampleA from raw epiData ---")
    print_via_message(d4)
    n_below_floor <- sum(d4$medianCov_raw < MIN_COV_INCLUDE, na.rm = TRUE)
    n_cov_na      <- sum(is.na(d4$medianCov_raw))
    message(sprintf("  Below %dx coverage floor: %d", MIN_COV_INCLUDE, n_below_floor))
    message(sprintf("  Coverage NA in raw epiData: %d", n_cov_na))
  } else {
    message("\n--- D4: skipped (raw epiData not loadable from PATHS$data_raw) ---")
  }

  # ---- D5: would the *other* sample in the pair have given us patientId? --
  d5 <- pairs_with_na %>%
    mutate(sampleB_in_meta = sampleB %in% gh_meta$Sample) %>%
    left_join(gh_meta %>% select(Sample, pid_from_B = patientId),
              by = c("sampleB" = "Sample"))
  message("\n--- D5: Could sampleB rescue patientId for these pairs? ---")
  print_via_message(d5 %>% select(pair_id, sampleA, sampleB,
                                   sampleB_in_meta, pid_from_B))
  message(sprintf("Pairs rescuable from sampleB: %d / %d",
                  sum(d5$sampleB_in_meta), nrow(d5)))

  # ---- D6: parse-fallback feasibility -------------------------------------
  parse_pid <- function(Sample) str_extract(as.character(Sample), "^KBTH\\d+")
  d6 <- pairs_with_na %>%
    mutate(parsed_from_A = parse_pid(sampleA),
           parsed_from_B = parse_pid(sampleB),
           agree         = parsed_from_A == parsed_from_B)
  message("\n--- D6: patientId parsable directly from Sample ID? ---")
  print_via_message(d6 %>% select(pair_id, sampleA, sampleB,
                                   parsed_from_A, parsed_from_B, agree))
  n_parsable <- sum(!is.na(d6$parsed_from_A))
  n_agreeing <- sum(d6$agree, na.rm = TRUE)
  message(sprintf("Parsable from sampleA: %d / %d", n_parsable, nrow(d6)))
  message(sprintf("Parsed sampleA == sampleB: %d / %d", n_agreeing, nrow(d6)))

  # ---- D7: if sampleA missing from gh_meta, what about HIV? ---------------
  # If the underlying issue is coverage filter, HIV from sampleB would
  # rescue HIV the same way patientId could. Confirm.
  d7 <- pairs_with_na %>%
    left_join(gh_meta %>% select(Sample, hiv_from_B = HIV),
              by = c("sampleB" = "Sample"))
  message("\n--- D7: HIV rescuable from sampleB for the same pairs? ---")
  print_via_message(d7 %>% select(pair_id, sampleA, sampleB, hiv_from_B))
}

message("\n[diagnose_patientId_missing] Done.")
