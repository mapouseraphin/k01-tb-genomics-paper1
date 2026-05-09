# ============================================================================
# Helper: build_sample_metadata()
# ----------------------------------------------------------------------------
# Cleans and parses `Sample` IDs, then enumerates same-visit replicate pairs.
# This is kept separate from the calibration function so we can validate and
# tweak parsing independently.
#
# Expected Sample pattern (after cleaning):
#   ^[^-]+-[SE]\d+(?:\.\d+)?$
# Examples:
#   KBTH001-E0   (early morning, baseline)
#   KBTH001-S0.1 (spot 1, baseline)
#   KBTH001-S0.2 (spot 2, baseline)
#   KBTH001-E1   (early morning, month 1)
#   KBTH001-S1.1 (spot 1, month 1)
#
# Special cases handled:
# - Drop trailing technical suffices: .a / .b (e.g., KBTH001-S0.2.a)
# - Fix early-morning wrongly suffixed with .1: -E0.1 -> -E0
# - Rare stray leading hyphen: -E0.1 -> E0 (will likely fail pattern and be flagged)
#
# Returns a list:
# - samples_clean: cleaned unique Sample IDs + original
# - meta: parsed fields (person_id, month_num, sample_type, spot, labels)
# - pairs_same_visit: all within-person, within-visit pairs (sampleA, sampleB)
# - variants_annot: input variants joined to `meta`
# - bad_rows: rows that failed parsing (for manual review)
# - dups: duplicates introduced by cleaning (for manual review)
#
# Dependencies: dplyr, stringr, tidyr, purrr, tibble
# ============================================================================
build_sample_metadata <- function(
    variants,
    sample_col      = "Sample",
    enforce_pattern = FALSE,
    pattern         = "^[^-]+-[SE]\\d+(?:\\.\\d+)?$",
    verbose         = TRUE
) {
  stopifnot(sample_col %in% names(variants))
  
  suppressPackageStartupMessages({
    require(dplyr, quietly = TRUE)
    require(stringr, quietly = TRUE)
    require(tidyr, quietly = TRUE)
    require(purrr, quietly = TRUE)
    require(tibble, quietly = TRUE)
  })
  `%>%` <- dplyr::`%>%`
  
  # 0) Distinct sample IDs
  samples_raw <- variants %>%
    dplyr::distinct(.data[[sample_col]]) %>%
    dplyr::filter(!is.na(.data[[sample_col]])) %>%
    dplyr::transmute(Sample = stringr::str_trim(.data[[sample_col]]))
  
  # 1) Sanitize
  samples_clean <- samples_raw %>%
    dplyr::mutate(Sample_original = Sample) %>%
    dplyr::filter(!stringr::str_detect(Sample, "\\.(?:a|b)$")) %>%                 # drop .a/.b
    dplyr::mutate(Sample = stringr::str_replace(Sample, "-E([0-9]+)\\.\\d+$", "-E\\1")) %>%  # E0.1 -> E0
    dplyr::mutate(Sample = stringr::str_replace(Sample, "^\\-E([0-9]+)\\.\\d+$", "E\\1")) %>% # "-E0.1" -> "E0"
    dplyr::distinct(Sample, Sample_original, .keep_all = TRUE)
  
  dups <- samples_clean %>% dplyr::count(Sample) %>% dplyr::filter(.data$n > 1)
  if (nrow(dups) > 0 && verbose) {
    warning("Duplicate Sample IDs created by cleaning. See `dups` in the return value.")
  }
  
  # 2) Parse fields
  meta <- samples_clean %>%
    dplyr::mutate(
      person_id   = stringr::str_extract(Sample, "^[^-]+"),
      tag         = stringr::str_extract(Sample, "(?<=-)[A-Za-z]\\d+(?:\\.\\d+)?"),
      sample_type = stringr::str_to_upper(stringr::str_sub(tag, 1, 1)),
      month_num   = suppressWarnings(as.integer(stringr::str_extract(tag, "\\d+"))),
      spot        = dplyr::if_else(sample_type == "S",
                                   suppressWarnings(as.integer(stringr::str_extract(tag, "(?<=\\.)\\d+"))),
                                   as.integer(NA_integer_))
    )
  
  meta <- if (enforce_pattern) {
    meta %>% dplyr::mutate(pattern_ok = stringr::str_detect(Sample, pattern))
  } else {
    meta %>% dplyr::mutate(pattern_ok = NA)
  }
  
  bad_rows <- meta %>% dplyr::filter(is.na(person_id) | is.na(tag) | is.na(month_num))
  if (nrow(bad_rows) > 0 && verbose) {
    warning("Some rows failed parsing. Inspect `bad_rows` in the return value.")
  }
  
  # Human-friendly labels
  meta <- meta %>%
    dplyr::mutate(
      visit_label = dplyr::case_when(
        !is.na(month_num) & month_num == 0 ~ "baseline",
        !is.na(month_num) & month_num == 1 ~ "month1",
        TRUE ~ paste0("visit", month_num)
      ),
      replicate_label = dplyr::case_when(
        sample_type == "E" ~ "E",
        sample_type == "S" & spot == 1 ~ "S1",
        sample_type == "S" & spot == 2 ~ "S2",
        TRUE ~ NA_character_
      )
    )
  
  # (unchanged) same-visit replicate pairs (biological replicates)
  pairs_same_visit <- meta %>%
    dplyr::group_by(person_id, month_num) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::summarise(
      pair_df = list(as.data.frame(t(combn(Sample, 2)), stringsAsFactors = FALSE)),
      .groups = "drop"
    ) %>%
    tidyr::unnest(pair_df) %>%
    dplyr::rename(sampleA = V1, sampleB = V2)
  
  # 4) Join annotations back to the variants
  variants_annot <- variants %>%
    dplyr::inner_join(
      meta %>% dplyr::select(Sample, person_id, month_num, visit_label,
                             sample_type, spot, replicate_label, Sample_original),
      by = setNames("Sample", sample_col)
    )
  
  list(
    samples_clean    = samples_clean,
    meta             = meta,     # now includes Sample_original
    pairs_same_visit = pairs_same_visit,
    variants_annot   = variants_annot,
    bad_rows         = bad_rows,
    dups             = dups
  )
}
