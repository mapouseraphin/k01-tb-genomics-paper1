suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(tibble); library(purrr) })

v_primary <- readRDS(file.path(PATHS$thresholded, "gh_variants_thr_primary.rds"))
gh_meta   <- readRDS(file.path(PATHS$meta,        "gh_meta.rds"))
cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))

pt_lkp <- gh_meta %>% select(Sample, patientId) %>% distinct()

# Cohort-wide site recurrence: distinct patients per (CHROM, POS, REF, ALT)
site_recur <- v_primary %>%
  inner_join(pt_lkp, by = "Sample") %>%
  distinct(patientId, CHROM, POS, REF, ALT) %>%
  count(CHROM, POS, REF, ALT, name = "n_patients_with_call")

cat("\n--- Background: cohort-wide site recurrence (all primary calls) ---\n")
print(site_recur %>% count(n_patients_with_call) %>% arrange(desc(n_patients_with_call)))

# Restrict to evaluable pairs at Primary
prim <- list(DP_min = 40, AD1_min = 6, MAF_min = 0.02)
eval_pairs <- cal_pairs %>%
  filter(DP_min == prim$DP_min, AD1_min == prim$AD1_min, MAF_min == prim$MAF_min,
         nA > 0, nB > 0) %>%
  select(sampleA, sampleB)

# For each evaluable pair, get the intersection (the calls driving Jaccard)
key <- function(df) paste(df$CHROM, df$POS, df$REF, df$ALT, sep = "|")

inter_sites <- map_dfr(seq_len(nrow(eval_pairs)), function(i) {
  sa <- eval_pairs$sampleA[i]; sb <- eval_pairs$sampleB[i]
  shared <- intersect(
    key(filter(v_primary, Sample == sa)),
    key(filter(v_primary, Sample == sb))
  )
  if (!length(shared)) return(tibble())
  do.call(rbind, strsplit(shared, "\\|")) %>%
    as.data.frame() %>%
    setNames(c("CHROM", "POS", "REF", "ALT")) %>%
    mutate(POS = as.integer(POS), pair_idx = i, sampleA = sa, sampleB = sb)
}) %>%
  left_join(pt_lkp, by = c("sampleA" = "Sample")) %>%
  left_join(site_recur, by = c("CHROM", "POS", "REF", "ALT"))

cat("\n--- Recurrence of sites in within-pair INTERSECTIONS ---\n")
print(inter_sites %>% count(n_patients_with_call) %>% arrange(desc(n_patients_with_call)))

cat("\n--- Top 15 sites recurring across evaluable pairs' intersections ---\n")
print(inter_sites %>%
        count(CHROM, POS, REF, ALT, n_patients_with_call, name = "n_eval_pairs") %>%
        arrange(desc(n_eval_pairs), desc(n_patients_with_call)) %>%
        head(15))
