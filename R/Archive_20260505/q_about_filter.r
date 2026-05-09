# Clean session at project root.
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
setwd("~/Library/CloudStorage/OneDrive-UniversityofFlorida/Research/manuscripts/Submitted/within_between_host_isnv_calibration")

SELECTION_METHOD <- "lex"
DROP_PPE_CAL     <- FALSE
DROP_PPE_APPLY   <- FALSE
SNV_ONLY         <- TRUE
M0_ONLY          <- TRUE
source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(ggplot2)
})

# Load shared inputs once
cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))
boot_picks <- readRDS(file.path(PATHS$bootstrap, "bootstrap_picks.rds"))
gh_meta    <- readRDS(file.path(PATHS$meta, "gh_meta.rds"))

# Modal cell from the lex bootstrap
MODAL <- list(DP = 70, AD1 = 6, MAF = 0.05)

# Q2 — Is the observed cohort an outlier from the bootstrap distribution at the modal cell?
# Per-rep metrics at (70, 6, 0.05). Recomputes the rep-level metrics directly
# from cal_pairs to build the empirical bootstrap distribution at MODAL.
n_pairs <- nrow(distinct(cal_pairs, sampleA, sampleB))
B <- 1000
set.seed(0xcafe)

cp_at_modal <- cal_pairs %>%
  filter(DP_min == MODAL$DP, AD1_min == MODAL$AD1, MAF_min == MODAL$MAF) %>%
  inner_join(distinct(cal_pairs, sampleA, sampleB) %>%
               mutate(pair_id = row_number()),
             by = c("sampleA", "sampleB"))

cp_split <- split(cp_at_modal, cp_at_modal$pair_id)

per_rep <- vapply(seq_len(B), function(b) {
  ids <- sample(seq_len(n_pairs), size = n_pairs, replace = TRUE)
  bb <- bind_rows(cp_split[as.character(ids)])
  if (nrow(bb) == 0) return(c(NA, NA, NA))
  evals <- bb$nA > 0 & bb$nB > 0
  c(median(bb$jaccard[evals & bb$n_union > 0], na.rm = TRUE),
    median(bb$confirm_both[evals & bb$n_union > 0], na.rm = TRUE),
    sum(bb$nA == 0 & bb$nB == 0) / nrow(bb))
}, numeric(3))

per_rep <- t(per_rep) %>% as_tibble(.name_repair = ~c("J", "C", "p_both0")) %>%
  filter(!is.na(J))

# Observed values at MODAL
obs <- cal_pairs %>%
  filter(DP_min == MODAL$DP, AD1_min == MODAL$AD1, MAF_min == MODAL$MAF) %>%
  summarise(J = median(jaccard[nA > 0 & nB > 0 & n_union > 0], na.rm = TRUE),
            C = median(confirm_both[nA > 0 & nB > 0 & n_union > 0], na.rm = TRUE),
            p_both0 = mean(nA == 0 & nB == 0))

cat("\n--- Q2: observed vs bootstrap at MODAL (70, 6, 0.05) ---\n")
cat(sprintf("Observed J = %.3f   bootstrap median = %.3f   pctile = %.0f%%\n",
            obs$J, median(per_rep$J), 100 * mean(per_rep$J <= obs$J)))
cat(sprintf("Observed C = %.3f   bootstrap median = %.3f   pctile = %.0f%%\n",
            obs$C, median(per_rep$C), 100 * mean(per_rep$C <= obs$C)))
cat(sprintf("Observed p_both0 = %.3f   bootstrap median = %.3f   pctile = %.0f%%\n",
            obs$p_both0, median(per_rep$p_both0),
            100 * mean(per_rep$p_both0 <= obs$p_both0)))

#Q3 — What does (70, 6, 0.05) actually contain on the observed cohort?
modal_cell <- cal_pairs %>%
  filter(DP_min == MODAL$DP, AD1_min == MODAL$AD1, MAF_min == MODAL$MAF)

cat("\n--- Q3: observed-cohort detail at MODAL ---\n")
cat(sprintf("Total pairs:           %d\n", nrow(modal_cell)))
cat(sprintf("Evaluable (nA>0 & nB>0): %d\n",
            sum(modal_cell$nA > 0 & modal_cell$nB > 0)))
cat(sprintf("Both-zero pairs:       %d (%.1f%%)\n",
            sum(modal_cell$nA == 0 & modal_cell$nB == 0),
            100 * mean(modal_cell$nA == 0 & modal_cell$nB == 0)))

evals <- modal_cell %>% filter(nA > 0, nB > 0)
cat(sprintf("\nAmong evaluable (n = %d):\n", nrow(evals)))
cat(sprintf("  Jaccard      median = %.3f   IQR = [%.3f, %.3f]\n",
            median(evals$jaccard), quantile(evals$jaccard, 0.25),
            quantile(evals$jaccard, 0.75)))
cat(sprintf("  Confirm      median = %.3f   IQR = [%.3f, %.3f]\n",
            median(evals$confirm_both), quantile(evals$confirm_both, 0.25),
            quantile(evals$confirm_both, 0.75)))
cat(sprintf("  Pair n_union median = %.0f     IQR = [%.0f, %.0f]\n",
            median(evals$n_union), quantile(evals$n_union, 0.25),
            quantile(evals$n_union, 0.75)))

# Slot composition of evaluable pairs
slot_a <- parse_slot(evals$sampleA)
slot_b <- parse_slot(evals$sampleB)
cat("\nEvaluable-pair slot composition:\n")
print(table(slot_a, slot_b))
  

# Q4 — Does Sens A (PE/PPE excluded) keep AD1 = 6?
# In a fresh R session — this is a pipeline rerun, not a diagnostic.
rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
SELECTION_METHOD <- "lex"
DROP_PPE_CAL     <- TRUE     # the only difference from canonical
DROP_PPE_APPLY   <- TRUE
SNV_ONLY <- TRUE; M0_ONLY <- TRUE
source("R/00_pipeline_config.R")
source("R/01_build_cal_pairs.R")
source("R/02_calibrate_lexicographic.R")
source("R/02b_calibration_bootstrap.R")
# ~2 minutes total; only need 01, 02, 02b for this question

# Q5 — Are failed reps concentrated at specific cohort substructures?
# The 25% failure rate is invariant across methods; instrument it directly.
n_pairs <- nrow(distinct(cal_pairs, sampleA, sampleB))
all_pairs <- distinct(cal_pairs, sampleA, sampleB) %>% mutate(pair_id = row_number())
all_pairs$slot_a <- parse_slot(all_pairs$sampleA)
all_pairs$slot_b <- parse_slot(all_pairs$sampleB)
all_pairs$slot_pair <- pmin(all_pairs$slot_a, all_pairs$slot_b) %>%
  paste(pmax(all_pairs$slot_a, all_pairs$slot_b), sep = "+")

# Reuse the boot_picks failed flag
failed_reps <- which(boot_picks$failed)
ok_reps     <- which(!boot_picks$failed)

# Resample using the SAME seed to recover which pair_ids landed in each rep
set.seed(0xcafe)
rep_pair_ids <- replicate(1000,
                          sample(seq_len(n_pairs), size = n_pairs, replace = TRUE),
                          simplify = FALSE)

# Slot composition per rep
slot_dist_per_rep <- function(rep_idx) {
  ids <- rep_pair_ids[[rep_idx]]
  table(all_pairs$slot_pair[ids])
}

# Aggregate over failed vs ok
agg_slot_dist <- function(rep_indices) {
  Reduce(`+`, lapply(rep_indices, slot_dist_per_rep))
}

cat("\n--- Q5: slot composition: failed reps vs ok reps ---\n")
failed_counts <- agg_slot_dist(failed_reps)
ok_counts     <- agg_slot_dist(ok_reps)

# Normalize to per-rep average for direct comparison
cat("\nFailed reps (avg per rep, n =", length(failed_reps), "):\n")
print(round(failed_counts / length(failed_reps), 1))
cat("\nOK reps (avg per rep, n =", length(ok_reps), "):\n")
print(round(ok_counts / length(ok_reps), 1))

# Test which slot pair is over- or under-represented in failures
cmp <- tibble(
  slot_pair  = names(failed_counts),
  per_failed = as.numeric(failed_counts) / length(failed_reps),
  per_ok     = as.numeric(ok_counts)     / length(ok_reps)
) %>%
  mutate(diff_per_rep = per_failed - per_ok,
         pct_diff = 100 * (per_failed - per_ok) / per_ok) %>%
  arrange(desc(abs(pct_diff)))

cat("\nSlot pair over/under-representation in failures (sorted by |%diff|):\n")
print(cmp)

# Which guard binds in failed reps?
cp_split_full <- split(cal_pairs, list(cal_pairs$DP_min,
                                       cal_pairs$AD1_min,
                                       cal_pairs$MAF_min), drop = TRUE)

# For 50 failed reps (sampling for speed), recompute grid and see which guard binds
set.seed(123)
sample_failed <- sample(failed_reps, min(50, length(failed_reps)))

guard_binds <- lapply(sample_failed, function(r) {
  ids <- rep_pair_ids[[r]]
  cp_b <- cal_pairs %>%
    inner_join(all_pairs[, c("sampleA", "sampleB", "pair_id")],
               by = c("sampleA", "sampleB")) %>%
    filter(pair_id %in% ids)
  grid_b <- cp_b %>%
    group_by(DP_min, AD1_min, MAF_min) %>%
    summarise(prop_both0      = mean(nA == 0 & nB == 0),
              pairs_evaluable = sum(nA > 0 & nB > 0),
              .groups = "drop")
  tibble(rep = r,
         n_pass_prop_both0   = sum(grid_b$prop_both0 <= MAX_PROP_BOTH0),
         n_pass_evaluable    = sum(grid_b$pairs_evaluable >= MIN_EVAL_PAIRS),
         n_pass_both         = sum(grid_b$prop_both0 <= MAX_PROP_BOTH0 &
                                     grid_b$pairs_evaluable >= MIN_EVAL_PAIRS))
})
guard_binds <- bind_rows(guard_binds)

cat("\nIn 50 sampled failed reps:\n")
cat(sprintf("  Mean cells passing prop_both0 guard:   %.1f / %d\n",
            mean(guard_binds$n_pass_prop_both0), nrow(cal_pairs %>% distinct(DP_min, AD1_min, MAF_min))))
cat(sprintf("  Mean cells passing pairs_evaluable guard: %.1f\n",
            mean(guard_binds$n_pass_evaluable)))
cat(sprintf("  Mean cells passing both guards:        %.1f\n",
            mean(guard_binds$n_pass_both)))

