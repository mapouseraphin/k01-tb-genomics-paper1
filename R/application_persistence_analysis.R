# =============================================================================
# application_persistence_analysis.R   (NEW 2026-05-03)
#
# Application: longitudinal iSNV persistence in the Ghana cohort.
#
# DESIGN INTENT (per handoff §6.6, §6.8)
# --------------------------------------
# This is a DESCRIPTIVE application of the calibrated rule (Primary tier)
# to longitudinal samples. It is NOT an LCA-based inference. Outputs report
# observed within-patient persistence of iSNV calls across visits (M0, M1,
# M2), stratified by HIV status. No causal claims about HIV are made; the
# stratification is descriptive, the discussion text in Methods/Results
# treats it as a clinical descriptor co-occurring with persistence patterns,
# not as an exposure under any identifying model.
#
# PRIMARY DELIVERABLES
# --------------------
#   Table 5      : patient-level persistence counts and rates by HIV status
#   Supp Table S4: trajectory-category breakdown (M0-only, M0+M1 only,
#                   M0+M2 only, M0+M1+M2, etc.)
#   Supp Table S5: MAF-stratified persistence (low / mid / high MAF bins)
#   Figure 5     : per-patient longitudinal trajectory line plot
#
# UNIT OF ANALYSIS
# ----------------
#   - For tables: patient (one row per patientId).
#   - For trajectories: (patientId, CHROM, POS, REF, ALT) tuple presence
#     indicator across visits.
#
# A "persistent" iSNV is defined as a (CHROM, POS, REF, ALT) tuple that is
# called (under Primary thresholds) in M0 AND in at least one of {M1, M2}
# from the SAME patient. Persistence rates are reported as:
#   - patient-level rate: proportion of patients with >= 1 persistent iSNV
#   - per-iSNV rate     : among M0 iSNVs in patients with >= 1 post-M0
#                          sample, proportion that re-appear at M1 or M2
#
# CAVEATS (Methods text; reproduced here for transparency)
# --------------------------------------------------------
#   - Sample-set restriction: only patients with >= 1 post-M0 sample
#     contribute to per-iSNV rates. Patient-level rates use the full cohort
#     in the denominator with explicit no-post-M0 noted.
#   - Coverage heterogeneity: a "non-persistent" call may reflect coverage
#     drop-out at the same site in the post-M0 sample, not biological loss.
#     Coverage at the site is reported in Supp Table S5 footnotes.
#   - Slot heterogeneity: M1/M2 slot composition differs from M0; persistence
#     comparisons are within-patient but not within-slot.
#
# Inputs:
#   PATHS$thresholded/gh_variants_thr_primary.rds
#   PATHS$meta/gh_meta.rds
#
# Outputs (under PATHS$application; created if missing):
#   table5_persistence_by_hiv.csv
#   table5_persistence_by_hiv.tex
#   supp_table_S4_trajectory_categories.csv
#   supp_table_S5_persistence_by_maf_bin.csv
#   figure5_persistence_trajectories.pdf
#   persistence_log.txt
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(purrr)
  library(ggplot2); library(stringr)
})

# Application output dir is not in PATHS by default; resolve under outputs/
APP_PATH <- if (!is.null(PATHS$application)) {
  PATHS$application
} else {
  file.path("outputs", "application", TAG$apply_tag)
}
dir.create(APP_PATH, recursive = TRUE, showWarnings = FALSE)

LOG <- character()
log_msg <- function(...) {
  args <- list(...)
  s <- if (length(args) == 1) as.character(args[[1]]) else do.call(sprintf, args)
  message(s); LOG <<- c(LOG, s); invisible(s)
}

log_msg("[application_persistence_analysis] Spec tag: %s", TAG$apply_tag)
log_msg("Output dir: %s", APP_PATH)

# ---- Inputs ----------------------------------------------------------------
v_path <- file.path(PATHS$thresholded, "gh_variants_thr_primary.rds")
m_path <- file.path(PATHS$meta,        "gh_meta.rds")
if (!file.exists(v_path)) stop("Missing primary variants: ", v_path)
if (!file.exists(m_path)) stop("Missing gh_meta:           ", m_path)

variants <- readRDS(v_path)
meta     <- readRDS(m_path)

req_v_cols <- c("Sample", "CHROM", "POS", "REF", "ALT", "DP", "MAF")
req_m_cols <- c("Sample", "patientId", "visit_label", "HIV")
miss_v <- setdiff(req_v_cols, names(variants))
miss_m <- setdiff(req_m_cols, names(meta))
if (length(miss_v))
  stop("variants missing columns: ", paste(miss_v, collapse = ", "))
if (length(miss_m))
  stop("meta missing columns: ", paste(miss_m, collapse = ", "))

# ---- Cohort scoping --------------------------------------------------------
# All Ghana samples assayed under M0/M1/M2 (visit_label in {M0, M1, M2}).
# The persistence definition requires M0 plus >= 1 post-M0 sample per patient.
meta_ml <- meta %>%
  filter(visit_label %in% c("M0", "M1", "M2")) %>%
  mutate(visit_label = factor(visit_label, levels = c("M0", "M1", "M2")))

n_samples_total   <- nrow(meta_ml)
n_patients_total  <- n_distinct(meta_ml$patientId)
log_msg("\nCohort (visit_label in M0/M1/M2):")
log_msg("  Samples : %d", n_samples_total)
log_msg("  Patients: %d", n_patients_total)

visit_composition <- meta_ml %>% count(visit_label, name = "n_samples")
log_msg("\nSamples per visit:")
log_msg(paste(capture.output(print(as.data.frame(visit_composition))), collapse = "\n"))

# Patient-level visit coverage matrix (one row per patient, cols = M0/M1/M2 boolean)
pt_visits <- meta_ml %>%
  distinct(patientId, visit_label) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = visit_label, values_from = present,
              values_fill = list(present = FALSE)) %>%
  mutate(across(any_of(c("M0", "M1", "M2")), ~replace_na(., FALSE)))

# Ensure all three columns exist (some specs may have no M2 samples at all)
for (v in c("M0", "M1", "M2")) {
  if (!v %in% names(pt_visits)) pt_visits[[v]] <- FALSE
}

# Patients eligible for persistence: have M0 and >= 1 of {M1, M2}
pt_visits <- pt_visits %>%
  mutate(has_postM0 = M1 | M2,
         eligible_persistence = M0 & has_postM0)

# Attach HIV status (one stable status per patient; warn if it varies)
hiv_per_pt <- meta_ml %>%
  distinct(patientId, HIV) %>%
  count(patientId, HIV) %>%
  group_by(patientId) %>%
  arrange(desc(n)) %>%
  slice(1) %>%
  ungroup() %>%
  select(patientId, HIV)

n_var_hiv <- meta_ml %>%
  distinct(patientId, HIV) %>%
  count(patientId) %>%
  filter(n > 1) %>% nrow()
if (n_var_hiv > 0) {
  log_msg("WARNING: %d patient(s) have varying HIV across visits; using modal status.",
          n_var_hiv)
}

pt_visits <- pt_visits %>% left_join(hiv_per_pt, by = "patientId")

n_eligible <- sum(pt_visits$eligible_persistence)
log_msg("\nEligibility for persistence analysis:")
log_msg("  Patients with M0:                    %d",
        sum(pt_visits$M0))
log_msg("  Patients with M0 AND >= 1 post-M0:   %d (eligible)", n_eligible)
log_msg("  Patients without post-M0 sample:     %d (excluded from per-iSNV rates)",
        sum(pt_visits$M0 & !pt_visits$has_postM0))

if (n_eligible == 0) {
  log_msg("\nNO PATIENTS ELIGIBLE for persistence analysis. Writing log only.")
  writeLines(LOG, file.path(APP_PATH, "persistence_log.txt"))
  invisible(return(NULL))
}

# ---- Build per-(patient, position) presence matrix -------------------------
# Join variants to meta to attach patientId and visit_label
v_with_visit <- variants %>% select(-visit_label)%>%
  inner_join(meta_ml %>% select(Sample, patientId, visit_label, HIV),
             by = "Sample") %>%
  mutate(site_key = paste(CHROM, POS, REF, ALT, sep = ":"))

# Per-(patient, site) presence indicator across visits
presence <- v_with_visit %>%
  distinct(patientId, HIV, site_key, CHROM, POS, REF, ALT, visit_label) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = visit_label, values_from = present,
              values_fill = list(present = FALSE),
              names_prefix = "at_") %>%
  mutate(across(any_of(c("at_M0", "at_M1", "at_M2")),
                ~replace_na(., FALSE)))
for (v in c("at_M0", "at_M1", "at_M2")) {
  if (!v %in% names(presence)) presence[[v]] <- FALSE
}

# Restrict per-iSNV analysis to eligible patients (M0 + >=1 post-M0)
elig_pts <- pt_visits %>% filter(eligible_persistence) %>% pull(patientId)
presence_elig <- presence %>% filter(patientId %in% elig_pts)

# Trajectory category per (patient, site). For sites called at M0 only, we
# include all eligible-patient sites (so M0-only is a real category).
m0_sites_elig <- presence_elig %>% filter(at_M0)

m0_sites_elig <- m0_sites_elig %>%
  mutate(persistent_M1     = at_M1,
         persistent_M2     = at_M2,
         persistent_either = at_M1 | at_M2,
         category = case_when(
           at_M0 &  at_M1 &  at_M2 ~ "M0 + M1 + M2",
           at_M0 &  at_M1 & !at_M2 ~ "M0 + M1 only",
           at_M0 & !at_M1 &  at_M2 ~ "M0 + M2 only",
           at_M0 & !at_M1 & !at_M2 ~ "M0 only",
           TRUE                    ~ "other"
         ))

# ---- Patient-level summary -------------------------------------------------
pt_persist <- m0_sites_elig %>%
  group_by(patientId, HIV) %>%
  summarise(
    n_M0_isnv             = n(),
    n_persistent_isnv     = sum(persistent_either),
    persistence_rate_isnv = ifelse(n_M0_isnv > 0,
                                    n_persistent_isnv / n_M0_isnv, NA_real_),
    has_persistent        = n_persistent_isnv > 0,
    .groups = "drop"
  )

# Patients eligible but with zero M0 iSNVs at all (no calls passing primary
# at M0) need to be represented too.
no_m0_isnv_pts <- setdiff(elig_pts, pt_persist$patientId)
if (length(no_m0_isnv_pts) > 0) {
  add_rows <- tibble(patientId = no_m0_isnv_pts) %>%
    left_join(hiv_per_pt, by = "patientId") %>%
    mutate(n_M0_isnv = 0L, n_persistent_isnv = 0L,
           persistence_rate_isnv = NA_real_, has_persistent = FALSE)
  pt_persist <- bind_rows(pt_persist, add_rows)
}

# ---- Table 5: persistence by HIV ------------------------------------------
build_strat_row <- function(df, label) {
  if (nrow(df) == 0) {
    return(tibble(stratum = label,
                  n_patients = 0L,
                  n_with_persistent = 0L,
                  pct_with_persistent = NA_real_,
                  median_persistence_rate = NA_real_,
                  total_M0_isnv = 0L,
                  total_persistent_isnv = 0L,
                  pooled_per_isnv_rate = NA_real_))
  }
  tibble(
    stratum                 = label,
    n_patients              = nrow(df),
    n_with_persistent       = sum(df$has_persistent),
    pct_with_persistent     = mean(df$has_persistent),
    median_persistence_rate = median(df$persistence_rate_isnv, na.rm = TRUE),
    total_M0_isnv           = sum(df$n_M0_isnv),
    total_persistent_isnv   = sum(df$n_persistent_isnv),
    pooled_per_isnv_rate    = ifelse(sum(df$n_M0_isnv) > 0,
                                     sum(df$n_persistent_isnv) / sum(df$n_M0_isnv),
                                     NA_real_)
  )
}

table5 <- bind_rows(
  build_strat_row(pt_persist,                                   "All eligible"),
  build_strat_row(pt_persist %>% filter(HIV == "HIV-"),         "HIV-"),
  build_strat_row(pt_persist %>% filter(HIV == "HIV+"),         "HIV+"),
  build_strat_row(pt_persist %>% filter(is.na(HIV)),            "HIV unknown")
) %>%
  filter(n_patients > 0) %>%
  mutate(
    pct_with_persistent     = round(100 * pct_with_persistent, 1),
    median_persistence_rate = round(median_persistence_rate, 3),
    pooled_per_isnv_rate    = round(pooled_per_isnv_rate, 3)
  )

write_csv(table5, file.path(APP_PATH, "table5_persistence_by_hiv.csv"))

# LaTeX rendering
build_latex_table5 <- function(t5) {
  hdr <- c("Stratum", "$n$", "$n$ persistent", "\\% persistent",
           "Median per-iSNV rate", "Total M0 iSNVs",
           "Total persistent", "Pooled per-iSNV rate")
  body <- apply(t5, 1, function(r) paste(r, collapse = " & "))
  paste(c(
    "\\begin{tabular}{lrrrrrrr}",
    "\\hline",
    paste(paste(hdr, collapse = " & "), "\\\\"),
    "\\hline",
    paste(body, "\\\\"),
    "\\hline",
    "\\end{tabular}"
  ), collapse = "\n")
}
writeLines(build_latex_table5(table5),
           file.path(APP_PATH, "table5_persistence_by_hiv.tex"))

log_msg("\nTable 5 (persistence by HIV):")
log_msg(paste(capture.output(print(as.data.frame(table5))), collapse = "\n"))

# ---- Supp Table S4: trajectory categories ---------------------------------
trajectory_summary <- m0_sites_elig %>%
  count(HIV, category) %>%
  group_by(HIV) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(HIV, factor(category, levels = c("M0 only", "M0 + M1 only",
                                            "M0 + M2 only", "M0 + M1 + M2",
                                            "other"))) %>%
  mutate(prop = round(prop, 3))

# Add an "All" row block
trajectory_all <- m0_sites_elig %>%
  count(category) %>%
  mutate(HIV = "All eligible",
         prop = round(n / sum(n), 3)) %>%
  select(HIV, category, n, prop)

trajectory_out <- bind_rows(trajectory_all, trajectory_summary)
write_csv(trajectory_out,
          file.path(APP_PATH, "supp_table_S4_trajectory_categories.csv"))

log_msg("\nSupp Table S4 (trajectory categories):")
log_msg(paste(capture.output(print(as.data.frame(trajectory_out))), collapse = "\n"))

# ---- Supp Table S5: MAF-stratified persistence -----------------------------
# Anchor MAF at the M0 call (where the persistence question begins).
maf_at_M0 <- v_with_visit %>%
  filter(visit_label == "M0", patientId %in% elig_pts) %>%
  mutate(site_key = paste(CHROM, POS, REF, ALT, sep = ":")) %>%
  group_by(patientId, site_key) %>%
  summarise(MAF_M0 = max(MAF, na.rm = TRUE), .groups = "drop")

m0_sites_with_maf <- m0_sites_elig %>%
  left_join(maf_at_M0, by = c("patientId", "site_key")) %>%
  mutate(maf_bin = case_when(
    is.na(MAF_M0)         ~ NA_character_,
    MAF_M0 < 0.05         ~ "[0.02, 0.05)",
    MAF_M0 < 0.10         ~ "[0.05, 0.10)",
    MAF_M0 < 0.25         ~ "[0.10, 0.25)",
    MAF_M0 <= 0.45        ~ "[0.25, 0.45]",
    TRUE                  ~ NA_character_
  ),
  maf_bin = factor(maf_bin,
                   levels = c("[0.02, 0.05)", "[0.05, 0.10)",
                              "[0.10, 0.25)", "[0.25, 0.45]")))

table_S5 <- m0_sites_with_maf %>%
  filter(!is.na(maf_bin)) %>%
  group_by(maf_bin, HIV) %>%
  summarise(
    n_M0_isnv         = n(),
    n_persistent      = sum(persistent_either),
    per_isnv_rate     = ifelse(n_M0_isnv > 0, n_persistent / n_M0_isnv, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(per_isnv_rate = round(per_isnv_rate, 3)) %>%
  arrange(HIV, maf_bin)

# Add overall (HIV-pooled) rows per MAF bin
table_S5_all <- m0_sites_with_maf %>%
  filter(!is.na(maf_bin)) %>%
  group_by(maf_bin) %>%
  summarise(HIV = "All eligible",
            n_M0_isnv = n(),
            n_persistent = sum(persistent_either),
            per_isnv_rate = round(ifelse(n() > 0, sum(persistent_either) / n(), NA_real_), 3),
            .groups = "drop") %>%
  select(maf_bin, HIV, n_M0_isnv, n_persistent, per_isnv_rate)

table_S5_out <- bind_rows(table_S5_all, table_S5)
write_csv(table_S5_out,
          file.path(APP_PATH, "supp_table_S5_persistence_by_maf_bin.csv"))

log_msg("\nSupp Table S5 (persistence by MAF bin):")
log_msg(paste(capture.output(print(as.data.frame(table_S5_out))), collapse = "\n"))

# ---- Figure 5: per-patient trajectory line plot ----------------------------
# For each eligible patient, count iSNVs per visit. One line per patient,
# colored by HIV. Lines connect counts across M0 / M1 / M2.
counts_per_visit <- v_with_visit %>%
  filter(patientId %in% elig_pts) %>%
  count(patientId, HIV, visit_label, name = "n_isnv") %>%
  complete(patientId, visit_label = c("M0", "M1", "M2"),
           fill = list(n_isnv = 0L)) %>%
  left_join(hiv_per_pt, by = "patientId", suffix = c("", ".pt")) %>%
  mutate(HIV = if (!"HIV" %in% names(.)) HIV.pt else coalesce(HIV, HIV.pt)) %>%
  select(-any_of("HIV.pt")) %>%
  # Patients only contribute a visit point if they actually have a sample
  # at that visit (avoid spurious 0s for missed visits).
  semi_join(meta_ml %>% select(patientId, visit_label),
            by = c("patientId", "visit_label")) %>%
  mutate(visit_label = factor(visit_label, levels = c("M0", "M1", "M2")))

fig5 <- ggplot(counts_per_visit,
               aes(x = visit_label, y = n_isnv,
                   group = patientId, color = HIV)) +
  geom_line(alpha = 0.5, linewidth = 0.5) +
  geom_point(alpha = 0.7, size = 1.6) +
  scale_color_manual(values = c("HIV-" = "#0072B2", "HIV+" = "#D55E00"),
                     na.value = "gray60") +
  labs(x = "Visit", y = "Number of iSNVs (Primary)",
       title = "Per-patient iSNV trajectories",
       subtitle = sprintf("Eligible patients (M0 + >= 1 post-M0 sample); n = %d",
                          n_eligible),
       caption = paste(
         "Each line is one patient. Calls passing the Primary calibrated rule.",
         "Stratification by HIV is descriptive; no causal claims attached.",
         "Slot heterogeneity across visits not controlled for.",
         sep = "\n")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 8, color = "gray30"))

ggsave(file.path(APP_PATH, "figure5_persistence_trajectories.pdf"),
       fig5, width = 7, height = 5)

log_msg("\n[application_persistence_analysis] Outputs written to %s/", APP_PATH)
log_msg("  table5_persistence_by_hiv.csv  + .tex")
log_msg("  supp_table_S4_trajectory_categories.csv")
log_msg("  supp_table_S5_persistence_by_maf_bin.csv")
log_msg("  figure5_persistence_trajectories.pdf")

writeLines(LOG, file.path(APP_PATH, "persistence_log.txt"))
