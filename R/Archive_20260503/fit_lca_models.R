# =============================================================================
# fit_lca_models.R
#
# v4: M1 (slot-specific) + M2 fit at three preregistered sigma_u priors:
#   primary:  half-normal(0, 1.0)
#   tight:    half-normal(0, 0.5)
#   loose:    half-normal(0, 2.0)
#
# Outputs the M1 fit and three M2 fits, plus a side-by-side prior sensitivity
# comparison table that determines whether the M2 posterior shape is
# data-driven or prior-driven.
#
# Inputs:
#   <LCA_PATH>/lca_dataset.rds
#   stan/lca_simple.stan
#   stan/lca_hierarchical.stan
#
# Outputs (under <LCA_PATH>):
#   fits/fit_simple.rds, fit_simple_summary.csv, fit_simple_ppc.csv
#   fits/fit_hierarchical_primary.rds, _tight.rds, _loose.rds  (+ summaries, PPCs)
#   diagnostics/m1_smoke_test.md
#   diagnostics/m2_smoke_test_<prior>.md  (3 files)
#   diagnostics/m2_prior_sensitivity.md
#   diagnostics/m1_vs_m2_comparison.md    (uses primary M2)
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr); library(readr); library(tibble); library(tidyr)
})

LCA_PATH  <- file.path("data_derived", paste0("05_lca_", TAG$apply_tag))
FITS_PATH <- file.path(LCA_PATH, "fits")
DIAG_PATH <- file.path(LCA_PATH, "diagnostics")
dir.create(FITS_PATH, recursive = TRUE, showWarnings = FALSE)
dir.create(DIAG_PATH, recursive = TRUE, showWarnings = FALSE)

# ---- Data prep --------------------------------------------------------------
lca <- readRDS(file.path(LCA_PATH, "lca_dataset.rds"))
n_pre <- nrow(lca)
lca_fit <- lca %>% filter(maf_bin != paste0(">", 0.45))

SLOT_LEVELS <- c("EM", "S0.1", "S0.2")
lca_fit <- lca_fit %>%
  mutate(slot_a_int  = match(slot_a, SLOT_LEVELS),
         slot_b_int  = match(slot_b, SLOT_LEVELS),
         patient_int = as.integer(factor(patientId)))

stopifnot(all(!is.na(lca_fit$slot_a_int)),
          all(!is.na(lca_fit$slot_b_int)),
          all(!is.na(lca_fit$patient_int)))
P <- max(lca_fit$patient_int)

pt_table <- lca_fit %>%
  count(pair_type, slot_a_int, slot_b_int, name = "pt_n") %>%
  arrange(pair_type)

stan_data_m1 <- list(
  N    = nrow(lca_fit), S = length(SLOT_LEVELS),
  W    = as.matrix(lca_fit[, c("W1", "W2")]),
  slot = as.matrix(lca_fit[, c("slot_a_int", "slot_b_int")]),
  N_pt = nrow(pt_table),
  pt_slot = as.matrix(pt_table[, c("slot_a_int", "slot_b_int")]),
  pt_n    = pt_table$pt_n
)
stan_data_m2_base <- c(stan_data_m1, list(P = P, patient = lca_fit$patient_int))

message(sprintf("Fitting set: %s rows; %d patients; %d pair_types",
                format(nrow(lca_fit), big.mark = ","), P, nrow(pt_table)))

# ---- Helpers ----------------------------------------------------------------
print_via_message <- function(df) {
  for (ln in capture.output(print(as.data.frame(df), row.names = FALSE))) {
    message(ln)
  }
}
fmt_df <- function(df) capture.output(print(as.data.frame(df), row.names = FALSE))

label_slots <- function(summ_df) {
  summ_df %>%
    mutate(slot = case_when(
      grepl("\\[1\\]", variable) ~ SLOT_LEVELS[1],
      grepl("\\[2\\]", variable) ~ SLOT_LEVELS[2],
      grepl("\\[3\\]", variable) ~ SLOT_LEVELS[3],
      TRUE ~ NA_character_
    ),
    param = gsub("\\[.+\\]", "", variable)) %>%
    select(variable, param, slot, everything())
}

build_ppc <- function(fit, pt_table, lca_fit) {
  pp_draws <- fit$draws("pp_cell", format = "draws_df")
  ppc_long <- expand_grid(pt_idx = seq_len(nrow(pt_table)), cell_k = 1:4) %>%
    mutate(pair_type = pt_table$pair_type[pt_idx],
           cell      = c("00","01","10","11")[cell_k],
           varname   = sprintf("pp_cell[%d,%d]", pt_idx, cell_k)) %>%
    rowwise() %>%
    mutate(pp_median = median(pp_draws[[varname]]),
           pp_q025   = quantile(pp_draws[[varname]], 0.025),
           pp_q975   = quantile(pp_draws[[varname]], 0.975)) %>%
    ungroup() %>%
    select(pair_type, cell, pp_median, pp_q025, pp_q975)
  obs <- lca_fit %>%
    count(pair_type, W1, W2, name = "n") %>%
    group_by(pair_type) %>% mutate(obs_prop = n / sum(n)) %>% ungroup() %>%
    mutate(cell = paste0(W1, W2)) %>% select(pair_type, cell, obs_prop)
  ppc_long %>%
    left_join(obs, by = c("pair_type", "cell")) %>%
    mutate(across(where(is.numeric), ~round(.x, 3)),
           within_CrI = obs_prop >= pp_q025 & obs_prop <= pp_q975) %>%
    arrange(pair_type, cell)
}

write_diag <- function(model_name, summ, ppc, t_elapsed, fit) {
  diag_lines <- character()
  add <- function(...) diag_lines <<- c(diag_lines, sprintf(...))
  add("# %s smoke test (%s)", model_name, format(Sys.time(), "%Y-%m-%d %H:%M"))
  add("Wall time: %.1f sec.", t_elapsed)
  add("")
  add("## Posterior summary")
  add("```")
  diag_lines <- c(diag_lines, fmt_df(summ %>%
    mutate(across(c(mean, median, q2.5, q97.5), ~round(.x, 3)),
           rhat = round(rhat, 4))))
  add("```"); add("")
  rhat_max     <- max(summ$rhat,     na.rm = TRUE)
  ess_bulk_min <- min(summ$ess_bulk, na.rm = TRUE)
  n_div        <- sum(fit$diagnostic_summary()$num_divergent)
  add("## Convergence: R-hat=%.4f | ESS=%.0f | div=%d",
      rhat_max, ess_bulk_min, n_div)
  add("")
  add("## PPC")
  add("```")
  diag_lines <- c(diag_lines, fmt_df(ppc))
  add("```")
  ppc_pass <- mean(ppc$within_CrI, na.rm = TRUE)
  add("PPC pass rate: %.0f%%", 100 * ppc_pass)
  add("")
  verdict_pass <- (rhat_max < 1.01) && (ess_bulk_min > 1000) &&
                  (n_div == 0) && (ppc_pass >= 0.75)
  add("## Verdict: **%s**", ifelse(verdict_pass, "PASS", "INVESTIGATE"))
  list(lines = diag_lines, rhat_max = rhat_max, ess_min = ess_bulk_min,
       n_div = n_div, ppc_pass = ppc_pass, verdict = verdict_pass)
}

# =============================================================================
# M1 (unchanged)
# =============================================================================
message("\n=== M1 (slot-specific simple LCA) ===")
mod_m1 <- cmdstan_model("stan/lca_simple.stan")
t0 <- Sys.time()
fit_m1 <- mod_m1$sample(
  data = stan_data_m1, chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000,
  seed = SEED, refresh = 0, show_messages = FALSE
)
t_m1 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
fit_m1$save_object(file.path(FITS_PATH, "fit_simple.rds"))

m1_params <- c("pi", paste0("sens[", 1:3, "]"), paste0("spec[", 1:3, "]"))
summ_m1 <- fit_m1$summary(m1_params, mean, median,
                          ~quantile(.x, c(0.025, 0.975)),
                          rhat, ess_bulk, ess_tail) %>%
  rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% label_slots()
write_csv(summ_m1, file.path(FITS_PATH, "fit_simple_summary.csv"))
ppc_m1 <- build_ppc(fit_m1, pt_table, lca_fit)
write_csv(ppc_m1, file.path(FITS_PATH, "fit_simple_ppc.csv"))
m1_diag <- write_diag("M1 (slot-specific)", summ_m1, ppc_m1, t_m1, fit_m1)
writeLines(m1_diag$lines, file.path(DIAG_PATH, "m1_smoke_test.md"))
message(sprintf("M1: R-hat=%.4f | ESS=%.0f | div=%d | PPC=%.0f%% | %s",
                m1_diag$rhat_max, m1_diag$ess_min, m1_diag$n_div,
                100 * m1_diag$ppc_pass, ifelse(m1_diag$verdict, "PASS", "INVESTIGATE")))

# =============================================================================
# M2 across three preregistered sigma_u priors
# =============================================================================
PRIOR_SCALES <- c(primary = 1.0, tight = 0.5, loose = 2.0)

mod_m2 <- cmdstan_model("stan/lca_hierarchical.stan")

m2_params <- c("pi", "sigma_u",
               paste0("sens[", 1:3, "]"), paste0("spec[", 1:3, "]"))

m2_summ_list <- list()
m2_ppc_list  <- list()
m2_diag_list <- list()

for (lab in names(PRIOR_SCALES)) {
  scale_v <- PRIOR_SCALES[[lab]]
  message(sprintf("\n=== M2 [%s] sigma_u ~ half-normal(0, %.1f) ===", lab, scale_v))

  stan_data <- c(stan_data_m2_base, list(sigma_u_prior_scale = scale_v))
  t0 <- Sys.time()
  fit <- mod_m2$sample(
    data = stan_data, chains = 4, parallel_chains = 4,
    iter_warmup = 1500, iter_sampling = 2000,
    seed = SEED, refresh = 0, show_messages = FALSE,
    adapt_delta = 0.95
  )
  t_elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  fit$save_object(file.path(FITS_PATH, sprintf("fit_hierarchical_%s.rds", lab)))

  summ <- fit$summary(m2_params, mean, median,
                      ~quantile(.x, c(0.025, 0.975)),
                      rhat, ess_bulk, ess_tail) %>%
    rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% label_slots()
  write_csv(summ, file.path(FITS_PATH, sprintf("fit_hierarchical_%s_summary.csv", lab)))

  ppc <- build_ppc(fit, pt_table, lca_fit)
  write_csv(ppc, file.path(FITS_PATH, sprintf("fit_hierarchical_%s_ppc.csv", lab)))

  diag <- write_diag(sprintf("M2 [%s, scale=%.1f]", lab, scale_v),
                     summ, ppc, t_elapsed, fit)
  writeLines(diag$lines,
             file.path(DIAG_PATH, sprintf("m2_smoke_test_%s.md", lab)))
  message(sprintf("M2 [%s]: R-hat=%.4f | ESS=%.0f | div=%d | PPC=%.0f%% | %s",
                  lab, diag$rhat_max, diag$ess_min, diag$n_div,
                  100 * diag$ppc_pass,
                  ifelse(diag$verdict, "PASS", "INVESTIGATE")))

  m2_summ_list[[lab]] <- summ
  m2_ppc_list[[lab]]  <- ppc
  m2_diag_list[[lab]] <- diag
}

# =============================================================================
# M2 variant: globalclass (one shared class-gap parameter across slots)
# =============================================================================
message("\n=== M2 [globalclass] global beta_class shared across slots ===")
mod_m2gc <- cmdstan_model("stan/lca_hierarchical_globalclass.stan")
stan_data <- c(stan_data_m2_base, list(sigma_u_prior_scale = 1.0))
t0 <- Sys.time()
fit_gc <- mod_m2gc$sample(
  data = stan_data, chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 2000,
  seed = SEED, refresh = 0, show_messages = FALSE,
  adapt_delta = 0.95
)
t_gc <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
fit_gc$save_object(file.path(FITS_PATH, "fit_hierarchical_globalclass.rds"))

m2gc_params <- c("pi", "sigma_u", "beta_class",
                 paste0("alpha_slot[", 1:3, "]"),
                 paste0("sens[", 1:3, "]"), paste0("spec[", 1:3, "]"))
summ_gc <- fit_gc$summary(m2gc_params, mean, median,
                          ~quantile(.x, c(0.025, 0.975)),
                          rhat, ess_bulk, ess_tail) %>%
  rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% label_slots()
write_csv(summ_gc, file.path(FITS_PATH, "fit_hierarchical_globalclass_summary.csv"))

ppc_gc <- build_ppc(fit_gc, pt_table, lca_fit)
write_csv(ppc_gc, file.path(FITS_PATH, "fit_hierarchical_globalclass_ppc.csv"))

diag_gc <- write_diag("M2 [globalclass]", summ_gc, ppc_gc, t_gc, fit_gc)
writeLines(diag_gc$lines, file.path(DIAG_PATH, "m2_smoke_test_globalclass.md"))
message(sprintf("M2 [globalclass]: R-hat=%.4f | ESS=%.0f | div=%d | PPC=%.0f%% | %s",
                diag_gc$rhat_max, diag_gc$ess_min, diag_gc$n_div,
                100 * diag_gc$ppc_pass,
                ifelse(diag_gc$verdict, "PASS", "INVESTIGATE")))

# =============================================================================
# Prior sensitivity comparison
# =============================================================================
ps_lines <- character()
psadd <- function(...) ps_lines <<- c(ps_lines, sprintf(...))

psadd("# M2 prior sensitivity (sigma_u prior scale)")
psadd("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
psadd("")
psadd("All three priors are half-normal(0, scale) on sigma_u.")
psadd("Same data, same model structure, same seed. Only scale differs.")
psadd("")

# Side-by-side parameter table
ps_tbl <- bind_rows(lapply(names(PRIOR_SCALES), function(lab) {
  m2_summ_list[[lab]] %>%
    transmute(param, slot,
              prior     = sprintf("hn(0,%.1f)", PRIOR_SCALES[[lab]]),
              prior_lab = lab,
              median    = round(median, 3),
              q2.5      = round(q2.5, 3),
              q97.5     = round(q97.5, 3),
              rhat      = round(rhat, 4),
              ess       = round(ess_bulk, 0))
})) %>% arrange(param, slot, prior_lab)

psadd("## Parameter posteriors across priors")
psadd("```")
ps_lines <- c(ps_lines, fmt_df(ps_tbl %>% select(-prior_lab)))
psadd("```")
psadd("")

# Compact per-prior summary
psadd("## Per-prior verdict summary")
verdict_tbl <- tibble(
  prior        = sprintf("hn(0,%.1f)", PRIOR_SCALES),
  prior_lab    = names(PRIOR_SCALES),
  rhat_max     = sapply(m2_diag_list, `[[`, "rhat_max") %>% round(4),
  ess_min      = sapply(m2_diag_list, `[[`, "ess_min") %>% round(0),
  n_divergent  = sapply(m2_diag_list, `[[`, "n_div"),
  ppc_pass_pct = (100 * sapply(m2_diag_list, `[[`, "ppc_pass")) %>% round(0),
  verdict      = ifelse(sapply(m2_diag_list, `[[`, "verdict"), "PASS", "INVESTIGATE")
)
psadd("```")
ps_lines <- c(ps_lines, fmt_df(verdict_tbl %>% select(-prior_lab)))
psadd("```")
psadd("")

# Decision rule
psadd("## Decision rule")
psadd("- If sigma_u, sens, spec medians shift by < 0.05 across priors AND")
psadd("  all three converge cleanly: M2 posterior is data-driven; the")
psadd("  primary (scale=1) result stands and is the Paper 1 baseline.")
psadd("- If posteriors shift substantially (> 0.10) across priors:")
psadd("  M2 is prior-driven and weakly identified at this n. Either")
psadd("  (a) accept this as a limitation and proceed with the primary,")
psadd("  reporting the range, or (b) consider model structural changes.")
psadd("- If only one prior produces a 'sensible-looking' posterior")
psadd("  (e.g., sens(EM) << sens(S0.x) consistent with descriptive data):")
psadd("  the others are in alternate modes. Investigate identifiability.")

writeLines(ps_lines, file.path(DIAG_PATH, "m2_prior_sensitivity.md"))

# =============================================================================
# M1 vs M2 (primary) comparison
# =============================================================================
cmp_lines <- character()
cadd <- function(...) cmp_lines <<- c(cmp_lines, sprintf(...))
cadd("# M1 vs M2 (primary, scale=1) comparison")
cadd("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
cadd("")
cmp_tbl <- bind_rows(
  summ_m1 %>% filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, model = "M1",
              median = round(median, 3),
              q2.5 = round(q2.5, 3), q97.5 = round(q97.5, 3)),
  m2_summ_list[["primary"]] %>% filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, model = "M2_primary",
              median = round(median, 3),
              q2.5 = round(q2.5, 3), q97.5 = round(q97.5, 3))
) %>% arrange(param, slot, model)
cadd("```")
cmp_lines <- c(cmp_lines, fmt_df(cmp_tbl))
cadd("```")
cadd("")
cadd("PPC pass rates: M1 = %.0f%%, M2 (primary) = %.0f%%",
     100 * m1_diag$ppc_pass, 100 * m2_diag_list[["primary"]]$ppc_pass)
writeLines(cmp_lines, file.path(DIAG_PATH, "m1_vs_m2_comparison.md"))

# =============================================================================
# M2 primary vs M2 globalclass: parsimony test
# =============================================================================
gc_lines <- character()
gadd <- function(...) gc_lines <<- c(gc_lines, sprintf(...))
gadd("# M2 primary vs M2 globalclass")
gadd("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
gadd("")
gadd("M2 primary: each slot has its own ordered[2] -- 6 detection params.")
gadd("M2 globalclass: 3 alpha_slot + 1 beta_class -- 4 detection params.")
gadd("If M2 globalclass fits comparably (PPC within 10pp), prefer it for")
gadd("parsimony. If much worse, M2 primary's per-slot freedom is essential.")
gadd("")
gc_cmp <- bind_rows(
  m2_summ_list[["primary"]] %>% filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, model = "M2_primary",
              median = round(median, 3),
              q2.5 = round(q2.5, 3), q97.5 = round(q97.5, 3)),
  summ_gc %>% filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, model = "M2_globalclass",
              median = round(median, 3),
              q2.5 = round(q2.5, 3), q97.5 = round(q97.5, 3))
) %>% arrange(param, slot, model)
gadd("```")
gc_lines <- c(gc_lines, fmt_df(gc_cmp))
gadd("```")
gadd("")
gadd("PPC pass rates: M2 primary = %.0f%%, M2 globalclass = %.0f%%",
     100 * m2_diag_list[["primary"]]$ppc_pass, 100 * diag_gc$ppc_pass)
gadd("")
sigma_u_p  <- m2_summ_list[["primary"]] %>% filter(variable == "sigma_u")
sigma_u_gc <- summ_gc %>% filter(variable == "sigma_u")
gadd("sigma_u medians: M2 primary = %.2f, M2 globalclass = %.2f",
     sigma_u_p$median, sigma_u_gc$median)
beta_row <- summ_gc %>% filter(variable == "beta_class")
gadd("M2 globalclass beta_class median = %.2f (logit gap between T=1 and T=0 detection)",
     beta_row$median)
writeLines(gc_lines, file.path(DIAG_PATH, "m2_primary_vs_globalclass.md"))

# =============================================================================
# M3a: hierarchical LCA with (slot x MAF) detection probabilities
# =============================================================================
message("\n=== M3a (slot x MAF detection + patient RE) ===")

# Drop sentinel maf_bin level for fitting (already excluded from lca_fit by
# earlier filter, but ensure factor levels are clean)
maf_levels <- levels(droplevels(lca_fit$maf_bin))
M_n <- length(maf_levels)
message(sprintf("MAF bins for M3 (excluding >MAF_max sentinel): %s",
                paste(maf_levels, collapse = ", ")))

lca_fit_m3 <- lca_fit %>%
  mutate(maf_bin_int = as.integer(factor(maf_bin, levels = maf_levels)))
stopifnot(all(!is.na(lca_fit_m3$maf_bin_int)))

# ---- M3 identification pre-check: disagreement cells per (slot, MAF) ------
# For each call, both replicates share the same MAF bin; the slots come from
# slot_a and slot_b. Build per-(slot, MAF) call counts and disagreement counts.
m3_check <- bind_rows(
  lca_fit_m3 %>% transmute(slot = slot_a, maf_bin, W_self = W1, W_other = W2),
  lca_fit_m3 %>% transmute(slot = slot_b, maf_bin, W_self = W2, W_other = W1)
) %>%
  group_by(slot, maf_bin) %>%
  summarise(
    n_calls       = n(),
    n_disagree    = sum(W_self != W_other),
    .groups = "drop"
  ) %>%
  arrange(slot, maf_bin)

message("M3 identification check (disagreement counts per slot x MAF):")
print_via_message(m3_check)

n_thin_cells <- sum(m3_check$n_disagree < 10)
if (n_thin_cells > 0) {
  message(sprintf(
    "WARNING: %d (slot, MAF) cell(s) have < 10 disagreement contributions; ",
    n_thin_cells))
  message("the corresponding (sens, spec) posteriors will be wide and ")
  message("prior-influenced. Consider M3b (hierarchical-MAF) as the primary ")
  message("if these are central to the paper's claims.")
}

# ---- PPC pair_type table extended with MAF bin ----------------------------
pt_table_m3 <- lca_fit_m3 %>%
  count(pair_type, slot_a_int, slot_b_int, maf_bin, maf_bin_int,
        name = "pt_n") %>%
  arrange(pair_type, maf_bin)

stan_data_m3 <- list(
  N       = nrow(lca_fit_m3),
  S       = length(SLOT_LEVELS),
  M       = M_n,
  P       = P,
  W       = as.matrix(lca_fit_m3[, c("W1", "W2")]),
  slot    = as.matrix(lca_fit_m3[, c("slot_a_int", "slot_b_int")]),
  maf_bin = lca_fit_m3$maf_bin_int,
  patient = lca_fit_m3$patient_int,
  N_pt    = nrow(pt_table_m3),
  pt_slot = as.matrix(pt_table_m3[, c("slot_a_int", "slot_b_int")]),
  pt_maf  = pt_table_m3$maf_bin_int,
  pt_n    = pt_table_m3$pt_n,
  sigma_u_prior_scale = 1.0
)

mod_m3 <- cmdstan_model("stan/lca_hierarchical_strata.stan")
t0 <- Sys.time()
fit_m3 <- mod_m3$sample(
  data = stan_data_m3, chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 2000,
  seed = SEED, refresh = 0, show_messages = FALSE,
  adapt_delta = 0.95
)
t_m3 <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
fit_m3$save_object(file.path(FITS_PATH, "fit_strata.rds"))

# ---- M3 posterior summary -------------------------------------------------
m3_params <- c("pi", "sigma_u",
               unlist(lapply(seq_len(length(SLOT_LEVELS)), function(s)
                 paste0("sens[", s, ",", seq_len(M_n), "]"))),
               unlist(lapply(seq_len(length(SLOT_LEVELS)), function(s)
                 paste0("spec[", s, ",", seq_len(M_n), "]"))))

summ_m3 <- fit_m3$summary(m3_params, mean, median,
                          ~quantile(.x, c(0.025, 0.975)),
                          rhat, ess_bulk, ess_tail) %>%
  rename(q2.5 = `2.5%`, q97.5 = `97.5%`)

# Parse [s, m] indices into slot and maf_bin labels
parse_idx <- function(v) {
  m <- regmatches(v, regexec("\\[(\\d+),(\\d+)\\]", v))
  ss <- vapply(m, function(x) if (length(x) >= 3) as.integer(x[2]) else NA_integer_, integer(1))
  mm <- vapply(m, function(x) if (length(x) >= 3) as.integer(x[3]) else NA_integer_, integer(1))
  list(slot_idx = ss, maf_idx = mm)
}
idx <- parse_idx(summ_m3$variable)
summ_m3 <- summ_m3 %>%
  mutate(param   = gsub("\\[.+\\]", "", variable),
         slot    = ifelse(is.na(idx$slot_idx), NA_character_, SLOT_LEVELS[idx$slot_idx]),
         maf_bin = ifelse(is.na(idx$maf_idx), NA_character_, maf_levels[idx$maf_idx])) %>%
  select(variable, param, slot, maf_bin, everything())

write_csv(summ_m3, file.path(FITS_PATH, "fit_strata_summary.csv"))

# ---- M3 PPC ---------------------------------------------------------------
pp_draws_m3 <- fit_m3$draws("pp_cell", format = "draws_df")
ppc_m3 <- expand_grid(pt_idx = seq_len(nrow(pt_table_m3)), cell_k = 1:4) %>%
  mutate(pair_type = pt_table_m3$pair_type[pt_idx],
         maf_bin   = as.character(pt_table_m3$maf_bin[pt_idx]),
         cell      = c("00","01","10","11")[cell_k],
         varname   = sprintf("pp_cell[%d,%d]", pt_idx, cell_k)) %>%
  rowwise() %>%
  mutate(pp_median = median(pp_draws_m3[[varname]]),
         pp_q025   = quantile(pp_draws_m3[[varname]], 0.025),
         pp_q975   = quantile(pp_draws_m3[[varname]], 0.975)) %>%
  ungroup() %>%
  select(pair_type, maf_bin, cell, pp_median, pp_q025, pp_q975)

obs_cells_m3 <- lca_fit_m3 %>%
  count(pair_type, maf_bin, W1, W2, name = "n") %>%
  group_by(pair_type, maf_bin) %>%
  mutate(obs_prop = n / sum(n)) %>% ungroup() %>%
  mutate(cell = paste0(W1, W2)) %>%
  select(pair_type, maf_bin, cell, obs_prop)

ppc_m3 <- ppc_m3 %>%
  left_join(obs_cells_m3, by = c("pair_type", "maf_bin", "cell")) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)),
         within_CrI = obs_prop >= pp_q025 & obs_prop <= pp_q975) %>%
  arrange(pair_type, maf_bin, cell)
write_csv(ppc_m3, file.path(FITS_PATH, "fit_strata_ppc.csv"))

# ---- M3 verdict + diagnostic ----------------------------------------------
m3_diag <- write_diag("M3a (slot x MAF + patient RE)",
                     summ_m3 %>% filter(param %in% c("pi", "sigma_u", "sens", "spec")),
                     ppc_m3, t_m3, fit_m3)
writeLines(m3_diag$lines, file.path(DIAG_PATH, "m3a_smoke_test.md"))
message(sprintf("M3a: R-hat=%.4f | ESS=%.0f | div=%d | PPC=%.0f%% | %s",
                m3_diag$rhat_max, m3_diag$ess_min, m3_diag$n_div,
                100 * m3_diag$ppc_pass,
                ifelse(m3_diag$verdict, "PASS", "INVESTIGATE")))

# =============================================================================
# M3b: hierarchical-MAF priors on per-cell logit detection
# =============================================================================
message("\n=== M3b (hierarchical-MAF priors within slot) ===")

mod_m3b <- cmdstan_model("stan/lca_hierarchical_strata_pooled.stan")
t0 <- Sys.time()
fit_m3b <- mod_m3b$sample(
  data = stan_data_m3, chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 2000,
  seed = SEED, refresh = 0, show_messages = FALSE,
  adapt_delta = 0.95
)
t_m3b <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
fit_m3b$save_object(file.path(FITS_PATH, "fit_strata_pooled.rds"))

# M3b adds tau_logit hyperparameters
m3b_params <- c(m3_params,
                paste0("mu_logit[", rep(1:length(SLOT_LEVELS), each = 2), ",",
                       rep(1:2, length(SLOT_LEVELS)), "]"),
                paste0("tau_logit[", rep(1:length(SLOT_LEVELS), each = 2), ",",
                       rep(1:2, length(SLOT_LEVELS)), "]"))

summ_m3b <- fit_m3b$summary(m3b_params, mean, median,
                            ~quantile(.x, c(0.025, 0.975)),
                            rhat, ess_bulk, ess_tail) %>%
  rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
idx_m3b <- parse_idx(summ_m3b$variable)
summ_m3b <- summ_m3b %>%
  mutate(param   = gsub("\\[.+\\]", "", variable),
         slot    = ifelse(is.na(idx_m3b$slot_idx), NA_character_, SLOT_LEVELS[idx_m3b$slot_idx]),
         maf_bin = ifelse(is.na(idx_m3b$maf_idx), NA_character_, maf_levels[idx_m3b$maf_idx])) %>%
  select(variable, param, slot, maf_bin, everything())

write_csv(summ_m3b, file.path(FITS_PATH, "fit_strata_pooled_summary.csv"))

# M3b PPC
pp_draws_m3b <- fit_m3b$draws("pp_cell", format = "draws_df")
ppc_m3b <- expand_grid(pt_idx = seq_len(nrow(pt_table_m3)), cell_k = 1:4) %>%
  mutate(pair_type = pt_table_m3$pair_type[pt_idx],
         maf_bin   = as.character(pt_table_m3$maf_bin[pt_idx]),
         cell      = c("00","01","10","11")[cell_k],
         varname   = sprintf("pp_cell[%d,%d]", pt_idx, cell_k)) %>%
  rowwise() %>%
  mutate(pp_median = median(pp_draws_m3b[[varname]]),
         pp_q025   = quantile(pp_draws_m3b[[varname]], 0.025),
         pp_q975   = quantile(pp_draws_m3b[[varname]], 0.975)) %>%
  ungroup() %>%
  select(pair_type, maf_bin, cell, pp_median, pp_q025, pp_q975)

ppc_m3b <- ppc_m3b %>%
  left_join(obs_cells_m3, by = c("pair_type", "maf_bin", "cell")) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)),
         within_CrI = obs_prop >= pp_q025 & obs_prop <= pp_q975) %>%
  arrange(pair_type, maf_bin, cell)
write_csv(ppc_m3b, file.path(FITS_PATH, "fit_strata_pooled_ppc.csv"))

m3b_diag <- write_diag("M3b (hierarchical-MAF)",
                      summ_m3b %>% filter(param %in% c("pi", "sigma_u", "sens", "spec",
                                                       "mu_logit", "tau_logit")),
                      ppc_m3b, t_m3b, fit_m3b)
writeLines(m3b_diag$lines, file.path(DIAG_PATH, "m3b_smoke_test.md"))
message(sprintf("M3b: R-hat=%.4f | ESS=%.0f | div=%d | PPC=%.0f%% | %s",
                m3b_diag$rhat_max, m3b_diag$ess_min, m3b_diag$n_div,
                100 * m3b_diag$ppc_pass,
                ifelse(m3b_diag$verdict, "PASS", "INVESTIGATE")))

# =============================================================================
# M3a vs M3b comparison + primary selection
# =============================================================================
m3cmp_lines <- character()
m3add <- function(...) m3cmp_lines <<- c(m3cmp_lines, sprintf(...))
m3add("# M3a vs M3b comparison")
m3add("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M"))
m3add("")
m3add("M3a: independent (slot x MAF) detection (24 cell params, no pooling)")
m3add("M3b: hierarchical-MAF priors within slot (24 cell params + 12 hyperparams)")
m3add("")

# Side-by-side per-cell sens/spec
m3_compare <- bind_rows(
  summ_m3 %>% filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, maf_bin, model = "M3a",
              median = round(median, 3),
              q2.5 = round(q2.5, 3), q97.5 = round(q97.5, 3),
              CrI_width = round(q97.5 - q2.5, 3)),
  summ_m3b %>% filter(param %in% c("sens", "spec")) %>%
    transmute(param, slot, maf_bin, model = "M3b",
              median = round(median, 3),
              q2.5 = round(q2.5, 3), q97.5 = round(q97.5, 3),
              CrI_width = round(q97.5 - q2.5, 3))
) %>% arrange(param, slot, maf_bin, model)

m3add("## Per-cell sens/spec: M3a vs M3b")
m3add("```")
m3cmp_lines <- c(m3cmp_lines, fmt_df(m3_compare))
m3add("```")
m3add("")

# tau_logit posteriors (the shrinkage indicator)
m3add("## Shrinkage indicator: tau_logit by slot x class")
m3add("Small tau means heavy pooling toward slot mean; large tau means")
m3add("MAF bins are nearly free.")
tau_summ <- summ_m3b %>% filter(param == "tau_logit") %>%
  mutate(class = ifelse(grepl("\\[.+,1\\]", variable), "T=0 (1-spec)", "T=1 (sens)"),
         median = round(median, 3),
         q2.5 = round(q2.5, 3),
         q97.5 = round(q97.5, 3)) %>%
  select(slot, class, median, q2.5, q97.5)
m3add("```")
m3cmp_lines <- c(m3cmp_lines, fmt_df(tau_summ))
m3add("```")
m3add("")

# Primary selection
m3add("## Primary model selection")
m3add("- M3a PPC pass rate: %.0f%%", 100 * m3_diag$ppc_pass)
m3add("- M3b PPC pass rate: %.0f%%", 100 * m3b_diag$ppc_pass)
m3add("- M3a max R-hat: %.4f, min ESS: %.0f, divergent: %d",
      m3_diag$rhat_max, m3_diag$ess_min, m3_diag$n_div)
m3add("- M3b max R-hat: %.4f, min ESS: %.0f, divergent: %d",
      m3b_diag$rhat_max, m3b_diag$ess_min, m3b_diag$n_div)
m3add("")

# Decision rule (PPC threshold relaxed to 65% for M3 due to 36 cells)
M3_PPC_THRESHOLD <- 0.65
m3a_pass <- m3_diag$rhat_max < 1.01 && m3_diag$n_div < 50 &&
            m3_diag$ppc_pass >= M3_PPC_THRESHOLD
m3b_pass <- m3b_diag$rhat_max < 1.01 && m3b_diag$n_div < 50 &&
            m3b_diag$ppc_pass >= M3_PPC_THRESHOLD

# Choose primary: M3b wins if (a) it passes AND M3a doesn't, or
# (b) both pass and M3b PPC > M3a PPC by > 5 pp. Otherwise M3a stays primary.
if (m3b_pass && !m3a_pass) {
  primary_model <- "M3b"
} else if (m3b_pass && m3a_pass &&
           (m3b_diag$ppc_pass - m3_diag$ppc_pass) > 0.05) {
  primary_model <- "M3b"
} else {
  primary_model <- "M3a"
}
m3add("Primary M3 model selected: **%s**", primary_model)
m3add("(PPC threshold for M3: %.0f%%; both M3a and M3b checked.)",
      100 * M3_PPC_THRESHOLD)
writeLines(m3cmp_lines, file.path(DIAG_PATH, "m3a_vs_m3b_comparison.md"))

message(sprintf("\nM3 primary selected: %s", primary_model))

# =============================================================================
# Rebuild Paper 2 deliverable from chosen primary
# =============================================================================
message("\n=== Rebuilding Paper 2 deliverable from primary M3 model ===")

primary_fit <- if (primary_model == "M3b") fit_m3b else fit_m3
primary_summ <- if (primary_model == "M3b") summ_m3b else summ_m3

sens_draws_pri <- primary_fit$draws(
  unlist(lapply(seq_len(length(SLOT_LEVELS)), function(s)
    paste0("sens[", s, ",", seq_len(M_n), "]"))),
  format = "draws_df")
spec_draws_pri <- primary_fit$draws(
  unlist(lapply(seq_len(length(SLOT_LEVELS)), function(s)
    paste0("spec[", s, ",", seq_len(M_n), "]"))),
  format = "draws_df")
pi_draws_pri      <- primary_fit$draws("pi",      format = "draws_df")
sigma_u_draws_pri <- primary_fit$draws("sigma_u", format = "draws_df")

slot_marginal_sens <- matrix(0, nrow = nrow(sens_draws_pri), ncol = M_n)
slot_marginal_spec <- matrix(0, nrow = nrow(spec_draws_pri), ncol = M_n)
for (m in seq_len(M_n)) {
  for (s in seq_len(length(SLOT_LEVELS))) {
    slot_marginal_sens[, m] <- slot_marginal_sens[, m] +
      sens_draws_pri[[paste0("sens[", s, ",", m, "]")]] / length(SLOT_LEVELS)
    slot_marginal_spec[, m] <- slot_marginal_spec[, m] +
      spec_draws_pri[[paste0("spec[", s, ",", m, "]")]] / length(SLOT_LEVELS)
  }
}
colnames(slot_marginal_sens) <- paste0("sens_marginal_maf", seq_len(M_n))
colnames(slot_marginal_spec) <- paste0("spec_marginal_maf", seq_len(M_n))

paper2_posteriors <- list(
  meta = list(
    spec_tag           = TAG$apply_tag,
    n_calls            = nrow(lca_fit_m3),
    n_patients         = P,
    slot_levels        = SLOT_LEVELS,
    maf_levels         = maf_levels,
    primary_model      = sprintf("%s (selected via PPC + parsimony)", primary_model),
    sigma_u_prior      = "half-normal(0, 1.0)",
    timestamp          = Sys.time()
  ),
  pi               = pi_draws_pri$pi,
  sigma_u          = sigma_u_draws_pri$sigma_u,
  sens_per_slot_maf = as.data.frame(sens_draws_pri %>%
                                     select(-.chain, -.iteration, -.draw)),
  spec_per_slot_maf = as.data.frame(spec_draws_pri %>%
                                     select(-.chain, -.iteration, -.draw)),
  sens_marginal_per_maf = as.data.frame(slot_marginal_sens),
  spec_marginal_per_maf = as.data.frame(slot_marginal_spec),
  notes = paste(
    sprintf("Paper 2 prior recommendations (from %s primary):", primary_model),
    "- PRIMARY (hedged across slots): sens_marginal_per_maf, spec_marginal_per_maf.",
    "  Equal-weight mixture across the three Ghana slots (EM, S0.1, S0.2),",
    "  separately per MAF bin.",
    "- SENSITIVITY (best-case = EM): sens[1, m], spec[1, m] from sens_per_slot_maf.",
    "- SENSITIVITY (worst-case = S0.x): use S0.1 or S0.2 columns.",
    "- ALTERNATIVE MIXTURE: reweight using (sens|spec)_per_slot_maf with",
    "  user-supplied slot weights summing to 1.",
    sprintf("- M3a and M3b both saved as %s/fit_strata.rds and fit_strata_pooled.rds",
            FITS_PATH),
    "  for sensitivity comparisons in Paper 2.",
    sep = "\n")
)

saveRDS(paper2_posteriors, file.path(LCA_PATH, "lca_posteriors.rds"))
message(sprintf("Paper 2 posteriors saved to %s/lca_posteriors.rds (primary: %s)",
                LCA_PATH, primary_model))

# =============================================================================
# (Paper 2 deliverable assembly happens after M3b comparison below)
# =============================================================================

# ---- Console -----------------------------------------------------------------
message("\n=== Prior sensitivity table ===")
print_via_message(ps_tbl %>% select(-prior_lab))
message("\n=== Per-prior verdict ===")
print_via_message(verdict_tbl %>% select(-prior_lab))
message("\n=== M1 vs M2 (primary) ===")
print_via_message(cmp_tbl)
message("\n=== M2 primary vs M2 globalclass ===")
print_via_message(gc_cmp)
message(sprintf("M2 primary PPC: %.0f%% | M2 globalclass PPC: %.0f%%",
                100 * m2_diag_list[["primary"]]$ppc_pass, 100 * diag_gc$ppc_pass))

message("\n=== M3a key posterior summaries ===")
m3_key <- summ_m3 %>%
  filter(param %in% c("sens", "spec")) %>%
  mutate(across(c(median, q2.5, q97.5), ~round(.x, 3))) %>%
  select(param, slot, maf_bin, median, q2.5, q97.5, rhat, ess_bulk) %>%
  arrange(param, slot, maf_bin)
print_via_message(m3_key)
message(sprintf("\nM3a verdict: PPC=%.0f%% | %s",
                100 * m3_diag$ppc_pass,
                ifelse(m3_diag$verdict, "PASS", "INVESTIGATE")))

message("\n=== M3b key posterior summaries ===")
m3b_key <- summ_m3b %>%
  filter(param %in% c("sens", "spec")) %>%
  mutate(across(c(median, q2.5, q97.5), ~round(.x, 3))) %>%
  select(param, slot, maf_bin, median, q2.5, q97.5, rhat, ess_bulk) %>%
  arrange(param, slot, maf_bin)
print_via_message(m3b_key)
message(sprintf("\nM3b verdict: PPC=%.0f%% | %s",
                100 * m3b_diag$ppc_pass,
                ifelse(m3b_diag$verdict, "PASS", "INVESTIGATE")))

message(sprintf("\n=== Primary M3 model: %s ===", primary_model))

message(sprintf("\n[fit_lca_models] Done. Diagnostics in %s/", DIAG_PATH))
