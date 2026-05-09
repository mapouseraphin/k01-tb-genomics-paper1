# =============================================================================
# 02c_depth_concordance_diagnostic.R
#
# Construct-validity diagnostic for the calibrated rule. At each tier
# (Looser, Primary, Tighter), regress two outcomes on pair minimum coverage:
#
#   (1) P(both_zero)  -- both specimens silent at this cell. Logistic.
#   (2) Jaccard       -- among evaluable pairs (nA > 0 AND nB > 0). Linear.
#
# Patient-clustered standard errors (sandwich::vcovCL) -- pairs nest in
# patients (169 pairs / 67 patients).
#
# Interpretation (per tier):
#   Measurement-like:   P(both_zero) DECREASES with depth (slope < 0)
#                       Jaccard      INCREASES with depth (slope > 0)
#   Stringency-cliff:   Both slopes ~ 0 (CIs cross zero)
#   Pathological:       P(both_zero) increases with depth, OR Jaccard falls
#
# Tighter tier may be NULL; script tolerates.
#
# Inputs:
#   PATHS$calibration/cal_pairs.rds
#   PATHS$calibration/thr_{looser,primary,tighter}.rds
#   PATHS$meta/gh_meta.rds
#
# Outputs:
#   PATHS$calibration/diagnostics/depth_concordance_table.csv
#   PATHS$calibration/diagnostics/depth_concordance_interpretation.csv
#   PATHS$calibration/diagnostics/depth_concordance_log.txt
#   PATHS$figures/figure_S_depth_concordance.pdf
# =============================================================================

source("R/00_pipeline_config.R")
source("R/isnv_helpers.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(tibble); library(purrr)
  library(ggplot2); library(sandwich); library(lmtest)
})

banner()
message("\n[02c_depth_concordance_diagnostic] depth-vs-concordance diagnostic")

DIAG_DIR <- file.path(PATHS$calibration, "diagnostics")
dir.create(DIAG_DIR,        recursive = TRUE, showWarnings = FALSE)
dir.create(PATHS$figures,   recursive = TRUE, showWarnings = FALSE)

LOG <- character()
log_msg <- function(...) {
  args <- list(...)
  s <- if (length(args) == 1) as.character(args[[1]]) else do.call(sprintf, args)
  message(s); LOG <<- c(LOG, s); invisible(s)
}

# ---- Load -------------------------------------------------------------------
cal_pairs <- readRDS(file.path(PATHS$calibration, "cal_pairs.rds"))
gh_meta   <- readRDS(file.path(PATHS$meta,        "gh_meta.rds"))

library(tibble)

coerce_thr <- function(x) {
  if (is.null(x)) return(NULL)
  # already a frame: take first row
  if (is.data.frame(x)) {
    if (nrow(x) < 1) return(NULL)
    return(x[1, , drop = FALSE])
  }
  # named list or vector: pluck the three fields
  nms <- names(x)
  need <- c("DP_min", "AD1_min", "MAF_min")
  if (is.null(nms) || !all(need %in% nms)) return(NULL)
  vals <- suppressWarnings(as.numeric(c(x[["DP_min"]],
                                        x[["AD1_min"]],
                                        x[["MAF_min"]])))
  if (any(is.na(vals))) return(NULL)
  tibble(DP_min = vals[1], AD1_min = vals[2], MAF_min = vals[3])
}

tiers <- list()
for (nm in c("looser", "primary", "tighter")) {
  fp <- file.path(PATHS$calibration, paste0("thr_", nm, ".rds"))
  if (!file.exists(fp)) next
  thr_raw <- tryCatch(readRDS(fp), error = function(e) NULL)
  thr <- coerce_thr(thr_raw)
  if (is.null(thr)) {
    log_msg("[load] %s: present but unusable; skipping.", fp)
    next
  }
  tiers[[nm]] <- thr
}
log_msg("Tiers loaded: %s", paste(names(tiers), collapse = ", "))
if (length(tiers) == 0) stop("No tier thresholds loadable from PATHS$calibration.")


# Resolve depth field defensively
depth_field <- intersect(c("depth_med", "coverage_median"), names(gh_meta))[1]
if (is.na(depth_field))
  stop("No depth field in gh_meta (expected depth_med or coverage_median).")
log_msg("Using depth field: %s", depth_field)

cov_lkp <- gh_meta %>%
  transmute(Sample, patientId,
            cov = suppressWarnings(as.numeric(.data[[depth_field]])))

# ---- Per-tier analysis -----------------------------------------------------
analyze_tier <- function(tier_name, thr) {
  cell <- cal_pairs %>%
    filter(DP_min  == thr$DP_min,
           AD1_min == thr$AD1_min,
           MAF_min == thr$MAF_min)
  
  if (nrow(cell) == 0) {
    log_msg("[%s] No pair rows for cell DP=%g AD1=%g MAF=%g; skipping.",
            tier_name, thr$DP_min, thr$AD1_min, thr$MAF_min)
    return(NULL)
  }
  
  d <- cell %>%
    left_join(cov_lkp %>% rename(sampleA = Sample, cov_a = cov, pid_a = patientId),
              by = "sampleA") %>%
    left_join(cov_lkp %>% rename(sampleB = Sample, cov_b = cov, pid_b = patientId),
              by = "sampleB") %>%
    mutate(
      pair_min_cov  = pmin(cov_a, cov_b),
      pair_mean_cov = (cov_a + cov_b) / 2,
      both_zero     = as.integer(nA == 0 & nB == 0),
      evaluable     = nA > 0 & nB > 0,
      patientId     = ifelse(pid_a == pid_b, pid_a, paste(pid_a, pid_b, sep = "|"))
    ) %>%
    filter(!is.na(pair_min_cov), !is.na(patientId))
  
  n_pairs <- nrow(d); n_pts <- n_distinct(d$patientId)
  log_msg("[%s] DP=%g AD1=%g MAF=%g: %d pairs / %d patients",
          tier_name, thr$DP_min, thr$AD1_min, thr$MAF_min, n_pairs, n_pts)
  
  # Slope expressed per +10x coverage so it is interpretable
  d$cov10 <- d$pair_min_cov / 10
  
  # (1) P(both_zero) ~ depth, logistic, patient-clustered SE
  fit_bz <- glm(both_zero ~ cov10, family = binomial(), data = d)
  vc_bz  <- sandwich::vcovCL(fit_bz, cluster = d$patientId, type = "HC0")
  ct_bz  <- lmtest::coeftest(fit_bz, vcov. = vc_bz)
  slope_bz <- unname(coef(fit_bz)[2])
  se_bz    <- unname(sqrt(diag(vc_bz))[2])
  p_bz     <- unname(ct_bz[2, 4])
  
  # (2) Jaccard ~ depth, linear, evaluable only, patient-clustered SE
  d_eval <- d %>% filter(evaluable, !is.na(jaccard))
  if (nrow(d_eval) >= 5 && n_distinct(d_eval$patientId) >= 2) {
    fit_jc <- lm(jaccard ~ cov10, data = d_eval)
    vc_jc  <- sandwich::vcovCL(fit_jc, cluster = d_eval$patientId, type = "HC0")
    ct_jc  <- lmtest::coeftest(fit_jc, vcov. = vc_jc)
    slope_jc <- unname(coef(fit_jc)[2])
    se_jc    <- unname(sqrt(diag(vc_jc))[2])
    p_jc     <- unname(ct_jc[2, 4])
  } else {
    slope_jc <- se_jc <- p_jc <- NA_real_
  }
  
  list(
    summary = tibble(
      tier            = tier_name,
      cell            = sprintf("DP=%g, AD1=%g, MAF=%g",
                                thr$DP_min, thr$AD1_min, thr$MAF_min),
      n_pairs         = n_pairs,
      n_patients      = n_pts,
      n_both_zero     = sum(d$both_zero == 1),
      n_evaluable     = sum(d$evaluable),
      slope_bz_per10x = slope_bz,
      se_bz_per10x    = se_bz,
      p_bz            = p_bz,
      slope_jc_per10x = slope_jc,
      se_jc_per10x    = se_jc,
      p_jc            = p_jc
    ),
    pair_data = d %>% mutate(tier = tier_name)
  )
}

results <- imap(tiers, ~ analyze_tier(.y, .x)) %>% compact()
if (length(results) == 0) stop("No tiers analyzable.")

# ---- Summary + plain-language interpretation -------------------------------
summary_tbl <- bind_rows(map(results, "summary"))
write_csv(summary_tbl, file.path(DIAG_DIR, "depth_concordance_table.csv"))
log_msg("\nSummary:")
log_msg(paste(capture.output(print(as.data.frame(summary_tbl))), collapse = "\n"))

interp <- function(slope, p, kind) {
  if (is.na(slope)) return("INSUFFICIENT DATA")
  sig <- !is.na(p) && p < 0.05
  if (kind == "bz") {
    if (!sig)         return("FLAT (cliff-like)")
    if (slope < 0)    return("OK (measurement-like): both-zero falls with depth")
    return("PATHOLOGICAL: both-zero rises with depth")
  } else {
    if (!sig)         return("FLAT (cliff-like or noisy)")
    if (slope > 0)    return("OK (measurement-like): Jaccard rises with depth")
    return("PATHOLOGICAL: Jaccard falls with depth")
  }
}
interp_tbl <- summary_tbl %>%
  rowwise() %>%
  mutate(
    interp_both_zero = interp(slope_bz_per10x, p_bz, "bz"),
    interp_jaccard   = interp(slope_jc_per10x, p_jc, "jc")
  ) %>% ungroup() %>%
  select(tier, cell, n_pairs, interp_both_zero, interp_jaccard)
write_csv(interp_tbl, file.path(DIAG_DIR, "depth_concordance_interpretation.csv"))
log_msg("\nInterpretation:")
log_msg(paste(capture.output(print(as.data.frame(interp_tbl))), collapse = "\n"))

# ---- Diagnostic figure -----------------------------------------------------
all_pairs <- bind_rows(map(results, "pair_data")) %>%
  mutate(tier = factor(tier, levels = c("looser", "primary", "tighter")))

pdf_path <- file.path(PATHS$figures, "figure_S_depth_concordance.pdf")
pdf(pdf_path, width = 10, height = 8)

p1 <- ggplot(all_pairs, aes(pair_min_cov, both_zero)) +
  geom_jitter(height = 0.04, width = 0, alpha = 0.4, size = 1) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"),
              se = TRUE, formula = y ~ x, color = "steelblue") +
  facet_wrap(~ tier, nrow = 1, labeller = label_both) +
  labs(x = "Pair min coverage (×)",
       y = "P(both specimens silent)",
       title = "(A) Both-zero probability vs depth",
       subtitle = "Measurement-like: slope < 0 | Cliff-like: slope ≈ 0") +
  theme_bw(base_size = 11)
print(p1)

p2 <- all_pairs %>%
  filter(evaluable, !is.na(jaccard)) %>%
  ggplot(aes(pair_min_cov, jaccard)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, color = "darkred") +
  facet_wrap(~ tier, nrow = 1, labeller = label_both) +
  labs(x = "Pair min coverage (×)",
       y = "Within-pair Jaccard (evaluable pairs)",
       title = "(B) Jaccard among evaluable pairs vs depth",
       subtitle = "Measurement-like: slope > 0 | Cliff-like: slope ≈ 0") +
  theme_bw(base_size = 11)
print(p2)

dev.off()
log_msg("\nFigure -> %s", pdf_path)

# ---- Log -------------------------------------------------------------------
writeLines(LOG, file.path(DIAG_DIR, "depth_concordance_log.txt"))

invisible(NULL)
