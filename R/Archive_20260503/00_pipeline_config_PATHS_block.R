# ---- Output paths (Paper 1 active set) -------------------------------------
# All per-spec outputs (data_derived AND outputs/...) carry the active spec's
# tag. Only paths used by the Paper 1 pipeline are auto-created. Legacy paths
# from earlier paper drafts are commented at the bottom; uncomment if you need
# to re-run those scripts.
PATHS <- list(
  data_raw      = "data_raw",
  meta          = "data_derived/00_metadata",
  variants      = "data_derived/00_variants",
  calibration   = file.path("data_derived",
                            paste0("01_calibration_",  TAG$cal_tag)),
  thresholded   = file.path("data_derived",
                            paste0("03_thresholded_", TAG$apply_tag)),
  lca           = file.path("data_derived",
                            paste0("05_lca_",          TAG$apply_tag)),
  lca_fits      = file.path("data_derived",
                            paste0("05_lca_",          TAG$apply_tag), "fits"),
  lca_diag      = file.path("data_derived",
                            paste0("05_lca_",          TAG$apply_tag),
                            "diagnostics"),
  tables        = file.path("outputs", "tables",       TAG$apply_tag),
  figures       = file.path("outputs", "figures",      TAG$apply_tag),
  supplemental  = file.path("outputs", "supplemental", TAG$apply_tag)
  # ---- Legacy (uncomment to re-enable) ----
  # external         = file.path("data_derived",
  #                              paste0("04_external_",    TAG$apply_tag)),
  # diagnostics_inv  = file.path("data_derived",
  #                              paste0("06_diagnostics_", TAG$apply_tag)),
  # contrast_tables  = file.path("outputs", "contrast", "tables"),
  # contrast_figures = file.path("outputs", "contrast", "figures")
)

# Create output directories ---------------------------------------------------
for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)
