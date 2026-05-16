# =============================================================================
# Makefile — Paper 1 (Microbial Genomics) reproducibility entry point
#
# Convenience wrapper around run_paper1.R. The manuscript primary spec is
# hardcoded in run_paper1.R; this file provides named targets for stage
# subsets and archival tasks.
#
# Usage:
#   make help        # list targets
#   make all         # run full pipeline end-to-end
#   make verify      # confirm expected output directories exist
#   make manifest    # write SHA256 manifest of derived files
#   make clean       # remove all derived data and outputs
# =============================================================================

R   := Rscript
RUN := $(R) run_paper1.R

# Tag produced by the manuscript primary spec
TAG := snv_only_ppe_excluded_lex

.PHONY: help all prep calibrate apply figures application dry-run \
        clean verify manifest renv-snapshot

.DEFAULT_GOAL := help

help:
	@echo "Paper 1 (Microbial Genomics) reproducibility targets:"
	@echo ""
	@echo "Full reproduction:"
	@echo "  make all           Run entire pipeline (prep -> application)"
	@echo ""
	@echo "Partial reruns:"
	@echo "  make prep          Metadata preparation only"
	@echo "  make calibrate     Calibration only (Stage 1)"
	@echo "  make apply         Apply thresholds only"
	@echo "  make figures       Rebuild figures/tables only (Fig 2, Table 2)"
	@echo "  make application   Application analyses (Tables 4, 5, S4, S5, Fig 5)"
	@echo ""
	@echo "Verification & archival:"
	@echo "  make dry-run       Print execution plan, no side effects"
	@echo "  make verify        Check expected output directories exist"
	@echo "  make manifest      Write SHA256 manifest of derived files"
	@echo "  make renv-snapshot Snapshot R package versions to renv.lock"
	@echo ""
	@echo "  make clean         Remove data_derived/, outputs/, logs/"

all:
	$(RUN) --from prep --to application

prep:
	$(RUN) --from prep --to prep

calibrate:
	$(RUN) --from calibrate --to calibrate

apply:
	$(RUN) --from apply --to apply

figures:
	$(RUN) --from figures --to figures

application:
	$(RUN) --from application --to application

dry-run:
	$(RUN) --from prep --to application --dry-run

clean:
	@echo "Removing data_derived/, outputs/, logs/..."
	rm -rf data_derived/ outputs/ logs/
	@echo "Clean complete."

verify:
	@echo "=== $(TAG) ==="
	@if [ -d "data_derived/01_calibration_$(TAG)" ]; then \
	  n=$$(find "data_derived/01_calibration_$(TAG)" -type f | wc -l); \
	  echo "  calibration:  $$n files"; \
	else \
	  echo "  calibration:  MISSING"; \
	fi
	@if [ -d "data_derived/03_thresholded_$(TAG)" ]; then \
	  n=$$(find "data_derived/03_thresholded_$(TAG)" -type f | wc -l); \
	  echo "  thresholded:  $$n files"; \
	else \
	  echo "  thresholded:  MISSING"; \
	fi
	@if [ -d "outputs/tables/$(TAG)" ]; then \
	  n=$$(find "outputs/tables/$(TAG)" -type f | wc -l); \
	  echo "  tables:       $$n files"; \
	else \
	  echo "  tables:       MISSING"; \
	fi
	@if [ -d "outputs/figures/$(TAG)" ]; then \
	  n=$$(find "outputs/figures/$(TAG)" -type f | wc -l); \
	  echo "  figures:      $$n files"; \
	else \
	  echo "  figures:      MISSING"; \
	fi

manifest:
	@echo "Writing SHA256 manifest..."
	@find data_derived outputs -type f \
	  \( -name "*.csv" -o -name "*.tsv" -o -name "*.rds" \
	     -o -name "*.pdf" -o -name "*.png" -o -name "*.svg" \) \
	  2>/dev/null | sort | xargs sha256sum > paper1_manifest.sha256 || \
	  echo "(No derived files found.)"
	@echo "Manifest: paper1_manifest.sha256"
	@wc -l paper1_manifest.sha256 2>/dev/null || true

renv-snapshot:
	$(R) -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv"); renv::snapshot()'
