# Paper 1 — Replicate-Anchored iSNV Detection Calibration in *Mycobacterium tuberculosis*

This repository accompanies the manuscript *"Replicate-anchored calibration of within-host single nucleotide variant detection in Mycobacterium tuberculosis whole genome sequencing"* (Séraphin et al., submitted to *Microbial Genomics* as a Short Communication, May 2026).

The work develops and validates a multi-criteria, lexicographically ranked, bootstrap-stabilized calibration framework for per-specimen intra-host single-nucleotide variant (iSNV) detection in *M. tuberculosis* whole-genome sequencing. Within-patient replicate sputum pairs from a pre-treatment TB cohort in Accra, Ghana, are scored across the joint (DP, AD, MAF) grid and ranked by reproducibility. Resulting thresholds are reported as a three-tier sensitivity ladder and applied to the full 97-patient cohort.

## Submission snapshot

The repository state at manuscript submission is tagged **`submission-microbial-genomics`** (commit d75fb5077fdd3c7c4ba55ef38134ebfabe008a66, dated 20260509). All results, figures, and tables in the submitted manuscript can be reproduced from the code at this snapshot.

GitHub release: https://github.com/MapouSeraphin/paper1_isnv_calibration/releases/tag/submission-microbial-genomics 

Zenodo archive (post-acceptance): `[DOI]`

## What was done

We defined the iSNV detection rule by three thresholds on read depth (DP), alternate-allele support (AD₁), and minor allele frequency (MAF), with MAF_max fixed at 0.50. Calibration used **169 within-visit replicate pairs collected at the pre-treatment timepoint (M0) from 67 patients**. A 612-cell grid (DP ∈ {40…200}, AD₁ ∈ {3…8}, MAF_min ∈ {0.02, 0.03, 0.05, 0.10, 0.15, 0.20}) was scored on six pair-level concordance metrics and ranked under a lexicographic selector minimizing the proportion of pairs with both specimens silent and maximizing reproducibility. Selection stability was quantified by **B = 1,000 nonparametric pair-level bootstrap resamples**, defining a Looser/Primary/Tighter sensitivity ladder.

The calibrated rule was applied to **282 cultured sputum specimens from 97 pre-treatment TB patients** (Korle-Bu Teaching Hospital, Greater Accra Region, July 2022 – June 2023). Construct validity was assessed by regressing pair-level concordance outcomes on minimum pair coverage with patient-clustered standard errors. Patient-level detection prevalence is reported under the calibrated rule, both overall and stratified by HIV status, sex, age, hemoptysis, and *M. tuberculosis* lineage.

## Headline results

- **Primary calibration cell:** DP ≥ 60×, AD₁ ≥ 3, MAF ∈ [0.02, 0.50].
- **Bracketing tiers:** Looser (DP ≥ 40×, AD₁ ≥ 3, MAF ∈ [0.02, 0.50]); Tighter (DP ≥ 100×, AD₁ ≥ 6, MAF ∈ [0.02, 0.50]).
- **Bootstrap stability:** Stage 1 primary cell coincided with the bootstrap modal cell in 45.7% of 703 successful replications; MAF = 0.02 selected in 100%, AD₁ = 3 in 90.2%.
- **Per-patient detectable iSNV prevalence at the Primary tier:** 16.5% (16/97, Wilson 95% CI 10.4%–25.1%).
- **Construct validity:** flat both-zero slope on minimum pair coverage (slope 0.087 per +10×, p = 0.15) and positive Jaccard slopes on coverage at all three tiers (Looser 0.030, p = 0.04; Primary 0.029, p = 0.05; Tighter 0.034, p = 0.02; n = 6 evaluable pairs, reported as exploratory).

## Reproducibility

- **R version:** 4.5.1
- **Random seed:** `20260428` (used for the calibration bootstrap and any other stochastic step).

To reproduce: clone the repository, install the listed package versions and run `run_pipeline.R` from the repository root. The pipeline regenerates derived data tables, main figures, and supplementary materials in the output directory.

## Data availability

- **Raw sequencing reads** (282 paired-end Illumina WGS samples) deposited at NCBI SRA under BioProject `[PRJNA######]`. Per-sample accessions and metadata in `Supplementary_Table_S1.csv` (and listed in the manuscript supplement).
- **Reference genome:** *M. tuberculosis* H37Rv, GenBank accession NC_000962.3.
- **Processed data tables** (calibration grid output, full variant call tables, derived analytic frames) included in this repository under `data/processed/`.

## Repository structure

```
run_paper1.R                         Master orchestrator (runs all 3 specs)
R/                                   Pipeline source
  ├── 00_pipeline_config.R           Canonical flags, tag construction, paths
  ├── 00_prep_metadata.R             Cohort cleaning, metadata harmonization
  ├── 01_build_cal_pairs.R           M0 within-visit replicate-pair construction
  ├── 02_calibrate_lexicographic.R   Stage 1: lex grid search (canonical)
  ├── 02b_calibration_bootstrap.R    Stage 2: B=1000 selection-stability bootstrap
  ├── 03_apply_thresholds.R          Apply calibrated rule to full cohort
  ├── application_per_patient_prevalence.R   Per-patient detection prevalence
  ├── application_persistence_analysis.R     Longitudinal persistence (M0→M1/M2)
  ├── make_paper1_tables_1to3.R      Tables 1, 2, 3
  └── 09_make_calibration_selection_figure_table.R   Figure 1, Table 4

manuscript/
  ├── calibration_validation_manuscript_20260506.docx
  └── supplement/

data_derived/                        Per-spec derived data (auto-generated)
  ├── 00_metadata/                   Spec-independent
  ├── 00_variants/                   Spec-independent (canonical universe)
  ├── 00_variants_ppe_excluded/      Spec-independent (Sens A universe)
  ├── 01_calibration_<cal_tag>/      Spec-tagged: thresholds, grid, bootstrap/
  └── 03_thresholded_<apply_tag>/    Spec-tagged: gh_variants_thr_*.rds

outputs/                              Final tables, figures, logs
  ├── tables/<lca_tag>/
  ├── figures/<lca_tag>/
  ├── supplemental/<lca_tag>/
  └── orchestrator_logs/

README.md                             This file
LICENSE                               MIT
```

*(Update file listing if directory structure differs at commit time.)*

## Citation

If you use this work, please cite:

> Séraphin MN, Afriyie-Mensah JS, Asare-Baah M, Chariker J, Domotey C, Kwarteng E, Zoungrana M, Mireku Appah S, Ganu H, Amo Omari M. Replicate-anchored calibration of within-host single nucleotide variant detection in Mycobacterium tuberculosis whole genome sequencing. *Microbial Genomics* (submitted, 2026). DOI: `[DOI when assigned]`

This calibration is also applied in the companion methods paper:

> Séraphin MN, et al. (in preparation). Operationalizing covariate-differential outcome misclassification for depth-dependent within-host pathogen sequencing outcomes: HIV and *M. tuberculosis* iSNV detection across cohorts. *American Journal of Epidemiology*. Repository: `github.com/MapouSeraphin/paper2-isnv-transport`.

## Funding

National Institutes of Health (K01AI153544); Gatorade Trust through the University of Florida Department of Medicine. The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

## Ethics

Approved by the Institutional Review Board of the University of Florida (IRB202003042) and the Institutional Review Board of Korle-Bu Teaching Hospital, Accra, Ghana (IRB/000135/2020). All participants provided written informed consent prior to enrollment.

## Contact

Corresponding author: **Marie Nancy Séraphin**, Department of Epidemiology, College of Public Health and Health Professions, College of Medicine, University of Florida, Gainesville, FL, USA.
Email:nseraphin@ufl.edu

## License

MIT, CC-BY-4.0 — see `LICENSE` file. Recommended: MIT or BSD-3-Clause for code; CC-BY-4.0 for the manuscript draft and figures.

---

## Placeholder summary (delete this section before publishing)

Fill in before tagging the submission snapshot:

- `[TAG_NAME]` — submission tag (e.g., `paper1-submission-2026-05-22` or `submission-microbial-genomics`).
- `[FULL_40_CHAR_SHA]` — commit hash, fill in after first commit.
- `[YYYY-MM-DD]` — commit date.
- `[https://...releases/tag/TAG_NAME]` — release URL.
- `[X.Y.Z]`, `[renv lockfile location]`, `[OS_VERSION]` — reproducibility metadata.
- `[PRJNA######]` — NCBI SRA BioProject accession.
- `[DOI]` (Zenodo) — once Zenodo snapshot created.
- `[DOI]` (manuscript) — once accepted/published.
- `[EMAIL]` — corresponding author email.
- `[LICENSE_NAME]` — choose license, add `LICENSE` file.
