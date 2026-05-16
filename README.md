# Paper 1 — Replicate-Anchored iSNV Detection Calibration in *Mycobacterium tuberculosis*

This repository accompanies the manuscript in preparation *"Replicate-anchored calibration of within-host single nucleotide variant detection in Mycobacterium tuberculosis whole genome sequencing"* (Séraphin et al., May 2026).

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

- **R version:** 4.6.0
- **Random seed:** `20260428` (used for the calibration bootstrap and any other stochastic step).

To reproduce: clone the repository, install the listed package versions and run `run_pipeline.R` from the repository root. The pipeline regenerates derived data tables, main figures, and supplementary materials in the output directory.

## Data availability

- **Raw sequencing reads** (282 paired-end Illumina WGS samples) deposited at NCBI SRA under BioProject `PRJNA1466981`. Per-sample accessions and metadata in `Supplementary_Data_S1.csv`.
- **Reference genome:** *M. tuberculosis* H37Rv, GenBank accession NC_000962.3.
- **Processed data tables** (calibration grid output, full variant call tables, derived analytic frames) included in this repository under `data/processed/`.

## Repository structure

```
.
├── run_paper1.R                                       Orchestrator entry point
├── Makefile                                           Named targets (make all, verify, etc.)
├── paper1_README.md                                   This file
├── LICENSE                                            MIT
├── CITATION.cff                                       Machine-readable citation metadata
├── .gitignore
│
├── R/                                                 Pipeline source
│   ├── 00_pipeline_config.R                           Flags, tag construction, paths
│   ├── build_sample_metadata_v2.R                     Sample metadata harmonization
│   ├── 00_prep_metadata.R                             Cohort cleaning (Table 1)
│   ├── isnv_helpers.R                                 Shared helper functions
│   ├── 01_build_cal_pairs.R                           M0 within-visit replicate-pair construction
│   ├── 02_calibrate_lexicographic.R                   Lexicographic threshold selection (Stage 1)
│   ├── 02b_calibration_bootstrap.R                    B=1000 selection-stability bootstrap (Stage 2)
│   ├── 02c_depth_concordance_diagnostic.R             Depth-concordance diagnostic (Figure S1)
│   ├── 03_apply_thresholds.R                          Apply calibrated rule to full cohort
│   ├── 09_make_calibration_selection_figure_table.R   Figure 1, Table 1
│   ├── application_per_patient_prevalence.R           Table 2
│   ├── application_logistic_regression.R              Table 3 (logistic regression)
│   └── application_persistence_analysis.R             Table 5 (persistence), Tables S4–S5, Figure 5
│
├── manuscript/                                        Manuscript draft and supplement
│   ├── calibration_validation_manuscript.docx
│   └── supplement/
│
├── data_derived/                                      Auto-generated by `make all` (gitignored)
│   ├── 00_metadata/                                   Cohort & sample metadata
│   ├── 00_variants/                                   Filtered variant universe
│   ├── 01_calibration_<tag>/                          Calibration grid, thresholds
│   │   └── bootstrap/                                 Selection-stability bootstrap
│   └── 03_thresholded_<tag>/                          Per-sample detection indicators
│
├── outputs/                                           Final figures, tables (gitignored)
|   ├── application/<tag>/                             Figure 2
│   ├── tables/<tag>/                                  Tables 2, 4, 5
│   ├── figures/<tag>/                                 Figures 2, 5, S
│   └── supplemental/<tag>/                            Supplementary tables S4, S5
│
└── logs/                                              Per-run timestamped logs (gitignored)
```

## Citation

If you use this work, please cite:

> Séraphin MN, Afriyie-Mensah JS, Asare-Baah M, Chariker J, Domotey C, Kwarteng E, Zoungrana M, Mireku Appah S, Ganu H, Amo Omari M. Replicate-anchored calibration of within-host single nucleotide variant detection in Mycobacterium tuberculosis whole genome sequencing. DOI: `[DOI when assigned]`


## Funding

National Institutes of Health (K01AI153544); Gatorade Trust through the University of Florida Department of Medicine. The funders had no role in study design, data collection and analysis, decision to publish, or preparation of the manuscript.

## Ethics

Approved by the Institutional Review Board of the University of Florida (IRB202003042) and the Institutional Review Board of Korle-Bu Teaching Hospital, Accra, Ghana (IRB/000135/2020). All participants provided written informed consent prior to enrollment.

## Contact

Corresponding author: **Marie Nancy Séraphin**, Department of Epidemiology, College of Public Health and Health Professions, College of Medicine, University of Florida, Gainesville, FL, USA.
Email:nseraphin@ufl.edu

## License

MIT, CC-BY-4.0 — see `LICENSE` file. 
---
