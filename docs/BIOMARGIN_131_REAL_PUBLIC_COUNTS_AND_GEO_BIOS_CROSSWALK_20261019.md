# Executed source audit: real BIOMARGIN/ROCKET 131 → GEO175718 (2026-10-19)

## Actual public files downloaded, verified and processed

Source exact author code: `InsermRightLab/MOFA2-Tinel-et-al` commit `afad658913857955abc9dbac160acecda22e92c5`, training notebook Git blob `616915d949797900e211d1595b07401d64c6c93c`. `scripts/qualify_public_biomargin_131_original_source.py` downloaded original author notebook, verified **Git blob identity**, extracted **131 UNIQUE** original `BIOS` people in **4 centre strata** and EXACT six named preprocessed MOFA2 input files. **Six views are RNA/miRNA measurements in blood, kidney biopsy and urine**, not six separate omics classes with measured proteins/metabolites. The source code references processed normalized view CSVs; it does **not** itself distribute their contents. As per original 2025 publication, **35 of the 131** have complete six-view availability, not 131. [Published article](https://www.sciencedirect.com/science/article/abs/pii/S0085253825005125); [source repository](https://github.com/InsermRightLab/MOFA2-Tinel-et-al).

Actual original GEO [GSE175718](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175718) **downloaded ephemerally directly from NCBI FTP by GitHub Actions**, compressed source `GSE175718_counts.csv.gz`, recorded SHA256:
`68936911cd82a03fa15e9e935dfb3074238398e475bffd1705a7d145dd392924`.

Source counts verified as exactly:
- **58,828** gene rows, **384** numeric RNAseq sample columns.
- Each original count field checked numerically at least on sampled source gene rows for nonnegative integral raw counts; this is **not** normalized VST expression.
- The **384 source sample column labels exactly match the 384 distinct titles** from the original separately downloaded `GSE175718_series_matrix.txt.gz` metadata.
- Source count labels do **NOT** directly match the public GEO `GSM` accessions, but can be linked via exact unique GEO sample title. No biological patient IDs or raw source expressions printed in the public logs.

**Critical unqualified join:** Neither raw GEO count column labels nor any inspected GEO Series Matrix fields (sample title, description, source, characteristics) contain the author original `BIOS.xx.xxxx` IDs. Aggregate exact overlap:
- `count_matrix_col ↔ original_GEO_title` **384/384 verified**
- `GEO_metadata ↔ original_MOFA_131_BIOS` **0 verified**
- `raw_GEO ↔ original_MOFA_131_BIOS` **0 verified**

This is an unresolved **original alias crosswalk**, not a missing numerical omics dataset. **Never assume sample order corresponds between the 384 GEO and 131 MOFA2 cohort, or that the published original MOFA2 VST training inputs can be exactly reconstructed from raw GEO counts without the original preprocessing parameters.**

Execution original source:
- [Actual public author notebook + GEO counts runner](https://github.com/rberrah/Omics2PK/actions/runs/39180780817)
- [Actual GEO Series Matrix crosswalk runner](https://github.com/rberrah/Omics2PK/actions/runs/39180958200)

Both successfully completed; all inputs were processed in-memory, no original molecular matrices or individual source IDs committed, logged or published as downloadable CI artifacts. The installed code is entirely reproducible using the original public URLs and exact author notebook blob; synthetic-only unit tests exercise ambiguity rejection.

## What this PROVES versus does NOT

**PROVES** the author original cross-compartment sample count and metadata-level correspondence, plus independent raw GEO expression source dimensions and identity QC. Enables research software pre-qualification for MultiOmics RNA import, but NOT exact MOFA2 numerical reproduction until missing author-processed VST input/CEL processing and authenticated `GEO ↔ BIOS` crosswalk are available. It is NOT a genetic genotype dataset, human CYP3A4/5 tissue activity study, tacrolimus drug AUC dataset or individual patient-level PK prediction demonstration.

**DOES NOT PROVE** that blood RNA represents hepatic/intestinal CYP3A function, that the 384 GEO patients have all 6 views, that 131 BIOMARGIN samples have tacrolimus PK, that clinical ABIS records can be joined without governance, or that omic addition improves tacrolimus prediction.

## Direct tacrolimus PK benchmark access: true original sources checked

### Vanhove et al. 2016 n147, stable kidney transplant

[Original free paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5099550/), DOI `10.1111/bcp.13083`. All 147 have actual tacrolimus AUC-based `CL/F`, CYP3A4/5 PGx and same-person plasma 4beta-OHC, oral midazolam functional phenotype. 70 also had intravenous probe. Original publisher/PMC supplementary inventory offers **only `BCP13083-sup-001-t1.docx`, a genotype frequency table** (not n147 individual molecular + tac PK source). The paper reports a **positive in-sample relationship in CYP3A5 non-expressers**, but neither molecular predictor explained tac clearance in the 29 CYP3A5 expressers. No independent sample-level external superiority validation is published in the cited supplement. Person-level original input needs authorized author transfer.

### Zhu et al. 2022 n137 with 693 serial observations, liver transplant

[Original open paper](https://www.frontiersin.org/journals/pharmacology/articles/10.3389/fphar.2022.1045843/full), DOI `10.3389/fphar.2022.1045843`. Metabolomics correlates with contemporaneous tac `C0/D` measured AFTER treatment began. Authors explicitly describe repeated cross-sectional design and state that prospective **pre-dose endogenous biomarker prediction** must be verified in longitudinal study. No published matched **donor and recipient CYP3A genotype** in original model. Article/Supplementary and author inquiry offered, but no accessible full patient × metabolite × clinical tac TDM source independently retrieved or fitted here. Do NOT claim independent prospective phenotype benefit from article-only test results.

## What the BIOMARGIN–ABIS request MUST include (privacy-preserving first)

1. **Source-key table:** author `BIOS.xx.xxxx` ↔ GEO (`GSE175718` original sample **title** and `GSM`) ↔ ABIS/BIOMARGIN pseudonymous participant keys ↔ biopsy visit & blood collection time. The request must specify whether an authorized linkage even exists and test the count of exact matches, not request all raw source files again.
2. **Scientific availability table by ORIGINAL whole-person count**: `clinical+PGx+TAC C0+dose`, plus separately matched 1 blood mRNA, biopsy mRNA, blood miRNA, biopsy miRNA, urine targeted mRNA; **35 of 131 have six molecular views, but actual n matched to TAC PGx/PK is UNKNOWN**.
3. Tacrolimus dose/formulation, actual administration times, TDM sample times, concentrations (whole blood and units), earlier trough and AUC if measured, hematocrit, steroids, interactions, kidney/liver function, stability stage. Pre-landmark omic timing must be verifiable; study is biopsy-timed and treatment-exposed, so the primary question likely concerns **future PK/changes** rather than untreated pretreatment.
4. Complete CYP3A5 variants `*3/*6/*7` appropriate to ancestry and CYP3A4 `*22`, allele call/diplotype provenance, possible missingness, donor vs recipient tissue biology where relevant.
5. Ethics, consent, data use restrictions, legal permission to connect research molecular and clinical rows. **Avoid any actual patient-level data in public GitHub or CI.**
6. First release only aggregated overlap by centre, stage, modality, molecular sampling date and number of authentic concentration–time PK records. Thereafter original source on approved institutional infrastructure.

### Locked analysis once patient-level source is authorized

Primary null `M_0 = clinical+real PGx`. Primary proposed effect: held-person-out incremental error of `M_0 + train-fold-only pre-landmark molecular score`, the **same** patients, same visits and target, with external cohort validation. Secondary gene-expression/RNA-SWEET, XGBoost limited dimensionality, true factor fusion MOFA and Pillar A restricted evidence, with consistent n and handling missing modalities. Tissue-specific caution: BIOMARGIN peripheral blood expression is not hepatocyte CYP3A4 abundance; mechanistically qualified scores may measure inflammation or systemic state. No interpretation of a positive retrospective association as prospective CYP-mediated clearance without timestamp and confounder work.

## Run from real public external sources

```bash
python -m scripts.qualify_public_biomargin_131_original_source --include-geo-counts
python -m scripts.qualify_public_biomargin_geo_gsm_bios_crosswalk
```

No dependency outside Python 3.11 standard library. The scripts print only source hashes, count dimensions and aggregate match counts. Original GEO sources are not checked into this repository or retained as a CI artifact.
