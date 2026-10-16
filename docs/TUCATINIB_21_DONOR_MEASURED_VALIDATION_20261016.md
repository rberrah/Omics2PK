# First measured Omics→in-vitro PK proof: tucatinib vs CYP protein and CYP activities

**Date:** 2026-10-16. **Status:** REAL published donor-level numerical exploratory reanalysis, NOT clinical PopPK and NOT an independent external validation.

## Provenance and source qualification

- Original publication: *Impact of Variation in CYP3A and CYP2C8 on Tucatinib Metabolic Clearance in Human Liver Microsomes*, 2025, DOI [10.1016/j.dmd.2025.100061](https://doi.org/10.1016/j.dmd.2025.100061); PMCID [PMC13095627](https://pmc.ncbi.nlm.nih.gov/articles/PMC13095627/), **original Table 1**, accessible as a rendered table.
- The **21 individual human liver microsome donors** have matched CYP2C8 protein abundance, total CYP3A4+CYP3A5 protein abundance, CYP2C8 enzyme activity, CYP3A enzyme activity and tucatinib clearance (all original Table 1 measurements). The **22nd row** is a composite **150-donor pooled HLM preparation**, excluded from donor-level validation.
- Researchers analyzed the 21 published rows **locally**, manually transcribing numerical columns from the accessible primary table and checking four independent Spearman correlations reported in the paper. Original individual rows and any reconstructed Table 1 CSV have **not** been committed to the public repository or CI artifact; only aggregate mathematical results and method code are published. This is a reanalysis of reported factual measurements; source reuse/automated extraction rights remain subject to publisher terms. Independent reproduction requires direct consultation of the primary publication.
- Original outcome: apparent **in-vitro microsomal metabolic clearance** in µL/min/mg microsomal protein; NOT measured clinical plasma CL/F, systemic IV clearance, or a compartmental PK curve.
- Original CYP protein predictors: CYP2C8 (pmol/mg protein) and total CYP3A (CYP3A4+CYP3A5 pmol/mg protein). Orthogonal enzyme activity comparator: CYP2C8 amodiaquine conversion and CYP3A midazolam hydroxylation (both pmol/min/mg protein). No population clinical covariates, no PGx or transcriptome were available as comparable complete measured predictors.

### Original Table 1 identity check

Published Spearman (donor CYP phenotype versus measured tucatinib clearance) was compared with the independently transcribed 21 numerical observations. All four match to **three decimal places**:

| Source phenotype | Reported Spearman rho | Recomputed rho |
|---|---:|---:|
| CYP2C8 protein abundance | 0.904 | 0.904190 |
| CYP3A total protein abundance | 0.882 | 0.881455 |
| CYP2C8 enzymatic activity | 0.781 | 0.781423 |
| CYP3A enzymatic activity | 0.907 | 0.907438 |

**Source-level inconsistency found:** the average of the 21 **rounded individual Table 1 clearance values** is **9.624 µL/min/mg**, whereas the article prose states **9.24 ± 7.94**. The source also presents a pooled-lot clearance in Table 1 distinct from the narrative description. This report uses only the 21 original printed **individual** Table 1 values, does not silently alter the source or infer missing precision, and excludes the inconsistent pooled reference entirely. These discrepancies should be queried with the article authors before using the dataset for independent calibrated IVIVE.

## Predeclared (before results) model and held-out design

Validated by existing source-locked implementation `src/omics2pk/tucatinib_local_validation.py`:
1. **Training-only mean** of observed clearance (same 20 training donors).
2. Two-protein **CYP2C8 + total CYP3A**, ridge regularization alpha=1 (fixed, no hyperparameter search).
3. Two measured enzyme activities **CYP2C8 + CYP3A**, same fixed ridge alpha=1, **distinct comparator** not merged into the protein model.
4. For **each** donor, hold out all their observations, learn centering/scaling on the **20** other biological donors, predict excluded donor, clip non-physical negative predictions to zero.
5. MAE in original clearance units; paired 1,500-draw donor bootstrap for protein minus mean MAE with reproducible Python seed 20261016. This simple interval is exploratory, does not account for new-donor population transport, and doesn't replace a model-refitting nested bootstrap.
6. Prespecified sensitivity analysis excludes the two highest-clearance biological donors **and refits 19 donor-held-out models**, without optimizing the model or reporting these as a second external cohort.

## Actual measured-source held-out results

| Reference/model | Donor-LODO MAE (µL/min/mg) | Donor-LODO RMSE | Interpretation |
|---|---:|---:|---|
| Mean clearance from 20 training donors | **6.246** | 8.263 | No omic covariate |
| Two measured CYP protein concentrations | **2.041** | 2.934 | Added protein phenotype |
| Two measured CYP enzyme activities | **1.940** | 2.430 | Functional enzymatic reference |

- Incremental protein MAE compared with the no-omic reference: **−4.204 µL/min/mg**, a **67.3% relative MAE reduction** for this small cohort.
- Protein-predictor individual absolute errors improved over the baseline in **17/21** donor-held-out predictions.
- Paired donor-resampled 95% percentile interval for **protein minus baseline** difference in held-out MAE: **[−5.985, −2.403]** (1,500 repetitions, fixed seed).
- Enzyme activity minus protein two-predictor MAE: −0.101; exploratory paired bootstrap interval **[−0.783, +0.573]** (20,000 resamples, independent fixed seed). The enzyme/protein difference **cannot be claimed significant** in these 21 donors.
- Robustness excluding the two highest-clearance donors (remaining **n=19**, redoing all held-out fits): baseline MAE **4.269**, protein MAE **2.073**, enzyme activity MAE **1.846**. This exploratory sensitivity is not an external test set.

**Interpretation:** Within this original single-source in-vitro dataset, donor-specific CYP phenotypes markedly reduce donor-held-out clearance prediction error relative to using only the donor cohort mean. The difference is biologically plausible, but **not new independent validation of the original publication's enzyme relationships**, nor an omics advantage over a qualified clinical+PGx PopPK reference, nor an RNA/SWEET result. CYP3A and CYP2C8 levels are correlated, n=21 is small, original data are rounded to one decimal, and there are no external donors. Claims remain restricted to this particular 2025 experiment.

## Reproducible engine and license boundary

Reproduction on a locally curated original source Table 1 CSV:

```bash
python -m scripts.validate_local_tucatinib /local/authorized/original_table1.csv
```

The CSV requires the 21 biological donor rows and one explicitly marked pooled row with the exact columns validated by the script. The validation now independently checks that all **four** source-published CYP–clearance Spearman correlations match before emitting the source-labeled measured-model report. Neither the donor data nor their source values are included in this GitHub repository. Synthetic-only CI demonstrates algorithm and input validation, not fresh original-donor fitting.

### Next empirical gates

1. Add truly **independent** organ preparation batches/drugs/donors, including Wegler 2021 only after original donor CLint numeric link is qualified.
2. For **transcriptomics**, Peng 2025 six PHH functional outcomes are matched by donor, but GEO supplementary files are Seurat **cell metadata**, not RNA expression values; recover gene count matrices from SRA or authors first. The deposited BioProject [PRJNA1223864](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1223864) reports six experiments and **approximately 0.19 TB raw sequencing**; full raw reconstruction is a separate controlled compute/storage task, not an appropriate GitHub Actions job.
3. Test whether measured transcriptome adds information over PGx/clinical reference only in a **genuinely paired** multiomic–PK study with concentrations, dose history, sampling dates and patient IDs.
4. Keep rberrah.github.io a public educational/research interface only after compatible scientific and rights gates; no MIPD clinical dosing claim from this in-vitro result.

**Scientific classification:** E2 preclinical actual measured **protein→in-vitro apparent clearance** proof; **NOT** E3 clinical Omics→PopPK; no patient-level source rows redistributed.
