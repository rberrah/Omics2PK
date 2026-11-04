# ORIGINAL NCI-60 three-omics versus single modality — real public GI50 results
**Date:** 4 November 2026. **Status:** `SOURCE_VERIFIED_ORIGINAL_REAL_CV_SUCCESS_EXPLORATORY`, in-vitro drug-response / PD, **NOT clinical PK, clearance, dose or proof of omic clinical benefit**.

## Research question and independently retrieved molecular sources
Does adding real measured protein and metabolite abundances improve the donor-disjoint (cell-line-disjoint) prediction of an **observed doxorubicin GI50 processed score** beyond a one-omic RNA representation?

- **Drug response:** original `mixOmics-org/mixOmics@9f9fd66b69105fdb618b54588b638b567e85c0d6`, `data/multidrug.rda`, doxorubicin fixed BEFORE modeling.
- **Transcriptome and metabolome:** original `Mathelab/NCI60_GeneMetabolite_Data@4c140708c18c41910458414a626b1a51d3dd2a85` (IntLIM dataset); original source after upstream curation **17,987 genes and 280 named metabolites**. NOT predose human metabolomics.
- **NEW original proteome:** Guo et al., *iScience* 2019, DOI **[10.1016/j.isci.2019.10.059](https://doi.org/10.1016/j.isci.2019.10.059)**, PMID 31733513, [PMC6889472](https://pmc.ncbi.nlm.nih.gov/articles/PMC6889472/); original Table S1E ("Averaged protein intensity") downloaded LIVE directly from [publisher supplemental mmc2.xlsx](https://ars.els-cdn.com/content/image/1-s2.0-S2589004219304407-mmc2.xlsx), **60 NCI-60 cell lines × 3,171 distinct SwissProt accession proteins**. Official original XLSX SHA-256: `e9741865cfb12cced0a471f534966e2bcd2f9c6979f591f454daf9c0206779cb`. Table source schema and file bytes separately verified in CI. The original proteome is based on SWATH measurements with **technical duplicates**, aggregated upstream in its original Table S1E.
- **Intersection:** **47 unique original cancer cell-line identities** shared by four source layers after **case and punctuation normalization only**, excluding source-discordant `MDA-MB-435` and nonmatching names. The **prefix-labeled original SWATH tissue class** (BR/CO/CNS/LC/LE/ME/OV/PR/RE) matches original independent mixOmics/IntLIM cancer lineage **for all 47**. No fuzzy or predicted cell-line identities, no joins by row order.

**Crucial:** measurements from the three molecular assays were obtained from **different studies, laboratories, collection periods and culture batches**, not from simultaneous aliquots of one specimen. A shared NCI-60 cell-line name is *not* proof of same-time multi-omics, and this is a cross-study predictive benchmark, not omics mechanism proof.

## Prespecified test

All **eleven arms** train on **identical 47 cell lines** and **identical 3 repetitions × 5 original cell-line-disjoint folds**. Each molecular modality individually receives:
1. source-only RNA / curated metabolite signal / original SWATH protein intensity (SWATH intensity transformed log10(1+x), as Table S1E is a raw-averaged intensity, not a log expression table);
2. median imputation (only where needed), feature filtering and train-fold-only variance ranking, caps **RNA 300 / metabolites 120 / proteins 300**;
3. train-fold-only standardization and **PCA 3 components** per layer;
4. Ridge **alpha=10**, fixed rather than selected based on outer-test results, tumor lineage one-hot encoded on **training** folds when used.

The common endpoint is the originally processed **in-vitro doxorubicin GI50 score**. GI50 values are not plasma concentrations, therapeutic drug exposure, microsomal CLint or PD patient survival. Each of 47 biological cell lines contributes a mean heldout loss across the three correlated validation repeats; 141 repeated test scores are **not 141 independent cell lines**.

## Actual first-run source-verified heldout results

**MAE in original processed GI50 score** (lower is better):

| Adjustments | Omics | Real 47-line heldout MAE |
|---|---|---:|
| None | **RNA only** | **0.308916** |
| None | Proteomics only | 0.364261 |
| None | Metabolomics only | 0.382598 |
| Cancer class | No omic (lineage baseline) | 0.342488 |
| Cancer class | **RNA only** | **0.334666** |
| Cancer class | Proteomics only | 0.348451 |
| Cancer class | Metabolomics only | 0.355840 |
| Cancer class | RNA + proteomics | 0.341053 |
| Cancer class | RNA + metabolomics | 0.351537 |
| Cancer class | Proteomics + metabolomics | 0.364030 |
| Cancer class | RNA + proteomics + metabolomics | **0.358124** |

**The numeric result is negative for automatic fusion**. Among *class-adjusted* arms, three layers give **+0.023458 MAE** versus class+RNA alone (roughly **7.0% higher error relative to class+RNA**), despite training all models on the **same real 47 cell lines**.

Exploratory paired **1,999-resample cell-line loss bootstrap** of three-omics minus comparator MAE (positive means omics fusion is worse):

| Three-omics versus | Delta MAE | Descriptive bootstrap 95% interval |
|---|---:|---|
| Cancer class only | +0.015636 | [−0.003562, +0.035026] |
| Cancer class + **RNA only** | **+0.023458** | **[+0.009865, +0.037185]** |
| Cancer class + proteins | +0.009673 | [−0.008200, +0.027536] |
| Cancer class + metabolites | +0.002284 | [−0.014027, +0.020005] |
| Cancer class + RNA + proteins | +0.017071 | [+0.004287, +0.029630] |
| Cancer class + RNA + metabolites | +0.006587 | [−0.003195, +0.017272] |
| Cancer class + proteins + metabolites | −0.005906 | [−0.020057, +0.008646] |

These **descriptive** intervals bootstrap already computed heldout errors, with no complete refit per bootstrap; neither family-wise inference nor independent external validation is supported. The choice of the best observed single modality is also post-hoc descriptive; do not turn its result into a multiplicity-corrected p-value.

**One nuanced observation:** RNA alone without class covariates has lower numeric MAE (0.308916) than class+RNA (0.334666), demonstrating the model and representation depend on covariate treatment, not a universal independent molecular mechanism.

## Scientifically defensible answer for Pierre Marquet

> On these 47 public NCI-60 cell lines, we independently reproduced a three-omics benchmark that compares every layer alone and each combination on identical heldout identities. Under a fixed, transparent Ridge-PCA protocol, adding proteomics and metabolomics to transcriptomics **does not improve** prediction of doxorubicin sensitivity; indeed the full fusion worsens class-adjusted MAE from **0.335 to 0.358**. This is a useful negative technical result, but not proof that multi-omics is generally inferior, because molecular samples came from different laboratories and the response is cellular **PD**, not clinical PK.

This directly motivates a **new scientific hypothesis** rather than more gratuitous feature fusion: if omics are to improve PK, source-specific measurement targets, transcript isoforms, intact functional proteins, metabolic substrate context and covariate time anchors must be qualified BEFORE integration.

## Complementary biological evidence, not a new NCI60 analysis

Published Guo 2019 same panel found wide cross-protein RNA/protein correlations, e.g. **MARCKSL1 0.93, LGALS3 0.90, ITGB1 0.88, TP53 0.14, TOP2B 0.09**. These are **the original authors' cross-sectional correlations**, NOT new independently refitted Omics2PK numbers. They demonstrate that RNA abundance is not generally interchangeable with protein abundance across genes.

For **CYP3A5**, CPIC [Birdwell et al. 2015 DOI 10.1002/cpt.113](https://pmc.ncbi.nlm.nih.gov/articles/PMC4481158/) confirms variant *3/*6/*7 can produce aberrant/truncated mRNA and lack of functional full enzyme. **Gene-level RNA counts** and **peptide-specific MS intensities** must not automatically be mapped to a functional CYP3A5-specific activity. The previously obtained original **HLC n171** CYP3A midazolam activity comparison is a partial single-SNP rs776746 vs gene-level RNA analysis, **not full PGx** and **not** isoform-resolved CYP3A5 function. Separate RNA splice junction and isoform-specific data are needed before promising GTEx CYP3A5 functional validation.

## Provenance / reproducibility

- Original complete code: `scripts/benchmark_original_nci60_threeomics_swath.py`.
- Synthetic identity/cv tests: `tests/test_nci60_threeomics_original.py`.
- 3-omics original source gate and aggregate-only CI: [run #41847348137](https://github.com/rberrah/Omics2PK/actions/runs/41847348137), **SUCCESS**; source-only aggregate artifact 13332087714.
- Structural original Table S1E audit: [run #41846796121](https://github.com/rberrah/Omics2PK/actions/runs/41846796121), **SUCCESS**.
- No original source donor/cell-line-level assay matrices, drug responses, original protein intensity vectors or individual predictions committed or uploaded; only research code and aggregate scores.

## Research next step requiring no wet-lab budget

Reproduce **on another drug** (drug fixed *a priori*) and independently characterize portability to other published NCI60 source preparations. Report the same fixed arms and account for cross-study domain shift. For **PK-specific E2**, keep original **human liver CYP→activity** HLC/EPA and Wegler/Latham preclinical drug metabolism as independent source experiments; their results **cannot be merged at a fictitious subject level**.

**Closing statement:** multiomics must earn its additional cost by improving a **measured biological or drug endpoint** relative to the best properly adjusted simple assay on identical biological units; this NCI60 test does not show that benefit under this protocol.
