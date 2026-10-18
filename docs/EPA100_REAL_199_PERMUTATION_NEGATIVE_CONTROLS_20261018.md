# Original real EPA100 source validation: sex-stratified negative-control RNA permutations (18 October 2026)

## Status and source

**REAL SOURCE ANALYSIS EXECUTED SUCCESSFULLY.** Full original-data source run: https://github.com/rberrah/Omics2PK/actions/runs/39141994686 ; US EPA official deposited source DOI https://doi.org/10.23719/1519197 ; original Liu et al. 2021 study DOI https://doi.org/10.1080/00498254.2020.1867929 . Original 60,345-byte workbook authentic SHA256 **cfc9dc334064ea61e99ca6fd39c2c8ef7c621d000a37c4da2c95f26ca85e36dd**. Individual donor qPCR RNA values and measured CYP in-vitro microsomal probe enzyme activities are independently crosswalked via **both original biological donor identifiers** in the two original worksheets; n100 distinct original human liver donors. Endpoint-specific complete n94–100.

This is **E2 original human tissue gene expression → quantitative functional metabolism**. It is **not** human subject administered dose and plasma concentrations, AUC, clearance, clinical PK or PopPK. Original observations remain ephemeral in GitHub Actions; no donor IDs or row values are committed or printed.

## Prespecified benchmark and negative control

For each of eight biologically specified CYP activity / matched CYP qPCR transcripts:

1. Comparator fixed ridge(alpha=1), age + sex only; experimental fixed ridge(alpha=1), age + sex + log1p matched gene expression. Enzyme activity response log1p transformed, reverse transformed to original **pmol/min/mg microsomal protein**. Fixed leave-one-biological-donor-out folds, training-only standardization, no learned tuning, same biological human per comparison.
2. Negative control: **199 independently pseudorandomly permuted RNA-to-donor matches within male and female strata** per CYP; keep original outcome/age/sex/subject keys fixed. Permute solely molecular donor linkage, without using any CYP activity outcomes. For each permutation, refit original ridge with each heldout donor excluded; train-only standardization. RNG seed 20261018; gene-specific stable offsets.
3. Statistic = delta MAE in original activity units: experimental RNA model minus age/sex comparator. **Lower is better**. Empirical one-sided Monte Carlo p = (1 + number of null deltas no larger than true delta)/(199+1). Lowest possible p=0.005, not 0 or less than 0.005. Bonferroni criterion is p ≤ 0.05/8 = 0.00625; eight endpoints are correlated within the same n100 humans. This is an exploratory negative-control association check, not a causal randomization trial or exact conditional exchangeability proof. Stratum-preserving shuffles do not preserve RNA-age dependence.

### Original n100 source results, not synthetic examples

| Matched CYP outcome | n donors | Baseline age/sex MAE | Age/sex + matched RNA MAE | Δ MAE RNA − baseline | RNA MAE reduction | Empirical p, 199 RNA mismatches | Below Bonferroni 0.00625? |
|---|---:|---:|---:|---:|---:|---:|---|
| CYP1A2 | 98 | 31.248882 | 25.659710 | −5.589172 | 17.89% | **0.005** | **Yes** |
| CYP2B6 | 100 | 72.296722 | 60.539867 | −11.756855 | 16.26% | **0.005** | **Yes** |
| CYP2C9 | 100 | 648.412453 | 628.739876 | −19.672577 | 3.03% | 0.015 | No |
| CYP2C19 | 94 | 72.396968 | 70.867054 | −1.529914 | 2.11% | 0.015 | No |
| CYP2D6 | 100 | 144.329993 | 126.326747 | −18.003246 | 12.47% | **0.005** | **Yes** |
| CYP2E1 | 100 | 953.020284 | 957.697013 | **+4.676728 (worse)** | −0.49% | 0.345 | No |
| CYP3A4 | 99 | 2297.026758 | 1959.749043 | −337.277715 | 14.68% | **0.005** | **Yes** |
| CYP4A11 | 98 | 611.976309 | 595.806794 | −16.169515 | 2.64% | 0.035 | No |

All MAE in **original assay activity** units pmol/min/mg protein. Do **not** average absolute MAE across incomparable enzyme activities. The eight outcomes are related observations from the same donors, not independent cohorts.

**Result:** seven of eight same-donor measured activity prediction MAE point estimates improve by adding the corresponding RNA, and **four CYPs (1A2, 2B6, 2D6, 3A4)** outperform all 199 tested within-sex donor-mismatched RNA permutations (empirical p=0.005 each) and meet the prespecified eight-endpoint Bonferroni p threshold. For CYP2E1 adding RNA slightly worsens out-of-donor prediction and its permutation p=0.345.

**Important complementary evidence:** previous paired bootstrap of out-of-fold errors with Bonferroni confidence intervals supported only three endpoints (CYP1A2, CYP2B6, CYP3A4): https://github.com/rberrah/Omics2PK/blob/main/docs/EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md . The methods test DIFFERENT inferential targets; CYP2D6's matching permutation p does **not** cancel the previous Bonferroni bootstrap CI that includes 0. Neither is proof of external transportability.

## Reproducibility and next study

- Negative-control code: src/omics2pk/epa_rna_negative_control.py
- Original source script: scripts/validate_epa100_rna_negative_controls.py
- Synthetic-only tests: tests/test_epa_rna_negative_control.py
- Original source CI: .github/workflows/epa100-rna-negative-control.yml
- Pinned study source original DOI and SHA; exactly 100 source donors and exactly eight prespecified CYP assays; strict fail/abstain if mismatch.
- The negative control is sex-stratified, not age/sex+batch+ancestry stratified. A new biological cohort with standardized target activity and independent assay/omics release is still needed.
- No patient concentration-time PK exists in this EPA original workbook. Silveira n40 author processed CSV are private and Wegler n54 publicly linked proteomics release lacks quantitative person-specific AUC/Cmax, so do not claim E3.
- Scientific implementation progress is **not** SWEET, global omics network M4, genotype-adjusted clinical PopPK or a medical device demonstration. This is a narrow, independently inspectable E2 mechanistic proxy feasibility result.

No original raw donor values or credentials appear anywhere in this public report.
