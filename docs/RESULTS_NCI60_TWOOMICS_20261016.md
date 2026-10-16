# First matched-lineage cross-study NCI-60 two-omics GI50 benchmark

**2026-10-16, exploratory only.**
Verified GitHub Actions: https://github.com/rberrah/Omics2PK/actions/runs/38917337764
Additional successful run including the aggregate JSON:
https://github.com/rberrah/Omics2PK/actions/runs/38917266859

## What was measured
- **Observed** in-vitro doxorubicin GI50, not clearance, concentration or clinical PK.
- 60 mixOmics source lines; 48 punctuation/case matching keys against IntLIM
  independent microarray RNA and metabolite matrices.
- Source tumor-type labels agree for **47** lines after transparent label-only
  LEUK→LEUKEMIA, MELAN→MELANOMA, OVAR→OVARIAN abbreviations.
  **MDA-MB-435 excluded** because mixOmics records BREAST while IntLIM
  records Melanoma. 12 other source lines have no syntax-only identity match.
- **17,987 RNA features, 280 metabolites, zero missing source matrix values**
  in the 47 eligible lines, as reported by the checked exporter. Features
  were limited via TRAINING-fold-only variance ranking to 300 and 120
  respectively. The primary endpoint compound was frozen before testing.
- Five outer subject-disjoint folds repeated with **three** predefined
  shuffles; 141 OOF predictions correspond to **47** biological cell-line
  identities, NOT 141 independent specimens.

## Actual GI50 performance

The unit is the original mixOmics processed drug-response score. Results
are pooled OOF metrics across the three correlated repeat folds, not
independent external validation.

| Prespecified comparator | MAE | MAE change vs tumor-only |
|---|---:|---:|
| Training-fold mean | 0.372543 | +0.030055 |
| Tumor type alone, Ridge | 0.342488 | 0 |
| RNA 3 PCs, Ridge | **0.308916** | -0.033572 |
| Metabolites 3 PCs, Ridge | 0.382598 | +0.040111 |
| RNA + metabolite PCs, Ridge | 0.318915 | -0.023573 |
| Tumor + RNA PCs, Ridge | 0.334666 | -0.007821 |
| Tumor + metabolite PCs, Ridge | 0.355840 | +0.013352 |
| Tumor + RNA + metabolite PCs, Ridge | 0.351537 | +0.009049 |
| Tumor + top-300 RNA raw, Ridge | 0.339704 | -0.002784 |
| Tumor + top-120 metabolite raw, Ridge | 0.392615 | +0.050128 |
| Tumor + both selected raw views, Ridge | 0.341443 | -0.001044 |
| Tumor + RNA PCs + shuffled metabolite PCs, Ridge | 0.343017 | +0.000529 |
| RNA PCs, fixed XGBoost | 0.334750 | -0.007737 |
| RNA + metabolite PCs, fixed XGBoost | 0.347879 | +0.005392 |
| Tumor + RNA PCs, fixed XGBoost | 0.337037 | -0.005450 |
| Tumor + RNA + metabolite PCs, fixed XGBoost | 0.351138 | +0.008651 |

**Main observation:** adding metabolite PCs to RNA PCs increases MAE
from 0.308916 to 0.318915 (+0.009999); adjusted tumor+RNA to
tumor+RNA+metabolite MAE rises from 0.334666 to 0.351537 (+0.016871).
The RNA-only PCA result beats tumor-only numerically by -0.033572, but
adding RNA to tumor covariates produces only a small MAE reduction of
0.007821. These are exploratory point estimates, not established
incremental value or statistical significance.

## Paired incremental metabolomic contrasts (same 47 lines)

Additional qualified CI: https://github.com/rberrah/Omics2PK/actions/runs/38917922590

Positive values below mean that **adding metabolites increases MAE**.
The intervals are descriptive 95% bootstrap percentiles: 1,000 resamples
of the **47 original line identities**, after averaging the errors across
the three correlated repetitions within each line. These are not
independent external-validation intervals or formal multiplicity-adjusted
tests.

| Comparison, same CV folds | ΔMAE joint minus RNA | Descriptive 95% interval |
|---|---:|---:|
| PCA/Ridge, without tumor covariates | +0.009999 | -0.001900 to +0.022202 |
| PCA/Ridge, with tumor covariates | **+0.016871** | +0.006110 to +0.028577 |
| Fixed XGBoost, without tumor covariates | +0.013129 | -0.000456 to +0.024515 |
| Fixed XGBoost, with tumor covariates | +0.014101 | +0.000686 to +0.025326 |

These exploratory intervals favor no added metabolomic benefit under
this design; the two unadjusted descriptive intervals include zero.
They do not establish a biological disadvantage of metabolomics in
general and cannot be transported to clinical PK or independent tissues.

## Scientific interpretation and gates
1. **No observed benefit of adding metabolomics** to the matched RNA
   representation for this response/drug/limited cohort in these
   fixed predeclared models. The negative result is retained.
2. The best numerical predictor (RNA-only PCA) must not be portrayed as
   independent RNA biology or a new drug mechanism: it can reflect
   cancer-lineage, platform/batch or other confounding. Cross-study
   molecular profiling was performed under different conditions.
3. Source IntLIM arrays were previously processed using all available
   material; that cannot be undone with outer-fold preprocessing.
   Training-only feature selection and PCA mitigate downstream
   leakage but do not create independent specimen-level pairing.
4. The bootstrap interval for paired tumor-only MAE differences is
   included in the aggregate CI artifact, but 47 clustered lineages
   and one drug are too limited to establish generalization.
5. Do not extrapolate this result to original SWEET, biochemical pathway
   networks, organ/tissue mechanisms, PopPK parameters, or clinical
   dosing. No original SWEET edges or actual NLME were fitted here.

## Follow-up prior to further P3 work
- Review the **paired descriptive uncertainty intervals** and sensitivity
  to exclusion of MDA-MB-435, feature caps, tumor-stratified splits,
  pre-existing IntLIM processing and domain shift, with hypotheses
  explicitly labeled as exploratory.
- Replicate the results on independent response dataset(s) and use a
  truly paired transcriptome/metabolome/response cohort before biological
  interpretation.
- Separate the public **Yin/Battail** gene-network/observed treatment
  response study and **Wegler** matched donor CLint/proteomic task;
  neither automatically establishes an Omics2PK PopPK improvement.
- Do not integrate anything into MultiOmics/PMx Explain as a validated
  prediction option at this stage.

Public upstream code and original data are not republished here.
No patient data have been uploaded. Aggregate JSON only.
