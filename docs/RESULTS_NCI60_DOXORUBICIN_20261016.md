# First measured-response results — NCI-60 / Doxorubicin

2026-10-16; pinned source: mixOmics-org/mixOmics at 9f9fd66b69105fdb618b54588b638b567e85c0d6.
CI evidence: https://github.com/rberrah/Omics2PK/actions/runs/38909940147

**Real experimental source:** 60 NCI-60 cancer cell lines; 48 ABC transporters measured by RT-qPCR and 1,429 compound activity measurements. Doxorubicin (index 170) was selected BEFORE running any predictions. All 60 outcomes finite, sample SD 0.56235. These are processed *in vitro* GI50 activity values, NOT patient PK clearance, AUC or a clinically validated exposure estimate.

## Five outer held-out folds, pre-fixed algorithms

| Model | Out-of-fold MAE | RMSE |
|---|---:|---:|
| Training-fold mean | **0.427104** | 0.568506 |
| 48-gene expression + Ridge | 0.492244 | 0.613293 |
| PCA-3 + Ridge | 0.439414 | 0.574503 |
| SWEET-inspired network strengths + Ridge | 0.538276 | 0.676192 |
| Raw expression + network + Ridge | 0.429252 | **0.522368** |
| Shuffled training-network negative control | 0.616112 | 0.808492 |
| M1: raw expression + fixed XGBoost | 0.441228 | 0.569714 |
| M1bis: network only + fixed XGBoost | 0.439839 | 0.598061 |
| M2: raw + network + fixed XGBoost | 0.448226 | 0.587711 |

**Result:** none improves MAE relative to the training-only population mean. Raw-plus-network Ridge improves RMSE only, an unconfirmed result requiring further sensitivity analysis. The negative control is worse than its network comparator, but that does not independently establish biological network utility.

## What was protected

- Same outer fivefold subject-level split (seed 20261016) for all models.
- Training-fold-only median imputation and nonconstant-gene filtering.
- Every training network computed on a leave-one-out reference; held-out networks computed using the frozen complete outer training cohort.
- Fixed Ridge alpha=10 and XGBoost hyperparameters; no post-hoc fit selection.
- Experimental SWEET-inspired method, NOT numerical equivalence to original SWEET, NOT enzyme activity.
- 3 test network weights clipped to reference range; potential distribution-shift warning.
- No private human data and no individual molecular profiles in saved CI artifacts.

## Blocking limitations

1. Only 60 cell lines, not an independent clinical cohort; target is drug growth inhibition, not tacrolimus PK or dosing.
2. Tumor lineage/histology may confound ABC expression; introduce a histology-only baseline before interpretation.
3. No nested model tuning, no separate calibration dataset, and only a single split seed.
4. Independent drugs must be prespecified; never choose the best result out of 1,429 compounds. Cisplatin, paclitaxel and etoposide are candidates for independently declared exploratory secondary comparisons.
5. This data set has no paired DNA PGx, proteomics/metabolomics or human concentrations: genotype vs molecular phenotype vs full multiomics and genuine PopPK remain to be tested using other legally matched datasets.
6. Training features use LOO reference size n-1 and test features use full n; examine scale transportability.
7. The original SWEET algorithm is transductive. These results evaluate a *SWEET-inspired inductive network* only.

**Qualification statement:** first real PD pilot completed; there is NO demonstrated predictive improvement in primary MAE and NO evidence of PopPK benefit.

## Additional exploratory histology adjustment and split sensitivity

Ran **three** independent fivefold partition seeds (20261016–20261018)
on the same 60 cell lines; NOT 180 independent observations. The source
tumor class was one-hot encoded only on each outer training fold.
Source code and run: https://github.com/rberrah/Omics2PK/actions/runs/38910460262

| Model | MAE first 5-fold | MAE 3×5-fold |
|---|---:|---:|
| Training-fold mean | 0.427104 | 0.429393 |
| Tumor class + Ridge | **0.361946** | **0.365900** |
| Raw ABC expression + Ridge | 0.492244 | 0.503744 |
| PCA + Ridge | 0.439414 | 0.438126 |
| Network strength + Ridge | 0.538276 | 0.526232 |
| Raw + network + Ridge | 0.429252 | 0.444364 |
| M1 raw + XGBoost | 0.441228 | 0.447923 |
| M1bis network + XGBoost | 0.439839 | 0.448934 |
| M2 combined + XGBoost | 0.448226 | 0.452749 |
| Tumor class + raw + network + Ridge | 0.429805 | 0.435737 |

**Confounding interpretation:** tumor lineage alone predicts the observed
cell-line drug response markedly better than either unadjusted or
tumor-class-adjusted molecular/network predictors in this small
fixed-hyperparameter experiment. A successful global network score
could instead encode cancer lineage rather than any CYP/ABC-mediated
drug sensitivity. Model comparison must adjust appropriately for
tumor or physiological context before claiming network-specific signal.

**Robustness warning:** 3/60 held-out weights were clipped in the first
fivefold evaluation and 9/180 across 3 repetitions. These are not
independent patients. Pooling repeated out-of-fold predictions provides
a sensitivity summary only; no inferential confidence interval or
independent validation has yet been computed.

The table is observational exploratory evidence: negative for the primary
network-benefit hypothesis **in this particular assay/drug/dataset/
fixed-model design**, not definitive evidence about the general SWEET
algorithm or use in clinical PopPK.
