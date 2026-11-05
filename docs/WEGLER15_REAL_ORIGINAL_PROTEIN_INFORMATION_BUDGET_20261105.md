# Omics2PK: measured donor-matched CYP proteins to drug intrinsic clearance — targeted vs panel vs XGBoost

**Source verification:** original Wegler et al. *Molecular Pharmaceutics* 2021, DOI
[10.1021/acs.molpharmaceut.1c00053](https://doi.org/10.1021/acs.molpharmaceut.1c00053),
PMCID PMC8041379. Original PMC AWS PDF source **SHA256**
\`5beaab1117a32bb281734504dc19237bd264b5db0857837dfd794a9319892bab\`.

**New original public execution, 2026-11-05:**
[GitHub Actions #42105373917 — SUCCESS](https://github.com/rberrah/Omics2PK/actions/runs/42105373917).
Both synthetic source-shape identity gates and original source replay passed.
The original data include **15 independent human liver donors** with same-donor
primary hepatocytes (**HH**) and liver microsomes (**HLM**), five probe drugs,
measured drug-depletion-derived CLint and distinct, quantified CYP abundances.
Donor-level values were processed transiently and never output or committed.

**Original specimen:** source HH/HLM are liver-derived, not whole blood,
PBMC or plasma metabolites. **Measured outcome:** source *in vitro* intrinsic
clearance scaled to mL/min/kg by the paper's total-protein/cell published
method (Table S3 \`eq.4\` for HH, \`eq.1\` for HLM). These are
**not CYP abundance-normalized \`eq.8\` endpoints** and therefore are
not a mathematical ratio of the tested CYP protein features.
Do not describe these rates as directly observed patient CL/F/AUC,
or CYP3A4-specific function (midazolam can involve CYP3A5).

## Method and honest evidence grade

5 **distinct assayed CYP proteins**, indexed on the same actual donor ID:
CYP3A4, CYP2C19, CYP2B6, CYP2C9, CYP2D6. These are five variables
from **one proteomic modality**, NOT five omics. Each drug's *corresponding*
CYP is selected biologically a priori, with known probe-specificity caveats.
For each drug-system, all arms use exactly the same donor-complete subset,
donor-heldout leave-one-out folds. Original HLM bufuralol has two ND CLint
rows (n13); other groups n15. Original nonpositive/unknown protein is
not silently imputed as a measured biological zero.

All outcomes are evaluated using **absolute error on log CLint**, so
numbers are not physical-unit absolute CL errors. No donor rows, source
identifiers or held-out predictions are exported.

Model arms:
- Train-only mean of log CLint (**reference**)
- The one biologically targeted CYP protein, log scaled, train-only
  standardization, fixed Ridge **alpha=1**
- All five CYP proteins, same Ridge alpha1
- All five proteins, fixed stronger Ridge alpha20
- Four nontarget proteins without the target, Ridge alpha1
- All five proteins, a **small predeclared XGBoost** (60 depth-1
  trees, learning rate .05, min child weight 3, L2 12,
  all rows/features, fixed seed), no feature selection or hyperparameter
  tuning on held-out donors.

Descriptive 95% intervals use *paired* original donor LOODO
loss resampling (1999 draws). Not nested-refit uncertainty,
not multiplicity corrected across drugs/systems/models.
Protocol follows prior analyses of this cohort, so this is an
exploratory robustness/falsification test, **not pre-registered**
and **not an external validation**.

## Original human hepatocytes: source n15 per drug

LODO **log CLint MAE**; lower is better.

| Real measured drug CLint | No protein (mean) | CYP target only (Ridge1) | 5 CYPs (Ridge1) | 5 CYPs (Ridge20) | 5 CYPs XGB |
|---|---:|---:|---:|---:|---:|
| Midazolam | 0.435032 | **0.162640** | 0.192036 | 0.267922 | 0.354374 |
| Omeprazole | 0.901056 | **0.428937** | 0.444898 | 0.512207 | 0.691177 |
| Diclofenac | 0.234435 | **0.143915** | 0.189084 | 0.199349 | 0.215862 |
| Bupropion | 0.282640 | **0.270248** | 0.356699 | 0.303159 | 0.314747 |
| Bufuralol | **0.340486** | 0.375142 | 0.477316 | 0.369753 | 0.356212 |

For **all five measured drug outcomes**, target CYP Ridge alpha1
has lower MAE than all-five-CYP Ridge alpha1. Differences
(five minus target) positive for each drug:
midazolam +0.029396,
omeprazole +0.015961,
diclofenac +0.045169,
bupropion +0.086451,
bufuralol +0.102174 log-MAE.
Descriptive person-bootstrap intervals excluding zero
(with **no adjustment for multiple tests**) are:
diclofenac [+0.007040,+0.096555] and
bupropion [+0.027544,+0.150318].
For midazolam [−0.030403,+0.093147] and omeprazole
[−0.127165,+0.165127], intervals include zero.
The XGB also worsens relative to target-only on the four
drugs for which target-only improves on no-omics;
it slightly outperforms target-only on bufuralol, yet
the no-omics mean remains numerically best for bufuralol.

**Scientific interpretation:** direct within-donor measured target
protein can be more useful for predicting an independent pharmacological
function than additional unrelated assayed CYP proteins with n15.
The success is drug-dependent; CYP2D6 abundance is not a reliable
bufuralol CLint predictor here. The 5-CYP excess error can be
finite-sample variance, correlated proteins, different functional
pathways, scaling artifacts or insufficient model tuning, rather
than intrinsic molecular-information harm.

## Human liver microsomes: source n15, bufuralol n13

| Real measured drug CLint | No protein (mean) | CYP target only (Ridge1) | 5 CYPs (Ridge1) | 5 CYPs (Ridge20) | 5 CYPs XGB |
|---|---:|---:|---:|---:|---:|
| Midazolam | 0.339026 | **0.219887** | 0.308052 | 0.279034 | 0.298768 |
| Omeprazole | 0.793596 | 0.767259 | 0.683326 | 0.678150 | **0.673489** |
| Diclofenac | 0.481246 | 0.488238 | 0.559236 | 0.465889 | **0.443843** |
| Bupropion | **0.486641** | 0.518351 | 0.733943 | 0.538996 | 0.496860 |
| Bufuralol (n13) | 0.417173 | 0.493169 | 0.518185 | 0.434201 | **0.366214** |

For HLM, five-protein Ridge1 is worse than target
for **4/5 substrates**; for omeprazole it is better.
The effect of stronger fixed shrinkage is substantial, and the
regularized tree model can outperform target-CYP on several
substrates. This is **not a general victory of white-box over
black-box**: it is model-, assay-system- and substrate-dependent.
Do not combine these five outcomes as independent donors.

## Implication for study requests and a clinical pilot

**Current positive public-source E2 result:** independently
measured, drug-specific human liver molecular abundance can predict
biochemical drug clearance on heldout same donors; for HH CYP3A4
midazolam, targeted protein substantially reduces original log
CLint error vs training mean. The added *breadth* of a five-CYP
panel need not help and can worsen prediction.

**What remains to prove clinically:** same patients with
actual **prior sampled omics + genuine CYP3A5/3A4 genotype
+ a patient-level dose/concentration/time course or
independently measured AUC/CL/F**. This cohort has no
patient systemic PK or genotype. The mechanistic tissue
signal from human liver cannot be assumed measurable from
plasma or PBMC. XGB success here is not a general
out-of-study model comparison.

**Clinical-first pilot:** reuse authorized stored plasma
from a tacrolimus-treated, PK+PGx-qualified cohort,
targeted endogenous pre-dose metabolites and relevant
inflammation/hematocrit/interaction variables;
compare model clinical+PGx alone, + mechanism-informed
panel, + unrestricted panel, and training-only nested
feature selection. Independent AUC/CL/F or time-stamped
PK predictions, same subject folds, calibration and
external clinical validation. Investigate whether
assay effects change across sample types. A collection
of omics alone or a retrospective association to
contemporaneous drug concentrations does not establish
a dose-prediction benefit.

Original-source runner:
\`scripts/benchmark_original_wegler15_target_vs_five_cyps.py\`.
Synthetic source safety tests:
\`tests/test_wegler15_target_vs_five_cyps.py\`.
Original-data CI:
\`.github/workflows/wegler15-original-5cyp-information-budget.yml\`.
