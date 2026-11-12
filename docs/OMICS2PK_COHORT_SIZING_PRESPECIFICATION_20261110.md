# Cohort sizing — protocol, frozen before execution

**Version 1.0.0 — 2026-11-10.**
Module: `src/omics2pk/cohort_size_simulation.py`.
Runner: `scripts/run_cohort_size_simulation.py`.
Artefact: `reports/data/cohort_size_simulation.json`.

---

## 0. What is being sized, and what sizing can mean here

Three questions are sized **separately**, because they have different
endpoints, different baselines and different costs. Nothing sized for one
transfers to another.

| | Question | Evidence level |
|---|---|---|
| **A** | does an enzyme's own abundance predict its activity? | **1 — transferable literature** |
| **B** | does an omic layer add beyond clinical covariates and pharmacogenetics? | **2 — declared scenario** |
| **C** | does an omic layer add for a clinical pharmacokinetic endpoint? | **2 — declared scenario** |

**The two-level rule.** Where the public record contains an effect of the same
kind, measured on the same sort of endpoint in the same sort of material, it is
used and its source is named. Where it does not, a declared grid of effects is
used and every number derived from it is a **scenario**, never an expected
effect. Each question states its level in the artefact.

**What a required n from a simulation is.** A statement about the assumed
generating process, not about biology. It is exactly as good as the assumed
effect, redundancy and noise. The protocol therefore requires a one-at-a-time
sensitivity table for every sized question, so that the assumptions which
decide the answer are visible rather than buried.

---

## 1. Why a simulation and not a formula

A closed-form sample size needs an effect expressed as a standardised
difference on a directly measured quantity. The quantity this project reports
is the **paired difference in out-of-sample absolute error between two arms**,
after a penalty chosen by inner cross-validation and, for a wide block, after
feature selection inside the training folds.

Shrinkage and selection change the distribution of that difference in ways no
closed form captures. The Axis A run demonstrated it concretely: a wider arm
scored *better* under a permutation null, with no information, purely because
it shrank harder. A formula that ignores that effect would size the wrong
quantity.

So the simulation runs **the whole pipeline** — selection inside folds,
penalty by inner cross-validation, out-of-sample scoring, multiplicity
adjustment — and measures the power of the contrast that would actually be
reported.

---

## 2. What the simulation is allowed to borrow from real data

**Only declared scalar parameters.** A dispersion, a correlation, a baseline
R². Read off published aggregates and named with their source.

**Never rows.** No individual from any source is joined to any individual from
any other. The standing prohibition on cross-cohort joins is not relaxed for a
simulation.

**Never a squared rank correlation.** A published Spearman enters through one
explicitly declared route, `spearman_to_pearson`, which applies the
Gaussian-copula mapping ρ<sub>P</sub> = 2 sin(π ρ<sub>S</sub> / 6) and records
that the mapping is exact only under that copula. Squaring a Spearman and
calling the result a fraction of explained variance is not permitted anywhere
in this module.

---

## 3. The model, the validation and the test

- Arms are the project's existing ridge arms. **No new predictive
  architecture.**
- **5-fold outer, 3-fold inner.** The executed pilots use
  leave-one-donor-out because they have fifteen donors. Across an n range from
  10 to 2000, mixing two schemes would make the required-n curve an artefact of
  the scheme, so one scheme is used throughout and the difference is recorded
  in the artefact.
- Feature selection is by absolute correlation with the target, computed on
  **training rows only**, inside every outer fold.
- The penalty is chosen by inner cross-validation on **training rows only**.
- **Detection** in a replicate requires a paired two-sided test on the
  per-individual differences in absolute error to clear α / family size **and**
  the point estimate to have the improving sign, so a significant degradation
  never counts as a detection.
- **Power** is the fraction of 200 replicates that detect.
- **Required n** is found by bisection, and the artefact records every n
  evaluated with its power, so the curve can be inspected rather than the
  single number trusted.

### 3.1 Choices that are deliberately optimistic, and are declared

Two simplifications make the layer look **better** than it would be, so every
required n they affect is a **lower bound**:

1. missing cells in the layer are filled with the column mean computed on the
   whole column rather than fold-safely;
2. the informative features carry a clean shared signal, with no batch
   structure and no confounding with the baseline beyond the declared
   redundancy.

A required n that is already large under optimistic assumptions is a stronger
argument than one that is large under pessimistic ones. That is why the
direction of each simplification is stated.

---

## 4. Question A — enzyme abundance predicts its own activity (level 1)

**Transferable inputs.** *Pharmaceuticals* 2025;18:1789 (PMC12735975, CC BY),
Table 3: Spearman correlations between each CYP's abundance and its own probe
reaction, in individual human preparations — ten pairs in hepatocytes (n = 20),
eleven in genotyped microsomes (n = 16), seven in intestinal mucosa (n = 15).

These are **in-sample rank correlations**. They are used only to set the
generating correlation of a simulation whose own **out-of-sample** performance
is then measured. They are never read as performance, and they are not a result
of this project: they are published aggregates with no per-donor row, so
nothing can be recomputed from them.

**Why this is level 1.** Same quantity, same kind of endpoint, same kind of
material as the question being sized.

**Grid.** Rank correlation 0.20, 0.40, 0.53, 0.70, 0.80, 0.90. The range spans
the published values rather than flattering them: the published hepatocyte
figures run from 0.29 (CYP2E1) to 0.80 (CYP2D6), and CYP3A4/midazolam — the
pair this project has twice tried to predict — is the **0.53** of that set,
which is why 0.53 is the reference.

**Reference scenario.** ρ<sub>S</sub> = 0.53, analytical dispersion 0.15 on
log₁₀ (about a 40 % coefficient of variation, the order of magnitude targeted
proteomics reports for absolute CYP quantification), no missingness.

**Sensitivity.** Rank correlation 0.20 / 0.70 / 0.90; analytical dispersion
0.05 / 0.30 / 0.50; missingness 5 % / 20 %.

### 4.1 Confirmatory sizing, reported separately

Discovery sizing asks how many individuals are needed to **find** an effect.
Confirmatory sizing asks how many are needed to **reject the null for one
prespecified effect**. The second is always smaller and the two are never
conflated.

The confirmatory figure is computed from the paired difference and its
dispersion **as measured** in the Lee / Jackson cohort, and the artefact
records that both come from fifteen donors, that a dispersion estimated from
fifteen is biased low more often than high, and that the figure is therefore
optimistic.

---

## 5. Question B — a layer beyond clinical covariates and pharmacogenetics (level 2)

**Why level 2, stated plainly.** No public source gives this project an
out-of-sample incremental effect for an omic layer beyond clinical covariates
and pharmacogenetics on a pharmacokinetic endpoint. ARIANA is the only internal
estimate, and it is one drug, 47 patients, one trough concentration and two
metabolite ratios. It cannot stand for "an omic layer". **No in-vitro effect is
carried in.**

**The declared grid** is 5 %, 10 %, 20 % and 30 % reduction in mean absolute
error, as the user's note specifies, labelled scenarios throughout.

`target_mae_reduction` is the reduction the layer would deliver **if measured
without error and carrying no redundancy with the baseline**. The realised
reduction is lower once analytical noise, redundancy, selection among many
features and missingness are applied — which is the reason to simulate rather
than solve.

**Reference scenario** (every number a declared assumption, none an estimate):

| Parameter | Value |
|---|---|
| baseline R² (clinical + pharmacogenetic) | 0.30 |
| assumed MAE reduction from the layer | 0.10 |
| features in the layer | 200 |
| informative features | 10 |
| redundancy with the baseline | 0.30 |
| analytical noise on each feature | 0.30 |
| missing cells | 5 % |
| features selected inside each training fold | 10 |

**Sensitivity.** Reduction 5 % / 20 %; redundancy 0 / 0.60 / 0.85;
features 20 / 2 000 / 20 000; analytical noise 0.10 / 0.60; baseline R²
0.10 / 0.50; missingness 20 %; features selected 2 / 30.

---

## 6. Question C — a clinical pharmacokinetic endpoint (level 2)

**Three endpoints, not one, and they are not interchangeable.** They differ in
how much a clinical and pharmacogenetic baseline already explains and in how
noisy the measurement is:

| Endpoint | baseline R² | noise | why |
|---|---|---|---|
| single trough concentration | 0.20 | 0.40 | assay, sampling-time and adherence error on top of between-subject variability. ARIANA's endpoint is of this kind. |
| AUC from dense sampling | 0.40 | 0.20 | averaging over a profile removes much of the sampling-time error, at the cost of a far more expensive protocol. |
| clearance from a population model | 0.50 | 0.15 | least noisy and most explained by covariates, because the model has already absorbed covariate effects — so an omic layer must beat a baseline that is partly the same model, and the apparent gain depends on the structural model chosen. |

**Explicit refusal.** No effect measured on in-vitro intrinsic clearance or on
an in-vitro activity rate is transferred to a clinical pharmacokinetic
endpoint. The in-vitro work sizes the in-vitro question only.

---

## 7. What the output may and may not be used for

| Reading | Permitted? |
|---|---|
| "under these stated assumptions, about N individuals would give 80 % power" | **yes** |
| "the assumption that decides the answer is X" | **yes**, that is the sensitivity table's purpose |
| "we need N patients" with no assumption quoted | **no** |
| a question-B or question-C number presented as a data-driven estimate | **no**, both are scenarios |
| an in-vitro effect used to size a clinical study | **no** |
| a required n read as a prediction of what a future study will find | **no** |

---

## 8. Freeze table

| Item | Value |
|---|---|
| seed | 20261110 |
| replicates per power estimate | 200 |
| outer / inner folds | 5 / 3 |
| ridge grid | 0.001, 0.01, 0.1, 1.0, 10.0, 100.0 |
| target power | 0.80 |
| α | 0.05, two-sided |
| search range for n | 8 to 2000 |
| rank-correlation grid, question A | 0.20, 0.40, 0.53, 0.70, 0.80, 0.90 |
| MAE-reduction grid, questions B and C | 0.05, 0.10, 0.20, 0.30 |
| Spearman→Pearson route | Gaussian copula, declared |

---

## 9. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | Frozen before execution. |
