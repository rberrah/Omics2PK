# H1 on ARIANA — the key secondary (M5 versus M4) and the hybrid arm

**Date:** 2026-11-09
**Protocol:** 1.4.7 at execution, recorded as 1.4.8
**Artefact:** `reports/data/h1_ariana_mechanistic_hybrid.json`
**Source SHA-256:** `7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37`
**Status: `exploratory`, and not upgradeable.** This cohort has been analysed
repeatedly by this project, so no ARIANA result may be reported as confirmatory
however favourable it is.

---

## 0. What was predicted, and what it is not

The target is `c(Metoprolol)[ng/ml]`: the **observed pre-dose metoprolol
concentration** in 47 patients, one observation each.

It is **not** an AUC, **not** a clearance and **not** a steady-state exposure
metric, and it is not called one anywhere in this document or in the artefact.
Two further exclusions are structural rather than stylistic:

- `lnMR OH-metoprolol/metoprolol` is **absent from every arm**, recorded in the
  artefact as `prohibited_predictor_absent`. The metabolite ratio is not used
  to predict the concentration it is derived from.
- The MAP-Bayesian estimate is **excluded by construction**: with one
  observation per patient, conditioning on it and then predicting it is
  circular. No figure of that kind appears here or may be cited as a
  performance.
- No NLME claim is made: between-subject and residual variability are not
  *separately* identifiable from one observation per subject.

---

## 1. The headline, with its conditions attached

**The key secondary succeeds under the exploratory criterion. It does not
survive its own sensitivity analyses intact.** Both halves of that sentence are
the result.

| | M4 (mechanistic + PGx) | M5 (M4 + omic) |
|---|---|---|
| MAE (ng/mL) | 39.44 | **32.83** |
| 10 % trimmed MAE | 25.92 | **21.84** |
| log MAE | 0.607 | **0.498** |
| RMSE | 63.87 | 53.11 |
| bias | −15.20 | −14.63 |
| within 20 % | 31.9 % | 25.5 % |

- Relative MAE improvement **16.8 %**, against the frozen 10 % threshold.
- Descriptive 95 % paired patient interval of the MAE difference
  **[−12.28, −1.26]**, excluding zero.
- **No reversal** on either robustness co-criterion: trimmed MAE and log MAE
  both improve, so the status is not capped at `PARTIALLY_SUPPORTED` by the
  reversal rule.
- 34 of 47 patients improved, 13 worsened; the gain is present in all four
  endpoint quartiles (mean gain 6.1 / 6.3 / 5.5 / 8.5 ng/mL from Q1 to Q4), so
  it is not an artefact of a handful of patients or of one end of the
  concentration range.

Note the one place M5 is **worse**: the fraction of predictions within 20 % of
the observed value falls from 31.9 % to 25.5 %. M5 reduces the average and the
log error while making fewer predictions land inside a tight relative band.
That is reported because it is the number a clinician would ask for, and the
frozen primary loss is MAE, not p20.

---

## 2. Why the 16.8 % must not be quoted on its own

The arm is built from constants transferred from the literature, and the gain
is conditional on them. The sweeps below were frozen before the run, use a
reduced bootstrap, and are **descriptive only**: their purpose is to show
whether the direction survives, never to choose a value. The published values
were fixed in advance and are not changed after the sweep.

### 2.1 The structural backbone decides the magnitude

| Backbone | M4 MAE | M5 MAE | relative | interval excludes zero |
|---|---|---|---|---|
| `FDA2006_ER` (primary) | 39.44 | 32.83 | −16.8 % | yes |
| `Chai2026_IR` | 37.82 | 30.43 | −19.5 % | yes |
| `Eugene2016_IR` | 39.95 | 39.88 | **−0.2 %** | **no** |

On one of the three published backbones the key secondary's effect is
**indistinguishable from zero**. The *ordering* M5 ≤ M4 holds on all three; the
*magnitude* does not. An earlier external bundle claimed the result was
identical across backbones: the ordering claim survives, the magnitude claim
does not, and the difference matters because the Eugene backbone is a
legitimate published choice, not a stress test.

The Eugene null was checked for the obvious artefact and is genuine: zero of
ten folds sat at a search bound on any backbone, the locked baseline clearances
are interior, and the clearance-multiplier spread is the identical 9.52× on all
three. The null is a property of the backbone, not of the fitting.

### 2.2 The transferred constants decide whether the interval excludes zero

| Constant | Value | relative | interval excludes zero |
|---|---|---|---|
| f(CYP2D6) | 0.6 | −8.2 % | yes |
| | 0.7 | −10.7 % | yes |
| | **0.80 (published)** | **−16.8 %** | **yes** |
| | 0.9 | −6.7 % | no |
| ΔAS per inhibitor | 0 | −17.8 % | yes |
| | **0.54 (published)** | **−16.8 %** | **yes** |
| | 1 | −10.4 % | no |
| | 3 | −8.4 % | no |
| omic exponent γ | 0.25 | −4.6 % | no |
| | 0.50 | −9.2 % | no |
| | 0.75 | −13.7 % | yes |
| | **1.0 (frozen)** | **−16.8 %** | **yes** |
| blend weight w | 0 | 0.0 % | no |
| | 0.25 | −8.6 % | yes |
| | **0.50 (frozen)** | **−16.8 %** | **yes** |
| | 0.75 | −18.4 % | yes |
| | 1.0 | −17.7 % | no |
| salt-dose equivalence | 1 → 1.0526 | −16.8 % → −16.8 % | yes |

Read plainly: **the published constants sit at or near the most favourable
point of three of the five sweeps.** f_m = 0.80 gives the largest gain of the
four values tried; γ = 1 gives the largest of four; ΔAS = 0.54 is second of
four. That is not evidence of tuning — every one of these values was fixed from
the literature before the run, and the sweep ranges were frozen with them —
but it is exactly the configuration in which a 16.8 % figure should be
distrusted as a point estimate. The defensible statement is the direction, not
the number.

Two sweeps deserve their own line:

- **γ = 0.25 and γ = 0.50 lose the interval.** The omic enters as
  `(ratio / reference)^γ`; halving the exponent halves the gain. The frozen
  γ = 1 is the strongest of the four.
- **w = 0.75 beats the frozen w = 0.50** (18.4 % against 16.8 %). The frozen
  value is therefore not the best available, which is the right direction for a
  prespecified choice to be wrong in.

---

## 3. The finding that most deflates the mechanistic story

**M5c — the omic module substituting for the genotype entirely — matches M5.**

| Arm | MAE | log MAE | free parameters |
|---|---|---|---|
| M4 (PGx, no omic) | 39.44 | 0.607 | 1 |
| M5 (PGx + omic) | 32.83 | 0.498 | 1 |
| **M5c (omic, no PGx)** | **32.47** | **0.490** | 1 |
| M4c (backbone only, no CYP2D6 information at all) | 50.76 | 0.757 | 1 |

M5c is marginally *better* than M5 on both losses while using **no genotype at
all**. Its contrast against M4 is −17.7 %, but its interval **[−18.76, 2.95]
includes zero**, so the comparison is not established; what is established is
that adding the genotype on top of the omic buys nothing measurable here.

This is the opposite of the architecture the arm was built to defend. The
mechanistic chain was specified as genotype → activity score → clearance, with
the omic as a *correction*. On these 47 patients the omic carries the signal
and the genotype is redundant to it. The honest reading is that the solanidine
ratio is a better measurement of the same thing the activity score estimates —
which is an argument for phenotype over genotype, not an argument for the
mechanistic blend.

M4c is the one comparison that confirms the structure carries content: removing
CYP2D6 information entirely costs 11.3 ng/mL of MAE (M4 versus M4c, −22.3 %,
threshold met). But **that contrast loses its interval under Bonferroni
adjustment**, like every other secondary in its family.

---

## 4. The hybrid arm fails, and so do both PopPK estimation arms

### 4.1 HY — structural prediction plus machine learning on the residual

| Contrast | relative | threshold met | nominal interval |
|---|---|---|---|
| HY versus P1 | **+9.9 % (worse)** | no | [−13.91, 5.87] |
| HY versus M5 | **−8.2 % (HY worse than M5)** | no | [−3.80, 10.03] |
| HYc (residual on the omic block only) versus P1 | −4.0 % | no | [−12.28, 14.86] |

HY reaches MAE 35.52 against P1's 39.44 and M5's 32.83. A ridge model of P1's
log residual, fitted inside each training fold, therefore lands **between** the
published-coefficient model and the mechanistic omic arm, closer to the latter,
and beats neither by enough to meet the criterion.

One caveat is recorded in the artefact and must travel with the number: **HY
versus M5 is not a fair arm comparison.** HY's residual stage uses the primary
contrast's M3 design, which contains both solanidine ratios, while M5 uses the
primary panel of one ratio. The two arms see different amounts of omic
information, so this contrast compares designs, not methods.

### 4.2 Estimating the coefficients on 47 patients makes the model worse

| Arm | MAE | free parameters |
|---|---|---|
| P1 (published coefficients) | 39.44 | 1 |
| P2 (sex and CYP2D6 estimated) | 43.20 | 3 |
| P3 (P2 plus an estimated omic coefficient) | 42.95 | 4 |

- P2 versus P1: **+9.6 % worse**, nominal interval **[0.55, 7.60] excluding
  zero on the wrong side.**
- P3 versus P1: **+8.9 % worse**, nominal interval [−8.69, 19.66].

This is the clearest positive-direction result in the whole set, and it points
away from estimation rather than towards omics: at n = 47, **fixing the
covariate coefficients from the literature beats estimating them.** The
coefficient diagnostics say why. The estimated activity-score coefficient
collapses from 0.526 in P2 (CV 13.5 %) to 0.217 in P3 (CV **60.0 %**) once the
omic is admitted — the two predictors compete for the same variance and neither
is identified. The estimated omic coefficient, 0.249, is about one third of its
mechanistic counterpart of 0.8.

An earlier external bundle described this collapse as a coefficient with a CV
above 100 %. The collapse is real and in the stated direction; the magnitude is
not — the measured CV is 60 %.

---

## 5. Multiplicity, stated rather than implied

- One primary contrast exists for the whole programme (ARIANA, Ridge, M3 versus
  M1) and is reported elsewhere. Nothing here changes it.
- The key secondary M5 versus M4 is reported **unadjusted and labelled
  secondary**, as protocol §9 prescribes.
- Every other contrast carries a Bonferroni-adjusted interval over a family of
  seven alongside its nominal one.
- **No other secondary survives adjustment.** `M4_vs_backbone_only` and
  `P2_vs_P1` have nominal intervals excluding zero and lose them at the
  adjusted level; the remaining five never excluded zero.
- Fold repetition is never counted as additional independent subjects. A
  patient's absolute error is averaged across the two repeats before any loss or
  bootstrap.

---

## 6. Two decisions that were taken against the available result

Both are recorded in protocol change-log row 1.4.8, because a choice made after
seeing a loss is not admissible as a primary.

1. **The primary omic panel is 4-OH-solanidine alone.** A prior external bundle
   promoted the OH + SSDA pair, and the pair does perform better here
   (−19.6 % against −16.8 %). The bundle's justification was the lower observed
   loss, which is result-guided; the panel of one ratio is the mechanistically
   prespecified one, and SSDA is semi-quantitative without an authentic
   analytical standard. The pair stays as a declared sensitivity and nothing
   more. Note that the SSDA-only panel reaches −16.9 % with an interval that
   **includes** zero, so the pair's advantage is not attributable to SSDA
   carrying independent signal.
2. **The ER backbone stays primary on a pharmaceutical argument.** 47.5 and
   95 mg unit doses are the German succinate extended-release presentations.
   The argument is about the product, never about the observed loss — which
   matters, since this is also the backbone where the key secondary's interval
   excludes zero.

A third fact belongs here because the protocol's own arm table is affected:

3. **P1 and M4 are the same model at ARIANA's covariate set**, measured across
   two independent code paths with a maximum absolute log difference of
   **5.77 × 10⁻¹⁵**. Protocol §3 lists them as separate arms, which is correct
   in general and degenerate on this cohort. "HY versus P1" is therefore also
   "HY versus M4", and P3 versus P1 is the fix-versus-estimate comparison at
   identical structure. They are not reported as two pieces of evidence.

---

## 7. A silent alignment, declared

**24 of 47 patients (51.1 %) receive an interpolated relative activity that the
published table does not contain.** Their CYP2D6 activity scores — 0.25 and
1.0 — are absent from the Rüdesheim 2020 grid (0, 0.5, 1.25, 1.5, 2, 3), which
predates the post-2021 CPIC change that scores `*41` as 0.25. The runner
applies monotone linear interpolation between published grid points with
constant extrapolation beyond 3, and `activity_score_grid_audit` reports the
off-grid scores rather than aligning them silently.

An interpolated relative activity is **not a value Rüdesheim 2020 reports**. It
is a decision of this analysis, it affects half the cohort, and it is declared
rather than buried because the alternative — dropping 24 patients — would be a
larger and less visible intervention.

---

## 8. Execution record

| | |
|---|---|
| Folds | 5 outer × 2 repeats, 3 inner, seed `20261107` |
| Folds identical to the primary contrast | yes — the fold list comes from the primary runner's own `make_folds`, imported rather than reimplemented |
| Fold fingerprint | `0cd3597f0d5fc08d0e7e3746139dd324d367d020e57691a2fb211f1bc3cd2eea` |
| Bootstrap | 2,000 draws over patients |
| Folds where the inhibitor-free requirement was relaxed | 0 |
| PopPK fits failing the convergence cross-check | **0** |
| Free parameters fitted on this cohort | M4/M4c/M5/M5c/P1: 1 · P2: 3 · P3: 4 · HY: 1 structural plus the residual stage inside the training fold |

The convergence figure needed a correction during execution and it is recorded
here rather than quietly fixed. The first run reported **20 of 20 P2/P3 fits
failing**. Five optimisers — BFGS at two tolerances, L-BFGS-B, Nelder-Mead and
Powell — returned identical estimates to ten decimal places, and the solution
was a verified local minimum in 12 of 12 and 16 of 16 probed directions. The
cause was an unreachable gradient tolerance for a finite-difference objective,
not a failure to converge. The fit now establishes convergence by a
derivative-free cross-check rather than by a solver flag: maximum disagreement
between the two methods 2.13 × 10⁻⁷, losses unchanged, zero failures.

---

## 9. What this cannot establish

Carried verbatim from the artefact, because these are the limits that decide how
the result may be cited.

1. No replication cohort exists for the solanidine axis in public data, so no
   result here can be upgraded beyond `exploratory`.
2. The galenic form and the exact time since the last dose are not documented in
   the source; every prediction assumes τ = 12 h.
3. SSDA is semi-quantitative without an authentic analytical standard, so any
   panel containing it carries an unquantified identification uncertainty.
4. The omic module's validation against the OH-metoprolol phenotype is internal
   to these 47 patients, one laboratory and one assay.
5. 24 of 47 patients receive an interpolated relative activity that the
   published table does not contain.

To which this document adds three:

6. The key secondary's magnitude is backbone-dependent, and on one of three
   published backbones its effect is indistinguishable from zero.
7. The genotype is redundant to the omic on this cohort (M5c matches M5), which
   contradicts the mechanistic ordering the arm was built to express.
8. A pre-dose concentration is not an exposure. Nothing here licenses a
   statement about AUC, clearance or dose adjustment.

---

## 10. What this does establish, stated narrowly

On 47 patients, one cohort, one assay, exploratory and not upgradeable:

- Adding a solanidine-based omic correction to a literature-fixed mechanistic
  CYP2D6 model **reduces the absolute error of a predicted pre-dose metoprolol
  concentration by about a sixth**, with an interval excluding zero, improvement
  on both robustness co-criteria, and the gain spread across all four
  concentration quartiles.
- That reduction is **conditional on the structural backbone and on the
  transferred constants**, and disappears on one published backbone.
- **Fixing the covariate coefficients from the literature beats estimating them
  at this sample size**, by about 9–10 % of MAE in both estimation arms.
- **A machine-learning residual stage on the published-coefficient model does
  not beat the mechanistic omic arm.**
- The omic alone performs as well as the omic plus the genotype, so the
  mechanistic blend is not what produces the gain.
