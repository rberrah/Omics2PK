# Axis B2 results — the replication fails, and the non-refitting interval is measured to be half as wide

**Version 1.0.0 — 2026-11-10.** Evidence type: **executed analysis, exploratory**.
Protocol: `docs/OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_PRESPECIFICATION_20261110.md` v1.0.0.
Runner: `scripts/run_plasma_biomarker_validation.py`.
Artefact: `reports/data/plasma_biomarker_validation.json`.
Source: `docs/OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_20261110.md`.

---

## 0. The result

**Zero of six contrasts are SUPPORTED.**

| Endpoint | Contrast | Δ log₁₀ (refitting) | Refitting interval | Null median | Rank p | Verdict |
|---|---|---|---|---|---|---|
| midazolam 1′-OH | **incremental** V3 vs V1 | −0.0348 | [−0.333, **+0.736**] | +0.0139 | 0.0498 | `NOT_SUPPORTED` |
| midazolam 1′-OH | **replication** V1 vs V0 | −0.0445 | [−0.234, +0.345] | +0.0068 | 0.0448 | `NOT_SUPPORTED` |
| midazolam 1′-OH | substitution V2 vs V0 | −0.0505 | [−0.293, +0.316] | +0.0074 | 0.0299 | `NOT_SUPPORTED` |
| ibrutinib M37 | incremental V3 vs V1 | −0.0169 | [−0.240, +0.612] | +0.0180 | 0.1642 | `NOT_SUPPORTED` |
| ibrutinib M37 | **replication** V1 vs V0 | −0.0897 | [−0.285, +0.236] | +0.0058 | **0.0100** | `NOT_SUPPORTED` |
| ibrutinib M37 | substitution V2 vs V0 | −0.0414 | [−0.227, +0.265] | +0.0063 | 0.0348 | `NOT_SUPPORTED` |

Every point estimate is in the improving direction. Five of six clear the
permutation null. **None** has a refitting interval that excludes zero.

Two readings of that, and the difference between them matters:

- **What it is:** at fifteen donors, with the uncertainty of learning priced
  in, this cohort does not establish any of the three relationships.
- **What it is not:** a demonstration that the relationships are absent. Every
  estimate points the right way and the design's power is low by construction.

---

## 1. The project's one positive functional claim is not confirmed

The Wegler pilot reported that a target CYP's own abundance predicts its
substrate's intrinsic clearance out of sample in hepatocytes, for CYP3A4 and
CYP2C19, surviving Bonferroni. That is the project's strongest functional
result and the one the funding case leans on.

**In an independent cohort, with the same probe reaction and the same model
family, it does not reproduce at the required bar.**

| | Wegler, hepatocytes, midazolam | Lee / Jackson, hepatocytes, midazolam |
|---|---|---|
| endpoint | CL<sub>int</sub>, equation 7 | 1′-hydroxylation rate, single concentration |
| baseline fold error | 1.55 | **3.196** |
| with CYP3A4 abundance | **1.18** | **2.976** |
| verdict | SUPPORTED after Bonferroni | `NOT_SUPPORTED` |

Three things have to be said about that comparison, in order.

**First, the endpoints differ.** A fitted intrinsic clearance and a
single-concentration formation rate are not the same quantity, and the
denominators differ too (per mg protein against per million cells). This is a
replication of a *relationship and a direction*, never of a number, and a
failure to reproduce is partly a statement about endpoint comparability.

**Second, out-of-sample prediction here is poor in absolute terms.** A fold
error of 2.98 means the prediction is off by a factor of three. The in-sample
rank correlation between CYP3A4 abundance and this activity is **+0.73** — and
it buys almost nothing out of sample at n = 15. That is the clearest
demonstration this project has produced of why a published correlation is not
a performance, and it is the same lesson the brief insisted on from the start.

**Third, the direction does hold.** The abundance arm beats the mean on both
endpoints, by more on ibrutinib M37 (fold error 2.545 → 2.116, a 17 % reduction)
than on midazolam, and M37's permutation rank p is 0.0100. A reader is entitled
to conclude that something real is probably there and that fifteen donors
cannot establish it. A reader is not entitled to conclude it is established.

---

## 2. The primary contrast: the plasma metabolite does not add detectably

`incremental_value_V3_vs_V1` is the contrast no published analysis of this
source performs, and the structural question the funding case rests on.

| | midazolam 1′-OH | ibrutinib M37 |
|---|---|---|
| CYP3A4 alone, fold error | 2.976 | 2.116 |
| CYP3A4 + plasma ratio, fold error | **2.532** | **1.993** |
| Δ log₁₀, stored-error estimator | −0.0702 | −0.0259 |
| Δ log₁₀, **refitting** estimator | −0.0348 | −0.0169 |
| refitting interval | [−0.333, **+0.736**] | [−0.240, **+0.612**] |

The improvement is in the right direction on both endpoints and is not
detectable. Note how asymmetric the refitting intervals are: the upper bound
reaches **+0.74** on the primary endpoint, meaning the draws in which the
two-predictor arm was learned from a poor bootstrap sample did very badly
indeed. **Adding a second predictor to fourteen training donors is itself a
risk**, and the refitting interval is the only estimator here that shows it.

### 2.1 The one genuinely interesting number

| Arm | midazolam fold error | M37 fold error |
|---|---|---|
| mean only | 3.196 | 2.545 |
| CYP3A4 abundance alone | 2.976 | 2.116 |
| **plasma ratio alone** | **2.796** | 2.233 |
| both | 2.532 | 1.993 |

**On the primary endpoint, the plasma metabolite alone predicts hepatic
activity slightly better than the hepatic protein alone** — 2.796 against
2.976. On the secondary endpoint the ordering reverses (2.233 against 2.116),
which is itself the honest caveat: at n = 15 and with these intervals the
ordering is not established and does not even hold across the two substrates
of the same cohort. It is recorded because it is the
opposite of what a tissue-first assumption predicts, and because it is cheap to
test properly in a larger cohort — a blood sample against a liver biopsy.

It is **not** evidence that a plasma marker can replace a liver measurement.
That would need an equivalence design, which fifteen donors cannot carry, and
the protocol declared this reading forbidden before the run.

---

## 3. The reservation about the Wegler intervals, now measured

The objection raised against the Wegler pilot was that its intervals resampled
stored out-of-fold errors **without refitting**, so they price the variability
of the donors used for scoring and not the variability of learning on fifteen
individuals, and are therefore narrower than a nested interval would be.

This run computes both estimators on the same data. The ratio of widths:

| Endpoint | Contrast | Stored width | Refitting width | **Ratio** |
|---|---|---|---|---|
| midazolam | replication V1 vs V0 | 0.294 | 0.579 | **1.97** |
| midazolam | substitution V2 vs V0 | 0.342 | 0.609 | **1.78** |
| midazolam | incremental V3 vs V1 | 0.284 | 1.069 | **3.76** |
| M37 | replication V1 vs V0 | 0.249 | 0.521 | **2.09** |
| M37 | substitution V2 vs V0 | 0.252 | 0.493 | **1.96** |
| M37 | incremental V3 vs V1 | 0.135 | 0.852 | **6.32** |

**The objection is correct and its size is about a factor of two for a
one-predictor arm, and three to six for a two-predictor arm.**

The consequence has to be stated plainly rather than left for a reader to
infer: **the Wegler pilot's positive findings were declared on the narrower
estimator.** Nothing in that frozen result is modified here — it stands as
published, under the protocol it was run against. But a reader of this project
should now know that the interval type behind it is roughly half the width of a
refitting one on a comparable design, and should discount it accordingly.

The direction of the error is the unfavourable one: a narrower interval
overstates the certainty of a **positive** finding more than of a negative one.

### 3.1 And the refitting estimator is a demanding bar, which is also worth saying

Each out-of-bag draw scores only the donors the draw left out — about a third
of fifteen, so roughly five donors. An interval built that way is wide for a
reason that is partly statistical honesty and partly the price of the
estimator: with five scoring donors per draw, a real effect of the size seen
here would be hard to detect however it was estimated.

So the `NOT_SUPPORTED` verdicts should be read as **"fifteen donors cannot
settle this, on an estimator that refuses to pretend otherwise"** — and not as
"the data argues against it". The bar was declared before the run precisely so
that this reading could not be chosen after seeing the numbers. The alternative
— declaring the contrasts on the narrower estimator — would have produced a
more quotable result and a less trustworthy one.

---

## 4. Controls

### 4.1 The permutation null is positive again

| Contrast family | Null median Δ log₁₀ |
|---|---|
| midazolam, three contrasts | +0.0068 to +0.0139 |
| M37, three contrasts | +0.0058 to +0.0180 |

The null is **not** centred on zero, exactly as the Axis A run established: a
wider arm shrinks harder under an inner-CV penalty and so scores better with no
information at all. The null is computed for **all six** contrasts, not only the
ones that looked promising — which is the discipline that caused three Axis A
contrasts to be withdrawn.

Five of six observed effects clear their null. That is the part of this result
that is *not* negative: the improvements are larger than destroying the
information would produce. They are simply not larger than the uncertainty of
learning on fifteen donors.

### 4.2 Declared sensitivity: CYP3A4 + CYP3A5 instead of CYP3A4

| Arm | midazolam fold error |
|---|---|
| CYP3A4 + CYP3A5 total | 2.913 (against 2.976 for CYP3A4 alone) |
| total + plasma ratio | 2.626 |

| Contrast | Δ log₁₀ | Interval |
|---|---|---|
| total vs mean | −0.0402 | [−0.139, +0.083] |
| total + plasma vs total | −0.0451 | [−0.110, +0.021] |

Adding CYP3A5 helps marginally and changes no verdict. It carries no
multiplicity adjustment of its own and corroborates the primary reading rather
than replacing it.

---

## 5. What this run establishes, and what it does not

| Claim | Status |
|---|---|
| the enzyme→activity relationship reproduces in an independent cohort | **NOT established** |
| a plasma metabolite adds to a hepatic protein | **NOT established** |
| a plasma metabolite alone predicts hepatic activity | **NOT established**, and it ranked ahead of the protein |
| a plasma metabolite could replace a liver measurement | not addressable by this design |
| any of these relationships is absent | **NOT established either** — every estimate points the improving way |
| a non-refitting paired bootstrap is about half the width of a refitting one here | **established**, on six contrasts |
| anything about a concentration, an AUC or a clearance in a patient | out of scope: this source has no dose and no drug concentration |

**Status ceiling:** `exploratory_independent_replication`. One cohort, fifteen
donors, one enzyme family, two substrates that are not two cohorts, an in-vitro
activity rate, and procurement plasma from deceased donors.

---

## 6. What this changes for the project

1. **The evidence ledger's one positive functional row must be downgraded** from
   a result awaiting replication to a result whose **attempted replication
   failed to confirm it** — which is a different and weaker position.
2. **The sizing for question A changes.** The confirmatory figure implied by
   Wegler's effect is not the figure implied by this cohort's. Both are carried
   into `docs/OMICS2PK_COHORT_SIZING_RESULTS_20261110.md`, and the
   disagreement between them is the honest input to a funding request.
3. **A new and cheap hypothesis worth testing properly:** that an accessible
   plasma metabolite may predict hepatic enzyme function about as well as the
   hepatic protein itself. This cohort neither establishes it nor orders the
   two consistently — the plasma arm leads on midazolam and trails on M37 —
   but the question was not on the project's list before today, and a blood
   sample against a liver biopsy is the cheapest high-value contrast
   identified so far.

---

## 7. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Executed against the protocol frozen the same day. |
