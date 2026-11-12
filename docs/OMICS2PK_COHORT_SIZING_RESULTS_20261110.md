# Cohort sizing results — three questions, two evidence levels, and one method finding

**Version 1.0.0 — 2026-11-10.** Evidence type: **simulation under declared assumptions**.
Protocol: `docs/OMICS2PK_COHORT_SIZING_PRESPECIFICATION_20261110.md` v1.0.0.
Runner: `scripts/run_cohort_size_simulation.py`, 200 replicates per power estimate.
Artefacts: `reports/data/cohort_size_simulation.json`,
`reports/data/selection_criterion_diagnostic.json`.

> **Every number here is a statement about an assumed generating process, not
> about biology.** Question A's correlation inputs come from published data and
> are labelled level 1. Questions B and C run on **declared scenarios** and
> every figure derived from them is a scenario, not an expected effect.

---

## 0. The four results worth carrying into a funding document

1. **Question A is answerable at a realistic size.** At the published
   hepatocyte rank correlation for CYP3A4/midazolam — 0.53, the pair this
   project has twice tried to predict — **≈95 donors** give 80 % power. The
   assumption that decides this is the correlation itself and nothing else: the
   required size runs from **14 donors at ρ = 0.90 to ≈900 at ρ = 0.20**.
2. **Question B is expensive, and dimensionality is what makes it expensive.**
   Under the reference scenario a 10 % error reduction needs **≈976**
   individuals and a 5 % reduction is **not reachable** below 1,200. Holding
   everything else fixed and widening the layer from 200 to 20,000 features
   drops power at n = 800 from **0.58 to 0.03**.
3. **Question C's three endpoints barely differ in required size**, which is
   the opposite of what the project assumed, and the reason is instructive.
4. **A method finding that was not asked for:** selecting candidate features by
   their correlation with the outcome **rewards redundancy** with the baseline.
   A protocol that asks whether a layer *adds* must select on the **residual**.

---

## 1. Question A — enzyme abundance predicts its own activity (level 1)

**Transferable inputs**, named and not invented: *Pharmaceuticals* 2025;18:1789
(PMC12735975, CC BY) Table 3 publishes Spearman correlations between each CYP's
abundance and its own probe reaction in individual human preparations — ten
pairs in hepatocytes (n = 20), eleven in genotyped microsomes (n = 16), seven in
intestinal mucosa (n = 15). They enter only as the generating correlation of a
simulation whose own **out-of-sample** performance is then measured, through a
declared Gaussian-copula mapping. They are never read as performance and their
square is never used as a variance fraction.

| Rank correlation | Required n for 80 % power |
|---|---|
| 0.20 | **902** |
| 0.40 | 189 |
| **0.53** — the published hepatocyte value for CYP3A4/midazolam | **95** |
| 0.70 | 45 |
| 0.80 | 28 |
| 0.90 | 14 |

**The published range is wide and the sizing follows it, not the flattering
end.** Hepatocyte values run from 0.29 (CYP2E1) to 0.80 (CYP2D6); microsomal
values from 0.65 to 0.98; intestinal mucosa is near nil for CYP3A4 (0.13
midazolam, 0.23 testosterone). So "how many donors" has no single answer: it
depends which enzyme, in which preparation.

### 1.1 One assumption dominates, and it is not the one that is easiest to fix

| Varied from the reference (ρ = 0.53, noise 0.15, no missingness) | Required n |
|---|---|
| ρ = 0.20 | **902** |
| analytical dispersion 0.50 | 125 |
| 20 % missing | 123 |
| analytical dispersion 0.30 | 106 |
| 5 % missing | 100 |
| reference | **95** |
| ρ = 0.70 | 45 |
| ρ = 0.90 | 14 |

**The correlation moves the answer by a factor of 64; everything else moves it
by about a third.** Tripling the analytical dispersion costs 11 extra donors.
Doubling the assumed correlation saves 50. A funding case should therefore
argue about *which enzyme and which preparation*, not about assay precision.

### 1.2 The measured effect is weaker than the correlation implies — read this before quoting ≈95

The same run sizes the **confirmatory** question from effects this project has
actually measured, which is a different question from discovery and is reported
separately:

| Measured in the Lee / Jackson cohort | Paired effect (log₁₀) | Dispersion | Required n |
|---|---|---|---|
| midazolam 1′-hydroxylation | −0.0309 | 0.2222 | **406** |
| ibrutinib M37 | −0.0802 | 0.2075 | **53** |

Two readings follow, and both matter.

**First, a confirmatory n from one small cohort is unstable.** The same fifteen
donors give 406 and 53 depending on which of two endpoints is chosen. A single
figure quoted from a cohort of fifteen is not a sizing; it is one draw from a
wide distribution.

**Second, and more important: the correlation-calibrated simulation is
optimistic.** It says ≈95 donors at ρ = 0.53, while the *observed*
out-of-sample effect on the same relationship implies 406. The in-sample
correlation in that cohort was **+0.73** — above the published 0.53 — and
out-of-sample prediction was still off by a factor of three. **So ≈95 is a
floor, not an estimate**, and the honest range for question A is roughly
**100 to 400 donors**, with the width coming from the gap between what a
correlation predicts and what cross-validation delivers.

---

## 2. Question B — a layer beyond clinical covariates and pharmacogenetics (level 2)

**Why these are scenarios and not estimates.** No public source gives this
project an out-of-sample incremental effect for an omic layer beyond clinical
covariates and pharmacogenetics on a pharmacokinetic endpoint. ARIANA is the
only internal estimate and it is one drug, 47 patients, one trough and two
metabolite ratios; it cannot stand for "an omic layer". **No in-vitro effect is
carried in.**

| Assumed reduction in mean absolute error | Required n |
|---|---|
| 5 % | **not reached below 1,200** |
| 10 % | **976** |
| 20 % | 529 |
| 30 % | 365 |

Reference scenario: baseline R² 0.30, 200 features of which 10 informative,
redundancy 0.30 with the baseline, analytical noise 0.30, 5 % missing, 10
features selected inside each training fold.

### 2.1 Dimensionality is the binding constraint

Reported as power at fixed cohort sizes rather than as a required n, because a
20,000-feature block costs about fifty times the reference per replicate and a
bisection would have cost more than the rest of the run together.

| Layer width | n = 100 | n = 200 | n = 400 | n = 800 |
|---|---|---|---|---|
| 200 features | 0.00 | 0.02 | 0.04 | **0.58** |
| 2,000 features | 0.00 | 0.00 | 0.00 | **0.20** |
| 20,000 features | 0.00 | 0.00 | 0.00 | **0.03** |

**A transcriptome-scale block with ten informative features and a 10 % true
gain is undetectable at 800 individuals.** Not weakly detectable — 3 % power.

This is the strongest quantitative argument the project has for its own
mechanistic approach: **restricting the candidate set before looking at the
outcome is not a stylistic preference, it is the difference between 0.58 and
0.03.** It is also a warning about the project's own Axis A negative result:
at 79 individuals and 127,620 features, that design had no power to find a
10 % effect, so its negative result is uninformative about effects of that
size and should never be quoted as evidence that layers do not complement.

### 2.2 The required-n ordering of the other assumptions is not interpretable

The reference scenario needs ≈976 individuals against a search ceiling of
1,200, so most sensitivity variations land in a region where the power curve is
flat and the bisection returns values that differ by less than Monte Carlo
error. The artefact records them; **they should not be ranked.** What is
interpretable is qualitative: **redundancy 0.85 makes the question
unreachable at any size inside the search range.** A layer that mostly repeats
the baseline cannot be rescued by recruitment.

---

## 3. Question C — a clinical pharmacokinetic endpoint (level 2)

**Three endpoints, sized separately, because they are not interchangeable.**
No effect measured on in-vitro intrinsic clearance or on an in-vitro activity
rate is transferred here.

| Endpoint | baseline R² | noise | 10 % reduction | 20 % reduction |
|---|---|---|---|---|
| single trough concentration | 0.20 | 0.40 | **1,051** | 538 |
| AUC from dense sampling | 0.40 | 0.20 | **939** | 501 |
| clearance from a population model | 0.50 | 0.15 | **929** | 501 |

**The three barely differ, and that is the finding.** The project had assumed a
cleaner endpoint would buy a much smaller cohort. It does not: halving the
measurement noise and raising the baseline R² from 0.20 to 0.50 moves the
required size by about 12 %.

The reason is structural. A cleaner endpoint reduces the error the layer could
remove *and* the error it must beat, roughly in proportion, so the
**fractional** gain the layer is asked to deliver is unchanged. Sizing responds
to the fractional effect, not to the absolute noise.

**The practical consequence is the opposite of what the project expected.**
Dense sampling for an AUC, or a population model for a clearance, has to be
justified by what the endpoint *means* — exposure is the quantity clinicians
act on — and not by a cohort-size saving, because there is almost none. A
cheaper single-concentration design at 1,051 patients is not meaningfully worse
statistically than an expensive AUC design at 939.

**Observed anchor, not an assumption:** the Axis C2 contrast measured a fold
error of ≈2.5 on a dose-normalised venetoclax trough in twelve patients
(`OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_CONTRAST_20261110.md`). That is the
scale of unexplained variation a real clinical trough carries with no covariate
model at all, and it is consistent with the 0.20 baseline R² assumed above.

---

## 4. The method finding: selection on the outcome rewards redundancy

Post-hoc, written after the frozen run returned something that looked like
noise and was not. Diagnostic:
`scripts/diagnose_selection_criterion_and_redundancy.py`,
`reports/data/selection_criterion_diagnostic.json`.

**Power is non-monotone in redundancy**, peaking at moderate redundancy:

| Power at fixed n | redundancy 0.0 | 0.3 | 0.6 | 0.85 |
|---|---|---|---|---|
| n = 800 | 0.342 | **0.554** | 0.314 | 0.012 |
| n = 1000 | 0.614 | **0.816** | 0.516 | 0.022 |

A layer with *less* overlap with the baseline is harder to detect than one with
moderate overlap. That is backwards for an incremental-value question, and the
cause is the selection criterion:

| Fraction of selected features that are genuinely informative | redundancy 0.0 | 0.3 | 0.6 | 0.85 |
|---|---|---|---|---|
| selecting on the **outcome** (what the frozen pipeline does) | 0.822 | 0.970 | 0.999 | **1.000** |
| selecting on the **residual** after the baseline | **0.920** | 0.879 | 0.659 | 0.234 |

**Selecting on the outcome rewards redundancy.** A redundant feature correlates
with the outcome through the baseline's own signal, so it is easier to find
while carrying less incremental information — selection accuracy rises to 1.000
exactly where the layer is least useful. Residual-based selection shows the
honest picture: a highly redundant layer has little unique signal, and that
criterion says so.

**The recommendation, in one line:** a protocol that tests whether a layer
*adds* to a baseline should select the layer's features on the **residual after
the baseline**, inside the training folds. Otherwise the selection step is
confounded with the question being asked.

**What this is not.** Not a biological result. It is a property of a
feature-selection criterion and would apply to any wide block in any domain.
It is recorded because it is a concrete, cheap change to a future protocol and
an expensive thing to discover after the data is collected.

---

## 5. What may and may not be quoted from this document

| Reading | Permitted? |
|---|---|
| "at the published correlation for CYP3A4/midazolam, ≈95 donors give 80 % power, and the measured effect implies up to 406" | **yes** |
| "≈95 donors would settle question A" | **no** — ≈95 is a floor; §1.2 |
| "a 10 % incremental gain needs ≈976 individuals **under the stated scenario**" | **yes** |
| "we need ≈976 patients" | **no** — the 10 % is an assumption, not an estimate |
| "a transcriptome-scale block is undetectable at n = 800 for a 10 % gain with 10 informative features" | **yes** |
| "Axis A's negative result shows layers do not complement" | **no** — that design had no power for an effect of this size; §2.1 |
| "a cleaner PK endpoint does not buy a smaller cohort" | **yes**, within these scenarios |
| any question-B or question-C figure as a data-driven estimate | **no**, both are scenarios |
| an in-vitro effect used to size a clinical study | **no**, and none was |

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. §2.2 records that the question-B sensitivity ordering is not interpretable near the search ceiling, rather than ranking Monte Carlo noise. §4 is post-hoc and labelled so. |
