# ARIANA results — what carries the signal, and what only represents it

**Version 1.0.0 — 2026-11-10.** Protocol:
`docs/OMICS2PK_ARIANA_STATISTICAL_VERSUS_MECHANISTIC_ABLATION_PRESPECIFICATION_20261110.md`
(v1.0.0, frozen before this run).
Artefact: `reports/data/ariana_statistical_versus_mechanistic.json`.
Status: **exploratory**. Evidence type: **new exploratory analysis of an
already-analysed cohort** — not a reproduction, not an independent validation,
not a clinical-PK proof.

Endpoint: **observed pre-dose metoprolol concentration**, ng/mL, one
observation per patient, n = 47. It is not an AUC, not a clearance and not a
steady-state exposure metric. The prohibited predictor
`lnMR OH-metoprolol/metoprolol` appears in no arm.

Folds: the statistical runner's own `make_folds`, 5 outer × 2 repeats × 3
inner, seed 20261107, fingerprint `0cd3597f0d5fc08d…`. **Identical to the
frozen mechanistic run's fingerprint**, verified in the artefact, which is what
makes a paired statistical-against-mechanistic contrast legitimate at all.

Every interval below is a paired bootstrap over **patients** (2,000 draws); a
fold repeat is never a second patient. Every claim is read at the
**Bonferroni-adjusted 99.55 % level** across the declared family of eleven.

---

## 1. The answer to the question that was asked

> *Does pharmacogenetics add information that the metabolomics does not
> already carry?*

**No, and the point estimate goes the wrong way.**

`A5_pgx_on_top_of_metabolomics` — adding the CYP2D6 activity score to an arm
that already holds the two solanidine ratios:

| | value |
|---|---|
| mean absolute error, `M2_clinical_omic` | **34.33 ng/mL** |
| mean absolute error, `M3a_clinical_genotype_omic` | **36.27 ng/mL** |
| difference | **+1.95 ng/mL (+5.68 %, a degradation)** |
| unadjusted 95 % interval | [−1.07, +4.96] |
| **adjusted 99.55 % interval** | **[−2.33, +6.15]** |
| claim | **none** |

The interval includes zero, so no effect was *detected*; that is never a proof
that there is none. But the direction of the point estimate is a degradation,
and the adjusted interval's best case for the genotype is a 2.33 ng/mL
improvement on a 34 ng/mL baseline — under 7 %. **The design can exclude a
PGx benefit of the prespecified size (10 %) on top of the metabolomics.**

The reverse direction is not symmetric at all.

---

## 2. The three-way ablation

All arms share the frozen clinical block (`Age`, `male sex`, `weight`, `BMI`,
`Daily Dose Metoprolol`), because a concentration model without dose is not a
weaker model but a wrong one (protocol §2).

| Contrast | Δ MAE | relative | adjusted 99.55 % interval | outcome |
|---|---|---|---|---|
| `A3_metabolomics_vs_clinical` | **−21.31** | **−38.30 %** | **[−46.81, −5.64]** | **claim carried** |
| `A4_both_vs_clinical` | **−19.36** | **−34.80 %** | **[−48.35, −4.93]** | **claim carried** |
| `A6_metabolomics_on_top_of_genotype` | **−17.23** | **−32.20 %** | **[−43.20, −1.24]** | **claim carried** |
| `A1_genotype_vs_clinical` | −2.14 | −3.84 % | [−6.88, +3.14] | nothing detected |
| `A2_comedication_vs_clinical` | −0.45 | −0.80 % | [−2.53, +1.54] | nothing detected |
| `A5_pgx_on_top_of_metabolomics` | +1.95 | +5.68 % | [−2.33, +6.15] | nothing detected (degradation) |

Arm-level mean absolute errors, ridge, same folds:

| Arm | MAE |
|---|---|
| `M0_clinical` | 55.64 |
| `M1b_clinical_comedication` | 55.19 |
| `M1a_clinical_genotype` | 53.50 |
| `M1_clinical_pgx` (genotype + comedication) | 52.81 |
| `M3a_clinical_genotype_omic` | 36.27 |
| `M3_clinical_pgx_omic` | 35.84 |
| `M3b_clinical_comedication_omic` | 34.45 |
| **`M2_clinical_omic`** | **34.33** |

**Read as a whole: the measurement carries the signal and the genotype does
not.** Three contrasts survive Bonferroni over eleven, and all three are
contrasts in which the two solanidine ratios enter. No contrast in which only
the genotype or only the comedication flag enters comes close. And the best
single arm in the statistical family is the one *without* the genotype.

This is the point the funding case rests on, and it is worth stating in the
narrow form the data supports: **in these 47 patients, a measured metabolite
ratio pair reduced the error of an observed trough concentration by about a
third, and the genotype that is already reimbursed in routine care did not
add to it.**

---

## 3. Statistical against mechanistic, on identical folds

Nothing in this section isolates "the representation". A mechanistic arm fits
**one** free parameter on these 47 patients and borrows every other constant
from the literature; a statistical arm fits a ridge penalty and six to eight
coefficients on the same 47. Family B therefore compares *transferred
structure plus one free parameter* against *estimation from 47 patients*. Any
sentence that drops one of those two ingredients is forbidden by the protocol.

Mechanistic arms at the **matched** panel (`oh_and_ssda`, two ratios, the same
two the statistical arms use):

| Arm | MAE |
|---|---|
| `M4c_backbone_only` | 50.76 |
| `M4_mechanistic_pgx` | 39.44 |
| `M5c_mechanistic_omic_only` | 32.69 |
| **`M5_mechanistic_pgx_omic`** | **31.69** |

| Contrast | Δ MAE | relative | adjusted 99.55 % interval | outcome |
|---|---|---|---|---|
| `B3_statistical_vs_mechanistic_genotype_only` | **+14.06** | **+35.65 %** | **[+3.08, +37.17]** | **difference detected, in favour of the mechanistic arm** |
| `B2_statistical_vs_mechanistic_with_genotype_matched` | +4.58 | +14.45 % | [−2.11, +12.66] | nothing detected |
| `B4_the_two_baselines` | +4.88 | +9.61 % | [−2.71, +13.78] | nothing detected |
| `B5_the_unmatched_panel` | +1.86 | +5.72 % | [−7.21, +18.43] | nothing detected |
| `B1_statistical_vs_mechanistic_omic_matched` | +1.64 | +5.02 % | [−11.63, +19.15] | nothing detected |

Two readings, and the difference between them matters.

**`B3` is a detected difference, and it is about the genotype.** The same
CYP2D6 activity score reduces error by 35.65 % more when it enters a
structural model with literature-transferred constants than when a ridge
estimates its coefficient from 47 patients (53.50 against 39.44 ng/mL,
adjusted interval [+3.08, +37.17], excluding zero). Under the protocol's
decision rule this carries **no claim**, because the rule only awards claims to
improvements of the *new* arm; it is nonetheless a finding, and it is filed as
`differences_detected_in_favour_of_the_baseline` rather than among the
contrasts where nothing was detected. The honest interpretation is not that
"mechanistic is better": it is that **at n = 47 borrowed knowledge substitutes
for data the cohort does not contain** — which is an argument for collecting
the data, not against it.

**`B1` detects nothing, and that is the more interesting negative.** At a
matched two-ratio panel the statistical ridge (34.33) and the mechanistic
representation (32.69) are indistinguishable: +5.02 %, adjusted interval
[−11.63, +19.15]. So the mechanistic family's advantage over its own PGx arm
(the frozen key secondary, M5 against M4, −16.8 %) is **driven by the omic
measurement entering at all**, not by how the omic is represented.

Put the two together and the pattern is coherent: representation matters for the
**genotype**, where the cohort carries almost no information to estimate from,
and does not detectably matter for the **omic**, where the measurement itself
is informative enough that two different models extract the same thing.

### 3.1 What the panel choice costs

The mechanistic protocol's primary panel is `oh_only` (one ratio). At that
panel `M5c` is 32.47 and `M5` is 32.83; at the matched `oh_and_ssda` panel they
are 32.69 and 31.69. Adding the SSDA ratio helps `M5` and very slightly hurts
`M5c`. `B5` exists so this is visible rather than buried: comparing the
statistical arms' two ratios against a one-ratio mechanistic arm would have
confounded representation with feature count, and it changes the contrast by
only +0.22 ng/mL in any case.

---

## 4. Analytical validity

- **Folds.** One fold list, the statistical runner's own, fingerprint
  `0cd3597f0d5fc08d…`, verified equal to the frozen mechanistic run's. Every
  arm in both families is scored on the same partitions.
- **No new model, no new arm, no new panel.** Every arm is imported from the
  two frozen runners; the runner refuses to start if their frozen constants,
  arm names or declared panels have moved.
- **The frozen primary contrast and key secondary are untouched.** This
  protocol restates neither.
- **Family A was not blind.** Its contrasts were already computed and visible
  in the committed `reports/data/h1_ariana_trough.json` before the protocol was
  written. Family A is an *adjusted reading of an existing computation*. This
  is recorded in the protocol (§1), in the artefact, and here.
- **Family B's direction was foreseeable.** The arm-level mean absolute errors
  of both families were already published before the protocol was written, so
  Family B is prespecified as to membership, estimator, adjustment and
  decision rule — not as to direction.
- **Multiplicity.** Bonferroni across eleven, applied to the same bootstrap
  percentiles. Bonferroni rather than the Holm the earlier note asked for,
  because it is uniformly more conservative and reuses tested code.
- **n = 47, one study, one drug, one observation per patient.** The three
  carried claims have adjusted intervals whose near end is −4.93, −5.64 and
  −1.24 ng/mL: the effects are large but their lower bounds are not, and the
  cohort cannot narrow them.
- **Repeated analysis.** ARIANA has now been analysed five times by this
  project. The status ceiling is `exploratory` and no result here may be
  reported as confirmatory or as an independent validation.
- **No patient identifier, value or prediction is in the artefact.**

---

## 5. What this does and does not support for the funding case

**Supports (exploratory, one cohort):**

1. A measured metabolite ratio pair reduced the error of an observed trough
   concentration by **38 %** relative to a clinical baseline, surviving
   Bonferroni over eleven contrasts.
2. The CYP2D6 genotype added **nothing detectable** on top of that
   measurement, and the design excludes a benefit of the prespecified 10 %
   size.
3. Where the cohort cannot support estimation — the genotype — **a structural
   model with transferred constants beat a fitted one by 36 %**, with an
   adjusted interval excluding zero.
4. Where the measurement is informative — the omic — **representation made no
   detectable difference**.

**Does not support:**

- Any statement about AUC, clearance, steady-state exposure or dose
  individualisation. The endpoint is a single observed trough.
- Any transfer to another drug, another CYP, or another cohort.
- Any claim that the mechanistic approach is superior in general: `B3`
  compares two things at once (structure *and* borrowed constants) and the
  protocol forbids collapsing them.
- Any confirmatory reading whatsoever.

**What it argues for measuring.** The one quantity that would turn point 3
from a statement about borrowed constants into a statement about biology is
**more patients with the same paired measurements**, so that the structural
model's constants can be estimated rather than transferred and the comparison
becomes like-for-like. The present result sizes that need: at n = 47 the
genotype's coefficient cannot be estimated well enough to match a literature
value, while the omic's can.

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First write-up. Reports all eleven declared contrasts with both unadjusted and Bonferroni-adjusted intervals, the three carried claims, the one difference detected in favour of a baseline, and the seven contrasts where nothing was detected. Records that Family A was not a blind prespecification, that Family B cannot isolate the representation, and that `B3` carries no claim under the decision rule while still being a finding. |
