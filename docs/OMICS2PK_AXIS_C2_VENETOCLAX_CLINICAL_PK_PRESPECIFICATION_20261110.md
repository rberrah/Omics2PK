# Axis C2 — a plasma metabolite against a clinical venetoclax concentration, frozen before execution

**Version 1.0.0 — 2026-11-10.** Frozen before execution.
Module: `src/omics2pk/venetoclax_clinical_pk_contrast.py`.
Runner: `scripts/run_venetoclax_clinical_pk_contrast.py`.
Artefact: `reports/data/venetoclax_clinical_pk_contrast.json`.

---

## 0. Why this exists

The evidence ledger records that this project has **no clinical-PK row at all**
beyond ARIANA: 47 patients, one drug, one observed trough, analysed five times.

A second source was found on 2026-11-10. Lee et al., *Clin Transl Sci*
2024;17:e70106 (PMC11645475, CC BY-NC), publishes for **twelve patients with
chronic lymphocytic leukaemia**, per patient and per visit:

| | Where |
|---|---|
| venetoclax daily dose and plasma concentration | main-text table |
| plasma 4β-HC, cholesterol and their ratio | **supplement Table S3** |

That is a clinical drug concentration paired with a plasma metabolite in the
same individuals — the structure the funding case is about, at a size that
cannot establish anything.

### 0.1 Not blind, and the marginal answer is already published

The article states that the dose-corrected trough was **not** significantly
correlated with the plasma 4β-HC/cholesterol ratio, and that the ratio did not
correlate with estimated trough concentrations either. **The expected result
here is therefore negative, and this protocol was written knowing that.**

What the article does not provide, and what this run adds:

- an **out-of-sample** reading under leave-one-patient-out, rather than an
  in-sample correlation;
- a **permutation null that refits the model**, to price the shrinkage
  artefact the Axis A run identified;
- an interval that **refits inside every draw**, resampling patients and not
  observations;
- a declared handling of the two not-determined concentrations.

A negative here is a **confirmation with controls**, not a discovery. Recording
it matters because it is the project's first clinical-PK result on a second
drug, and because it bounds the ARIANA claim rather than extending it.

---

## 1. The cohort, and the three things wrong with it

Twelve patients, two visits each. Two concentrations are marked not determined
(one patient had completed treatment, one is recorded `ND`), and one biomarker
row is `NA`. The declared rule: **an observation missing either the endpoint or
the predictor is excluded and counted in the artefact**, never dropped
silently and never imputed.

**(i) The baseline is the mean, not a clinical model.** The article publishes
patient characteristics — age, weight, creatinine clearance, bilirubin,
transaminases — **by sex group only**. No per-patient covariate exists. So the
comparator cannot be "clinical plus pharmacogenetics", which is what the
funding question asks. A gain measured here would be a gain over **no
information at all** — a far easier bar. **This design is biased in favour of
the biomarker**, which makes a non-detection harder to dismiss and a detection
much weaker than it would look.

**(ii) No genotype.** The contrast that matters most to the funding case — an
omic layer against pharmacogenetics — cannot be formed in this cohort at all.

**(iii) Twelve patients, median age 71.8, one drug, and three dose reductions
for toxicity.** Nothing here can establish an effect size.

---

## 2. The endpoint

**log₁₀( venetoclax plasma concentration ÷ daily dose ).**

Dividing by the dose is not a modelling choice dressed as a transform: three
patients had their dose reduced between visits, so an unnormalised
concentration would carry the prescriber's response to toxicity as if it were
pharmacokinetics.

**What it is not.** Not an AUC. Not a clearance. Not a steady-state exposure
metric. It is one timed concentration normalised by dose, and it is the
**reciprocal** of the quantity the article calls "estimated venetoclax CL" — so
it is a concentration-scale reading of the same ratio and **not an independent
endpoint**.

---

## 3. Arms and the contrast

| Arm | Predictors |
|---|---|
| `C0_mean` | none — the training-fold mean |
| `C1_plasma_metabolite` | plasma 4β-HC / cholesterol ratio |

One contrast, `C1` against `C0`. One predictor. No selection.

**Primary analysis: visit 1 only**, one observation per patient, twelve
patients. This removes the repeated-measures problem entirely rather than
modelling around it.

**Declared sensitivity: both visits**, with every observation of a held-out
patient held out together, and the patient as the resampling and permutation
unit. Two visits of one patient are not two individuals.

**Family size 2** (the primary and the sensitivity), Bonferroni at that level.

---

## 4. Validation and controls

- **Leave-one-patient-out.** The penalty is chosen by an inner leave-one-out
  among the training patients only.
- **Permutation null, 200 replicates, model refitted.** The biomarker is
  permuted **across patients, carrying each patient's whole visit vector**.
  Permuting observations would destroy the within-patient structure as well as
  the between-patient association and make the null easier to beat than the
  design deserves.
- **Two intervals.** A stored-error paired bootstrap, kept for comparability
  with the other pilots and explicitly labelled as overstating independence
  when a patient contributes two visits; and a **refitting out-of-bag interval
  that resamples patients**. The refitting one governs.
- **`SUPPORTED` requires all three:** the refitting interval excludes zero, its
  estimate is an improvement, and the permutation rank p-value is ≤ 0.05.

---

## 5. Verification of the extraction

The article publishes three figures about these patients' exposure, and the
extraction must reproduce them before anything runs:

| Published | How it is checked |
|---|---|
| estimated venetoclax CL varied **28.8-fold** between patients | CL is proportional to dose ÷ concentration, so the fold range is reproducible **exactly** from the two extracted columns |
| median **4.0 L/h**, range **1.02 to 29.4 L/h** | the single proportionality constant is recovered from the published median, then the two bounds must fall within 0.15 L/h |
| patients **3 and 7 at visit 2** have the highest CL | the two largest dose ÷ concentration ratios must be exactly those two patient-visits |

**What this proves and does not prove.** It proves the dose and concentration
columns were read correctly — the main-text table carries the patient number
only on its first row per patient, which is the same shape of trap that
mis-assigned six of Wegler's fifteen donors. It does **not** validate the
article's clearance method, which lives in a Supporting Information file this
project does not hold: the constant is recovered from the published median
rather than derived.

---

## 6. Freeze table

| Item | Value |
|---|---|
| seed | 20261110 |
| primary analysis | visit 1 only, 12 patients |
| sensitivity | both visits, patient held together |
| validation | leave-one-patient-out |
| ridge grid | 0.001, 0.01, 0.1, 1.0, 10.0, 100.0 |
| endpoint | log₁₀ (concentration ÷ daily dose) |
| primary arm / comparator | `C1_plasma_metabolite` / `C0_mean` |
| permutation replicates | 200, refitted, permuted across patients |
| bootstrap draws | 2000 stored-error, 2000 refitting out-of-bag |
| declared family size | 2 |
| missing-data rule | exclude and count; never impute, never drop silently |
| status ceiling | `exploratory_clinical_pk` |

---

## 7. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | Frozen before execution. §0.1 records that the article's marginal answer was known at freeze and that this is therefore not a blind prespecification. |
