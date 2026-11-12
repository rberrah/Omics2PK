# Axis C2 results — the project's first clinical-PK row, and it is negative

**Version 1.0.0 — 2026-11-10.** Evidence type: **executed analysis, exploratory**.
Protocol: `docs/OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_PRESPECIFICATION_20261110.md` v1.0.0.
Runner: `scripts/run_venetoclax_clinical_pk_contrast.py`.
Artefact: `reports/data/venetoclax_clinical_pk_contrast.json`.
Source: Lee J, et al. *Clin Transl Sci* 2024;17:e70106, PMCID PMC11645475, CC BY-NC.

---

## 0. The result

**A plasma CYP3A metabolite does not predict a dose-normalised venetoclax
trough out of sample, in twelve patients — against a baseline of no
information at all.**

| | Primary: visit 1, 12 patients | Sensitivity: both visits, 22 observations |
|---|---|---|
| mean-only fold error | 2.609 | 2.330 |
| with the 4β-HC/cholesterol ratio | **2.541** | **2.402** (worse) |
| Δ log₁₀, stored-error | −0.0114 [−0.0471, +0.0284] | +0.0132 [−0.0068, +0.0455] |
| Δ log₁₀, **refitting, patients resampled** | **+0.0006** [−0.111, +0.190] | **+0.0107** [−0.038, +0.132] |
| permutation null median | +0.0091 | +0.0049 |
| rank p | 0.1493 | 0.8109 |
| verdict | `NOT_SUPPORTED` | `NOT_SUPPORTED` |

The primary estimate is **+0.0006 log₁₀** — not a small gain, but no gain at
all to four decimal places. In the both-visits analysis the biomarker arm is
**worse than the mean**.

**The bar here was the lowest the project has ever set.** No per-patient
covariate is published for these patients — age, weight and creatinine
clearance appear by sex group only — so the comparator is the training-fold
mean rather than a clinical model. The biomarker had to beat *nothing*, and did
not. That makes this non-detection harder to explain away than most.

---

## 1. Why this row matters more than its size suggests

The evidence ledger recorded that the project had **no clinical-PK row at
all** beyond ARIANA. It now has one, and it points the other way.

| | ARIANA | Here |
|---|---|---|
| drug | metoprolol | venetoclax |
| endpoint | observed pre-dose concentration | dose-normalised plasma concentration |
| predictor | two **solanidine** ratios, a CYP2D6 phenotype marker | **4β-HC/cholesterol**, a CYP3A phenotype marker |
| baseline | clinical covariates, then genotype | the mean — no covariate is published |
| n | 47 | 12 |
| result | error cut by **38 %**, surviving Bonferroni over 11 | **no gain**, and a degradation with both visits |

**Two drugs, two endogenous phenotype markers, two enzyme families, opposite
answers.** That is a boundary on the ARIANA claim, and boundaries are what a
funding case needs in order not to overreach. The honest reading is *not*
"endogenous metabolite markers do not work" — n = 12, and the two designs
differ in the drug, the enzyme, the marker and the baseline. It is: **the
ARIANA result does not generalise on the evidence available, and nothing in the
public record lets us say when it would.**

### 1.1 A difference that may matter, stated as a hypothesis and not a finding

Solanidine's metabolites are formed **by the enzyme whose activity is being
inferred**, on a substrate the patient eats. 4β-hydroxycholesterol is also
formed by CYP3A, but its plasma concentration depends on cholesterol turnover,
hepatic and systemic, over weeks. The first is closer to a direct phenotype
measurement than the second. Whether that explains the difference in these two
cohorts is **not testable here**, and it is recorded because it is cheap to
design for in a future cohort — measure both marker classes in the same
patients — and expensive to discover afterwards.

---

## 2. The extraction is verified against the article's own exposure figures

| Published | Recomputed |
|---|---|
| estimated venetoclax CL varied **28.8-fold** between patients | **28.76-fold** |
| range **1.02** to **29.4** L/h | **1.02** to **29.29** L/h |
| patients **3 and 7 at visit 2** have the highest CL | patients 3 and 7 at visit 2 |

The article's clearance is proportional to dose ÷ concentration; the single
proportionality constant recovers as **8.641** from the published median of
4.0 L/h, and both bounds then fall within 0.15 L/h.

This proves the two columns were read correctly. The main-text table carries
the patient number only on its **first row per patient** — the same shape of
trap that mis-assigned six of Wegler's fifteen donors — so the reader carries
it forward and the anchors confirm it did so.

It does **not** validate the article's clearance method, which lives in a
Supporting Information file this project does not hold.

---

## 3. Controls

**The permutation null is positive again** — median +0.0091 on the primary and
+0.0049 on the sensitivity. A wider arm shrinks harder under an inner-CV
penalty and so scores better with no information, exactly as the Axis A run
established. The primary observed effect (+0.0006) sits **below** its null
median, which is another way of saying there is nothing there.

**The biomarker is permuted across patients, carrying each patient's whole
visit vector**, so the null destroys only the between-patient association being
tested and not the within-patient structure. Permuting observations would have
made the null easier to beat than this design deserves.

**The interval resamples patients, not observations**, and refits inside every
draw. The stored-error interval is reported for comparability and is explicitly
labelled as overstating independence when a patient contributes two visits.

**Two observations were excluded by a rule declared before the run** — one
patient had completed treatment and one concentration is recorded `ND` — and
both are counted in the artefact. Neither was imputed.

---

## 4. What this cannot claim

| Claim | Status |
|---|---|
| a plasma CYP3A marker predicts venetoclax exposure | **not established**; estimate is zero to four decimals |
| a plasma CYP3A marker is useless for venetoclax exposure | **not established either** — n = 12 |
| anything about an **AUC** or a **clearance** | out of scope: the endpoint is one timed concentration ÷ dose, and is the reciprocal of the article's own CL proxy |
| an omic layer fails to add **beyond clinical covariates and pharmacogenetics** | **cannot be tested here**: no per-patient covariate and no genotype is published |
| this bounds ARIANA | **yes, weakly** — a different drug, enzyme, marker and baseline |
| the in-vitro donors in the same article inform this | **no.** Those are fifteen hepatocyte donors and twenty microsome lots — different people. No analysis joins them |

**Status ceiling:** `exploratory_clinical_pk`. Twelve patients, one drug, one
visit in the primary, a baseline of no information, and a published marginal
result this run confirms rather than discovers.

---

## 5. What it changes

1. **The ledger gains a C-level row, and it is negative.** The statement "no
   clinical-PK proof" becomes "one clinical-PK contrast executed on a second
   drug, with no gain detected".
2. **The question-C sizing gains a real dispersion.** The fold error of ≈2.5 on
   a dose-normalised trough is an observed quantity, not an assumption, and it
   belongs in the sizing of any future clinical cohort.
3. **A design requirement follows directly:** a future cohort should measure
   **both** a substrate-derived phenotype marker and a turnover-dependent one
   in the same patients, because the public record now contains one of each
   with opposite results and no way to tell which difference is responsible.

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Source found hours earlier, after a note in `OMICS2PK_INDEPENDENT_REPOSITORY_SEARCH_20261110.md` wrongly inferred that these patients' biomarker values were figure-only; that sentence is struck there. |
