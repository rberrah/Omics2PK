# Evidence ledger — what Omics2PK has established, by level of proof

**Version 1.0.0 — 2026-11-10.** Consolidates every claim the project currently
makes, each with its evidence type, its cohort, its number and the thing it
cannot support. Nothing new was computed for this document; every figure is
read from a committed artefact or results document, and each row names it.

This is **not** the October 2026 `docs/EVIDENCE_MATRIX.md`, which is a
dataset-admissibility matrix and a gate list. This document classifies
*results*.

---

## 0. The five levels, defined before the table

| Level | Means | What it licenses |
|---|---|---|
| **R — reproduction** | a published result or a published artefact's identity recomputed from its own source | that we are using the data the authors published, and arithmetic we can check |
| **V — independent validation** | a relationship established in a cohort **not** used to form it, with a prespecified criterion | a claim about the relationship, bounded by that cohort |
| **E — exploratory** | prespecified but on a cohort this project has inspected before, or with an operationalisation chosen at analysis time | a hypothesis worth testing elsewhere; never a confirmation |
| **F — functional proof** | an in-vitro molecular→function relationship, out-of-sample | a statement about enzyme function; **never** about exposure in a patient |
| **C — clinical-PK proof** | a molecular→exposure relationship in patients, prespecified, in an independent cohort | a claim about pharmacokinetics |

**The project holds no C-level result, and this document does not contain one.**
Its single clinical-PK cohort is 47 patients it has now analysed five times;
that is E, by construction.

---

## 1. The ledger

Statuses are the ones the source documents assign. A number without an
interval is not quoted.

### 1.1 Reproduction (R)

| Claim | Cohort | Number | Source |
|---|---|---|---|
| The ARIANA file in use **is** the file the publication deposited | ARIANA n=47 | SHA-256 `7397b13e…fe37`, identical to the public repository's; 47 × 16 | `ARIANA_PHASE1_AUDIT_AND_PUBLISHED_RESULT_RECONCILIATION_20261107.md` |
| Protocol arms `P1` and `M4` are the **same model** at ARIANA's covariate set, measured rather than asserted | ARIANA n=47 | max absolute log difference **5.77 × 10⁻¹⁵** across two independent code paths | `OMICS2PK_H1_ARIANA_MECHANISTIC_HYBRID_RESULTS_20261109.md` |
| Wegler's Table S3 as recovered reproduces the article's own published statistics | Wegler n=15 | **20 of 20** Table S4 correlations agree (18 within 0.02, worst 0.0355); median midazolam `CLint,mic` **374.7** against a published **375** | `OMICS2PK_WEGLER_PER_DONOR_CLINT_QUALIFICATION_20261110.md` |

The third row is the strongest R-level evidence the project holds, and it
exists only because a transcription was checked against the source rather than
trusted.

### 1.2 Independent validation (V)

| Claim | Cohort | Number | Status |
|---|---|---|---|
| The CYP transcript adds to age and sex in predicting **in-vitro CYP activity** | **EPA n=100 donors**, independent of ARIANA and of HLC | **2 of 8 enzymes** survive Bonferroni over eight: CYP1A2 ΔMAE **−5.589** [−10.310, −1.087], CYP2B6 **−11.757** [−20.195, −3.224]. CYP2E1 is a **degradation** (+4.677). Five others not detected. | V, bounded |

This is the project's only V-level row, and it is modest on purpose: two of
eight, with one enzyme going the wrong way. The eight endpoints are **the same
100 donors**, so they are not eight cohorts.

### 1.3 Functional proof (F)

| Claim | Cohort | Number | Status |
|---|---|---|---|
| A single enzyme's **absolute protein abundance** predicts its own substrate's **intrinsic clearance**, out-of-sample | Wegler n=15, **hepatocytes** | **2 of 5** substrates after Bonferroni over five: midazolam/CYP3A4 ΔMAE(log₁₀) **−0.1168** [−0.2276, −0.0303], omeprazole/CYP2C19 **−0.2012** [−0.3646, −0.0204]. Error falls from ×1.55 to **×1.18** and from ×2.46 to **×1.55**. | **F**, exploratory |
| The same relationship in **microsomes** | Wegler n=15 | **0 of 5.** Closest is midazolam −0.0518 [−0.1352, +0.0226] | F, negative |
| A two-component summary of the remaining 4,037 microsomal proteins adds to that abundance | Wegler n=15 | **0 of 10 cells.** One interval excludes zero and it is a **degradation** (hepatocyte diclofenac +0.0196 [+0.0069, +0.0326]) | F, negative |
| The benefit differs between the two preparations, paired on the same donors | Wegler n=15 | **3 of 5** exclude zero: omeprazole +0.2095 [+0.0637, +0.3533], midazolam +0.0650 [+0.0014, +0.1279], diclofenac +0.0397 [+0.0060, +0.0717] | F, **not causal** |
| Signed CollecTRI transcription-factor activities add to the target enzyme's own transcript | HLC n=397 | **No.** `F5` vs `F1` **+9.88 MAE (+1.36 % degradation)**, [−9.03, +28.23]; the interval **excludes** a 10 % gain. The rewired null brackets `F4p` (real 827.56 inside [815.68, 845.50]) | F, negative |

**Every interval in the Wegler rows resamples stored out-of-fold errors without
refitting**, so it prices the scoring donors' variability and not the variance
of learning on 14 donors; all of them are narrower than a fully nested interval
would be, which matters most for the two positive rows.

### 1.4 Exploratory (E)

| Claim | Cohort | Number | Status |
|---|---|---|---|
| Omics add to clinical + PGx in predicting an observed metoprolol trough | ARIANA n=47 | ΔMAE **−16.97 ng/mL (−32.1 %)**, [−33.76, −4.69]; frozen ≥10 % criterion **met** | `SUPPORTED` as an exploratory success |
| Metabolomics add to clinical alone | ARIANA n=47 | **−21.31 (−38.3 %)**, Bonferroni-adjusted over 11 **[−46.81, −5.64]** | claim carried |
| Metabolomics add on top of the genotype | ARIANA n=47 | **−17.23 (−32.2 %)**, adjusted **[−43.20, −1.24]** | claim carried |
| **PGx adds on top of the metabolomics** | ARIANA n=47 | **+1.95 (+5.68 %, a degradation)**, adjusted [−2.33, +6.15]; a ≥10 % PGx benefit is **excluded** | `NOT_SUPPORTED` at the margin |
| PGx adds to clinical alone | ARIANA n=47 | not detected; intervals reach −11.6 % | `INCONCLUSIVE`, **not** negative |
| A mechanistic structure uses the genotype better than a fitted model does | ARIANA n=47 | statistical `M1a` vs mechanistic `M4`: **+14.06 (+35.65 %)** in favour of the mechanistic arm, adjusted **[+3.08, +37.17]** | detected; **compares structure *and* borrowed constants**, not representation alone |
| At a matched two-ratio panel, statistical and mechanistic treatments of the omic differ | ARIANA n=47 | **+5.02 %**, adjusted [−11.63, +19.15] | nothing detected |
| The mechanistic omic module adds to the mechanistic PGx model | ARIANA n=47 | `M5` vs `M4` **−16.8 %** [−12.28, −1.26] | `SUPPORTED` but **backbone-dependent, constant-dependent**, and `M5c` (omic alone) matches `M5` |
| A multi-omic ensemble reconstructs a masked molecular layer beyond demographics | EATRIS-Plus n=79 / 55 / 114 | **`NOT_SUPPORTED` in all four cells.** Of **60** contrasts, **one** survives both the bootstrap and the permutation null — and it is **intra-family** (amino acids → acylcarnitines, +0.0194 [+0.0071, +0.0363]) | E, negative; **no pharmacological reading** |
| Solanidine marker robustness | Müller n=20 | quantified vs semi-quantitative ICC difference **+0.362** [+0.084, +0.596], **direction opposite to expectation**; headline ICC moves 0.942 → ≈0.40 on subject inclusion | `SUPPORTED` / `PARTIALLY_SUPPORTED`, operationalisation not frozen |

### 1.5 Methodological findings that constrain everything above

| Finding | Number | Where |
|---|---|---|
| **A permutation null is positive whenever the compared arms differ in width.** Three single-family contrasts cleared a 2,000-draw bootstrap interval while lying **entirely inside a wholly positive null** (methylome +0.0091 [+0.0044, +0.0118] against a null of [+0.0034, +0.0109]) | null median **+0.0021 to +0.0040** median ΔR² in every cell, reaching +0.012 | Axis A §1.2 |
| Same-donor readiness, executed on real counts | **`fit_allowed` on 0 of 10** registered sources; the HLC+Wegler fusion refused at 466 / 15 / 488 with intersection **0** | `reports/data/hepatic_multiview_readiness_verdict.json` |
| Candidate paired sources, qualified under a declared rule | **0 of 7** usable as datasets: 2 on access, 5 on the absence of a per-donor layer paired with a per-donor outcome | `OMICS2PK_CANDIDATE_PAIRED_LIVER_SOURCE_QUALIFICATION_20261110.md` |
| Two of three published clearance equations in Wegler are **scaled by the CYP abundance** that serves as the predictor | eq. 8 divides by the probed CYP and multiplies back by the homogenate CYP; eq. 9 adds Kpuu | Axis B §3 |
| The two halves of the Wegler source number their donors differently | donors 1–9 identity, **10–15 permuted**; assuming the identity mis-pairs **6 of 15** | Axis B §3 |

---

## 2. The register of withdrawn and narrowed claims

A ledger that only grows is not a ledger. Eight claims have been withdrawn,
reversed or narrowed, all by this project on its own evidence.

| Claim | What happened | Date |
|---|---|---|
| "Protein explains CYP activity far better than mRNA: published r = 0.944 against ρ = 0.508, above the interval on 9 of 9" | **Interpretation withdrawn.** It compared a figure with an interval against one without: Grangeon publishes no per-donor value. Wegler, recomputable, gives a microsomal median ρ of **0.412**; Grangeon falls outside Wegler's interval on 2 of 5 enzymes and the HLC transcript falls **inside** it on **5 of 5** | 2026-11-10 |
| "mRNA accounts for about 26 % of rank variance, protein for about 89 %" | **Struck as an invalid statistic.** It squared a **Spearman** coefficient and read the result as a variance fraction | 2026-11-10 |
| "The omic substitutes for the genotype" (ARIANA) | **Withdrawn**: every PGx contrast has an interval covering zero | 2026-11-07 |
| H5 sampling-time effect | **Downgraded** `SUPPORTED` → `PARTIALLY_SUPPORTED` | 2026-11-07 |
| H5 "nothing here is tuned on the results" | **Corrected**: objectives and criterion families were frozen, the statistics were not, and the 29 h reference was added after a result | 2026-11-07 |
| F-series: "the rewired null shows no CollecTRI-specific information" | **Narrowed**: for `F5` the real network beats all 20 nulls by under 1 %; 20 nulls are an exploratory control, not a proof of absence | 2026-11-09 |
| Axis B: "the rest of the proteome adds nothing" | **Rescoped** to the two prespecified principal components actually tested | 2026-11-10 |
| Al-Majdoub: `ABSTAIN_TOO_FEW_DONORS` | **Corrected** to `ABSTAIN_NO_PER_DONOR_VALUES` after opening the file: abundance is measured on pooled lots, so the pairing never existed | 2026-11-10 |

---

## 3. What the body of evidence supports as a funding case

Four statements, each traceable to a row above.

1. **The pipeline runs end-to-end on real individual-level multi-omic data**,
   under protocols frozen before the first fit, on a 12-layer cohort of 125
   individuals and 127,620 features, reproducing each configuration's
   commissioned size and checksum. (Axis A)
2. **A mechanistic molecular measurement predicts function where mechanism says
   it should, out-of-sample, and fails where it should not.** One enzyme's
   abundance cuts its substrate's hepatocyte clearance error from ×1.55 to
   ×1.18 for CYP3A4 and from ×2.46 to ×1.55 for CYP2C19; it does nothing for
   CYP2B6 or CYP2D6, and nothing in microsomes. An independent 100-donor cohort
   gives the transcript analogue: 2 of 8 enzymes after Bonferroni. **Two
   sources, two layers, the same modest and enzyme-specific answer.**
3. **Measurement beats representation where data is scarce, and representation
   beats nothing where it is not.** In ARIANA, two measured metabolite ratios
   cut an observed trough's error by 38 % while the reimbursed genotype adds
   nothing on top of them; a mechanistic structure uses the *genotype* 36 %
   better than a fitted model can at n = 47, and makes no detectable difference
   to the *omic* at a matched panel.
4. **The negative results are controlled, not absent.** The project has
   executed, and published in full, a permutation grid that explains away three
   of its own bootstrap-significant contrasts; a rewired-network null that
   brackets one of its own arms; a readiness gate that refuses all ten of its
   own registered sources; and eight withdrawals of its own claims.

---

## 4. What it does not support, stated so no reader has to infer it

- **No clinical-PK proof.** No row is C-level. The single clinical cohort is 47
  patients, one drug, one observed trough, analysed five times by this project.
- **No statement about AUC, clearance in patients, steady-state exposure or
  dose individualisation.** The ARIANA endpoint is an observed pre-dose
  concentration and is never called anything else.
- **No multi-omic complementarity.** Axis A is negative and its single robust
  contrast is intra-family. Axis B's proteome contrast is negative for the two
  components tested.
- **No transfer between drugs, enzymes, tissues or cohorts.** Every positive
  row is bounded by its own cohort.
- **No causal compartment effect**, despite three paired differences excluding
  zero.
- **No equivalence between RNA and protein.** Wide intervals at n = 15 and
  n = 16 mean the *size* of the difference is unestablished, which is not the
  same as the layers being interchangeable.
- **No fusion of individuals across cohorts**, ever: the readiness gate exists
  to refuse it and has.

---

## 5. The gap, as precisely as the evidence allows

Three requirements, and the one number the project can put on any of them.

| Requirement | Why the existing data cannot meet it | Status |
|---|---|---|
| **More donors with same-donor abundance and drug-specific function** | the only qualified source has **15** | **≈150–200 donors** would establish the CYP3A4/CYP2C19 relationship under the weakest effect the pilot's interval allows (80 % power, five-substrate correction). This sizes **that contrast only** |
| **A transcriptome on the *same* donors** | no qualified source has both layers per donor; six candidates measured both and five do not publish the individual values, two are paywalled | **not sized** — needs an effect estimate the project does not have |
| **A donor-matched in-vivo endpoint** | no public source pairs a molecular layer with an exposure in the same individuals beyond ARIANA's 47 | **not sized** |

**The sharpest version of the gap is not that the data has never been
collected.** Six of seven candidate cohorts measured both a molecular layer and
either a second layer or an activity in the *same* donors. Four of them are
open-access CC BY articles with supplements that publish medians, ranges and
correlation coefficients and no per-donor row; two are behind a subscription.
Wegler was usable only because its Table S3 happens to be a per-donor table.

That points at a second instrument the project has not used: **a data-sharing
request to five identified corresponding authors** is a different thing from a
grant, and it is cheaper. The standing constraint is to use only data
accessible without a request, so this ledger names the option and leaves the
decision where it belongs.

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First consolidation. Classifies every current claim at five levels defined before the table, with its cohort, its interval and its limit; records that the project holds **no** clinical-PK result and **one** independent validation; lists the eight claims it has withdrawn, reversed or narrowed on its own evidence; and states the data gap in the sharper form the candidate-source audit established — the measurements largely exist and the individual values are largely not published. No figure was recomputed for this document. |
