# H1–H5 — the scientific state, favourable, unfavourable and inconclusive

**Version 1.0.0 — 2026-11-10.** Evidence type: **consolidation, nothing recomputed**.
Hypotheses as frozen in `OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md` v1.4.10.
**No frozen protocol, threshold, contrast or result is modified by this document.**

Every figure here is traceable to a results document and an artefact. Where a
result is exploratory, it says so; where a cohort has been analysed repeatedly,
no result from it may be read as confirmatory however favourable.

---

## 0. The state in one table

| | Hypothesis | Favourable | Unfavourable | Verdict |
|---|---|---|---|---|
| **H1** | omics add beyond clinical + genotype | ARIANA: error cut **38 %**, Bonferroni-surviving | venetoclax n=12: **no gain** against a baseline of no information | **exploratory, one drug. Does not generalise on the evidence available** |
| **H2** | mechanism-guided selection is reliable | non-inferior in **45/80** cells; demotes a destroyed signal in **10/10** folds | better than a univariate filter in **1/80** cells; the hybrid is **worse** in 21/72 | **reliable, not superior.** The simplest statistical selector is still the arm to beat |
| **H3** | a small panel keeps the gain | two ratios carry ARIANA's whole gain; one protein carries Wegler's; 20,000 features have **3 % power** at n=800 | — | **the best-supported of the five** |
| **H4** | omics become harmful by estimation, not by information | **four** distinct degradation mechanisms now quantified | — | **supported, and mechanisms identified** |
| **H5** | layers are useful when origin, measurement and timing match the mechanism | compartment and marker class change the answer; no layer is intrinsically superior | the project's own strongest marker's ICC falls **0.942 → ≈0.40** on subject inclusion | **supported as a caution, not as a ranking** |

**The single most important sentence across all five:** the project's two
strongest positive results — ARIANA's 38 % and Wegler's protein→clearance —
each failed to reproduce when tested on a second cohort, and in both cases
every point estimate still ran the right way. That is the pattern to argue
from: **consistent direction, unestablishable magnitude.**

---

## 1. H1 — incremental value of omics

> *A mechanistically relevant molecular biomarker adds predictive information
> beyond clinical covariates and genotype.*
> Primary contrast, frozen: ARIANA **Ridge M3 vs Ridge M1**, MAE, identical folds.

### 1.1 Favourable

| Contrast | Cohort | Result |
|---|---|---|
| **Primary.** Ridge M3 vs M1 | ARIANA n=47 | ΔMAE **−16.97 ng/mL (−32.1 %)**, [−33.76, −4.69]; frozen ≥10 % criterion **met** |
| metabolomics vs clinical alone | ARIANA n=47 | **−21.31 (−38.3 %)**, Bonferroni-adjusted over 11 **[−46.81, −5.64]** |
| metabolomics on top of the genotype | ARIANA n=47 | **−17.23 (−32.2 %)**, adjusted **[−43.20, −1.24]** |
| key secondary, mechanistic M5 vs M4 | ARIANA n=47 | **−16.8 %** [−12.28, −1.26] |
| transcript adds to age and sex for in-vitro activity | **EPA n=100** | **2 of 8** enzymes after Bonferroni; CYP2E1 a degradation |

### 1.2 Unfavourable

| Contrast | Cohort | Result |
|---|---|---|
| **PGx on top of metabolomics** | ARIANA n=47 | **+1.95 (+5.68 %, a degradation)**, adjusted [−2.33, +6.15]; a ≥10 % PGx benefit is **excluded** |
| **plasma CYP3A marker → venetoclax trough** | **n=12 patients** | **`NOT_SUPPORTED`.** Effect +0.0006 log₁₀, refitting interval [−0.111, +0.190]. With both visits the arm is **worse than the mean** |
| plasma metabolite + hepatic protein → activity | Lee / Jackson n=15 | **`NOT_SUPPORTED`**, refitting interval [−0.333, +0.736] |

### 1.3 Inconclusive

- **PGx vs clinical alone** (ARIANA): not detected, intervals reach −11.6 %.
  `INCONCLUSIVE`, explicitly **not** negative.
- `M5c` (omic alone) matches `M5` (mechanistic PGx + omic), so the key
  secondary's gain is **not attributable to the mechanistic backbone**.

### 1.4 What H1 now licenses

**One drug, one cohort, exploratory, and bounded by a second drug.** ARIANA has
been analysed five times by this project; its status cannot rise above
exploratory. The venetoclax contrast is the first clinical-PK test on a second
drug and it found nothing — against the **easiest possible baseline**, because
no per-patient covariate is published for those patients.

The two differ in drug, enzyme, marker class and baseline, so neither refutes
the other. The honest statement is: **H1 is supported on one drug and does not
generalise on the evidence available.**

---

## 2. H2 — reliability of mechanism-guided selection

> *Mechanism-guided selection is reliable: it can work with few variables,
> resists a signal that is not there, and recognises where its knowledge is
> insufficient.*

### 2.1 Favourable

| Property | Evidence |
|---|---|
| **non-inferiority** | NetPK stays inside a 5 % margin of **all three** statistical arms in **45 of 80** cells (HLC, 10 endpoints × 8 panel sizes); at every panel size on CYP2B6, CYP2C8, CYP2D6 |
| **beats the heavier learners** | better than **LASSO entry order in 31/80** and **XGBoost gain in 27/80** cells, intervals excluding zero |
| **the hybrid detects a dead variable** | with the target gene's donor order permuted, the hybrid **demoted it in 10 folds of 10, on every endpoint** |
| **the hybrid repairs knowledge's failure** | better than its own blind ordering in **10 of 72** cells, repairing **4.3 % to 17.2 %** of the blind arm's loss on 7 of 8 assessable endpoints |
| **abstention works** | an empty retained set is an allowed answer and the model falls back to clinical covariates |
| **mechanistic structure uses a genotype better than a fitted model** | ARIANA `M1a` vs `M4`: **+35.65 %** in favour of the mechanistic arm, adjusted **[+3.08, +37.17]** |

### 2.2 Unfavourable

| Property | Evidence |
|---|---|
| **no superiority over the simplest selector** | NetPK beats the **univariate Spearman filter in 1 of 80** cells. The hybrid beats it in **1 of 72** and is **detectably worse in 21 of 72**, concentrated on CYP2C19 (8/8) and CYP2E1 (6/8) |
| **the repair criterion mostly fails** | met on **2 of 8** assessable endpoints (CYP2C8, CYP2C9); genuine inferiority on CYP2B6 and both CYP3A4 substrates |
| **a network ranking cannot choose the panel size** | recorded in protocol v1.4.0 as a negative finding before any activity value was read |
| **mechanistic relevance ≠ predictive relevance** | the blind arm kept choosing a gene whose information had been destroyed, and lost its entire advantage |

### 2.3 Inconclusive

- ARIANA at a **matched two-ratio panel**: statistical vs mechanistic treatment
  of the same omic differs by **+5.02 %**, adjusted [−11.63, +19.15] — nothing
  detected. So the `M1a`-vs-`M4` result above **confounds structure with
  borrowed literature constants** and is not evidence about representation
  alone.
- H2D's hybrid is `INCONCLUSIVE` at the margin on CYP2C19 and CYP2A6.

### 2.4 What H2 now licenses

**Reliable, not superior.** Mechanism-guided selection is defensible as a
*constraint* — it keeps panels small, resists decoys, abstains, and beats
LASSO and XGBoost — and it does **not** beat a univariate filter that reads the
endpoint. The hybrid rule fixes the specific failure it was designed for
(keeping a dead variable) without closing the gap to the simplest comparator.

**This is the central fact for any new method in this family**, and §2 of
`OMICS2PK_METHOD_AUDIT_AND_HYBRID_FUSING_ASSESSMENT_20261110.md` begins from it.

---

## 3. H3 — parsimony and minimal information

> *A small, mechanistically relevant panel retains most of the attainable gain.*

### 3.1 Favourable, and this is the project's best-supported hypothesis

| Evidence | Figure |
|---|---|
| ARIANA's entire gain is carried by **two** metabolite ratios | −38 % at k=2; the genotype adds nothing on top |
| Wegler's entire functional gain is carried by **one** protein | the target enzyme's own abundance; two PCA components of the other 4,037 proteins add nothing in 10 of 10 cells |
| HLC non-inferiority often holds at **k = 1** | every arm picks the same single transcript at k=1 |
| **width is actively harmful** | at a 10 % true gain with 10 informative features, power at n=800 falls **0.58 → 0.20 → 0.03** for 200 → 2,000 → 20,000 candidates |
| Axis A's one surviving contrast of 60 is **intra-family** | amino acids → acylcarnitines, +0.0194 [+0.0071, +0.0363] |

### 3.2 The reframing this produces

Parsimony is usually argued as interpretability or cost. The sizing result makes
it a **statistical necessity**: at the cohort sizes available, a wide candidate
set does not merely add noise, it removes the ability to detect a real effect at
all. **Restricting the candidate set before looking at the outcome is the
difference between 0.58 and 0.03 power.**

### 3.3 Unfavourable / inconclusive

- `k80`/`k90` are **undefined** wherever no panel improves on the baseline, and
  the protocol forbids manufacturing a parsimony score there. That applies to
  every Axis A cell and to the microsomal Wegler arms.
- No evidence yet distinguishes *mechanistic* parsimony from *statistical*
  parsimony: §2.2 shows the univariate filter reaches the same small panels.

---

## 4. H4 — when omics become harmful

> *Degradation comes not from the extra information itself but from estimating
> and exploiting it with finite samples.*
> **Frozen prohibition:** no universal `p/n ≈ 2` boundary may be asserted;
> synthetic experiments need ≥ 100 repetitions.

### 4.1 Four degradation mechanisms, each now quantified

| # | Mechanism | Quantified as |
|---|---|---|
| **1** | **Differential shrinkage.** A wider arm shrinks harder under an inner-CV penalty and therefore scores better *with no information at all* | permutation-null medians **+0.0021 to +0.0040** ΔR², reaching **+0.012**; three of the project's own bootstrap-significant contrasts withdrawn on this ground |
| **2** | **Dimensionality.** Candidate width destroys detectability before it destroys fit | power at n=800: **0.58 → 0.03** from 200 to 20,000 candidates |
| **3** | **Selection on the outcome rewards redundancy.** A redundant feature correlates with the outcome *through the baseline's own signal* | informative fraction of the selected set rises to **1.000** at redundancy 0.85 selecting on the outcome, and falls to **0.234** selecting on the residual. Power is **non-monotone** in redundancy, peaking at 0.3 |
| **4** | **Interval understatement.** A non-refitting paired bootstrap prices the scoring donors and not the learning | refitting intervals are **1.78–2.09×** wider for one predictor and **3.76–6.32×** for two |

All four are **method** failures, not data failures — which is exactly what H4
predicted. Mechanism 3 was not anticipated by the protocol and is the only one
with a one-line fix: **residualise before selecting.**

### 4.2 Observed degradations in real arms

- `F5` signed transcription-factor activities vs the target's own transcript:
  **+9.88 MAE (+1.36 %)**, [−9.03, +28.23]; the interval **excludes** a 10 % gain.
- PGx on top of metabolomics: **+5.68 %**.
- Axis B's one interval excluding zero is a **degradation** (hepatocyte
  diclofenac +0.0196 [+0.0069, +0.0326]).
- venetoclax both-visits: biomarker arm **worse** than the mean.

### 4.3 What is still missing for H4

The protocol requires an **empirical reliability region** from ≥100
repetitions. The cohort-sizing run supplies 200 replicates per power estimate
over cohort size, candidate width, redundancy, analytical noise and
missingness — which is the first evidence meeting that bar. It is a
**simulation** and is labelled so; the reliability region it describes is a
property of the assumed generating process.

---

## 5. H5 — biological, analytical and temporal robustness

> *Biomarkers close to the biological function are more useful when origin,
> measurement and timing are compatible with the PK mechanism.*

### 5.1 Favourable — and the finding is that context decides, not the layer

| Question | Evidence |
|---|---|
| **does the compartment matter?** | Wegler: the enzyme predicts its substrate's clearance in **hepatocytes** (2 of 5) and in **microsomes for none**. Published independently: abundance→activity ρ **0.54** in hepatocytes vs **0.82** in microsomes, and **0.13** in intestinal mucosa. A compartment can abolish the relationship |
| **is a layer intrinsically superior?** | **no.** In Lee / Jackson the **plasma** metabolite alone predicted hepatic activity better than the **hepatic protein** alone on midazolam (2.796 vs 2.976) and worse on M37 (2.233 vs 2.116) |
| **does the marker class matter?** | two endogenous phenotype markers, opposite answers: solanidine ratios **−38 %** on metoprolol; 4β-HC/cholesterol **no gain** on venetoclax |
| **does sampling time matter?** | H5 on Müller n=20: `PARTIALLY_SUPPORTED`; the ratio's dependence on sampling instant is quantified across 9 timepoints |

### 5.2 Unfavourable — about the project's own best marker

- **Quantified vs semi-quantitative identification differs by ICC +0.362**
  [+0.084, +0.596], **in the direction opposite to expectation**.
- The headline ICC moves **0.942 → ≈0.40** on subject inclusion, and the
  operationalisation was **not frozen** before that was seen.
- **Limit-of-quantification imputation**: participant 1, the cohort's only poor
  metaboliser — *the subject the marker most needs to identify* — is below the
  limit of quantification at all 9 timepoints for both metabolites.
- CYP2D6 activity-score conventions differ between the two cohorts
  (`*41` = 0.25 post-2021 CPIC vs the Rüdesheim 2020 grid) and are **not
  aligned**; no cross-cohort mapping is applied.

### 5.3 Inconclusive

- RNA vs protein: the *size* of the difference is unestablished. Wide intervals
  at n = 15 and n = 16 are not evidence of interchangeability, and the published
  record disagrees with this project's own out-of-sample result on the
  microsomal side.
- No source pairs RNA and protein **per donor** with a function, so the
  comparison remains cross-cohort and therefore not a comparison.

---

## 6. What the five hypotheses jointly imply for the next method

1. **H2 is the binding constraint.** A new selection method must beat a
   univariate filter that reads the endpoint. Nothing in this project does.
2. **H3 says the answer is small.** Any method whose advantage requires many
   variables is arguing against the project's own best-supported hypothesis.
3. **H4 names four ways a method can look good without being good**, and every
   one of them is a control a new method must carry: a refitting permutation
   null, a width-matched comparison, residualised selection, and a refitting
   interval.
4. **H5 says context decides.** A method that conditions on context is therefore
   aiming at the right target — and H1's two-drug disagreement is the clearest
   case where context-dependence would have had to be learned from n ≤ 47.
5. **H1 says the endpoint is the scarce resource.** Three usable public
   clinical-PK sources exist in total, of sizes 47, 12 and (functional) 100.

---

## 7. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Consolidation only; no frozen protocol or result is altered. |
