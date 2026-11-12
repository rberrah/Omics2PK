# Axis B prespecification — a parsimonious mechanistic pilot: does the proteome add to the target enzyme?

**Version 1.0.1 — frozen 2026-11-10, before any model was fitted on this cohort.**
Version 1.0.1 adds only a §10 entry recording a defect found by running the
freeze, and one clearly-labelled post-hoc projection in the artefact; no
design value changed and no result was recomputed.
Objective: Axis B (functional). Status ceiling: **`exploratory`**.
Endpoint class: `functional_in_vitro`. Outcome kind: `drug_specific_clint`.

Three arms, one primary contrast, five substrates, 15 donors,
leave-one-donor-out. **No free selection among 5,172 proteins**: the extra
proteomic information enters as **two prespecified principal components**,
computed inside training folds only.

Nothing here is a clinical-pharmacokinetic result. The endpoint is an
**in-vitro intrinsic clearance** and the conclusions are specific to it.

---

## 1. The question, in three parts that must not be merged

1. **Does the target CYP's own abundance predict its substrate's intrinsic
   clearance?** — `W1` against `W0`.
2. **Does the rest of the microsomal proteome add information the target
   enzyme does not already carry?** — `W2` against `W1`. **Primary.**
3. **Do these relationships differ between microsomes and hepatocytes?** — a
   **paired comparison of two experimental systems on the same donors**, and
   explicitly *not* evidence of a causal compartment effect. The two
   preparations differ in enzyme content, in membrane environment, in the
   presence of cytosolic cofactors and in how clearance is scaled; this design
   cannot attribute a difference to any one of those.

The prior is unfavourable and is recorded here. **n = 15**, with 14 training
donors per fold and a design of at most three columns. A null result is the
expected outcome, the effect size and its uncertainty are the deliverable, and
§7 turns them into a sample-size requirement rather than into a claim.

---

## 2. The circularity check, performed before this freeze

The commissioning note requires that the CLint definition used must not
itself contain the CYP abundance that serves as the predictor. **It was
checked, and two of the three published equations fail it.**

From the article's own Figure S11 caption, verbatim:

> "a–e) Unbound in vitro clearance was scaled to kg body weight (bw) from HH
> with hepatocellularity and HLM with literature MPPGL **(eq. 7)** […]
> f–j) Unbound in vitro clearance was scaled to kg bw **with specific CYP
> amount** (midazolam → CYP3A4, omeprazole → CYP2C19, diclofenac → CYP2C9,
> bupropion → CYP2B6, and bufuralol → CYP2D6) in HH and HLM, respectively and
> **CYP concentration in HL (eq. 8)**. k–o) HLM clearance **scaled with CYP
> amount** was adjusted with Kpuu determined in HH **(eq. 9)**."

The article's **Methods** are more specific still, and were read rather than
inferred from the caption:

- **eq. 7** — *"CL<sub>int,u</sub> values from HH and HLM (eq 1, mL/min/mg
  protein) were scaled to mL/min/kg bw with literature scaling factors."*
  The scaling factors are liver-content constants from the literature. **No
  donor CYP measurement enters.**
- **eq. 8** — *"CL<sub>int,u</sub> values from HH and HLM normalized against
  the specific probed CYP (eq 2, mL/min/pmol CYP)"*, then *"scaled to
  mL/min/kg bw using the amount of the corresponding CYP in the liver"*, where
  *"Hom<sub>CYP</sub> is the specific probe CYP concentration in HL (pmol
  CYP/mg total protein)"*. So the CYP amount is in the **denominator** and the
  **homogenate** CYP concentration is a **multiplier**: eq. 8 contains **two**
  of Table S3's three CYP columns.
- **eq. 9** — the eq. 2 CYP normalisation carried forward and *"adjusted for
  the intracellular unbound concentrations (Kp<sub>uu</sub>) from HH"*.

| Equation | Scaled by | Admissible as an endpoint here? |
|---|---|---|
| **eq. 7** | hepatocellularity (HH) or literature MPPGL (HLM) — literature constants, **no donor CYP measurement** | **YES — the only admissible endpoint** |
| eq. 8 | the specific CYP amount **in the denominator** and the homogenate CYP concentration **as a multiplier** | **NO — circular against two of the three CYP columns** |
| eq. 9 | eq. 8's CYP normalisation, further adjusted by Kpuu | **NO — circular, and Kpuu is itself a Table S3 column** |

**`clint_mic_eq8`, `clint_mic_eq9` and `clint_hep_eq8` are forbidden as
endpoints in this protocol**, in every arm and every sensitivity. The runner
asserts the endpoint name is an `eq7` column and stops otherwise. The
`cyp_pmol_in_incubation_*` columns are likewise never endpoints and never
predictors: they are the CYP concentration multiplied by the incubation
protein, so they are the baseline predictor in another unit.

This is also visible in the data, which is why it would have mattered: CYP2D6
abundance spans 587-fold across the 15 microsomal donors, and dividing by it
inflates the bufuralol fold-range from **3.5** under eq. 7 to **40.9** under
eq. 8 and **116.7** under eq. 9.

---

## 3. The donor mapping, proven before this freeze

**The two halves of this source do not number their donors the same way, and
assuming they did would have mis-paired six of fifteen donors.**

Table S3 labels donors `UU1`–`UU15`. Data S1 labels them `1`–`15`. The mapping
was recovered by matching the per-donor CYP concentrations that **both** files
publish — four genes (`CYP3A4`, `CYP2C19`, `CYP2C9`, `CYP2B6`) × three
preparations (microsomes, hepatocytes, homogenate) = **twelve independent
vectors**:

| Data S1 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | **10** | **11** | **12** | **13** | **14** | **15** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Table S3 | UU1 | UU2 | UU3 | UU4 | UU5 | UU6 | UU7 | UU8 | UU9 | **UU14** | **UU15** | **UU10** | **UU11** | **UU12** | **UU13** |

Every one of the fifteen indices resolved **unanimously, 12 votes of 12**, to
a single batch, and the result is a bijection. Donors 1–9 are the identity;
**donors 10–15 are permuted.**

The runner re-derives this mapping from the data at every run and **refuses to
start** unless it is a unanimous bijection. It is never hard-coded as a
constant to be trusted.

---

## 4. Data, as measured

### 4.1 The endpoint

`clint_mic_eq7` (primary) and `clint_hep_eq7` (paired sensitivity), from Table
S3, recovered and verified against twenty of the article's own published
correlations and the quoted median
(`scripts/audit_wegler_table_s3_per_donor_kinetics.py`; see
`docs/OMICS2PK_WEGLER_PER_DONOR_CLINT_QUALIFICATION_20261110.md`).

| Substrate | CYP | microsomal n | `ND` | median | fold range | hepatocyte n | median |
|---|---|---|---|---|---|---|---|
| Midazolam | CYP3A4 | 15 | 0 | 374.7 | 4.0 | 15 | 134.0 |
| Omeprazole | CYP2C19 | 15 | 0 | 39.8 | 64.2 | 15 | 16.6 |
| Diclofenac | CYP2C9 | 15 | 0 | 494.1 | 6.1 | 15 | 392.3 |
| Bupropion | CYP2B6 | 15 | 0 | 24.5 | 8.3 | 15 | 68.9 |
| **Bufuralol** | CYP2D6 | **13** | **2** | 37.0 | 3.5 | 15 | 65.7 |

**The two `ND` values are handled by exclusion, never by imputation.** The
bufuralol microsomal analysis therefore runs on **13 donors and 13 folds**,
and that is stated beside every bufuralol number rather than left to be
inferred. Its hepatocyte analysis has all 15. No other cell is affected.

### 4.2 The baseline predictor

The target CYP's own abundance, **from Table S3** (`cyp_pmol_per_mg_HLM` for
the microsomal analysis, `cyp_pmol_per_mg_HH` for the hepatocyte one).
Table S3 is used rather than Data S1 for one reason: **CYP2D6 is absent from
Data S1 entirely**, in both preparations, so only Table S3 covers all five
substrates. Where both files carry the value they agree to within 0.05
pmol/mg, which is how §3's mapping was proven.

### 4.3 The proteome block

Data S1, microsomes sheet (hepatocytes for the sensitivity), with three
filters applied before any fit and all three independent of the endpoint:

| Filter | Microsomes | Hepatocytes |
|---|---|---|
| rows deposited | 5,172 | 4,784 |
| **complete across all 15 donors** (no imputation) | **4,039** | 3,786 |
| named gene (one blank row dropped) | 4,038 | 3,785 |
| minus the target CYP itself | **4,037** (4,038 for bufuralol) | 3,784 (3,785 for bufuralol) |

**Removing the target CYP from the proteome block is mandatory, not
cosmetic.** Leaving it in would put the baseline predictor inside the block
that the permutation control scrambles, which would both inflate `W2` and
break the null. For bufuralol the question does not arise: CYP2D6 is not in
Data S1.

Values are strictly positive and span 0.000425 to 1,726 fmol/µg — a ratio of
4.06 × 10⁶ — so they are **`log10`-transformed** before standardisation. The
hepatocyte sheet spans 1.07 × 10⁶ and is treated identically. One duplicated
gene name exists in the hepatocyte sheet (`MARCH2`, multiplicity 2) and both
copies are kept, since neither is a target.

---

## 5. Arms

Every arm predicts **`log10` of CLint**. Nothing is fitted on the raw scale.

| Arm | Predictors | Columns |
|---|---|---|
| **`W0_mean`** | none — the training-fold mean of `log10` CLint | 0 |
| **`W1_target_cyp`** | `log10` of the target CYP's abundance | 1 |
| **`W2_target_cyp_plus_proteome`** | `W1` **+ k principal components of the rest of the microsomal proteome** | 1 + k |

- **Primary contrast: `W2_target_cyp_plus_proteome` versus `W1_target_cyp`.**
- **Secondary contrast: `W1_target_cyp` versus `W0_mean`** — does the direct
  measurement of the enzyme buy anything at all?
- `W2` versus `W0` is reported descriptively and carries no separate claim.

**k = 2 is primary; k = 1 is a declared sensitivity**, reported whatever it
shows. Two is the commissioning note's upper bound and still leaves the design
at three columns against fourteen training donors; one is the most
parsimonious reading and is reported beside it so that a null at k = 2 cannot
be blamed on the extra column, nor a null at k = 1 on having had only one
direction.

### 5.1 Declared sensitivity: is it just other cytochromes?

The proteome block retains the **other** cytochromes and the redox partners
(32 such genes in the complete microsomal set, including `POR`, `CYB5A`,
`CYP3A5`, `CYP2C8`). Co-regulated cytochromes could make `W2`'s apparent
addition a re-measurement of the baseline rather than new information. A
declared sensitivity therefore removes **all 32** and re-runs `W2`. Declared
here, before any fit, and reported whatever it shows.

---

## 6. Design, estimator and controls

### 6.1 Validation

| Item | Value |
|---|---|
| Scheme | **leave-one-donor-out**, one fold per donor |
| Folds | 15 (13 for bufuralol in microsomes) |
| Folds identical across | all arms, all k, within a substrate and system |
| Estimator | **ridge** on standardised predictors, intercept fitted |
| Penalty | chosen by an **inner leave-one-out inside each training fold**, over the fixed grid 0.001, 0.01, 0.1, 1, 10, 100 |
| Seed | `20261110` |

Every constant — the proteome's log transform is fixed, but its column means
and SDs, the principal-component loadings, the predictor means and SDs, the
target mean, and the penalty — is estimated **inside the training donors
only**. The held-out donor's proteome row is projected with the training
loadings and never contributes to any of them.

`W0` is the training-fold mean of the log endpoint: a model with no predictor,
which is the reference the other two must beat.

### 6.2 Losses, in units that can be read

Fitting is on `log10`. Three losses are reported for every arm:

1. **MAE on the `log10` scale** — the primary loss.
2. **Geometric fold error, `10^MAE`** — "predictions are off by a factor of
   *x* on average", which is the interpretable form.
3. **MAE in mL/min/kg bw**, after back-transformation, so the error is also
   expressed in the endpoint's own units.

The primary contrast is on loss 1. Losses 2 and 3 are the same quantity
re-expressed and carry no separate test.

### 6.3 Controls

| Control | Definition | Replicates |
|---|---|---|
| **Proteome permutation** | the **donor labels of the proteome block only** are permuted; the target CYP abundance and the endpoint stay paired. The **entire** leave-one-donor-out is re-run, including the per-fold standardisation, the principal components and the inner penalty selection. | **200**, seed `20261110` |
| **`W0` floor** | an arm that cannot beat the training-fold mean has not been shown to carry anything | every arm |
| **Degenerate-component guard** | a training block yielding fewer than k positive eigenvalues contributes what it has, and the shortfall is reported | every fold |

The permutation is the control the commissioning note specifies, and it is the
right one for this contrast: it destroys the link between the *extra*
proteomic information and the outcome while leaving the baseline arm
untouched. **Models are refitted under every replicate** — the components are
not reused, because permuting the donor labels changes the training block, and
the penalty is reselected by the inner leave-one-out inside each replicate's
training donors.

**Which cells get the control, fixed here:** the primary contrast
(`W2` versus `W1`) at k = 2, in **both** systems — microsomes and hepatocytes
— for **all five substrates**. That is ten nulls of 200 replicates. The
k = 1 sensitivity and the cytochrome-free sensitivity of §5.1 are reported
with their intervals and are read against the k = 2 null of the same cell,
which is stated rather than implied.

**The null will not be centred on zero, and that is the point.** `W2` has more
columns than `W1`, and under a penalty chosen by inner CV a wider arm shrinks
harder and lands closer to the training mean. The permutation measures the
size of that artefact in this exact design. A small apparent improvement that
sits inside the null is reported as an artefact, not a finding. (This was
measured in Axis A, where three contrasts cleared a bootstrap interval while
lying inside a wholly positive null.)

The null is reported by its **minimum, deciles, median and maximum**, and by
the rank p-value. The 200 replicate values are not published individually: an
array of 200 numbers in a committed artefact is how a donor-level vector would
escape, and the summary carries everything the decision rule uses.

### 6.3.1 The paired system comparison, defined

Question 3 of §1 is answered by a **difference of differences over donors**,
never by comparing two separately-reported numbers. For each substrate, and
for each of the two contrasts:

> `d_i = (error_new − error_baseline)` in **microsomes** for donor *i*
>       − `(error_new − error_baseline)` in **hepatocytes** for the same donor *i*

computed on the donors present in **both** systems (13 for bufuralol, 15
otherwise), with a donor bootstrap interval. A non-zero `d` means the arm's
benefit is not the same in the two preparations. It does **not** say why: the
two systems differ in enzyme content, membrane environment, cytosolic
cofactors and scaling constant simultaneously, and this design separates none
of them. Any sentence attributing the difference to "the compartment" as a
cause is forbidden by §7.4.

### 6.4 Uncertainty

Paired bootstrap over **donors**, 2,000 draws, seed `20261110`, resampling the
per-donor out-of-fold errors without refitting. A fold is a donor here, so the
resampling unit and the fold unit coincide.

**An interval is reported with its exploratory character stated, every time.**
At n = 15 the intervals will be wide; `how_to_read` is written into the
artefact.

### 6.5 Multiplicity

Five substrates. **They are not five independent cohorts**: the same 15
donors, the same preparation, the same laboratory, the same proteome. The five
primary contrasts are therefore reported with **Bonferroni-adjusted intervals
at 95 %/5**, beside the unadjusted ones, and a claim requires the adjusted
interval. The microsomal and hepatocyte readings of the same donors are
likewise **not two cohorts**, and the system comparison of §1.3 is paired.

---

## 7. What may not be concluded, and what the result is for instead

### 7.1 Absence of difference is not equivalence

An interval containing zero means **no difference was detected**. It is never
a demonstration that there is none, and no sentence in the results document may
read that way. At n = 15 the design's power is low by construction, and the
honest output is **the estimated effect, its uncertainty, and the feasibility
of the measurement** — not a significance verdict.

### 7.2 No variance fraction from a rank correlation

No percentage of explained variance may be derived from any Spearman
coefficient, in this protocol or anywhere in this project. The square of a
rank correlation is not a variance fraction. (A sentence violating this in
`docs/OMICS2PK_PUBLIC_SOURCE_QUALIFICATION_20261109.md` §4 was withdrawn on
2026-11-10.)

### 7.3 The sample-size projection, with its assumptions on the table

If the gain is small or uncertain, the results document must convert it into a
**required sample size for a future cohort**, and must do so by stating its
assumptions and exploring **several scenarios**, not one. Fixed here:

- The projection is for the **primary contrast** (`W2` versus `W1`) on the
  `log10` MAE, using the **observed per-donor paired differences** as the
  dispersion estimate.
- **Scenarios on the true effect**: the observed point estimate; half of it;
  and the most favourable value inside the bootstrap interval. Each is run
  separately and labelled.
- **Scenarios on the design**: 80 % and 90 % power, two-sided α = 0.05, with
  and without the Bonferroni correction over five substrates.
- The projection is reported as an **order of magnitude with its assumptions
  attached**, never as a single number, and it is labelled as resting on a
  dispersion estimated from 15 donors — which is itself imprecise.

### 7.4 Not to be done after seeing the results

- No new arm, no third component, no substrate added or dropped, no change of
  endpoint equation, no change of fold scheme, loss, penalty grid or
  replicate count.
- No use of `eq. 8` or `eq. 9` as an endpoint, in any arm or sensitivity.
- No free selection among the 4,037 proteins, and no reporting of individual
  protein coefficients as a biomarker panel.
- No causal statement about the compartment (§1.3).
- No clinical-pharmacokinetic statement of any kind.
- No reading of a non-detection as equivalence.
- A defect found after this freeze is recorded in §10, not repaired silently.

---

## 8. Reporting, in three separated sections

1. **Predictive performance** — arms, the three losses, the contrasts with
   unadjusted and Bonferroni-adjusted intervals, the permutation nulls.
2. **Biological reading** — the three questions of §1, each with its caveat,
   and the cytochrome sensitivity of §5.1.
3. **Analytical validity** — the circularity check, the donor mapping, the
   `ND` handling, component shortfalls, fold safety, and the sample-size
   projection with its assumptions.

---

## 9. Freeze

| Item | Value |
|---|---|
| Cohort | `wegler_15` |
| Donors | 15 (13 for bufuralol in microsomes) |
| Validation | leave-one-donor-out |
| Seed | `20261110` |
| Ridge grid | 0.001, 0.01, 0.1, 1, 10, 100 |
| Components (k) | 2 |
| Component sensitivity | 1 |
| Permutation replicates | 200 |
| Bootstrap draws | 2,000 |
| Primary endpoint | `clint_mic_eq7` |
| Paired sensitivity endpoint | `clint_hep_eq7` |
| Forbidden endpoints | `clint_mic_eq8`, `clint_mic_eq9`, `clint_hep_eq8` |
| Primary arm | `W2_target_cyp_plus_proteome` |
| Primary comparator | `W1_target_cyp` |
| Secondary comparator | `W0_mean` |
| Substrates | 5, Bonferroni over 5 |
| Target scale | `log10` |
| Status ceiling | `exploratory` |

`scripts/run_wegler_protein_clint_pilot.py` reads this table back and stops
before loading any data if a value differs.
`tests/test_wegler_protein_clint_pilot.py` fails if this document and the
module disagree.

---

## 10. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.1 | 2026-11-10 | **One defect of this freeze, found by running it, recorded here rather than repaired.** §7.3 points the sample-size projection at the **primary** contrast (`W2` versus `W1`) — which turned out to be the contrast with no effect. Its most favourable in-interval value is therefore near zero and the donors required run into the hundreds of thousands, which is a real answer (the proteome-beyond-the-enzyme question is not worth sizing a study for) but not the useful one. The contrast that **does** show an effect is the secondary `W1` versus `W0`, in hepatocytes, and that is what a future cohort would need sized. A projection for it was added to the artefact under `post_hoc_q1_projection_not_prespecified`, with `status: POST_HOC_ADDED_AFTER_SEEING_THE_RESULTS`, and its third scenario is the **least** favourable in-interval value rather than the most, because when the effect is real the conservative question is how many donors the weakest plausible version needs. **The prespecified projection is unaltered, no arm, endpoint, fold, loss, threshold or replicate count changed, and no result was recomputed.** This is a reporting addition and is labelled as one everywhere it appears. |
| 1.0.0 | 2026-11-10 | First freeze. Written after two checks that both changed the design and before any model was fitted. (1) **The circularity check fires:** the article's own Figure S11 caption and then its Methods text show that `eq. 8` divides clearance by the specific probed CYP amount and multiplies it back by the **homogenate** CYP concentration — two of Table S3's three CYP columns — and that `eq. 9` carries that normalisation and adds Kpuu, so both are circular against a CYP-abundance predictor and are forbidden as endpoints; `eq. 7` scales by hepatocellularity or literature MPPGL and is the only admissible one. (2) **The donor mapping is not the identity:** Table S3's `UU1`–`UU15` and Data S1's `1`–`15` agree for donors 1–9 and are permuted for 10–15, proven by twelve independent CYP-concentration vectors resolving unanimously to a bijection; assuming the identity would have mis-paired six of fifteen donors. Also records in advance: the unfavourable prior at n = 15 (§1); that the target CYP must leave the proteome block or the permutation control breaks (§4.3); that the two bufuralol `ND` values are excluded and never imputed (§4.1); that the permutation null will not be centred on zero and why (§6.3); that five substrates are not five cohorts (§6.5); that no variance fraction may come from a rank correlation (§7.2); and the sample-size projection's scenario structure (§7.3). |
