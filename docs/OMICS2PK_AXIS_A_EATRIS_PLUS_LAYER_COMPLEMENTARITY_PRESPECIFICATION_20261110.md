# Axis A prespecification — reconstructing a masked biological layer in EATRIS-Plus v4

**Version 1.0.1 — frozen 2026-11-10, before any model was fitted on this cohort.**
Version 1.0.1 adds only §14 entries recording defects found by running the
freeze; no design value changed.
Objective: Axis A (multi-omic method). Status ceiling: `exploratory_methodological`.

This document fixes, in advance, the target modality, the family grouping, the
individual sets, the arms, the estimator, the folds, the controls, the summary
statistic, the success criterion and the decision rules for the EATRIS-Plus
masked-modality reconstruction. It is committed before the first fit. Anything
decided after a result appears belongs in §13 and nowhere else.

---

## 1. The question, stated so that it can fail

**Does the ensemble of molecular layers carry information about one masked
biological layer that demographics do not already carry — and is that
information *inter-family*, rather than a re-measurement of the same molecules
on a second platform?**

The second clause is the whole experiment. A plasma proteome predicted from a
second plasma-proteome assay is an assay-agreement result, not multi-omic
complementarity. This design therefore separates the two by construction
(§6, §7.3) and prespecifies a rule under which the complementarity claim
**fails even when the prediction succeeds** (§8.3).

The prior is unfavourable and is recorded as such. The predictors are a
**leukocyte** transcriptome, a **leukocyte** methylome and **plasma** small
molecules; the target is a **plasma** protein panel. Cross-compartment
prediction of individual protein levels at n = 79 is not expected to be
strong. A null result is a valid outcome of this experiment and will be
reported as one.

---

## 2. What this experiment is not

- **Not pharmacological.** EATRIS-Plus holds no drug, no dose, no
  concentration, no time, no enzyme activity and no clinical endpoint. No
  number produced under this protocol may be presented as evidence about drug
  metabolism, enzyme activity, intrinsic clearance or pharmacokinetics, in any
  document, figure or summary. The supervised target is a **molecular layer**;
  where a person-level target is used it is **demographic**, declared as such.
- **Not a validation of Axis B or C.** Blood is not liver. This cohort is in
  Axis A and nowhere else.
- **Not a clinical biomarker study.** The cohort is a healthy reference set
  from a single biobank (all identifiers carry the `CZC` prefix: one site, so
  no centre confound and no internal replication either).
- **Not confirmatory.** Its container was enumerated on 2026-11-09 before this
  freeze (layer names, individual counts, feature counts, covariate
  completeness). No model was fitted and no target-versus-predictor
  association was examined. The ceiling is nevertheless
  `exploratory_methodological`.

---

## 3. Source, as measured

| Item | Value |
|---|---|
| Record | Zenodo 19136013, DOI `10.5281/zenodo.19136013`, published 2026-03-20 |
| Title | *EATRIS-Plus multi-omics data of a human reference cohort* |
| Licence | MIT, `access_right: open`, files unrestricted |
| File | `EATRIS_PLUS_mudata.h5mu`, 141,704,166 B |
| SHA-256 | `a764a6cd7246edbd5c3a4c4ce166cd2517dc93eaa11834f2d1e8530245581c7c` |
| Layers | 12 |
| Union of individuals | 125 |
| Person-level rows | 125 |

The deposited ISA workbook claims 128 rows (127 individuals) for every assay
sheet. **The data does not agree**, and the counts below are the container's,
measured by `scripts/audit_eatris_plus_layer_configurations.py`, which fits
nothing. Anyone sizing this analysis from the ISA sheet would over-count the
transcriptome by 37 individuals.

| Layer | individuals | features |
|---|---|---|
| `acylcarnitines` | 125 | 33 |
| `aminoacids` | 125 | 46 |
| `fattyacids` | 125 | 5 |
| `lipidomics-neg` | 125 | 164 |
| `lipidomics-pos` | 125 | 196 |
| `proteins` | 125 | 1,393 |
| `lipidomics-spline` | 124 | 257 |
| `proteomics-ipp` | 123 | 6,831 |
| `EM-seq` | 117 | 100,000 |
| `mRNA` | 90 | 16,386 |
| `miRNA-mature` | 89 | 2,130 |
| `miRNA-qRT-PCR` | 86 | 179 |

---

## 4. Families — fixed before any fit

Twelve layers, **six biological families**, an exact partition. The runner
fails if the partition is not exhaustive and disjoint.

| Family | Layers |
|---|---|
| `transcriptome` | `mRNA` |
| `small_rna` | `miRNA-mature`, `miRNA-qRT-PCR` |
| `methylome` | `EM-seq` |
| `proteome` | `proteins`, `proteomics-ipp` |
| `lipidome` | `lipidomics-neg`, `lipidomics-pos`, `lipidomics-spline`, `fattyacids` |
| `metabolome_targeted` | `acylcarnitines`, `aminoacids` |

Two groupings deserve their justification, because both could be argued
otherwise and both are fixed here rather than after a result:

- `fattyacids` (VLCFA, 5 features) sits in `lipidome` because it measures
  lipid species, not because of its platform.
- `acylcarnitines` and `aminoacids` sit together because they are the same
  targeted small-molecule panel. Acylcarnitines are nevertheless fatty-acid
  conjugates, so `lipidome` is declared **biochemically adjacent** to
  `metabolome_targeted` (§5.2). Adjacency is not the same as identity, and it
  changes which arm carries the headline for the secondary target only.

---

## 5. Target modality

### 5.1 Primary target: `proteins` (1,393 features, 125/125 individuals)

Chosen for four reasons, all stated before any fit:

1. **Proteins are the functional layer.** The hepatic counterpart of this
   layer — enzyme abundance — is what Axis B needs. A method demonstrated on
   the protein layer here is the method that would be applied there. (This is
   a statement about *method transfer*, not about biology: see §2.)
2. **It costs no individuals.** `proteins` is present for all 125, so the
   complete-case size of every configuration is set by the predictors, not by
   the target.
3. **It is measurable and reportable.** 1,393 features give a distribution of
   per-feature performance rather than a single number, and 1,393 is small
   enough to report in full.
4. **It makes the redundancy warning testable.** `proteins` has a same-family
   sibling, `proteomics-ipp`. That turns "beware of intra-family redundancy"
   from a caveat into a measured quantity (§7.3).

### 5.2 Secondary target: `acylcarnitines` (33 features, 125/125)

Declared now, reported regardless of outcome. It is the layer with the
strongest *a priori* biochemical coupling to another family, so it functions
as a positive control for the machinery: if cross-family reconstruction works
anywhere in this cohort, it works here. If **both** targets return a null, the
limiting factor is n or the method, not the particular target — a distinction
that cannot be drawn from one target alone.

Its headline arm excludes **both** its own family (`aminoacids`) **and** the
biochemically adjacent `lipidome`, because predicting acylcarnitines from
lipid species is chemical overlap, not cross-layer regulation.

### 5.3 Masking

For a held-out individual, the **entire target-layer row is removed**: it is
never a predictor, never contributes to a standardisation constant, never
contributes to a component basis, and is revealed only to score the
prediction. The target layer is also absent from the predictor set of every
arm, and its same-family sibling is absent from the primary arm.

---

## 6. Configurations and individual sets

The commissioning note froze **eleven layers on 79 individuals** as primary,
with two prespecified sensitivities. The counts below were measured from the
container and agree exactly.

| Configuration | Role | Layers | Complete-case n | Families | Features | Individual-set SHA-256 (first 24) |
|---|---|---|---|---|---|---|
| `C11_primary` | primary | 11 (drop `miRNA-qRT-PCR`) | **79** | 6 | 127,441 | `265ffd52c9ce1f9f25c5e26c` |
| `C12_all_layers` | secondary — what the qRT-PCR layer adds | 12 | **55** | 6 | 127,620 | `b10acdb92dd6df92cba10585` |
| `C9_no_rna` | secondary — wider coverage without the transcriptome | 9 (drop both miRNA and `mRNA`) | **114** | 4 | 108,925 | `3ba65c29b0c988fd6d9e4a66` |

Measured coverage profile: 55 individuals hold all twelve layers, 32 hold
eleven, 25 hold ten, 13 hold nine, **none holds fewer than nine**. Dropping
`miRNA-mature` in addition to `miRNA-qRT-PCR` buys zero further individuals,
which is why there is no ten-layer configuration.

**The sets are nested — measured, not assumed:** C12 ⊂ C11 ⊂ C9. This is what
makes the sensitivity questions answerable (§9).

**No identifier is published.** Each set is recorded by its size and by the
SHA-256 of its sorted identifiers, which is sufficient to prove later that a
run used the same people without stating who they are.

### 6.1 The baseline block, and why haematology is not in it

| Block | Columns | Completeness on all three sets | Imputation |
|---|---|---|---|
| `K0` — demographics | `Age`, `Sex` (0/1), `BMI`, `Height`, `Weight` | **complete, 0 missing** | none |
| `H` — haematology | `Leukocytes`, `Neutrophils`, `Lymphocytes`, `Hematocrit`, `Erythrocytes`, `Hemoglobin`, `MCV`, `MCH`, `MCHC`, `MXD`, `Platelets` | incomplete (C11: 4 missing for `MXD` and `Neutrophils`, 1 for nine others) | training-fold median, secondary arms only |

`K0` is the prespecified baseline. It is complete, so the primary contrast
needs no imputation at all and the commissioned n = 79 is preserved exactly.

`H` is excluded from the baseline and from every omic arm, and appears only in
one declared secondary reference arm (§7.1, `A4`). Haematology is itself a
measured assay; putting it in the baseline would make the baseline stronger
but would require imputation inside the primary contrast. Both readings are
reported, and which is the primary is decided here.

`Age.group` and `BMI.group` are excluded throughout: they are coarsenings of
`Age` and `BMI` and would enter the same information twice.

---

## 7. Arms

Every arm contains `K0`. The omic contribution is therefore always
**incremental over demographics**, which is this project's standing
requirement: the primary contrast is baseline-plus-feature against baseline,
never feature against baseline.

| Arm | Predictors | Role |
|---|---|---|
| `A0_covariates` | `K0` | the baseline |
| `A1_other_families` | `K0` + every layer outside the target's family | **primary** |
| `A2_same_family_sibling` | `K0` + the target's same-family sibling only | the redundancy arm |
| `A3_all_but_target` | `K0` + every layer except the target | upper reference; includes the sibling |
| `A4_covariates_haematology` | `K0` + `H` | secondary reference, imputed |
| `B_<family>` | `K0` + one non-target family | single-family ablation, one arm per family |
| `L_<family>` | `K0` + all non-target families except that one | leave-one-family-out, one arm per family |
| `S_<subset>` | `K0` + an arbitrary subset of non-target families | parsimony enumeration (§10) |

### 7.1 Counts, fixed here

For the primary target `proteins` in `C11_primary`: 9 non-proteome layers
(119,217 features) across 5 non-target families; 5 single-family arms, 5
leave-one-out arms, and 2⁵ = **32** subsets in the parsimony enumeration (the
empty subset is `A0`). In `C9_no_rna`: 7 layers, 3 families, 2³ = 8 subsets.
In `C12_all_layers`: 10 layers, 5 families, 32 subsets.

### 7.2 The arms exist to answer named questions

- `A1` vs `A0` — **does the multi-omic ensemble add to demographics?** Primary.
- `A2` vs `A0` — how much does a *second assay of the same molecules* add?
- `A1` vs `A2` — **is the ensemble's information inter-family, or is it the
  sibling's?** This is the complementarity question.
- `A3` vs `A1` — what does adding the sibling back buy on top of the ensemble?
- `B_<family>` — which single family carries information?
- `L_<family>` — which family is *necessary*, given the others?
- `A4` vs `A0` — does a cheap clinical blood count already carry it?

### 7.3 The redundancy trap, as a measured quantity

`A2` is the arm the commissioning note's warning requires. If `A2` alone
reaches or exceeds `A1`, then what looks like multi-omic reconstruction is a
second measurement of the same plasma proteins, and §8.3 forbids the
complementarity claim regardless of how `A1` compares to `A0`.

---

## 8. Design, losses and decision rules

### 8.1 Validation

| Item | Value |
|---|---|
| Scheme | repeated K-fold over **individuals** |
| Outer folds × repeats | 5 × 2 |
| Inner folds | 3 |
| Seed | `20261110` |
| Folds identical across | all arms, all target features, all subsets, within a configuration |
| Folds across configurations | **necessarily different** (different n). No paired contrast crosses a configuration. |

A repeat is never a second individual. All intervals are bootstrapped over
individuals.

### 8.2 Predictor construction, estimator, and the dimensionality rule

p = 127,441 features against n = 79 individuals. Nothing is fitted on p.
Inside **each outer training fold**, for each retained predictor layer:

1. standardise the layer's columns with **training** means and SDs;
2. extract principal directions from the **training** rows only, via the
   n × n Gram matrix (exact, and cheap even for the 100,000-feature methylome);
3. keep the first **k = 5** components;
4. project the held-out rows with the **training** loadings;
5. standardise the assembled design — the `K0` columns together with every
   retained component score — once more with **training** means and SDs, so
   that ridge penalises every column on the same footing and no layer is
   favoured by the scale of its leading component.

The design is `K0` concatenated with the retained layers' component scores:
at most 5 + 5 × (number of layers in the arm) columns. For the primary arm
`A1` on `proteins` in `C11_primary` that is 5 + 45 = **50 columns against
n_train ≈ 63**.

Estimator: **ridge**, one model per target feature, with the penalty chosen by
the 3-fold **inner** CV inside the training fold over the fixed grid
`0.01, 0.1, 1, 10, 100, 1000, 10000`. Target features are standardised with
training means and SDs and predictions are returned to the original scale.
Every constant — layer means, layer SDs, component loadings, target means,
target SDs, penalty — is estimated inside the training fold and nowhere else.

k is prespecified at **5**. Sensitivities at **k = 2** and **k = 10** are
declared here and reported regardless of outcome (§9.3).

### 8.2.1 Value scales and missing cells inside layers

Both of the following were measured in the container before this freeze, and
both need a rule, because neither is uniform across the twelve layers.

**Scales are not uniform.** Measured on the full deposited layers:
`lipidomics-neg` spans 1,926 to 2.128 × 10⁸ and `lipidomics-pos` spans 9,439
to 8.162 × 10⁹ — raw intensities over five to six orders of magnitude —
whereas `proteins` spans 11.9 to 33.6, `EM-seq` −0.96 to 0.86 and
`proteomics-ipp` −6.4 to 10.3, all of which are already on a compressed scale.
`mRNA` (−4.2 to 2.4 × 10⁴) and `miRNA-mature` (−302 to 2.3 × 10⁵) carry
negative values and therefore cannot be logged.

The rule, fixed here:

1. **Log rule.** A layer is `log10`-transformed if, on the **full deposited
   layer** (all its rows, independent of configuration, fold and target), its
   finite values are strictly positive **and** max/min > 1000. Evaluated in
   code, this fires on exactly **`lipidomics-neg` and `lipidomics-pos`** and on
   no other layer. It is a statement about the assay's units, uses no
   individual's relation to the target, and is identical in every fold.
2. **Standardisation.** Columns are then centred and scaled with **training**
   means and SDs. A column with zero training SD is dropped from that fold.
3. **Clip.** Standardised values are clipped to **±4** training-fold SD, in
   both the training and the held-out rows, using the training SD only. This
   is a declared robustness measure against the heavy tails that survive
   column standardisation; it is applied uniformly to every layer so that no
   layer receives special treatment after the fact.

**Missing cells exist, in two layers only.** Measured on the configurations'
individuals:

| Configuration | Layer | columns | columns with any missing | cells | individuals affected |
|---|---|---|---|---|---|
| `C11_primary` (n = 79) | `acylcarnitines` | 33 | 1 | 4 | 4 |
| `C11_primary` | `aminoacids` | 46 | 9 | 55 | 35 |
| `C12_all_layers` (n = 55) | `acylcarnitines` | 33 | 1 | 2 | 2 |
| `C12_all_layers` | `aminoacids` | 46 | 8 | 43 | 26 |
| `C9_no_rna` (n = 114) | `acylcarnitines` | 33 | 1 | 5 | 5 |
| `C9_no_rna` | `aminoacids` | 46 | 9 | 82 | 51 |

Every other layer, including the primary target `proteins`, has **zero**
missing cells.

**The rule: a column with any missing cell on the configuration's individuals
is dropped — there is no imputation inside a layer.** In `C11_primary` this
costs **10 columns of 127,441**, leaving 127,431. The screen reads only the
*pattern* of missingness, never a value and never the target, so it cannot
leak; and because it is evaluated once per configuration, the design is
identical in every fold.

Dropping rather than imputing is chosen because the cost is four parts in
50,000 and because an imputed predictor would make the no-imputation claim of
the primary contrast false. When `acylcarnitines` is the **secondary target**
(§5.2), the same rule removes 1 of its 33 target features, so that target is
scored on 32.

### 8.3 Summary statistic, interval, and the success criterion

Per target feature, the out-of-fold coefficient of determination R² is
computed against the **training-fold mean** of that feature, pooled over all
held-out predictions.

**The primary contrast is `A1_other_families` versus `A0_covariates`**, on the
primary target `proteins`, in configuration `C11_primary`.

**Primary summary: the median across the target layer's features of
ΔR² = R²(`A1_other_families`) − R²(`A0_covariates`).**

Interval: paired bootstrap over **individuals**, 2,000 draws, seed `20261110`.
Each draw resamples individuals, recomputes both arms' R² per feature from the
already-stored out-of-fold predictions, and takes the median difference. No
refit.

**The primary claim is supported only if all three hold:**

1. median ΔR² > 0;
2. the 95 % bootstrap interval on the median ΔR² excludes 0;
3. the observed median ΔR² exceeds **every** one of the 20 permutation
   replicates (§8.4), i.e. rank p ≤ 1/21.

**And the complementarity claim additionally requires:**

4. median ΔR²(`A1`) > median ΔR²(`A2_same_family_sibling`), with a paired
   bootstrap interval on that difference excluding 0.

If 1–3 hold but 4 fails, the reported conclusion is *"the ensemble predicts
the masked layer, and the prediction is attributable to a second assay of the
same molecules, not to inter-family information."* That sentence is written
here, before the result, so that it cannot be softened afterwards.

### 8.4 Controls

| Control | Definition | Replicates |
|---|---|---|
| **Individual permutation** | the individual ↔ target-layer pairing is permuted **before** the entire cross-validation; predictors, covariates and folds are untouched; the full pipeline including inner CV is re-run | **20**, seed `20261110` |
| **Baseline-only floor** | `A0` is itself a control: a layer that cannot beat age, sex, BMI, height and weight has not been shown to carry anything | every arm |
| **Degenerate-component guard** | a layer whose training-fold Gram matrix has fewer than k positive eigenvalues above 1e-10 contributes only the components it has, and the shortfall is reported | every fold |

The permutation control decouples the target from everything else and must
return a median ΔR² at or below 0. If it does not, the pipeline leaks, and the
run is reported as invalid rather than interpreted.

Predictor-side PCA is **unchanged** under target permutation, so the component
bases are computed once and reused across replicates. This is an exact
identity, not an approximation, and it is what makes 20 replicates affordable.

### 8.5 The fold-safety assertions the runner must pass before fitting

1. no held-out individual's target row enters any training statistic;
2. no layer mean, SD or loading is computed on a set containing a held-out
   individual;
3. the target layer's name is absent from every arm's predictor list;
4. the target's same-family sibling is absent from `A1`;
5. the declared family grouping partitions the deposited layers exactly;
6. the individual set's SHA-256 matches the value in §6 for the configuration;
7. `K0` has zero missing values on the configuration's individuals.

A failure is a hard stop, not a warning.

---

## 9. Secondary analyses — all declared here

### 9.1 What the miRNA qRT-PCR layer adds, at fixed n

**Inside `C12_all_layers` (the same 55 individuals, the same folds):** the arm
containing all non-target families *including* `miRNA-qRT-PCR` against the
same arm with that single layer removed. Because the individuals and folds are
identical, the difference is attributable to the layer.

### 9.2 Coverage against modalities — and the comparison that is forbidden

**Inside `C11_primary` (the same 79 individuals, the same folds):** the
11-layer arm against the 9-layer arm. This isolates the contribution of `mRNA`
and `miRNA-mature` at fixed n and fixed people.

**`C9_no_rna` (n = 114)** is reported for its own absolute performance.
The comparison between the 9-layer arm on 114 individuals and the 9-layer arm
on 79 individuals is **descriptive and unpaired**: different people, different
folds. No bootstrap interval is computed for it, and no sentence in any
results document may attribute such a difference to a modality. The
commissioning note's rule is adopted verbatim: *performance differences
between distinct sample sizes are not modality effects.*

### 9.3 Component-count sensitivity

k = 2 and k = 10, on the primary target, in `C11_primary`, for `A0`, `A1` and
`A2` only. Declared here and reported whatever they show.

### 9.4 Secondary target

`acylcarnitines`, full arm set, in `C11_primary`, with the adjacency rule of
§5.2.

---

## 10. Parsimony — produced by the algorithm, not by hand

All 2^(F−1) subsets of the non-target families are enumerated (32 in
`C11_primary`, 8 in `C9_no_rna`). For each, the median ΔR² against `A0` is
computed on the same folds.

**The parsimonious set is defined by this rule, executed in code:**

> Among subsets whose median ΔR² is at least (the best subset's median ΔR²
> **minus 0.01**), take the one with the fewest **families**; ties are broken
> by the fewest **features**, then by lexicographic order of family names.

The tolerance 0.01 and the tie-break order are fixed here.

**The selected subset's own performance is optimistically biased**, because it
was chosen on out-of-fold performance over 32 candidates. It is reported as a
**descriptive** selection, its performance number is labelled as biased, and
no unbiased estimate is claimed for it. The full enumeration table is
published so that the selection can be checked against every alternative.

`tests/` must contain a test that fails if the selected subset is written into
the code rather than returned by this rule.

---

## 11. Three separate reports, never merged

The results document reports, in three separate sections that are never
combined into a single verdict:

1. **Predictive performance** — median ΔR², intervals, the arm table, the
   permutation null.
2. **Biological relevance** — which families and which target features, with
   the compartment caveat attached to every statement, and with the explicit
   statement that no pharmacological interpretation is available.
3. **Analytical validity** — fold safety, degenerate components, the
   individual-set checksums, the covariate completeness, the optimism of the
   parsimony selection, and the single-site limitation.

---

## 12. What may not be done after seeing the results

- No new target modality, family grouping, arm, k, fold scheme, penalty grid,
  threshold, tolerance or summary statistic.
- No re-running of a configuration with different individuals.
- No paired contrast across configurations.
- No reinterpretation of a `C9`-versus-`C11` difference as a modality effect.
- No pharmacological, enzymatic or pharmacokinetic statement, in any register,
  including a hedged one.
- No presentation of the parsimony selection's performance as unbiased.
- A defect found after the freeze is **recorded in §13**, with its effect on
  every number, and not repaired silently.

---

## 13. Freeze

| Item | Value |
|---|---|
| Fold seed | `20261110` |
| Outer folds × repeats × inner | 5 × 2 × 3 |
| Bootstrap draws | 2,000 |
| Ridge grid | 0.01, 0.1, 1, 10, 100, 1000, 10000 |
| Components per layer (k) | 5 |
| Component sensitivities | 2, 10 |
| Permutation replicates | 20 |
| Log10 rule | strictly positive **and** max/min > 1000 on the full deposited layer |
| Standardised-value clip | ±4 training-fold SD |
| Within-layer missing columns | dropped; no imputation inside a layer |
| Parsimony tolerance (median ΔR²) | 0.01 |
| Primary target | `proteins` |
| Secondary target | `acylcarnitines` |
| Primary arm | `A1_other_families` |
| Primary comparator | `A0_covariates` |
| Complementarity comparator | `A2_same_family_sibling` |
| Primary configuration | `C11_primary` (11 layers, n = 79) |
| Families | 6, as partitioned in §4 |
| Baseline block | `K0` = Age, Sex, BMI, Height, Weight |
| Status ceiling | `exploratory_methodological` |

`scripts/run_eatris_plus_masked_modality_reconstruction.py` re-reads this table
from `omics2pk.multiomic_layer_complementarity.FROZEN` and stops before loading
any data if a value differs.
`tests/test_multiomic_layer_complementarity.py` fails if a value in this
document and the module disagree.

---

## 14. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.1 | 2026-11-10 | **Three defects found after the freeze and after the run, recorded here as §12 requires rather than repaired in place.** None changed a target, family, arm, configuration, individual set, fold, threshold or summary statistic. (1) **The parsimony tolerance 0.01 is wider than the entire effect range (≈0.013 in median ΔR² across the 32 subsets),** so 30 of 32 subsets cleared the floor and §10's rule reduced to its tie-breaks, selecting a one-family subset whose own ΔR² (+0.0043) lies inside the primary contrast's permutation null. The constant is **not** changed; the selection is published with that fact attached and its number labelled optimistically biased. (2) **The nine-layer-subset contrast of §9.2 is degenerate in `C9_no_rna`:** that configuration holds no RNA layer, so the subset is identical to `A1` and the first execution reported an arm against itself as a ΔR² of exactly 0 with interval [0, 0]. The guard compared arm names rather than layer content. Fixed in the runner; the comparison is now absent in that cell instead of reported as a zero; the primary contrast reads +0.00682 before and after. (3) **§8.4 does not state which arms the permutation control covers.** The clause defines the control on the individual-to-target pairing, not on an arm pair, and the first execution nevertheless applied it to the primary contrast alone. It was extended to **every** reported contrast, uniformly — an extension that can only weaken conclusions, and the one that revealed that three single-family contrasts with bootstrap intervals excluding zero sit entirely inside a wholly positive null. Recorded as a clarification of an under-specified clause rather than a change of design. Results in `docs/OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_RESULTS_20261110.md`. |
| 1.0.0 | 2026-11-10 | **One correction to this document's own first draft, made before the freeze and before any model was fitted, found by measuring the layers rather than reading the deposit description.** The first draft specified standardisation and said nothing about value scales or missing cells inside layers, both of which turn out to be non-uniform: `lipidomics-neg` and `lipidomics-pos` are raw intensities spanning five to six orders of magnitude while every other layer is already compressed, and `aminoacids` and `acylcarnitines` carry missing cells (55 cells over 9 columns and 4 cells over 1 column respectively, on the primary configuration's 79 individuals) while the other ten layers carry none. §8.2.1 now fixes a log rule evaluated on the full deposited layer (it fires on exactly those two lipidomics layers), a ±4 SD clip applied uniformly, and the no-imputation column drop that costs 10 columns of 127,441. No target, family, arm, configuration, individual set, fold, threshold, tolerance or summary statistic changed. |
| 1.0.0 | 2026-11-10 | First freeze. Written after enumerating the container (`scripts/audit_eatris_plus_layer_configurations.py`, which fits nothing) and before fitting any model. Records in advance: the unfavourable cross-compartment prior (§1); the family partition with its two arguable groupings justified (§4); `proteins` as primary target with its same-family sibling named and the redundancy arm `A2` built to measure it (§5.1, §7.3); the rule under which the complementarity claim fails even when the prediction succeeds (§8.3 criterion 4); haematology excluded from the baseline so that the primary contrast needs no imputation and n = 79 is preserved exactly (§6.1); the forbidden `C9`-versus-`C11` paired comparison (§9.2); and the parsimony rule with its declared optimism (§10). |
