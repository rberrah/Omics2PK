# Axis A results — the masked layer is not reconstructed, and the only robust signal is intra-family

**Version 1.0.0 — 2026-11-10.** Protocol:
`docs/OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_PRESPECIFICATION_20261110.md`
(v1.0.0, frozen before any model was fitted).
Artefact: `reports/data/axis_a_eatris_masked_modality_reconstruction.json`.
Status: **`exploratory_methodological`**. Evidence type: **new exploratory
analysis on independent public data**.

**No pharmacological reading is available anywhere in this document.**
EATRIS-Plus v4 holds no drug, no dose, no concentration, no time, no enzyme
activity and no clinical endpoint. The predictors are a leukocyte
transcriptome and methylome and plasma small molecules; the primary target is a
plasma protein panel. Blood is not liver. Nothing here is evidence about drug
metabolism, enzyme activity, intrinsic clearance or pharmacokinetics.

Reported in the three separate sections the protocol requires (§11), which are
never merged into a single verdict.

---

## 0. The result in one paragraph

Across **four cells and sixty prespecified contrasts**, the multi-omic ensemble
did **not** reconstruct the masked layer beyond demographics: every verdict is
`NOT_SUPPORTED`, and **exactly one** contrast in the whole experiment passes
both the patient-level bootstrap interval and the permutation null. That one
contrast is the one that adds the target's **same-family sibling** on top of
everything else. The design was built to tell inter-omic complementarity apart
from a second assay of the same molecules, and the only robust signal it found
is the second kind.

A second finding is methodological and matters more than the first. **The
permutation null is positive in every cell** — median +0.0021 to +0.0040 in
median ΔR² — because a wider arm shrinks harder under an inner-CV penalty and
lands closer to the training mean than a narrow arm does, with no information
at all. Three single-family contrasts have bootstrap intervals that exclude
zero while sitting **entirely inside a wholly positive null**. Without the
uniform permutation grid, those would have read as findings.

---

## 1. Predictive performance

### 1.1 The primary cell: `proteins`, eleven layers, 79 individuals

Median out-of-fold R² across the 1,393 target features, against each feature's
training-fold mean:

| Arm | predictors | design columns | median R² | features with R² > 0 |
|---|---|---|---|---|
| `A0_covariates` | age, sex, BMI, height, weight | 5 | **−0.0195** | 11.5 % |
| `A1_other_families` | + 9 layers, 5 non-protein families (119,217 features) | 50 | **−0.0110** | 19.9 % |
| `A2_same_family_sibling` | + `proteomics-ipp` only (6,831 features) | 10 | −0.0137 | 12.5 % |
| `A3_all_but_target` | + all 10 non-target layers | 55 | −0.0114 | 18.8 % |
| `A4_covariates_haematology` | + 11 haematology indices | 16 | −0.0225 | 13.0 % |

**Every arm's median R² is negative.** At n = 79, no arm — including
demographics alone — predicts the median plasma protein better than that
protein's own training-fold mean.

The primary contrast, `A1_other_families` versus `A0_covariates`:

| | value |
|---|---|
| median ΔR² | **+0.0068** |
| 95 % paired bootstrap over individuals | **[−0.0024, +0.0149]** |
| permutation null (20 replicates) | median **+0.0021**, range **[−0.0040, +0.0120]** |
| rank p | **0.190** |
| criterion 1, median ΔR² > 0 | met |
| criterion 2, interval excludes 0 | **not met** |
| criterion 3, exceeds every replicate | **not met** |
| **verdict** | **`NOT_SUPPORTED`** |

The observed +0.0068 lies **inside** the permutation null. Two of the three
prespecified criteria fail, and the failure of the third is the informative
one: the small positive point estimate is what this pipeline returns when the
target is decoupled from everything else.

The complementarity comparison never became answerable, because the primary
claim failed. For the record, `A1` against `A2_same_family_sibling` is
**+0.0001**, interval [−0.0063, +0.0084]: **nine layers across five families
and 119,217 features are indistinguishable from the single sibling assay's ten
components.**

### 1.2 Single-family arms, and why two "significant" results are not findings

| Contrast | median ΔR² | 95 % interval | excludes 0 | permutation null | rank p |
|---|---|---|---|---|---|
| `B_methylome` vs `A0` | +0.0091 | [+0.0044, +0.0118] | **yes** | **[+0.0034, +0.0109]** | 0.143 |
| `B_small_rna` vs `A0` | +0.0064 | [+0.0019, +0.0103] | **yes** | [−0.0001, +0.0066] | 0.143 |
| `B_metabolome_targeted` vs `A0` | +0.0040 | [−0.0030, +0.0109] | no | [−0.0015, +0.0074] | 0.286 |
| `B_lipidome` vs `A0` | +0.0036 | [−0.0040, +0.0122] | no | [−0.0057, +0.0110] | 0.190 |
| `B_transcriptome` vs `A0` | +0.0030 | [−0.0050, +0.0100] | no | [−0.0005, +0.0078] | 0.238 |

The methylome row is the one to read carefully. Its bootstrap interval excludes
zero and would, on its own, look like a positive finding. **Its permutation
null is entirely positive and brackets the observed value**: twenty replicates
in which the methylome cannot possibly carry information about any individual's
protein levels produced median ΔR² between +0.0034 and +0.0109, and the
observed +0.0091 sits inside that. The same pattern recurs in `C12`
(+0.0087, null [+0.0044, +0.0170], p = 0.667) and in `C9` (+0.0041, null
[+0.0015, +0.0102], p = 0.571).

Leaving any single family out of `A1` changes nothing: the five
leave-one-family-out deltas run from −0.0001 to +0.0014 and not one interval
excludes zero. The ensemble is not more than its parts; four of the five
single-family arms are at least as good as all nine layers together.

### 1.3 The only contrast that survives both controls

Secondary target `acylcarnitines`, 32 of 33 features scored (one dropped for
missingness), same 79 individuals, same folds:

| Arm | median R² | features with R² > 0 |
|---|---|---|
| `A0_covariates` | −0.0092 | 31.2 % |
| `A2_same_family_sibling` (`aminoacids`) | −0.0050 | 46.9 % |
| `A1_other_families` (5 families, no `aminoacids`) | +0.0177 | 62.5 % |
| **`A3_all_but_target`** (`A1` + `aminoacids`) | **+0.0391** | **81.2 %** |

| Contrast | median ΔR² | 95 % interval | excludes 0 | null range | rank p |
|---|---|---|---|---|---|
| **`A3_all_but_target` vs `A1`** | **+0.0194** | **[+0.0071, +0.0363]** | **yes** | **[−0.0024, +0.0029]** | **0.048** |
| `primary_A1_vs_A0` | +0.0383 | [−0.0058, +0.0672] | no | [−0.0154, +0.0127] | 0.048 |
| `sibling_A2_vs_A0` | +0.0323 | [−0.0041, +0.0642] | no | [−0.0027, +0.0151] | 0.048 |
| `complementarity_A1_vs_A2` | +0.0059 | [−0.0471, +0.0556] | no | [−0.0187, +0.0137] | 0.143 |

**`A3` against `A1` is the single contrast in the entire experiment — sixty
contrasts, four cells — whose bootstrap interval excludes zero *and* whose
observed value exceeds all twenty permutation replicates.** What it measures is
the contribution of **`aminoacids`**, the target's own same-family layer, added
on top of five other families.

That is intra-family, by the protocol's own definition, and §8.3 criterion 4
exists precisely to stop it being called complementarity. The amino-acid and
acylcarnitine panels are the same targeted small-molecule assay, and
acylcarnitines are fatty-acid conjugates whose levels are coupled to the
mitochondrial amino-acid pool. Whether the signal is biochemical coupling or
shared assay structure, **the public record cannot separate the two**, and
neither can this analysis.

Note also that `primary_A1_vs_A0` for this target exceeds every permutation
replicate (+0.0383 against a null topping out at +0.0127) while its bootstrap
interval still spans zero. The reading is exact: **the signal is not a pipeline
artefact, and 79 individuals are not enough to pin it down.** Those two
statements are not in tension; they are the two different questions the two
controls ask.

### 1.4 The two commissioned sensitivities

**What the miRNA qRT-PCR layer adds, at fixed n.** Inside `C12_all_layers` —
the same 55 individuals and the same folds — the arm with all non-target
families including `miRNA-qRT-PCR` against the same arm with that one layer
removed:

> median ΔR² **−0.00056**, interval [−0.0026, +0.0007], rank p 0.619.

**The layer adds nothing, and the point estimate is slightly negative.** The
24 individuals that the twelfth layer costs buy no information in exchange.

**Coverage against modalities, at fixed n.** Inside `C11_primary` — the same 79
individuals and the same folds — the eleven-layer arm against the nine-layer
arm (dropping `mRNA` and `miRNA-mature`):

> median ΔR² **+0.00005**, interval [−0.0023, +0.0023], rank p 0.238.

**The transcriptome and the mature-miRNA layer together add nothing** to the
reconstruction of the plasma protein panel, at fixed n and fixed people. This
is the clean comparison; the one below is not.

**`C9_no_rna`, n = 114, nine layers.** Primary ΔR² **+0.0050**, interval
[−0.0008, +0.0101], null median +0.0025, range [−0.0001, +0.0089], p = 0.095.
Verdict `NOT_SUPPORTED`.

> **The comparison between the nine-layer arm on 114 individuals and the
> nine-layer arm on 79 is descriptive and unpaired: different people,
> different folds. No interval is computed for it, and no difference between
> them is attributed to a modality. Performance differences between distinct
> sample sizes are not modality effects.**

That rule was written into the protocol before the run and is restated here
because the two numbers invite exactly the comparison it forbids.

### 1.5 Component-count sensitivity

Declared in advance at k = 2 and k = 10, reported whatever they showed:

| k | A0 | A1 | A2 | primary ΔR² | interval |
|---|---|---|---|---|---|
| 2 | −0.0195 | −0.0140 | −0.0146 | +0.0044 | [−0.0043, +0.0138] |
| **5 (frozen)** | −0.0195 | −0.0110 | −0.0137 | **+0.0068** | [−0.0024, +0.0149] |
| 10 | −0.0195 | −0.0101 | −0.0098 | +0.0067 | [−0.0012, +0.0145] |

The conclusion does not depend on k. No interval excludes zero at any k.

### 1.6 Parsimony — the rule ran, and it exposed its own tolerance

The enumeration covers all 2⁵ = 32 subsets of the five non-target families in
`C11_primary`. The best subset is `S_methylome+small_rna` at median ΔR²
**+0.01289**. The prespecified tolerance is 0.01, so the floor is +0.00289 and
**30 of the 32 subsets are eligible**. The rule then takes the fewest families,
breaking ties by fewest features, and selects:

> **`S_metabolome_targeted`** — one family, 69 features, median ΔR² +0.00431.

Three things must be said about that, and all three are consequences of the
frozen rule rather than reasons to change it:

1. **The selected subset's performance is optimistically biased** and is not an
   unbiased estimate of anything. It was chosen on out-of-fold performance
   across 32 candidates. The protocol said so in advance (§10).
2. **The tolerance was set wider than the entire effect range.** The whole
   enumeration spans about 0.013 in median ΔR², and the tolerance is 0.01, so
   almost every subset qualifies and the rule degenerates to "fewest families,
   fewest features". That is a defect in the prespecified constant, discovered
   by running it. It is recorded here, in §4 and in the protocol's change log —
   **not repaired**, because retuning a tolerance after seeing the landscape is
   exactly what §12 forbids.
3. **The selection is selecting among noise.** +0.00431 lies inside the
   primary contrast's permutation null, [−0.0040, +0.0120]. A parsimonious
   subset chosen inside a null band is a statement about the rule, not about
   biology.

`C12_all_layers` selects `S_lipidome` (24 of 32 eligible, best +0.01193);
`C9_no_rna` selects **`S_none`** — the empty subset, i.e. demographics alone —
with all 8 of 8 subsets eligible and the best at +0.00601. In the configuration
with the most individuals, the frozen rule declines to keep any omic family at
all.

---

## 2. Biological relevance

Every sentence in this section carries the compartment caveat, and none
carries a pharmacological reading.

- **There is no demonstrated inter-family complementarity in this cohort.**
  Nine molecular layers across five biological families, 119,217 features, add
  nothing detectable to five demographic variables in predicting individual
  plasma protein levels — and nothing detectable beyond a single second assay
  of those same proteins.
- **The one robust signal is within a family.** Amino acids improve the
  reconstruction of acylcarnitines on top of five other families, robustly to
  both controls. Acylcarnitines are acyl conjugates of carnitine and their
  levels are coupled to the mitochondrial amino-acid pool, so a biochemical
  reading is available; so is a shared-assay reading, because both panels are
  the same targeted small-molecule measurement. **The public record cannot
  separate them**, and this analysis does not claim to.
- **The cross-compartment prior was recorded in advance and was correct.** The
  protocol stated before the run that predicting individual plasma protein
  levels from a leukocyte transcriptome and methylome at n = 79 was not
  expected to be strong. It was not. That this was anticipated does not make it
  a weaker result; it makes it a result that was not shaped by its outcome.
- **The transcriptome specifically buys nothing here.** Its layer contributes
  median ΔR² +0.0030 alone (interval spanning zero) and +0.00005 when added to
  seven other layers. In a cohort where the transcriptome costs 35 individuals
  of complete-case coverage, that is a measured price for a measured return.
- **Haematology is not a cheap substitute.** `A4` is *worse* than demographics
  alone for the protein target (−0.0225 against −0.0195).

**What may not be concluded.** Nothing above is a statement about hepatic
enzyme abundance, drug-metabolising activity, intrinsic clearance or
pharmacokinetics. A negative result on leukocyte-to-plasma reconstruction in
healthy volunteers says nothing about whether a liver transcriptome predicts a
liver proteome in the same donors, which is a different question on different
tissue and remains unanswerable from any dataset this project can reach.

---

## 3. Analytical validity

- **Fold safety.** Layer means, SDs, component loadings, design means, design
  SDs, target means, target SDs and the ridge penalty are all estimated inside
  the training rows. A test replaces the entire held-out block and asserts the
  training scores are bit-for-bit identical; another reproduces the held-out
  projection by an independent route; another asserts the first j components of
  a k-component extraction equal the j-component extraction.
- **Target masking.** The target layer is absent from every arm's predictor
  list in the artefact, verified by test across all four cells. Its same-family
  sibling is absent from `A1`, also verified.
- **Individual sets.** Each configuration's complete case equals the
  commissioned count exactly — **79**, **55**, **114** — and each set's
  SHA-256 matches the value frozen in the protocol before the run
  (`265ffd52…`, `b10acdb9…`, `3ba65c29…`). The sets are nested,
  C12 ⊂ C11 ⊂ C9, which is what makes the two sensitivities answerable at
  fixed n. No identifier is published.
- **Baseline completeness.** `K0` has zero missing values on all three sets, so
  the primary contrast involves **no imputation at all**. Haematology, which
  does have missing values (17 cells in `C11`), is confined to the single
  declared `A4` arm and imputed with training-fold medians.
- **Within-layer missingness.** Columns with any missing cell were dropped, not
  imputed: 10 columns of 127,441 in `C11` (1 acylcarnitine, 9 amino acids).
  The screen reads the pattern of missingness only.
- **Scale rule.** The `log10` rule fired on exactly `lipidomics-neg` and
  `lipidomics-pos`, as the protocol predicted from the deposited ranges, and on
  no other layer. Standardised values were clipped at ±4 training-fold SD
  uniformly.
- **Component shortfalls.** `fattyacids` has 5 features and therefore supplied
  5 components rather than the 10 requested for the cache, in all 10 folds of
  all three configurations. This is the guard reporting a hard limit of the
  layer, not a defect; the frozen k is 5, so no arm was short.
- **A degenerate contrast was found and removed before publication.** In
  `C9_no_rna` the "nine-layer subset" is identical to `A1` — the configuration
  holds no RNA layer — so the comparison was an arm against itself and returned
  a tautological ΔR² of exactly 0 with interval [0, 0]. The guard compared arm
  *names* rather than layer *content*. Fixed, and the comparison is now simply
  absent in that cell rather than reported as a zero. All other numbers are
  unchanged: the primary contrast reads +0.00682 before and after.
- **The permutation grid is uniform, not targeted.** The null is computed for
  **every** reported contrast, in every cell, not only for the primary one. A
  null added afterwards to whichever contrast looked significant would have
  been selection; added to all sixty, it can only temper claims. This decision
  is the reason §1.2 exists.
- **Why the null is not centred on zero.** A wider arm shrinks harder under an
  inner-CV penalty, so it lands closer to the training mean than a narrow arm
  and scores a slightly higher R² with no information at all. The null measures
  that artefact in this exact design. Any study reporting an incremental R² of
  0.002–0.009 between arms of different width, without such a control, is at
  risk of reporting this artefact.
- **Single site, healthy reference cohort.** All individuals carry the `CZC`
  prefix: one biobank, so no centre confound and no internal replication.
- **Runtime** 123.4 s, 1,393 + 32 target features, 60 contrasts, 20
  permutation replicates per contrast. Wall clock on one machine; a suspended
  process inflates it without computing anything, so it is never read as a cost
  model.

---

## 4. Defects in the prespecification, found by running it

Recorded here as §12 requires, not repaired in place.

1. **The parsimony tolerance (0.01) is wider than the whole effect range
   (≈0.013).** Thirty of 32 subsets qualified, so the rule reduced to its
   tie-breaks. The selected subset is reported with that fact attached and its
   number labelled biased. The constant is **not** changed.
2. **The nine-layer subset contrast is degenerate in `C9_no_rna`** and was
   reported as a tautological zero in the first execution. Fixed in the runner
   before publication, with no effect on any other number (§3).
3. **§8.4 does not say which arms the permutation control covers.** The first
   execution applied it to the primary contrast only. Reading the clause as
   written — the control is defined on the target pairing, not on an arm pair —
   it was extended to every reported contrast uniformly. That extension can
   only weaken conclusions, and it is what revealed §1.2. Recorded as a
   clarification of an under-specified clause, not as a change of design.

---

## 5. What this contributes to the funding case

**Evidence type: new exploratory analysis on independent public data.** The
pipeline was executed end to end on a real 12-layer individual-level
multi-omic cohort, 127,620 features, 125 individuals, under a protocol frozen
before any fit. That demonstration is the deliverable; the negative result is
the finding.

**Supports:**

1. **The method handles real multi-omic data at individual level**, including a
   100,000-feature methylome, with every transform fold-safe and every
   configuration reproducing its commissioned size and checksum.
2. **Complementarity and redundancy are separable, and the separation
   matters.** Built to tell inter-family information from a second assay of the
   same molecules, the design found exactly one robust signal and correctly
   classified it as intra-family.
3. **A quantified warning about incremental-R² reporting.** In this regime a
   permutation null sits at +0.002 to +0.004 median ΔR² with an upper range to
   +0.012, purely from differential shrinkage. Three contrasts cleared a
   bootstrap interval while lying inside that band.
4. **Priced layer choices.** The twelfth layer (miRNA qRT-PCR) costs 24
   individuals and returns −0.0006. The transcriptome and mature-miRNA layers
   cost 35 individuals and return +0.00005.

**Does not support, and must not be cited for:**

- Any statement about drug metabolism, enzyme activity, intrinsic clearance or
  pharmacokinetics.
- Any statement about hepatic tissue.
- Any claim that multi-omic integration does not work in general: this is one
  cohort, one compartment, 79 complete cases, one target layer, and a
  cross-compartment question with an unfavourable prior recorded in advance.

**What it argues for measuring.** The acylcarnitine cell is the precise shape
of the data requirement: a signal that **survives the permutation null** while
its **patient-level interval still spans zero** is a signal that exists and is
under-powered. Either more individuals at the same layers, or — for the
question this project actually needs — the same design in the **right
compartment**, with a drug and a measured outcome. Neither is available in any
public dataset this project has qualified.

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First write-up. Reports all four cells and sixty contrasts, the single contrast passing both controls and its intra-family nature, the positive permutation null and the three bootstrap-significant contrasts it explains, both commissioned sensitivities, the component-count sensitivity, the parsimony enumeration with its three stated caveats, and the three defects of the prespecification found by running it. |
