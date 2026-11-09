# F-series prespecification — do signed regulatory activities add anything to the target transcript?

**Status: frozen 2026-11-09, before any F-series model was fitted.**
**Version 1.0.0.**

Every number in this document was measured from the sources before the
specification was written, and none of them is a model result. The distinction
matters because this document's only purpose is to make the experiment
unfalsifiable-by-hindsight: the arms, the folds, the losses, the thresholds,
the controls and the abstention rules are all fixed here, and the runner
refuses to start if they have moved.

The identifiers `F0`–`F5` belong to this functional benchmark alone. They do
**not** replace `M0`–`M5`, `P1`–`P3` or `HY`, which keep their meaning in
`OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md`. Nothing here modifies the
programme's frozen primary contrast (ARIANA, Ridge, M3 versus M1).

---

## 1. The question, stated so that it can fail

A transcript tells you how much message for one enzyme is present. A signed
regulatory network additionally claims to know *why*: which transcription
factors drive or repress that enzyme, and in which direction. If that claim
carries information, then a score summarising the state of a transcription
factor's regulon should predict enzyme activity **beyond** the enzyme's own
transcript.

**Primary question.** Does adding CollecTRI transcription-factor activities to
the target enzyme's own transcript reduce the out-of-fold absolute error of
microsomal CYP3A4 activity?

**The honest prior is unfavourable.** Three facts say so in advance:

1. In H2a on this same cohort, no selector beat the univariate filter in 0 of
   80 comparisons, and the target-transcript-alone arm was already hard to
   improve on. A functional layer is a *selector of a different kind*; it
   inherits that difficulty.
2. A transcription-factor activity computed by a univariate linear model over a
   regulon of 200-odd genes is, mechanically, a weighted average of a large
   fraction of the transcriptome. Fifty such scores span a space close to the
   leading principal components of the array. Ridge already had access to that
   space in H2a's no-selection arm, where it did not help.
3. Transcription-factor *protein* activity is post-translationally regulated.
   An mRNA-footprint proxy for it is a proxy for a proxy, measured in
   post-mortem liver.

A null result is therefore the expected outcome and is reported as a result,
not as a failure of execution. This document fixes what would have to be true
for the opposite conclusion, so that the conclusion cannot be chosen later.

---

## 2. What this experiment is not

- **It is not a pharmacokinetic validation.** The Human Liver Cohort has no
  clinical pharmacokinetics: no dose, no concentration, no time. Every
  F-series endpoint is `functional_in_vitro`, the same class as H2a, and is
  never pooled with, or reported in the same column as, a clinical result.
- **It will not be extended to clinical exposure by simulation.** Running a
  physiologically based model on an activity predicted here and presenting the
  resulting concentration as evidence of transfer is forbidden. The runner
  writes `clinical_pk_transfer: "NOT_EVALUATED"` with its reason into every
  report, and no F-series artefact contains a simulated concentration.
- **It is not confirmatory.** The Human Liver Cohort has been inspected
  repeatedly by this project (H2a, H2d, the pool construction, the
  destroyed-signal diagnostic). Every F-series result is capped at
  `exploratory` by prior exposure, whatever its direction, exactly as ARIANA
  is. The cap is not a hedge about the arithmetic; it is a statement about how
  many times these donors have been looked at.
- **It builds no other layer.** SIGNOR, PROGENy, COSMOS, causal-flow
  reconstruction and the physiologically based scaling step are out of scope
  and no file is created for them. `F3` is reserved for a pathway-footprint arm
  and is deliberately not implemented; see §11.

---

## 3. Sources, as measured on 2026-11-09

| Source | What was measured | Value |
|---|---|---|
| CollecTRI, served by OmniPath | interactions retrieved | 64,516 |
| | rows whose `sources` names CollecTRI | 64,516 of 64,516 |
| | `is_directed` false | 0 |
| | distinct complex regulators (`A_B` form) | 25, covering 16,970 edges |
| HLC `curatedExpressionLiver.zip` | probes × donors | 40,638 × 466 |
| HLC `curatedPhenotypesynapse.zip` | donors with a CYP3A4/midazolam value | 487 |
| | **donors with both expression and that value** | **397** |
| `QUALIFIED_PROBE_POOLS.json`, pool `A_qualified_human` | admissible probes on this array | 26,066 |
| | distinct symbols they carry | 21,634 |
| | **symbols left after the cytochrome guards of §4.3** | **21,580** |

The donor count that matters is **397**, not 487: the phenotype table covers
donors the expression table does not. Any report quoting 466 or 487 as the
analysis n would be wrong.

---

## 4. Building the signed regulons — every rule fixed here

### 4.1 Sign resolution

CollecTRI carries two flag pairs, `is_stimulation` / `is_inhibition` and the
two consensus columns. **Frozen rule: the consensus columns decide.** The
weight of an edge is `+1` when the consensus is stimulation alone, `-1` when it
is inhibition alone, and the edge is **dropped** when the consensus says both
or neither. A dropped edge is not given weight zero and silently kept in the
count; it leaves the regulon, because a zero-weight member would still count
towards the evaluability floor of §4.5.

Measured on the 64,516 rows:

| | edges |
|---|---|
| consensus determinate | 61,926 |
| — stimulation | 52,782 |
| — inhibition | 9,144 |
| consensus indeterminate, **dropped** | 2,590 (4.01 %) |
| — both columns true | 1,020 |
| — both columns false | 1,570 |

The consensus columns are **not** a restatement of the flags, and this is the
reason the rule names them explicitly: 4,776 edges carry both `is_*` flags, of
which 1,280 stay indeterminate under the consensus; and 1,310 edges carrying
exactly one `is_*` flag are nonetheless indeterminate under the consensus. An
implementation reading the flags instead of the consensus would retain 1,310
edges this one drops and would resolve 3,496 it keeps differently.

### 4.2 Duplicate edges

Measured: after sign resolution, 2,004 rows duplicate a `(source, target)` pair
already seen — `RELA → CYP3A4` among them, recorded twice with identical signs.
120 pairs resolve to *opposite* signs across their rows.

**Frozen rule.** Duplicates collapse to one edge. A pair whose rows resolve to
different signs is dropped: 120 such pairs leave.

**58,991 signed edges survive all of §4.1 to §4.3**, which is
`64,516 − 2,590 − 811 − 2,004 − 120`.

### 4.3 Removing every cytochrome, in three places

This is the defect that made the first draft of this experiment circular, and
it is removed in three independent places rather than one, because one place
can be bypassed by an alias.

Measured: each study enzyme sits inside the regulons of its own candidate
regulators — CYP3A4 is a member of 51 of them, CYP2D6 of 6. Scoring a
transcription factor from a regulon containing CYP3A4 and then using that score
to predict CYP3A4 activity predicts a variable from itself.

**Frozen rule, applied in all three places:**

1. **Edge level.** Every edge whose target symbol matches
   `^CYP[0-9]+[A-Z]*[0-9]*(P[0-9]+)?$` is removed from the network before any
   regulon is built. Measured: 822 of 64,516 edges (1.27 %) have a cytochrome
   target; 11 of them have already left as sign-indeterminate, so this guard
   removes the remaining **811**. All cytochromes go, not only the ten
   endpoints, so that an endpoint's own paralogue cannot re-enter.
2. **Matrix level.** Every probe whose array symbol matches that pattern is
   removed from the expression matrix the scores are computed from. Measured:
   100 admissible probes.
3. **Alias level.** Every probe whose *resolved* human gene matches that
   pattern is removed even when its array symbol does not. Measured: 4 such
   admissible probes — `FLJ20359 → CYP2W1`, `FLJ39501 → CYP4F22` (two probes),
   `FLJ16008 → CYP27C1`. None of these three symbols is a CollecTRI target, so
   on this network the alias path happens to be empty; the guard stays, and the
   audit reports the number found rather than asserting zero.

A test fails if any of the three guards is removed, and the runner stops if a
CYP-resolving probe reaches the scoring matrix.

### 4.4 Candidate regulators, fixed per endpoint

**Frozen rule.** The candidate regulators of endpoint `e` are the transcription
factors with a **sign-resolvable** CollecTRI edge to the enzyme of `e`, read
*before* the cytochrome removal of §4.3.

Two separate points, and they must stay separate. The enzyme's identity selects
the regulators — contextualisation the endpoint's own column name supplies, as
H2a already declares. The enzyme is then removed from every regulon used for
scoring. Selecting is a mechanistic prior; keeping would be circular.

The sign has to be resolvable, because the whole premise under test is a
*signed* network: a regulator whose edge to the enzyme carries no direction
states no mechanism, so it is not a candidate. This costs three enzymes a
regulator or two and is recorded here rather than discovered later.

Measured, with the number that survive the evaluability floor of §4.5:

| Enzyme | regulators in the raw table | sign-resolvable candidates | evaluable |
|---|---|---|---|
| CYP3A4 | 51 | 50 | 49 |
| CYP2B6 | 45 | 43 | 42 |
| CYP1A2 | 31 | 31 | 30 |
| CYP2C9 | 21 | 21 | 19 |
| CYP2E1 | 15 | 15 | 15 |
| CYP2A6 | 12 | 12 | 11 |
| CYP2C19 | 10 | 10 | 10 |
| CYP2C8 | 7 | 7 | 7 |
| CYP2D6 | 6 | 5 | 5 |

Of CYP3A4's 50 candidates, **26 are activating and 24 repressing**, and the
repressing set is the whole NF-κB family (`NFKB1`, `NFKB2`, `REL`, `RELA`,
`RELB` and their dimers), together with `NR0B2`, `NCOR1`, `NCOR2`, `NR1H3`,
`NR2F1`, `NR2F2`, `ESR1`, `MYC` and `NFIL3`. The signed
inflammation-to-CYP3A4 link that motivated this work is therefore already in
the public resource; this experiment tests whether it *predicts*, not whether
it exists.

### 4.5 Evaluability floor and `NOT_EVALUABLE`

**Frozen rules.** A regulator is evaluable for an endpoint only if its
cytochrome-free regulon retains **at least 5** members with an admissible
probe. An endpoint is `NOT_EVALUABLE`:

- for `F1`, `F2` and `F5` when the target enzyme has **no admissible probe**;
- for `F4p` and `F5` when **fewer than 5** evaluable regulators remain;
- for any arm when fewer than **50** donors carry both the endpoint value and
  expression.

Measured consequences, fixed in advance:

- **`CYP1A2_activity` is `NOT_EVALUABLE` for `F1`, `F2` and `F5`.** Its four
  probes are all graded `conflicting` and none is admissible. `F4p` does not
  need the target transcript, so CYP1A2 is evaluable for `F4p` against `F0`
  alone and is reported there, separately, as the one endpoint where the
  functional layer is the only available omic.
- Six regulator-endpoint pairs fall under the regulon floor, with the
  admissible member count that refused each: `NCOA6` 3 (CYP3A4 and CYP2C9),
  `SUPT20H` 4 (CYP2C9), `AHRR` 1 (CYP1A2), `UBP1` 0 (CYP2B6 and CYP2A6). `UBP1`
  has no surviving regulon at all once the sign rule has run.
- No endpoint falls under the regulator floor or the donor floor. The smallest
  regulator set is CYP2D6's five, which clears the floor exactly.
- Donors carrying both the endpoint and expression: 397 for CYP3A4M, CYP3A4T,
  CYP2A6 and CYP2D6; 398 for CYP1A2, CYP2B6, CYP2C9, CYP2C19 and CYP2E1; 387
  for CYP2C8.

`NOT_EVALUABLE` is reported as itself. It is never silently dropped, never
filled by a neighbouring endpoint, and never counted as a negative result.

*Naming note.* The arm written `F4′` in the commissioning text is written
`F4p` in code, filenames and JSON keys, because a prime is not safely
round-trippable through every tool that reads these artefacts. They are the
same arm.

---

## 5. Scoring — the univariate linear model, written out

For donor `i` and regulator `e`, the activity is the t-statistic of the slope in
the ordinary least squares fit of that donor's standardised expression profile
on the regulon's weight vector:

```
x_i  in R^G    donor i's standardised expression over the G scoring genes
w_e  in R^G    +1 / -1 for regulon members, 0 elsewhere
design = [1, w_e]
activity(e, i) = beta_1 / se(beta_1)
```

This is the statistic decoupleR implements as `run_ulm` with CollecTRI as the
network. It is implemented here in the repository, in about twenty lines of
linear algebra, rather than by importing decoupleR: the statistic is short
enough to write and test against its closed form, and a new runtime dependency
would have to be vetted, pinned and justified for a single number. A test
checks the vectorised implementation against a per-donor `numpy.linalg.lstsq`
fit and against the textbook t-statistic.

**Frozen rules.**

- `G` is the admissible, cytochrome-free gene universe of §4.3, one probe per
  gene, **not** the regulon. The regression is over the whole measured profile
  with zeros outside the regulon, as the published statistic requires.
- The probe for a gene is the one with the highest mean **in the training
  fold**, among probes missing in no more than 20 % of training donors; missing
  values take the training-fold mean. This is the rule H2a uses, so that
  F-series and H2a features mean the same thing.
- Standardisation is per gene, using the **training fold's** mean and standard
  deviation, applied unchanged to the test fold. The statistic is unsupervised,
  but it is still computed inside the fold: nothing about the test donors enters
  the transformation.
- A regulator whose weight vector has zero variance over `G` yields no score and
  is dropped with a recorded reason.

---

## 6. Arms

Every arm sees the same donors, the same folds and the same estimator — Ridge
with the penalty chosen by 3-fold inner cross-validation on the training fold,
over the frozen grid `(0.01, 0.1, 1, 10, 100, 1000)`. Arms differ only in which
columns they carry.

| Arm | Columns | What it isolates |
|---|---|---|
| `F0` | age, sex | the covariate floor |
| `F1` | `F0` + target enzyme transcript | **the primary comparator** |
| `F2` | `F1` + the evaluable regulators' own transcripts | the same genes, no network |
| `F4p` | `F0` + the regulators' activities (cytochrome-free regulons) | the functional layer alone |
| `F5` | `F1` + the regulators' activities | **the primary arm** |

`F2` expands a complex regulator (`NFKB1_RELA`) into its named members and
keeps each member's transcript if admissible; a complex with no admissible
member contributes no column. `F2` and `F5` therefore carry a comparable number
of columns derived from the same regulators, which is what makes `F5` versus
`F2` a test of the *aggregation* rather than of the gene list. Measured on the
primary endpoint: `F5` carries 49 activity columns and `F2` 38 transcript
columns, the difference being the complexes whose members collapse onto the
same few proteins — `F2` is therefore the narrower arm, and a `F5` advantage
over it is not automatically an advantage of the network.

`F3` is reserved for a pathway-footprint arm (PROGENy) and is not implemented.
Reason, measured: PROGENy's 14 footprints contain no xenobiotic-metabolism
pathway, so a PROGENy arm would carry no column with a stated mechanistic
relation to a cytochrome. Implementing it would add arms without adding a
hypothesis.

---

## 7. Design, losses and decision rules

- **Folds.** Donor-disjoint repeated nested cross-validation: 5 outer folds × 2
  repeats, 3 inner folds, seed `20261109`. Identical folds for every arm,
  endpoint and control. A donor's absolute error is averaged across repeats
  **before** any loss or bootstrap, so a repeat is never a second donor.
- **Primary loss.** Mean absolute error on the endpoint's own scale.
- **Mandatory robustness co-criteria.** 10 % one-sided trimmed MAE and log MAE.
  A reversal on **both** caps the status at `PARTIALLY_SUPPORTED`.
- **Intervals.** Paired donor bootstrap, 2,000 draws, seed `20261109`.

### 7.1 The one primary contrast

**`F5` versus `F1` on `CYP3A4M_activity`, by MAE.**

CYP3A4 is fixed as the primary enzyme because it has the largest candidate
regulator set (51), because the midazolam reaction is the one CYP3A4 substrate
the source documentation names, and because the proposal this experiment comes
from named CYP3A4 as its demonstrator. None of those three reasons is a result.

`CYP3A4T_activity` (testosterone) is reported as the **same-enzyme replicate**
of the primary endpoint and is never counted as an independent validation: two
reactions measured in the same microsomes from the same donors are one
biological test of one enzyme.

### 7.2 Success criterion, fixed

`F5` is declared better than `F1` on the primary contrast only if **all three**
hold:

1. relative MAE reduction **≥ 10 %**;
2. the 95 % paired bootstrap interval of the MAE difference **excludes zero**;
3. no reversal on **both** robustness co-criteria.

Any other outcome is reported as "no gain detected", which under this design
means exactly that and never "there is no gain". Status ceiling `exploratory`
regardless of direction.

### 7.3 Key secondary contrasts

On the primary endpoint, and indispensable to interpreting it:

- **`F4p` versus `F1`** — can the functional layer stand without the target
  transcript?
- **`F5` versus `F2`** — does weighting by the network beat simply measuring the
  same regulators?

### 7.4 Other secondary contrasts and multiplicity

The same three contrasts on the eight remaining enzymes, with
Bonferroni-adjusted bootstrap intervals inside each contrast family. The primary
contrast is one test and is not adjusted. CYP3A4T is reported but not counted in
the family size, being the same enzyme.

---

## 8. Controls — and why the obvious one is worthless

### 8.1 The full sign inversion is vacuous, and that is provable

Inverting every edge sign is the control a reader expects, and it cannot work.
Under the univariate linear model, replacing `w` by `-w` leaves the residuals
unchanged and negates the slope, so `activity → -activity` exactly. Ridge on
standardised columns responds by negating that column's coefficient, and the
predictions are **identical to floating-point error**. A full inversion
therefore tests nothing; reporting it as an unchanged result would be reporting
an algebraic identity as a biological finding.

It is nonetheless run, at perturbation fraction 1.0, **in order to show the
identity holds** — a prediction difference above `1e-8` would mean the
implementation is not the model this document describes. It is labelled
`algebraic_identity_check`, never `control`.

### 8.2 Partial sign perturbation — the control that bites

**Frozen rule.** For each fraction `f` in `{0.25, 0.50, 0.75}`, flip the sign of
a random `f` of the edges of every regulon, independently per regulon, 20
replicates per fraction, seeds `20261109 + replicate`. A partial flip is not an
overall sign change: it changes the *weighting* across regulon members, so the
score genuinely changes and is no longer the published statistic.

Expected if the network carries information: performance degrades monotonically
with `f`. Expected if it does not: performance is flat, and the reported `F5`
gain, if any, came from the gene list rather than from the signs.

### 8.3 Structurally comparable null networks

**Frozen rule.** Rewire each regulon by resampling its members uniformly from
the scoring universe, preserving exactly the regulon's size and its count of
`+1` and `-1` edges, 20 replicates, seeds `20261109 + 100 + replicate`. The null
network has the same shape and the same signed composition as the real one and
differs only in *which* genes it names.

This is the control that separates "the network knows which genes" from "a
weighted average of two hundred arbitrary transcripts predicts liver enzyme
activity". Given that CYP3A4's evaluable regulons hold a median of 195
admissible members out of 214 (coverage 0.879), ranging from 8 to 1,162, that
second possibility is the serious one.

### 8.4 Destroyed outcome link

**Frozen rule.** Permute the donor order of the scoring expression matrix before
computing activities, keeping the network and the endpoint untouched, 20
replicates, seeds `20261109 + 200 + replicate`. Any apparent `F5` gain surviving
this is an artefact of dimension, not of biology.

### 8.5 Reading the controls

Controls qualify the primary contrast; they cannot create it. If the primary
contrast shows no gain, the controls are reported as executed and the conclusion
stands on the primary contrast alone. If the primary contrast shows a gain that
the null networks of §8.3 reproduce, **the gain is reported as refuted by its
own control** and the conclusion is negative.

---

## 9. Three separate reports, never merged

The run writes three independent blocks, and no block's conclusion is allowed to
borrow another's.

1. **`predictive_performance`** — MAE, trimmed MAE, log MAE per arm and
   endpoint; the contrasts of §7; the controls of §8. This block alone decides
   the primary question.
2. **`biological_relevance`** — which regulators receive non-negligible Ridge
   weight; whether the sign of that weight agrees with CollecTRI's recorded sign
   for the regulator's edge to the target enzyme (the edge removed from scoring
   is kept for exactly this check); regulon sizes and admissible coverage. A
   sign agreement above chance is biologically interesting **and is not
   predictive evidence**; a model can be right about direction and useless about
   magnitude.
3. **`analytical_validity`** — probe admissibility counts; the three cytochrome
   guards and what each removed; the alias audit; `NOT_EVALUABLE` endpoints with
   reasons; missingness; fold fingerprints; the algebraic identity check.

A favourable number in block 2 or 3 never upgrades block 1.

---

## 10. What may not be done after seeing the results

- No arm is added, removed or redefined.
- The primary endpoint, contrast, loss and threshold do not move.
- The regulator sets, the evaluability floors, the sign rule and the fold seed do
  not move.
- No endpoint is promoted from secondary to primary, and `CYP3A4T` never becomes
  an independent confirmation.
- A control is not reinterpreted after its result is known.
- If an unforeseen defect is found, it is recorded as a dated amendment with the
  result it would change, and the affected number is withdrawn rather than
  replaced in place.

---

## 11. Deliberately not built

| Not built | Why, measured where possible |
|---|---|
| SIGNOR layer | 102 of 40,677 rows touch any cytochrome, nearly all steroidogenic; **zero** inflammation-to-CYP edges. It cannot carry this layer. |
| PROGENy arm (`F3`) | 14 footprints, none xenobiotic. No column would have a stated relation to a cytochrome. |
| COSMOS / causal flow | Needs paired transcriptome and metabolome in the same donors; HLC has no metabolome. |
| Scaling to clearance | Needs clinical pharmacokinetics in the donors who have the transcriptome. No registered public dataset has both. |
| Clinical transfer | Same reason. Declared `ABSTAIN` with a named data requirement, not scheduled. |

The data requirement, stated so it can be checked rather than hoped for:
individual-level liver or intestinal transcriptome **and** individual
concentration-time pharmacokinetics **in the same subjects**, publicly
downloadable without authorisation. The project's registry currently contains no
such dataset; ARIANA has clinical pharmacokinetics and no transcriptome, HLC has
a transcriptome and no pharmacokinetics.

---

## 12. Freeze

| Item | Value |
|---|---|
| Fold seed | `20261109` |
| Outer folds × repeats × inner | 5 × 2 × 3 |
| Bootstrap draws | 2,000 |
| Ridge grid | 0.01, 0.1, 1, 10, 100, 1000 |
| Relative MAE threshold | 0.10 |
| Trim fraction | 0.10 |
| Regulon floor (admissible members) | 5 |
| Regulator floor per endpoint | 5 |
| Donor floor per endpoint | 50 |
| Max probe missing fraction | 0.20 |
| Control replicates | 20 per setting |
| Perturbation fractions | 0.25, 0.50, 0.75 (+ 1.00 as identity check) |
| Status ceiling | `exploratory` |

`scripts/run_f_series_hlc_functional_omics.py` re-reads this table from
`omics2pk.functional_omics.FROZEN` and stops before loading any data if a value
differs. `tests/test_functional_omics.py` fails if a value in this document and
the module disagree.

---

## 13. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-09 | First freeze. Written after measuring the sources and before fitting any F-series model. Records in advance: the unfavourable prior (§1); that the full sign inversion is an algebraic identity rather than a control, with the partial perturbation and the size-and-sign-preserving null replacing it (§8); the three-place cytochrome removal with its measured cost (§4.3); `CYP1A2_activity` as `NOT_EVALUABLE` for the transcript arms because none of its four probes is admissible (§4.5); CYP3A4T as a same-enzyme replicate and not a second validation (§7.1); and the clinical-transfer `ABSTAIN` with its data requirement (§11). |
| 1.0.0 | 2026-11-09 | **Two corrections to this document's own first draft, made before the freeze and before any model was fitted, both found by running the implementation against the sources with `--dry-run`.** (1) The sign rule had been written from the `is_stimulation` / `is_inhibition` flags and the counts quoted accordingly; the implementation reads the consensus columns, and the two do not agree — 1,310 single-flag edges are indeterminate under the consensus and 3,496 both-flag edges resolve through it. §4.1 now states the rule as "the consensus decides" and carries the consensus counts (2,590 dropped, not 1,280). (2) The candidate-regulator rule said "whose regulon contains the enzyme"; the implementation additionally requires the regulator-to-enzyme edge to have a resolvable sign, which is the stricter and the correct reading of a *signed* network, and which costs CYP3A4 one candidate, CYP2B6 two and CYP2D6 one. §4.4 now states that requirement and gives both counts. No threshold, floor, arm, endpoint, loss or primary contrast changed. |
