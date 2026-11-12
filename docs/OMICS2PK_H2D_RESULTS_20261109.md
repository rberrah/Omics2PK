# H2d — hybrid mechanistic selection on the Human Liver Cohort

**Document version:** 1.1.0
**Date:** 2026-11-09
**Prespecification:** `OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md` v1.0.0
**Erratum carried:** `OMICS2PK_H2D_ERRATUM_20261108.md` v1.0.2
**Protocol:** `OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md`, run under v1.4.5, read under v1.4.7 (which changed no computation)
**Universes:** pool `A_qualified_human` (21,634 candidate symbols) and pool `B_qualified_human_on_network` (11,030), both under rule `probe-identity-1.1.0`. Sections 1–8 report pool A; §10 compares the two.
**Endpoint class:** `functional_in_vitro`. Nothing here can promote a clinical claim above `PARTIALLY_SUPPORTED`.

**What the nine endpoints are, and are not.** They are nine endpoint
measurements from a **single cohort**, none of which was used to design the
rule. Their pairwise Spearman correlations run from 0.144 to 0.841 with a median
of 0.446, the first principal component takes 49.6 % of their variance, and
their correlation matrix carries a **spectral dependence index of about 3.55** —
an index, never a count of independent replications. The two CYP3A4 substrates
correlate at 0.841 and count as **one** endpoint in any tally. No sentence here
may read "confirmed on nine endpoints".

**Why the generating endpoint is separate.** `CYP2D6_activity` produced the
hypothesis, so its result is capped at `exploratory` and reported in §6 alone.

---

## 1. The rule, in one paragraph

Knowledge proposes, data disposes. The outcome-blind mechanistic ordering of
protocol §5.1 supplies a shortlist of three times the panel budget; the training
donors then keep, demote or abstain, by partial rank correlation against the
clinical baseline, calibrated on a permutation null at a 5 % target false
discovery rate. Folds, donors, predictor, panel grid, loss hierarchy, bootstrap,
one-feature-per-gene universe, missing-value rule and baseline covariates are
**imported from the H2a runner unchanged**, so the only difference between the
hybrid arm and the blind mechanistic arm is the retention step.

---

## 2. The declared failure mode did not occur

The prespecification named, in advance, the outcome most likely to make the
hybrid worthless: **collapse onto the univariate selector**, which would mean
the knowledge contributed nothing beyond shortlisting.

**It does not happen on any endpoint.** Exact set equality with the univariate
top-`k` is **0.00 at every panel size from 5 upward, on all nine endpoints**
(0.0–1.0 at `k = 1`, where every arm picks the same single transcript and
equality means nothing). The hybrid selects genuinely different panels.

**This is a control that passed, not evidence of mechanistic utility.** A
different panel can predict as well, better or worse; difference alone says
nothing about whether the knowledge contributed information. The question "did
the mechanism add anything?" is answered in §3 and §4, where the comparison is
against error, not against overlap — and the honest summary of those sections is
that the knowledge adds something to *its own blind ordering* on a minority of
cells and nothing to the univariate filter.

---

## 3. Prediction: non-inferior in three cells of five, better than its own blind arm in one of seven

72 cells — nine endpoints by eight panel sizes.

| endpoint | non-inferior to **all three** statistical arms, at budget | …at matched dimension | better than the **blind** arm, interval excluding zero | **worse** than the univariate arm, interval excluding zero | repair criterion |
|---|---:|---:|---:|---:|---|
| CYP1A2 *(no admissible target feature)* | 1/8 | 1/8 | 0/8 | 4/8 | not assessable |
| CYP2A6 | 7/8 | 7/8 | 1/8 | 1/8 | not met |
| CYP2B6 | **8/8** | **8/8** | 0/8 | 0/8 | not met |
| CYP2C8 | 7/8 | 7/8 | 2/8 | 0/8 | **met** |
| CYP2C9 | 2/8 | 2/8 | 3/8 | 0/8 | **met** |
| CYP2C19 | 0/8 | 0/8 | 0/8 | **8/8** | not met |
| CYP2E1 | **8/8** | **8/8** | 1/8 | 6/8 | not met |
| CYP3A4 (midazolam) | 5/8 | 5/8 | 1/8 | 2/8 | not met |
| CYP3A4 (testosterone) | 6/8 | 6/8 | 2/8 | 0/8 | not met |
| **total** | **44/72** | **44/72** | **10/72** | **21/72** | 2 of 8 assessable |

Three readings, all needed:

- **Non-inferiority holds in 44 cells of 72**, and identically at matched
  dimension. The identity is **expected, not corroborating**: the hybrid filled
  its budget in every fold of every endpoint except a few of CYP2E1's at
  `k = 100`, so the two comparisons are made at the same dimension, and the
  matched-dimension loss is numerically identical to the budget loss in **71
  cells of 72**. The matched-dimension column exists to catch a hidden advantage
  that did not arise here, because the hybrid never became smaller (§5).
- **Against its own blind arm the hybrid is better in 10 cells and worse in 3.**
  The retention step is therefore a real improvement on knowledge alone in a
  minority of cells and a real cost in a smaller one — not a uniform gain.
- **Against the univariate filter it is detectably worse in 21 cells and
  detectably better in 1** (CYP2C8 at `k = 1`, −2.3 %), the rest being
  undetermined. The worse cells concentrate on CYP2C19 (8 of 8) and CYP2E1 (6 of
  8). The simplest outcome-reading selector remains the arm to beat, and the
  hybrid does not beat it.

### 3.1 How these counts may be read, and how they may not

**72 cells are not 72 experiments.** The eight panel sizes are nested, every
cell uses the same 387 to 398 donors, and the endpoints are correlated — the
dependence index recorded in the header is about 3.55 for ten endpoint
measurements. Each interval is a 95 % bootstrap interval on **one** paired
contrast, computed on 2,000 resamples of donor-level errors.

Protocol §9 places every H2 comparison in the **exploratory** class and requires
Bonferroni adjustment within a family sharing an endpoint, with the nominal
level stated. The counts in the table above are therefore **descriptive
summaries at a nominal 95 % level per contrast**, and no confirmatory claim is
made from any of them: not from 44 of 72, not from 10 of 72, not from 21 of 72.
The frozen run reports unadjusted 95 % intervals and stores no bootstrap draws,
so an adjusted interval cannot be recovered from it after the fact; that is a
limit of this run, stated rather than worked around.

What the counts legitimately support is a **direction and a spread**: the hybrid
is usually within the margin of the statistical arms, it improves on its own
blind ordering in a minority of cells, and it loses to the univariate filter on
two endpoints in most of their cells. Each individual endpoint's repair contrast
in §4 is a single prespecified test and is read as such.

---

## 4. The destroyed-signal repair: met on two endpoints of eight assessable

This is the experiment the rule was designed for. The target gene's donor order
is permuted within each cross-validation block, so the gene keeps its rank in
the mechanistic ordering while carrying nothing. The blind arm must keep
choosing it; the hybrid must demote it.

**The hybrid demoted the dead gene in 10 folds of 10, on every endpoint.** The
detection half of the rule works without exception.

| endpoint | hybrid | blind | univariate | against blind | against univariate | criterion |
|---|---:|---:|---:|---|---|---|
| CYP2C8 | **52.27** | 63.09 | 54.10 | −17.2 %, interval [−21.9, −12.2] | **−3.4 %, interval [−5.8, −1.0]: better** | **met** |
| CYP2C9 | **247.62** | 274.49 | 245.88 | −9.8 %, interval [−14.1, −5.5] | +0.7 %, interval [−3.4, +4.9]: within margin | **met** |
| CYP2C19 | 51.29 | 55.51 | 49.07 | −7.6 %, interval [−11.6, −3.8] | +4.5 %, interval [+0.7, +8.7] | inconclusive at the margin |
| CYP2A6 | 600.78 | 644.28 | 577.80 | −6.8 %, interval [−10.6, −3.1] | +4.0 %, interval [+1.3, +6.7] | inconclusive at the margin |
| CYP2B6 | 257.59 | 269.27 | 244.33 | −4.3 %, interval [−7.7, −1.1] | +5.4 %, interval [+2.2, +9.0] | not met |
| CYP3A4 (testosterone) | 3359.73 | 3638.86 | 3133.52 | −7.7 %, interval [−11.3, −4.0] | +7.2 %, interval [+3.6, +11.1] | not met |
| CYP3A4 (midazolam) | 801.73 | 876.38 | 714.26 | −8.5 %, interval [−11.9, −4.9] | +12.2 %, interval [+8.3, +16.6] | not met |
| CYP2E1 | 1334.20 | 1338.37 | 1316.25 | −0.3 %, **not detectable** | +1.4 %, within margin | undiscriminating |
| CYP1A2 | — | — | — | — | — | **not assessable** |

The criterion needs **both** halves: within the frozen margin of the univariate
arm, *and* an interval against the blind arm excluding zero on the favourable
side. The second half holds on **seven endpoints of eight assessable** — the
hybrid repairs 4.3 % to 17.2 % of the blind arm's loss, with intervals
excluding zero, everywhere except CYP2E1. The first half is what decides, and the outcomes
are of four kinds, which must not be merged:

1. **Met, on CYP2C8 and CYP2C9.** On CYP2C8 the hybrid is **better than the
   univariate arm** by 3.4 % with an interval excluding zero, while repairing
   17.2 % against its own blind arm; on CYP2C9 it is 0.7 % above with an
   interval inside the margin.
2. **Inconclusive at the margin, on CYP2C19 and CYP2A6.** The point estimate
   sits inside the 5 % margin (+4.5 % and +4.0 %) while the interval runs past
   it (+8.7 % and +6.7 %). By the protocol's own rule that is `INCONCLUSIVE`,
   not a failure, and not a success either.
3. **Genuine inferiority, on CYP2B6 and both CYP3A4 substrates.** The point
   estimate is beyond the margin: +5.4 %, +7.2 % and +12.2 %. The hybrid repairs
   a large share of the blind arm's loss and remains clearly worse than the
   univariate filter.
4. **Undiscriminating, on CYP2E1.** Every arm sits at the clinical baseline
   (1337.24), so no interval can exclude zero. The experiment did not
   discriminate; it did not refute the rule.

**Where the univariate arm's own number comes from matters.** Protocol §5.5 as
amended in 1.4.7 records that this comparator's error on a wrecked matrix is
only interpretable once the feature it falls back on is named. On CYP2A6 it
falls back on `CYP2A7`, a human paralogue sharing about 95 % of its sequence
with the target; on CYP2B6 on `CYP2B7P1`, a pseudogene of the same family. Both
are admissible and **no annotation can expose them as duplicates of the
destroyed transcript**. So on two of the three endpoints where the hybrid loses
this comparison — CYP2A6 inconclusively, CYP2B6 beyond the margin — the
comparator is plausibly reading the same molecule under another name rather than
finding a replacement for it.

That is **not** a defence of the hybrid. It is a limit of the comparison, it
cuts in the hybrid's favour, and it is stated because the opposite reading would
be unearned. It does not apply to the CYP3A4 pair, where the comparator's
fallback is `IL12B-AS1`, a different gene; the hybrid's inferiority there is not
explained away.

---

## 5. Neither abstention nor parsimony was demonstrated

One of the three capabilities the hybrid was meant to demonstrate is recognising
when its mechanistic knowledge is insufficient. **It never fired.** Across nine
endpoints, eight panel sizes and ten folds, the abstention count is **zero**,
and the retained count equals the budget in every fold but a few of CYP2E1's at
`k = 100` (minimum 92 of 100).

The reason is structural, not accidental: the shortlist is three times the
budget, so a demoted candidate is replaced by the next one down the mechanistic
order, and the panel refills. The retention step is doing work — on CYP2C9 at
`k = 10` it demoted CYP2B6, CYP1A1, CYP2D6, CYP3A5, CYP4A11, UGT1A6, AOX1,
UGT1A4, PPIG and UGT1A1, that is the blind arm's **entire** panel, and replaced
all ten — but it never runs out of candidates that pass its null.

**Two consequences, and the second was missing from version 1.0.0.**

- **No evidence at all about abstention.** No text may claim the rule abstains
  appropriately on the basis of this run. Testing abstention needs a design
  where the mechanistic shortlist genuinely contains nothing informative, which
  this cohort's endpoints do not provide.
- **No evidence about parsimony either.** The hybrid does not reduce the number
  of variables: it fills `k` almost always, so what it changes is *which* `k`
  genes are used, not *how many*. A claim that hybrid selection buys parsimony
  would be unsupported here, and the matched-dimension column cannot support it
  either, for the reason given in §3.

**A measurement this run cannot give.** How much of the mechanistic ordering
survives retention is not recoverable from the committed report: the runner
stores at most ten demoted examples for the first fold of each panel size, not
the count of demotions per fold. The one visible case is striking — on CYP2C9 at
`k = 10` the blind arm's entire top ten was demoted and replaced — and it raises
the right question, which is how much of the initial knowledge the hybrid
actually keeps. The run as frozen cannot answer it, and no answer is implied.

---

## 6. The generating endpoint, capped at exploratory

`CYP2D6_activity` produced the hypothesis and cannot confirm it. Reported for
completeness.

| k | retained | hybrid | blind | univariate | LASSO | XGBoost |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 91.43 | 91.43 | 91.43 | 91.43 | 91.43 |
| 10 | 10 | 92.77 | 92.78 | 92.07 | 95.15 | 94.66 |
| 50 | 50 | 94.90 | 96.65 | 96.51 | 110.36 | 100.29 |
| 100 | 100 | **92.58** | 97.23 | 99.40 | 110.29 | 104.92 |

On the destroyed-signal test the hybrid demoted the dead gene in 10 folds of 10
and reached 107.56 against the blind arm's 109.53 and the univariate arm's
108.88, with an interval of [−4.10, +0.25] against blind: **no difference
detected**, criterion not met. Under this pool every arm is at the clinical
baseline of 109.64, because the duplicate the univariate arm used to exploit
(`CYP2D17`) is refused and nothing admissible replaces it. **The experiment has
nothing to find on this endpoint, so it cannot discriminate between arms** — the
reading protocol §5.5 makes mandatory.

---

## 7. Analytical traceability of every selection

Because the probe standing for a gene is chosen per training fold, a selection
recorded as a gene symbol is not traceable to a measurement. Every feature every
arm actually used was therefore checked against the pool, fold by fold:

- about **9,560 selections checked per endpoint**, and **zero** outside the
  qualified pool. The restriction held; it was not assumed to hold.
- **1,106 symbols traced** across the nine endpoints for the hybrid's own panels
  and the target, of which **28 used a different probe in different folds** —
  direct evidence that the per-fold probe choice is not constant, and that a
  symbol-level audit would have missed it.
- the hybrid selected a decoy **zero** times. That is not a merit: a decoy
  carries the symbol of a randomly drawn network gene, almost always far from
  the seed, so it was never a candidate for the top of the shortlist. The
  statistical arms made 986 decoy selections in total, which is what the decoy
  block exists to measure.

---

## 8. What is not claimed

- Not that the hybrid is better than statistical selection. Across 72 cells it
  is detectably worse than the univariate filter in 21 and detectably better in
  1.
- Not that the repair generalises. The criterion is met on two endpoints of
  eight assessable, both in the CYP2C family, and both are endpoints where the
  univariate comparator had no stable replacement to exploit. What does hold
  more widely is the weaker statement: the retention step repairs a detectable
  share of the blind arm's loss on seven endpoints of eight.
- Not that the rule abstains appropriately: it never abstained (§5).
- Not nine independent replications. Nine endpoint measurements from one cohort
  with a spectral dependence index of about 3.55, the two CYP3A4 substrates
  counting as one.
- Not a clinical pharmacokinetic result. The endpoint is microsomal activity.
- Not a demonstrated calibration of the retention null. The prespecification's
  120 draws across four regimes are a sanity check, not a calibration proof, and
  the erratum says so.
- Not parsimony, and not abstention (§5).
- Not a confirmatory reading of any cell count (§3.1).

---

## 9. Provisional verdict

**What is established.** The hybrid recognises, systematically and on every
endpoint, that a transcript the mechanistic ordering insists on has lost its
information, where the blind arm keeps privileging it. That is a targeted
feasibility proof of **mechanistic revision**, and the complete repair criterion
is met on two of eight assessable endpoints.

**What is not.** No general predictive superiority: the hybrid is non-inferior
to the statistical arms in most cells and detectably worse than the univariate
filter on two endpoints in most of theirs. No parsimony. No abstention.

**The admissible form of words**, which this document adopts: *hybrid selection
mitigates a weakness of the blind mechanistic ordering, without guaranteeing
performance equivalent to statistical selectors.*

**What follows for the programme.** This is not a justification for making
hybrid selection NetPK's default selection method. It is a candidate
**experimental module** for the question it actually answers — when should a
mechanistic hypothesis be kept, and when should it be contested — and that is
where it is proposed. §10 narrows the feasibility claim further: of the two
endpoints meeting the complete criterion, only one meets it in both universes.
The selector, the threshold and the abstention rule are **not** modified in the
light of these results.

---

## 10. Pool A against pool B, and what moves between them

Pool B is pool A intersected with the genes the public network contains: 11,030
candidate symbols instead of 21,634. Both runs are frozen, use the same folds,
seeds and rule, and differ only in the candidate universe.

### 10.1 The arms under test are unaffected; only the comparators move

| arm | cells where pool A and pool B give the same MAE | largest relative change |
|---|---:|---|
| hybrid, at budget | **72 of 72** | 0.000 |
| blind mechanistic | **72 of 72** | 0.000 |
| univariate Spearman | 9 of 72 | 2.3 % |
| LASSO entry order | 1 of 72 | 11.1 % |
| XGBoost gain | 0 of 72 | 6.1 % |

The hybrid and the blind arm are **numerically identical** in every one of the
72 cells, and so are their contrasts against each other on the destroyed-signal
test. That is not a coincidence and it is worth stating plainly: both draw their
candidates from the mechanistic ordering, whose genes are network genes by
construction, and pool B removes only genes the network does not contain. Those
arms never selected such a gene, so nothing changes for them.

**Every verdict that moves between pool A and pool B therefore moves because the
comparator changed, not because the method under test did.** The A-to-B
comparison measures how much each statistical selector relied on features
outside the network; it is not a sensitivity analysis of the hybrid.

### 10.2 The tallies

| | pool A | pool B |
|---|---:|---:|
| non-inferior to all three statistical arms, at budget | 44/72 | 49/72 |
| better than the blind arm, interval excluding zero | 10/72 | 10/72 |
| worse than the univariate filter, interval excluding zero | 21/72 | 22/72 |
| collapse onto the univariate selector | none (0.00) | none (worst 0.10, on CYP2B6) |
| complete repair criterion met | CYP2C8, CYP2C9 | CYP2C8, **CYP2B6** |

Descriptive counts at a nominal 95 % level per contrast, read under §3.1.

### 10.3 The repair criterion is robust on one endpoint, not two

| endpoint | against blind, both pools | against univariate, pool A | against univariate, pool B | criterion |
|---|---|---|---|---|
| **CYP2C8** | −17.2 % | −3.4 % [−5.8, −1.0] | −3.4 % [−5.8, −1.0] | **met in both** |
| CYP2C9 | −9.8 % | +0.7 % [−3.4, +4.9] | +2.4 % [−1.7, +6.9] | met in A, inconclusive in B |
| CYP2B6 | −4.3 % | +5.4 % [+2.2, +9.0] | +0.7 % [−1.3, +2.9] | not met in A, met in B |
| CYP2A6 | −6.8 % | +4.0 % [+1.3, +6.7] | +4.0 % [+1.3, +6.7] | inconclusive in both |
| CYP2C19 | −7.6 % | +4.5 % [+0.7, +8.7] | +4.5 % [+0.7, +8.7] | inconclusive in both |
| CYP3A4 (midazolam) | −8.5 % | +12.2 % [+8.3, +16.6] | +9.0 % [+5.6, +12.6] | not met in both |
| CYP3A4 (testosterone) | −7.7 % | +7.2 % [+3.6, +11.1] | +5.0 % [+1.8, +8.2] | not met in both |
| CYP2E1 | −0.3 %, not detectable | +1.4 % | +1.0 % | undiscriminating in both |

**Only CYP2C8 meets the complete criterion in both universes.** The column
against the blind arm is identical in both, as §10.1 predicts, so every movement
in the last column is a movement of the comparator:

- **CYP2C9 loses the criterion in pool B** because the univariate arm got
  *better* when its universe shrank — its destroyed-signal error falls from
  245.88 to 241.75 — and the hybrid, unchanged at 247.62, slips from +0.7 % to
  +2.4 %. The point estimate stays well inside the 5 % margin; it is the
  interval that now runs past it, so the outcome is `INCONCLUSIVE`, not a
  reversal.
- **CYP2B6 gains the criterion in pool B** because the univariate arm got
  *worse*: it loses `CYP2B7P1`, the same-family pseudogene it was falling back
  on, which is admissible but absent from the network, and its error rises from
  244.33 to 255.82. The hybrid did nothing differently.

That second case is the one to keep in mind. A criterion that a method passes
because its comparator lost a feature is not evidence about the method, and the
pass is reported here only so that the asymmetry is visible. **The defensible
feasibility claim is CYP2C8 alone**, where the hybrid is better than the
univariate arm by 3.4 % with an interval excluding zero in both universes while
repairing 17.2 % of the blind arm's loss.

### 10.4 Stability of the rest

Non-inferiority rises from 44 to 49 cells, concentrated on the two CYP3A4
substrates (5→7 and 6→7) and CYP1A2 (1→2), for the same reason: the comparators
lost reach. The count of cells where the hybrid beats its own blind arm is
**10 of 72 in both pools**. Collapse does not occur in either. Abstention never
fires in either.

### 10.5 The generating endpoint under pool B

`CYP2D6_activity`, capped at `exploratory`: the hybrid reaches 92.58 at `k = 100`
against the blind arm's 97.23 and the univariate arm's 99.52, repairs nothing
detectable on the destroyed-signal test (hybrid 107.56, blind 109.53, univariate
108.20, all at the clinical baseline of 109.64), and does not collapse. Identical
to pool A for the two mechanistic arms, as everywhere else.
