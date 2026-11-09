# H2d — hybrid mechanistic selection on the Human Liver Cohort

**Document version:** 1.0.0
**Date:** 2026-11-09
**Prespecification:** `OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md` v1.0.0
**Erratum carried:** `OMICS2PK_H2D_ERRATUM_20261108.md` v1.0.2
**Protocol:** `OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md`, run under v1.4.5, read under v1.4.7 (which changed no computation)
**Universe:** pool `A_qualified_human`, rule `probe-identity-1.1.0`, 21,634 candidate symbols
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
  dimension, which matters because the hybrid is allowed to retain fewer
  variables than its budget (see §5) and crediting selection quality for
  dimension reduction would be a hidden advantage.
- **Against its own blind arm the hybrid is better in 10 cells and worse in 3.**
  The retention step is therefore a real improvement on knowledge alone in a
  minority of cells and a real cost in a smaller one — not a uniform gain.
- **Against the univariate filter it is detectably worse in 21 cells and
  detectably better in 1** (CYP2C8 at `k = 1`, −2.3 %), the rest being
  undetermined. The worse cells concentrate on CYP2C19 (8 of 8) and CYP2E1 (6 of
  8). The simplest outcome-reading selector remains the arm to beat, and the
  hybrid does not beat it.

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

## 5. Abstention was never exercised, so H2d says nothing about it

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

**Consequence.** H2d supplies evidence about parsimony and about resistance to a
destroyed signal. It supplies **no evidence at all** about abstention, and no
text may claim the rule abstains appropriately on the basis of this run. Testing
abstention needs a design where the mechanistic shortlist genuinely contains
nothing informative, which this cohort's endpoints do not provide.

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
