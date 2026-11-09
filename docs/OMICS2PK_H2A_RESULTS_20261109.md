# H2a — mechanistic versus statistical selection on the Human Liver Cohort

**Document version:** 1.1.0
**Date:** 2026-11-09
**Protocol:** `OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md` v1.4.7
**Erratum carried:** `OMICS2PK_H2A_ERRATUM_20261109.md` v1.0.3
**Endpoint class:** `functional_in_vitro` — measured microsomal activity, never
pooled with a clinical pharmacokinetic result
**Status of the H2 claim:** `PARTIALLY_SUPPORTED` at best from this evidence, by
the status cap in §5.3 of the protocol. Nothing here can promote the clinical
form of the claim.

---

## 1. What was asked, and what was run

The question is not whether a mechanistic ordering obtains the lowest error. It
is whether mechanism-guided selection is **reliable**: whether it can work with
few variables, whether it resists a signal that is not there, and whether it
recognises the cases where its own knowledge is insufficient.

H2a is the **contextualised** form of that question, and the label matters. The
endpoint's own name supplies the enzyme that performs the measured reaction, so
the chain is told the mechanism and asked only which *further* genes to retain.
This is not autonomous discovery of a mechanism; that is H2b, which is not run
here.

Ten enzyme-activity endpoints, 387 to 398 donors each, donor-disjoint repeated
nested cross-validation (5 outer folds × 2 repeats, 3 inner folds, seed
`20261108`), a donor's absolute errors averaged across repeats before any loss
or bootstrap. Panel sizes 1, 2, 3, 5, 10, 20, 50, 100. Baseline covariates age
and sex. Eight arms, every one selecting the same *number* of variables:

1. **NetPK**, the executable mechanistic ranking: public network, random walk
   with restart, seeds obtained automatically from the endpoint's own column
   name; no gene symbol appears anywhere in the selector;
2. NetPK seeded with a **deliberately wrong** enzyme, by a rotation fixed before
   running;
3. NetPK seeded from a **public pathway annotation**;
4. **univariate** Spearman filtering;
5. **LASSO** entry order;
6. **XGBoost** gain;
7. the **target transcript alone** — the information the endpoint's label gives
   away for free;
8. **no selection**, and the **baseline** and a **random panel** (≥ 100
   replications) as references.

---

## 2. Primary result, on the unrestricted universe

Run under protocol v1.4.3, 31,862 candidate genes, of which 11,182 are reachable
by the mechanistic ordering. Frozen and not recomputed; the pool-restricted runs
in §4 are additions, not replacements.

| endpoint | n | non-inferiority cells (of 8) | NetPK best MAE | univariate best MAE | target alone | baseline | no selection |
|---|---:|---:|---:|---:|---:|---:|---:|
| CYP1A2 *(withdrawn, §5)* | 398 | 0 | 434.62 | 419.85 | 492.58 | 492.43 | 451.17 |
| CYP2A6 | 397 | 5 | 572.28 | 568.93 | 572.28 | 645.41 | 751.24 |
| CYP2B6 | 398 | **8** | 236.66 | 240.10 | 241.59 | 270.37 | 286.82 |
| CYP2C8 | 387 | **8** | 51.13 | 50.91 | 52.97 | 62.66 | 59.96 |
| CYP2C9 | 398 | 1 | 238.53 | 230.49 | 247.96 | 273.65 | 234.02 |
| CYP2C19 | 398 | 0 | 48.60 | 46.06 | 51.18 | 55.52 | 57.67 |
| CYP2D6 | 397 | **8** | 91.43 | 91.43 | 91.43 | 109.64 | 116.74 |
| CYP2E1 | 398 | 7 | 1331.31 | 1302.14 | 1346.03 | 1337.24 | 1725.20 |
| CYP3A4 (midazolam) | 397 | 3 | 726.61 | 712.78 | 726.61 | 872.13 | 826.30 |
| CYP3A4 (testosterone) | 397 | 5 | 3181.27 | 3130.26 | 3181.27 | 3628.37 | 3321.41 |

A **non-inferiority cell** is one endpoint at one panel size where NetPK stays
within the frozen 5 % relative margin of **every one of the three** statistical
methods separately. Choosing the comparator by its out-of-fold error would use
the test folds twice, so no "best statistical arm" exists in this report.

**45 of 80 cells.** Non-inferiority holds at every panel size on three endpoints
(CYP2B6, CYP2C8, CYP2D6) and at none on two (CYP1A2, CYP2C19).

**How this count may be read, and how it may not.** 80 cells are not 80
experiments: the eight panel sizes are nested, every cell uses the same 387 to
398 donors, and the endpoints are correlated, with a spectral dependence index
of about 3.55 for the ten measurements. Each interval is a 95 % bootstrap
interval on one paired contrast. Protocol §9 places every H2 comparison in the
**exploratory** class and requires Bonferroni adjustment within a family sharing
an endpoint, with the nominal level stated; the counts here are therefore
**descriptive summaries at a nominal 95 % level per contrast**, and no
confirmatory claim is made from any of them. What they support is a direction
and a spread, not a pass rate.

### 2.1 Superiority, stated per comparator rather than in one word

| comparator | cells where NetPK is better with an interval excluding zero (of 80) |
|---|---:|
| univariate Spearman filter | **1** |
| LASSO entry order | 31 |
| XGBoost gain | 27 |
| **all three at once** | **1** (CYP2D6 at k = 100) |

This is the result, and it needs both halves. NetPK **is** better than LASSO
entry order and XGBoost gain in a substantial minority of cells — roughly a
third of them — and that is not noise: the intervals exclude zero. Against the
**univariate filter** it is better in one cell out of eighty. A simple
outcome-reading filter is the comparator that matters, and mechanistic ordering
does not beat it here.

An earlier internal summary of this run said "NetPK is never superior". That was
too coarse: it is accurate against the univariate arm and against the three
jointly, and wrong about LASSO and XGBoost.

### 2.2 The wrong-mechanism control is void on one endpoint of ten

The rotation that assigns each endpoint a deliberately wrong enzyme was fixed
before running, and it does not guarantee a *different* gene: this cohort
measures one enzyme with two substrates, so for CYP3A4 (midazolam) the next
label resolves to the same gene. For that endpoint the control compares an arm
against itself and says nothing. **9 of 10 endpoints carry a valid control; 1 is
void**, declared by the merge script rather than left to be noticed in a table
of near-identical numbers. The frozen rotation is not re-run to repair it.

### 2.3 No selection helps as often as it harms

Against the age-and-sex baseline, using every admissible probe without selection
is **worse on 5 endpoints and better on 5** — by +388 on CYP2E1 and by −307 on
CYP3A4 (testosterone). High-dimensional molecular data is therefore not uniformly
harmful nor uniformly useful at this sample size, which is the direct evidence
for H4 and is reported there as well.

---

## 3. The destroyed-signal diagnostic, and what it actually measures

The target gene's donor order is permuted. The gene keeps its place in the
mechanistic ordering while carrying no information at all, so an outcome-blind
arm must keep choosing it and an arm that reads the endpoint should choose
something else.

| endpoint | blind mechanistic arm, k = 1 | univariate arm | folds where the univariate arm still picked the dead gene |
|---|---:|---:|---:|
| CYP1A2 | 491.93 | 450.89 | 0/10 |
| CYP2A6 | 645.68 | 571.97 | 0/10 |
| CYP2B6 | 269.44 | 244.33 | 0/10 |
| CYP2C8 | 62.56 | 54.01 | 0/10 |
| CYP2C9 | 273.81 | 245.88 | 0/10 |
| CYP2C19 | 55.54 | 49.07 | 0/10 |
| CYP2D6 | 110.15 | 96.86 | 0/10 |
| CYP2E1 | 1338.27 | 1319.89 | 0/10 |
| CYP3A4 (midazolam) | 877.86 | 720.88 | 0/10 |
| CYP3A4 (testosterone) | 3662.00 | 3133.52 | 0/10 |

On every endpoint the blind arm falls to about the baseline, and the univariate
arm was **never once fooled, in 100 folds**: it never selected the dead gene.

### 3.1 What the lower error of the univariate arm actually was

The table above was first read as the price of being outcome-blind: the arm that
reads the endpoint notices the destruction *and recovers*, the blind arm cannot.
The first half is right. The second cannot be read from the number at all, and
establishing that took naming the feature each arm fell back on, fold by fold,
on every endpoint and in both universes.

| endpoint | replacement, unrestricted | replacement, pool A | what it is |
|---|---|---|---|
| CYP2A6 | `CYP2A24` 6/10, `CYP2A23` 4/10 — refused | **`CYP2A7` 10/10 — admissible** | a human paralogue sharing about 95 % of its sequence with the target |
| CYP2B6 | `CYP2B7P1` 10/10 — admissible | `CYP2B7P1` 10/10 — admissible | a pseudogene of the same family |
| CYP2D6 | `CYP2D17` 10/10 — refused | *no stable replacement*: `F11` 2/10 | a macaque CYP2D name; nothing admissible takes its place |
| CYP3A4 (midazolam) | `LOC285626` 6/10, `CYP3A64` 4/10 | `LOC285626` 10/10 — admissible | `IL12B-AS1`, an antisense RNA of a pro-inflammatory cytokine |
| CYP3A4 (testosterone) | `LOC285626` 10/10 | `LOC285626` 10/10 | the same |
| CYP2C19 | `APOM` 9/10 — admissible | `APOM` 9/10 | a different gene, stable across folds |
| CYP2C9 | `AMT` 6/10 — admissible | `AMT` 6/10 | a different gene, semi-stable |
| CYP1A2 | `LAMB3` 4/10, `PON3` 3/10, `F5` 1/10 | **not assessable** | fold-specific; no admissible target feature under the pool |
| CYP2C8 | `CYP3A4` 4/10, `SLC27A5` 3/10 | `CYP3A4` 5/10 | fold-specific |
| CYP2E1 | `ACAT2` 3/10, … | `HSD17B4` 3/10 | fold-specific, \|ρ\| ≤ 0.38 |

Three regimes, not one, and they do not mean the same thing:

1. **A plausible duplicate of the destroyed transcript.** Permuting a gene's
   column does not permute a *differently named* feature that measures the same
   molecule. On CYP2D6 that feature is `CYP2D17`, refused by every qualified
   pool, and once it is gone nothing admissible replaces it: the arm drops to
   108.88 against a baseline of 109.64, beside the blind arm's 110.15. On CYP2A6
   and CYP2B6 the duplicate is **admissible** — a human paralogue and a
   pseudogene of the same family — and **no annotation-based qualification can
   expose it**, because establishing cross-hybridisation needs the
   oligonucleotide sequence, which no source distributes. That limit is
   structural.
2. **A different gene with a stable association.** On CYP2C19 the replacement is
   `APOM` in 9 folds of 10; on both CYP3A4 substrates it is `IL12B-AS1` in 10 of
   10, at |ρ| ≈ 0.6 against midazolam hydroxylation. The suppression of CYP3A4
   by inflammation is documented in the literature, which makes the observation
   coherent; **this analysis establishes no mechanism** and the finding is
   hypothesis-generating only. It is not an annotation artefact.
3. **A fold-specific pick at noise level.** On CYP1A2, CYP2C8 and CYP2E1 the arm
   finds a different weak association in almost every fold.

What follows, kept apart:

- **The detection claim stands.** In 100 folds across ten endpoints the
  univariate arm never selected the dead gene. Reading the endpoint lets a
  selector notice that a column has stopped carrying information.
- **The robustness reading is withdrawn.** The number is compatible with all
  three regimes, so on its own it supports none of them. Any earlier statement
  in this project that the diagnostic exposes the mechanistic arm's main
  weakness was reporting a mixture of three phenomena as one.
- **A qualified pool is necessary and not sufficient.** It removes the
  duplicates annotation can see, and only those. An intermediate version of this
  document claimed the pool repairs the diagnostic; that was a generalisation
  from CYP2D6 alone and is withdrawn — the same error, from the same cause, as
  the inflated conflict count: one case taken for the population.
- **The diagnostic is reportable only with its replacement named**, which
  `scripts/name_destroyed_signal_replacement.py` does and the two committed
  artefacts record. It changes no arm, loss, margin or contrast.

**Consequence for H2d.** Where no admissible feature carries the target's
signal — CYP2D6 under pool A is exactly that case — this experiment cannot
discriminate between arms: all of them sit at the baseline, so a repair
criterion requiring an interval to exclude zero fails for want of anything to
find, not for want of a method that could find it.

This diagnostic is what generated the H2d hybrid hypothesis, on one endpoint,
which is why that endpoint cannot confirm it — and the duplicate behind its
original reading has since been refused.

---

## 4. The same comparison on analytically qualified universes

*(§4 is completed as the pool-restricted runs finish; see §5 for why they exist.
Pool definitions and counts are in protocol §5.6 and in
`registry/QUALIFIED_PROBE_POOLS.json`.)*

### 4.1 Restricting every arm to the genes the network can reach

Run before the probe qualification, on the unrestricted annotation: 11,182
candidates instead of 31,862, every arm selecting from the same pool.

**51 of 80 non-inferiority cells**, against 45 in the primary run, and
superiority against all three comparators falls from 1 cell to **0**. The
restriction moves the comparison in both directions at once — it removes
unannotated features a statistical arm could exploit, and it removes candidates
from every arm — so neither number is the honest answer alone. Both are
reported.

### 4.2 Pool A, the qualified human features

21,634 candidate symbols and 26,066 admissible features, every arm selecting
from the same pool, under rule `probe-identity-1.1.0`.

**47 of 80 non-inferiority cells**, against 45 unrestricted and 51 restricted to
the network alone. The ranking of endpoints barely moves: non-inferiority at
every panel size on CYP2B6 and CYP2D6, at none on CYP2C19.

| endpoint | NI cells (of 8) | NetPK best MAE | univariate best MAE | baseline |
|---|---:|---:|---:|---:|
| CYP1A2 *(no admissible target)* | 1 | 434.61 | 419.55 | 492.43 |
| CYP2A6 | 7 | 572.28 | 572.28 | 645.41 |
| CYP2B6 | **8** | 236.75 | 237.99 | 270.37 |
| CYP2C8 | 7 | 51.13 | 50.85 | 62.66 |
| CYP2C9 | 1 | 237.44 | 232.47 | 273.65 |
| CYP2C19 | 0 | 48.31 | 45.96 | 55.52 |
| CYP2D6 | **8** | 91.43 | 91.43 | 109.64 |
| CYP2E1 | 7 | 1330.46 | 1302.16 | 1337.24 |
| CYP3A4 (midazolam) | 3 | 726.61 | 714.52 | 872.13 |
| CYP3A4 (testosterone) | 5 | 3181.27 | 3133.52 | 3628.37 |

**A second wrong-mechanism control goes void under the pool, and the reason is
worth recording.** Unrestricted, 9 of 10 controls are valid and 1 is void
because the frozen rotation sends CYP3A4 (midazolam) to the other CYP3A4
substrate. Under pool A, **8 are valid and 2 void**: the rotation sends CYP3A4
(testosterone) to CYP1A2, which no longer has an admissible feature, so the
wrong-seeded walk falls back on a neighbour that turns out to be CYP3A4 itself —
the target arm's own gene. The control then compares an arm against itself and
says nothing. It is declared by the merge script rather than noticed by eye in a
table of near-identical numbers, and the frozen rotation is not re-run to repair
it.

**CYP1A2's target-only arm equals the baseline exactly** (492.43), because under
the pool it has no feature to select: the arm fits age and sex alone. Its row is
kept for completeness and excluded from any primary tally.

### 4.3 Pool B: the comparison at identical candidates

11,030 candidate symbols and 14,230 features — the intersection of the
analytically admissible features with the genes the public network contains. It
is the only universe that is **both** analytically admissible and reachable by
the ordering — the network-only run of §4.1 is reachable but admits features
whose identity is unestablished — so it is the fair test of **ordering quality**
rather than of pool size.

**49 of 80 non-inferiority cells.** And the answer to the question the pool was
built to settle is negative:

| universe | candidates | NetPK better than the univariate filter | worse | better than LASSO entry order | better than XGBoost gain | NI cells | superior to all three |
|---|---:|---:|---:|---:|---:|---:|---:|
| unrestricted | 31,862 | **1**/80 | 24/80 | 31/80 | 27/80 | 45/80 | 1 |
| network only | 11,182 | 1/80 | 21/80 | 31/80 | 22/80 | 51/80 | 0 |
| pool A | 21,634 | **0**/80 | 21/80 | 31/80 | 20/80 | 47/80 | 0 |
| **pool B** | 11,030 | **0**/80 | 21/80 | 31/80 | 24/80 | 49/80 | 0 |

Counts of cells with a 95 % interval excluding zero, descriptive and
unadjusted, read under the rule stated in §2.

**Restricting the universe narrows the gap and never creates an advantage.**
Non-inferiority rises from 45 to 49 cells, and NetPK remains at **zero** cells
better than the univariate filter against 21 worse. Where the mechanistic
ordering does have a real and strikingly stable advantage is against the other
two statistical selectors: **31 of 80 cells against LASSO entry order, the same
count in all four universes**, and 20 to 27 against XGBoost gain. The honest
summary is therefore that mechanism-guided ordering beats two commonly used
data-driven selectors and does not beat the simplest one.

### 4.4 Stability of the verdict across the four universes

| endpoint | unrestricted | network only | pool A | pool B |
|---|---:|---:|---:|---:|
| CYP1A2 | 0 | 0 | 1 | 1 |
| CYP2A6 | 5 | 7 | 7 | 7 |
| CYP2B6 | **8** | **8** | **8** | **8** |
| CYP2C8 | 8 | 7 | 7 | 6 |
| CYP2C9 | 1 | 0 | 1 | 1 |
| CYP2C19 | **0** | **0** | **0** | **0** |
| CYP2D6 | **8** | **8** | **8** | **8** |
| CYP2E1 | 7 | 8 | 7 | 8 |
| CYP3A4 (midazolam) | 3 | 6 | 3 | 6 |
| CYP3A4 (testosterone) | 5 | 7 | 5 | 4 |

**The endpoint-level verdict is stable on eight endpoints of ten.** CYP2B6 and
CYP2D6 hold at 8 of 8 and CYP2C19 at 0 of 8 in every universe; CYP1A2, CYP2A6,
CYP2C8, CYP2C9 and CYP2E1 move by at most one cell once any restriction is
applied.

**The two CYP3A4 substrates are where the universe matters, and they disagree
with each other.** Under pool B midazolam rises to 6 cells while testosterone
falls to 4. Two measurements of one enzyme, correlated at 0.841 on these donors,
moving in opposite directions under the same restriction: that is evidence that
these two estimates are unstable, not that the restriction has an effect on
them. No conclusion is drawn from either.

### 4.5 What the destroyed-signal numbers do across universes

The univariate arm's error on a wrecked matrix changes with the universe exactly
as the identity of its fallback predicts (§3.1), which is the clearest
confirmation that the diagnostic measures what feature survives rather than
what method is robust:

| endpoint | unrestricted | network only | pool A | pool B | the fallback, and why it moves |
|---|---:|---:|---:|---:|---|
| CYP2D6 | 96.86 | 108.20 | 108.88 | 108.20 | `CYP2D17` is refused by both restrictions; nothing admissible replaces it |
| CYP3A4 (midazolam) | 720.88 | 807.18 | 714.26 | 735.85 | `IL12B-AS1` is admissible but **absent from the network**, so pool A keeps it and the network restrictions do not |
| CYP2B6 | 244.33 | 257.29 | 244.33 | 255.82 | `CYP2B7P1`, a pseudogene: admissible, not in the network |
| CYP2A6 | 571.97 | 577.80 | 577.80 | 577.80 | the macaque names go, `CYP2A7` takes over, and it is in both |
| CYP2C19 | 49.07 | 49.07 | 49.07 | 49.07 | `APOM`, admissible and on the network: unaffected |

The blind mechanistic arm's error is **identical in all four universes** on every
endpoint, which is expected: it selects the target transcript, and no
restriction removes it except on CYP1A2.

---

## 5. What the measurements are, before what they predict

Three questions are kept apart, and a variable can pass the third while failing
the first:

| level | question | finding |
|---|---|---|
| analytical quality | is the announced biomarker the molecule actually measured? | about a third of candidate symbols carry no feature whose human identity can be established |
| biological validity | does the measured molecule belong to the mechanism? | the ordering is produced by a public network, not asserted |
| predictive usefulness | does it improve prediction out of fold? | §2 and §4 |

**CYP1A2 is withdrawn.** All four of its features are `conflicting`: neither
accession resolves to a human gene and the platform record assigns a non-human
one. The defensible statement is not that the probes are macaque — a probe
designed from a non-human sequence can hybridise the human orthologue, and no
source distributes the oligonucleotide sequence — but that **the identity cannot
be established from any public source**. Its activity measurement remains usable
as a secondary endpoint for mechanisms other than its own transcript.

**CYP2E1 is mixed and keeps what the evidence allows.** Two of its three
features are concordant and one is conflicting; the gene is not dropped whole,
because excluding it would discard a real measurement on account of a different,
badly annotated feature.

**Neither CYP1A2 nor CYP2E1 refutes the biological mechanism.** What they refute
is the assumption that the measured transcript is a sufficient predictive
biomarker of the enzyme's activity in these data. A pseudogene or an
undetermined feature is not a false variable; an undetermined identity is not a
wrong measurement, only one whose biological meaning the available annotation
does not pin down.

---

## 6. What is not claimed

- No clinical pharmacokinetic claim. The endpoint is microsomal activity.
- No frontier in `p/n`. A single realisation cannot establish one.
- No autonomous mechanism discovery. The enzyme is given by the endpoint's name;
  the substrate-seeded variant reaches only the two endpoints whose substrate the
  source documents.
- No claim that the mechanistic arm's zero decoy retention is a merit: a decoy
  carries the symbol of a randomly drawn network gene, almost always far from the
  seed, so it was never a candidate for the top of the ordering.
- No stability claim from NetPK's Jaccard of 1.0 at fixed `k`: that is a
  structural consequence of never reading the outcome.
- The two CYP3A4 endpoints are two substrates of one enzyme, correlated at 0.841
  on these donors. They count as **one** endpoint in any tally.
