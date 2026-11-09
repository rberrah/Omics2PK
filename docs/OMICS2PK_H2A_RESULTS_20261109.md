# H2a — mechanistic versus statistical selection on the Human Liver Cohort

**Document version:** 1.0.0 (draft; the pool-restricted sections are filled as
their runs complete)
**Date:** 2026-11-09
**Protocol:** `OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md` v1.4.5
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

## 3. The destroyed-signal diagnostic: the robust finding

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
arm was **never once fooled, in 100 folds**. This is the price of being
outcome-blind, measured rather than argued, and it is what generated the H2d
hybrid hypothesis — on one endpoint, which is why that endpoint cannot confirm
it.

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

*(pending)*

### 4.3 Pool B, qualified and reachable

*(pending)*

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
