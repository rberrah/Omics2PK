# H2d — erratum to the prespecification, before any execution

**Erratum version:** 1.0.1 to `OMICS2PK_H2D_HYBRID_SELECTION_PRESPECIFICATION_20261108.md` v1.0.0
**Date:** 2026-11-08
**Status:** written **before the hybrid rule has been run on any cohort endpoint**. No success threshold is changed. Three definitions are made precise, one reporting requirement is added, and one claim about the validation set is weakened by a measurement.

---

## E1 — The convergence measure was a containment ratio, not set agreement

**What was wrong.** The prespecification declared the failure mode as the hybrid "collapsing onto the univariate selector", and the code measured it as

> share of the hybrid's retained genes that the statistical arm also chose.

That is containment, not agreement. Demonstrated on synthetic data: a hybrid retaining a single informative gene, compared against a statistical arm's top five that happens to contain it, scored a perfect **1.0** while the two panels were plainly different sets. Read as convergence, that number would have condemned a method that had in fact selected one gene where the comparator selected five.

**What replaces it.** Three measures, reported together at every panel size, because they answer different questions:

| Measure | Question |
|---|---|
| **containment** | Is the hybrid's panel a *subset* of the statistical panel? High containment with a small hybrid panel is a **parsimony** claim, not a collapse. |
| **Jaccard** | Graded set agreement; falls when either panel holds genes the other does not. |
| **exact set equality** | Are the panels *the same*? This is the one the failure mode is about. |

**The 90 % threshold is unchanged in value and is now attached to its real definition:** the collapse verdict is reached when the hybrid's retained set is **exactly equal** to the statistical arm's top-`k` in at least 90 % of folds.

**One qualification, prespecified here rather than discovered later.** At `k = 1` the H2a run already showed every arm choosing the same single transcript, so exact equality at `k = 1` is expected and carries no information about collapse. The collapse verdict is therefore assessed on **`k ≥ 5`**. All three measures are still reported at every `k`, including `k = 1`, so nothing is hidden by this restriction.

**An abstention has no panel to compare.** Every measure is undefined for it, not zero. Scoring an abstention as disagreement would flatter the method.

---

## E2 — Equal dimension is no longer guaranteed, so three numbers are required

**What was wrong.** H2a could compare arms at identical dimension because every arm retained exactly `k`. The hybrid may retain fewer, or nothing. Comparing it against statistical arms that always retain `k` would credit selection quality for what is partly dimension reduction — and in the direction that flatters the hybrid, since a smaller panel is often the better model here: H2a's no-selection arm, carrying every candidate gene, was worse than age and sex alone.

**What is required.** Every hybrid result reports, jointly and never one without the others:

1. **performance at budget `k`** — the hybrid against comparators each holding `k` genes;
2. **the number actually retained**, per fold: its distribution, its median, and the abstention count;
3. **performance at matched dimension** — for each fold, after the hybrid has retained `m ≤ k` genes, the comparators are re-evaluated at their own top-`m`, so both sides hold `m`.

Item 3 is the contrast that isolates selection quality. Items 1 and 2 are what make item 3 interpretable. A favourable reading of item 1 alone is not reportable as a selection result.

---

## E3 — Calibration, extended across regimes

**What was wrong.** The prespecification recorded a single null design, 12 of 200 draws against a nominal 0.05, and that is compatible with correct calibration without demonstrating it.

**What was done.** Four further synthetic regimes, none touching the cohort, each 120 draws except the last at 40, budget 2 except the last at 1:

| Regime | Selection fired | Mean retained |
|---|---:|---:|
| pure noise, 6 null candidates | 0.083 | 0.09 |
| correlated nulls, pairwise *r* ≈ 0.8 | 0.025 | 0.03 |
| endpoint driven entirely by the baseline, with candidate echoes of age and sex | 0.075 | 0.10 |
| one weak real signal, partial *r* ≈ 0.4 | **1.000** | 1.11 |
| one strong real signal | **1.000** | 1.00 |

Reading: specificity sits at or slightly above the nominal rate in all three null regimes, correlated candidates do not inflate it, a candidate that merely echoes a baseline covariate is not mistaken for information, and sensitivity is complete at a partial correlation of 0.4. The rule is calibrated on these designs. It is **not** thereby calibrated on real transcriptome dependence structure, which is far richer than a single shared factor, and no claim of general FDR control is made.

---

## E4 — The nine remaining endpoints are not nine independent replications

**What was overstated.** The prespecification called the other nine endpoints "the confirmatory set", which is correct about *design independence* — they did not generate the hypothesis — and misleading about *statistical independence*. They are measured in the same livers and the enzymes are co-regulated.

**Measured on the donors who carry expression data, before any hybrid run:**

- 45 endpoint pairs, Spearman **median 0.446**, mean 0.425, range 0.144 to 0.841; 13 pairs above 0.5.
- The first principal component of the endpoint correlation matrix explains **49.6 %** of the variance.
- Effective number of independent endpoints, from the eigenvalue spread: **3.55 of 10**.
- `CYP3A4M_activity` and `CYP3A4T_activity` correlate at **0.841**. They are two substrates of one enzyme and **count as one endpoint**, never two, in any tally of how many endpoints agreed.

**Consequences, frozen:**

1. No text may describe an H2d result as "confirmed on nine endpoints". The admissible form is: *nine endpoint measurements from a single cohort, of effective independent dimensionality about three and a half, none of which was used to design the rule.*
2. Agreement across endpoints is reported as a count **and** with this effective dimensionality beside it.
3. This is **validation on unused endpoints, not independent replication.** Independent replication would need another cohort, and §4.1 of the parent protocol already records that no admissible public cohort supplies one for this axis.

---

## What is unchanged

Every success and failure threshold of v1.0.0: the 5 % relative margin, the `q = 0.05` retention rate, the oversampling factor of 3, the 200 permutations, the seed, the parsimony requirement, the abstention requirement, the `exploratory` cap on the generating endpoint, and the statement that lowest MAE is not a criterion. Nothing in this erratum makes a criterion easier to satisfy; E1 and E2 each make a claim harder to make, and E4 weakens what the result may be called.

No H2a number is recomputed, reinterpreted or re-frozen.
