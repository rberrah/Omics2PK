# H2d — hybrid mechanistic selection: prespecification

**Version:** 1.0.0
**Written:** 2026-11-08
**Status:** frozen before any hybrid selection has been computed on any endpoint.
**Parent protocol:** `OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md` v1.4.4, §5.3.

---

## 0. Where this hypothesis came from, stated first

**This rule was designed after seeing a result, and that result was obtained on one endpoint.**

The H2a destroyed-signal diagnostic on `CYP2D6_activity` permuted the target transcript's donor order, leaving the gene in place in the mechanistic ordering but stripping its information. The outcome-blind mechanistic arm went on selecting it and lost its entire advantage — MAE 110.15 against a clinical baseline of 109.64 — while the univariate arm chose a different variable and reached 96.86, and never once picked the dead gene across ten folds.

The hybrid rule below is the response to that observation. It is therefore **hypothesis-generating with respect to `CYP2D6_activity`** and cannot be confirmatory there, whatever it produces. Three consequences, frozen:

1. `CYP2D6_activity` is the **generating endpoint**. Its hybrid result is reported separately, labelled as the case that produced the hypothesis, and **capped at `exploratory`**.
2. The **confirmatory set is the other nine endpoints**, which were running under the frozen H2a protocol when this document was written and whose destroyed-signal diagnostics had not been inspected.
3. No statement in any later document may describe this rule as having been frozen before the observation that motivated it. The sequence was: observe on one endpoint, design, freeze, then evaluate elsewhere.

This is weaker than a genuine prior prespecification and is labelled as such. It is the strongest honest position available, because the alternative — evaluating on the endpoint that generated the idea and calling it prespecified — would be a misrepresentation.

**Naming.** The external review that proposed this strategy called it "H2c". That label was already taken, two protocol revisions earlier, by the mechanistic-integration question. The hybrid strategy is therefore **H2d**, and the review's "H2c" means this document.

---

## 1. The question

Does hybrid selection keep the parsimony of mechanistic selection while shedding its vulnerability to mechanisms that are biologically plausible but carry no information in the data?

This is a sharper question than which of NetPK and a gradient-boosting selector attains the lower MAE, and it is the one on which the method's methodological contribution rests.

---

## 2. The rule, frozen

| Step | Specification |
|---|---|
| 1. Knowledge proposes | The top `oversample × k` genes of the outcome-blind mechanistic ordering of protocol §5.1. Nothing in this step reads the endpoint. |
| 2. Data disposes | For each shortlisted gene, the absolute **partial rank correlation** with the endpoint, after rank-transforming both and residualising each on the baseline covariates. Computed on **training donors only**. |
| 3. Calibration | A null built by permuting the **residualised** endpoint among training donors, so the baseline adjustment survives the permutation. Retention by the permutation false-discovery estimator on counts, the one already used and tested in this project. An empirical p-value is not used: from `B` permutations it cannot fall below `1/(B+1)`, which is the defect recorded in protocol v1.4.0. |
| 4. Keep, demote or abstain | Retain at most `k`. If nothing survives, retain nothing and let the model fall back to the clinical covariates. |

Frozen values: `oversample = 3`, `null_permutations = 200`, `target_fdr = 0.05`, `permutation_seed = 20261108`.

Everything else — folds, donors, predictor, panel-size grid, loss hierarchy, bootstrap, the one-feature-per-gene universe, the missing-value rule — is **identical to H2a** and inherited unchanged, so the only difference between the hybrid arm and the blind mechanistic arm is step 2 onward.

**The oversampling factor is deliberately small.** A wide shortlist would let the data step do all the work.

### 2.1 Calibration of the retention step, measured on synthetic data

Checked before any cohort endpoint was touched, on simulated blocks rather than
on the cohort, so this is a property of the rule and not a result:

- a shortlist of five variables that carry **no** information produces a
  selection in **12 of 200 draws (0.060)** against a nominal 0.05, and never
  retains more than one;
- a variable that merely repeats a baseline covariate scores **0.04**, because
  the adjustment is made in the same rank space as the variables it adjusts.
  Residualising ranks on raw covariates leaves curvature behind and gave that
  same variable **0.20**, a partial correlation it had not earned;
- the permuted endpoint is **re-residualised** inside each permutation. Without
  that step the permutation reintroduces a baseline component to which the
  residualised columns are orthogonal, shrinking the null and making the test
  anti-conservative.

Both of the last two points were defects in the first implementation of this
rule, found by its own tests before any endpoint was read.

---

## 3. What success and failure look like, decided now

### Success requires all four

1. **The stated weakness is repaired.** On the destroyed-signal diagnostic, the hybrid's MAE is within the frozen 5 % relative margin of the univariate arm's, *and* its interval against the blind mechanistic arm excludes zero on the favourable side. Equivalently: the hybrid must demote the dead gene.
2. **No cost at equal dimension.** At every panel size, the hybrid is non-inferior to the blind mechanistic arm within the 5 % margin.
3. **Parsimony is kept.** At matched MAE, the hybrid's median retained size does not exceed the statistical arms'.
4. **Abstention is real.** The hybrid retains nothing on at least one endpoint-and-budget where the mechanism carries no information, rather than always returning its best guess.

### Failure, in the specific form most likely to occur

**The hybrid collapses onto the univariate selector.** If the retained sets coincide with the univariate arm's top-`k` in at least 90 % of folds, then the knowledge contributed nothing beyond shortlisting, and the honest conclusion is that the mechanistic layer is a dimensionality-reduction convenience rather than a source of information. The overlap fraction is reported with every result, at every panel size, so this cannot be quietly omitted.

Two further failure modes: the hybrid keeps the dead gene (weakness not repaired), or it degrades beyond the margin against the blind arm (the data step is destroying real signal).

### What is not a criterion

Beating a statistical selector on MAE. A favourable MAE accompanied by a high overlap with the univariate choice, or by no abstention capability, is a **weaker** result than a modest MAE with low overlap and demonstrated abstention.

---

## 4. Independence of the comparison

The hybrid, the blind mechanistic arm and the three statistical arms are compared on **identical folds**, with every selection, residualisation, permutation and penalty choice confined to training donors. No test fold enters any selection step, and no comparator is chosen by its out-of-fold error — the defect corrected in protocol v1.4.3 applies here too: the 5 % margin must hold against **each** statistical arm separately.

---

## 5. What the hybrid gives up

Its selection varies between folds because it reads the endpoint. The blind ordering's Jaccard stability of 1.0 at fixed `k` was a structural property and never a merit, as protocol §5 records; the hybrid forfeits it and in exchange makes the stability comparison meaningful for the first time. Both facts are reported together.

---

## 6. Endpoint class and ceiling

The endpoint remains measured microsomal activity, `functional_in_vitro`. A result here is the **functional analogue** of H2 and is never pooled with a clinical pharmacokinetic result nor reported in the same effect column. It cannot promote the clinical H2 claim above `PARTIALLY_SUPPORTED`.

---

## 7. Order of execution, frozen

1. H2a completes on all ten endpoints under the protocol frozen before this document.
2. The hybrid is computed on the **nine confirmatory endpoints**.
3. The hybrid is computed on `CYP2D6_activity` and reported separately as the generating case, capped at `exploratory`.
4. H2b (mechanism recovery from the drug) and H2c (mechanistic integration) remain untouched by this document.

No H2a number is recomputed, reinterpreted or re-frozen on account of this document.
