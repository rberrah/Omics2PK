# F-series results — signed regulatory activities do not add to the target transcript

**Date:** 2026-11-09
**Prespecification:** `docs/OMICS2PK_F_SERIES_FUNCTIONAL_OMICS_PRESPECIFICATION_20261109.md` v1.0.0, frozen at commit `dc37476` **before** this run
**Artefacts:** `reports/data/f_series_hlc_functional_omics.json`, `reports/data/f_series_activity_collinearity.json`
**Network:** CollecTRI via OmniPath, SHA-256 `062529d997ef6ff2d599030f5d155810fb076d3c68d482366172b148553fb878`
**Status: `exploratory`**, as frozen, because this cohort has been inspected repeatedly by this project.

**`clinical_pk_transfer: NOT_EVALUATED`.** This cohort has no dose, no
concentration and no time. Nothing here was scaled to a clearance or an
exposure, and no simulated concentration appears in any F-series artefact.

---

## 1. The answer

**No. The primary contrast fails, and it fails in the strong form.**

| CYP3A4/midazolam activity, 397 donors | MAE | 10 % trimmed MAE |
|---|---|---|
| `F0` covariates only | 876.74 | — |
| **`F1` + target transcript (comparator)** | **728.82** | **518.25** |
| `F2` + regulator transcripts | 761.75 | — |
| `F4p` activities, no target transcript | 827.56 | — |
| **`F5` + activities (primary arm)** | **738.70** | **525.16** |

`F5` versus `F1`: **+9.88 MAE, a 1.36 % *degradation*.**

- 95 % paired donor bootstrap of the difference: **[−9.03, +28.23]**, includes
  zero.
- On the relative scale: **[−1.21 %, +3.99 %]**.
- Reversal on the trimmed MAE co-criterion: **yes** (525.16 against 518.25).

The relative interval is the sentence that matters. The frozen success
criterion asked for a reduction of at least 10 %, and **the interval excludes
−10 %**. This is therefore not the weak statement "no gain was detected"; the
data are **incompatible with a gain of the prespecified size**. At 397 donors
the experiment had the resolution to find a 10 % effect and did not.

The same-enzyme replicate agrees and is not a second validation: on
`CYP3A4T_activity` (testosterone, same microsomes, same donors) `F5` reaches
3138.89 against `F1`'s 3190.66, a 1.62 % improvement with interval
[−4.78 %, +1.88 %] — the opposite sign to the primary and equally far from the
threshold.

---

## 2. The two key secondaries, which say why

### 2.1 The functional layer cannot replace the transcript

`F4p` versus `F1` on the primary endpoint: 827.56 against 728.82, **+13.55 %
worse, interval [+7.86 %, +19.59 %] excluding zero.**

Transcription-factor activities alone are **significantly worse** than the
enzyme's own transcript. That is a clear, interval-backed result in the
unfavourable direction, and it holds on three further endpoints: CYP2D6
+17.08 % [+7.23, +27.33], CYP2A6 +11.38 % [+5.42, +19.04], CYP2B6 +8.57 %
[+2.16, +15.94], all excluding zero.

### 2.2 Weighting by the network does not beat measuring the regulators

`F5` versus `F2`: 738.70 against 761.75, −3.03 % with interval
**[−7.16 %, +0.85 %]**, including zero.

`F5` carries 49 activity columns, `F2` 38 transcript columns from the same
regulators. The regulon-weighted aggregation is not distinguishable from
simply putting the regulators' own transcripts in the model. No other endpoint
reaches an interval excluding zero on this contrast.

---

## 3. The eight other endpoints: no signal anywhere

`F5` versus `F1`, Bonferroni-adjusted over the family of seven evaluable
endpoints:

| Endpoint | `F1` | `F5` | relative | adjusted interval |
|---|---|---|---|---|
| CYP2C9 | 248.80 | 237.82 | **−4.41 %** | [−9.06, +0.51] |
| CYP2C8 | 53.18 | 52.68 | −0.94 % | [−5.04, +3.45] |
| CYP2E1 | 1349.81 | 1340.27 | −0.71 % | [−2.40, +0.75] |
| CYP2C19 | 51.29 | 51.09 | −0.38 % | [−4.05, +3.40] |
| CYP2D6 | 91.19 | 92.70 | +1.66 % | [+0.05, +3.32] |
| CYP2B6 | 240.89 | 246.14 | +2.18 % | [−2.63, +7.56] |
| CYP2A6 | 571.06 | 588.03 | +2.97 % | [−0.60, +6.61] |
| CYP1A2 | — | — | `NOT_EVALUABLE` | no admissible target probe |

Four endpoints marginally better, three marginally worse, **none within a
factor of two of the 10 % threshold**, and the single interval excluding zero
(CYP2D6) does so on the *wrong* side. The best result anywhere is CYP2C9's
−4.41 %, whose interval still contains zero.

`CYP1A2_activity` was declared `NOT_EVALUABLE` for `F1`, `F2` and `F5` in the
prespecification, before the run, because none of its four probes is
admissible. It is reported as itself, not as a negative result.

---

## 4. The control that refutes the premise

This is the most informative part of the run, and it is the one the
prespecification was rewritten to make possible.

| Setting on CYP3A4M | `F4p` MAE | `F5` MAE |
|---|---|---|
| **real CollecTRI network** | **827.56** | **738.70** |
| size-and-sign-preserving rewired null, 20 replicates | 830.38 mean, **[815.68, 845.50]** | 769.16 mean, [745.02, 786.45] |
| partial sign flip, 25 % | 821.91 | 745.56 |
| partial sign flip, 50 % | 823.66 | 764.76 |
| partial sign flip, 75 % | 828.19 | 748.05 |
| donor-permuted expression, 20 replicates | 882.08 [873.45, 891.30] | 752.09 [737.08, 760.54] |
| `F0` floor, for reference | 876.74 | — |

Read in order:

1. **The activities do carry donor-specific information.** `F4p` at 827.56 is
   far below its donor-permuted control at 882.08 and below the covariate floor
   at 876.74. Destroying the link between expression and donor costs 55 MAE, so
   the scores are not noise.
2. **But that information is not network-specific.** A null network with each
   regulon's *exact size* and *exact count of activating and repressing edges*,
   differing only in **which genes it names**, reaches 830.38 on average and
   spans [815.68, 845.50]. **The real network's 827.56 lies inside that
   range.** CollecTRI's gene membership does no better than a size-matched
   random draw from the same 21,580-gene universe.
3. **The signs do not order the result.** Flipping a quarter, a half or three
   quarters of each regulon's edges gives `F4p` 821.91, 823.66 and 828.19 — all
   within noise of the real 827.56, and the 25 % flip is nominally *better*.
   There is no monotone degradation with perturbation, which is what a network
   carrying sign information would have to show.
4. **One asymmetry, reported because it is real and does not rescue anything.**
   For `F5` the real network's 738.70 is better than all twenty rewired nulls
   (minimum 745.02). So when the target transcript is present, the genuine gene
   membership does produce slightly better columns than a random one — by about
   6 MAE on 739, under 1 %. That is a measurable network effect and it is too
   small to turn `F5` into an improvement over `F1`. It is reported as an
   observation, not as a success.

The prespecification fixed in advance how to read this: *"If the primary
contrast shows a gain that the null networks reproduce, the gain is reported as
refuted by its own control."* There is no gain to refute; the controls instead
explain the absence of one. The functional layer behaves like a weighted
average of two hundred liver transcripts, which is what it mechanically is.

### 4.1 The algebraic identity holds exactly

The full sign inversion, run at perturbation fraction 1.0 and labelled
`algebraic_identity_check` rather than `control`: **maximum absolute prediction
difference 0.0 across all ten folds**, against a tolerance of 1 × 10⁻⁸.

Inverting every edge sign negates every activity and ridge answers by negating
one coefficient, so the predictions are identical. Had this been reported as a
control, it would have read as "the result is robust to destroying the
network's signs" — a true sentence about an algebraic identity and a false one
about biology. It is reported here only as evidence that the implementation is
the model the prespecification describes.

---

## 5. Biological relevance: a significant number that means nothing about biology

The frozen reporting rule keeps this block separate, and it has to be read with
that separation enforced.

**Measured:** of the 49 evaluable CYP3A4 regulators, the fitted ridge
coefficient's sign agrees with CollecTRI's recorded sign for the
regulator-to-enzyme edge in **15 cases, 30.6 %**, two-sided exact binomial
**p = 0.0094** — significantly *below* chance.

Two readings are possible: the resource's signs are wrong, or the coefficient
signs are not identified. **The second is measurable and it is the right one.**
A post-hoc diagnostic (`scripts/diagnose_f_series_activity_collinearity.py`,
labelled `POST_HOC_DIAGNOSTIC`, changing no arm, loss, threshold or contrast):

| | |
|---|---|
| maximum absolute pairwise correlation among the 49 activity columns | **1.0000** |
| median absolute pairwise correlation | 0.424 |
| pairs above \|r\| = 0.9 | 9.0 % |
| condition number of the standardised design | **1.66 × 10³³** |
| effective rank (Shannon) of 49 columns | **7.16** |
| variance in the first component | 50.3 % |
| components for 95 % of the variance | 19 |
| maximum regulon Jaccard overlap | **1.000** |

At least one pair of regulons is **identical** after the sign and cytochrome
rules, so at least one pair of activity columns is a perfect duplicate; the
design is numerically singular. Ridge distributes weight among near-duplicate
columns and the sign of any individual one is arbitrary.

**The 30.6 % is therefore a statement about the design, not about CollecTRI.**
It is reported because the prespecification required it, and it is explicitly
**not** evidence that the resource's signs are wrong. The transferable
methodological point is the one worth keeping: *reading decoupleR-style
activity scores through individual regression coefficients is invalid at this
collinearity*, and anyone proposing to interpret such coefficients should
measure the effective rank first.

This also explains §1 mechanically. Forty-nine columns carrying about seven
effective dimensions, half of them in a single component, add almost nothing to
a model that already has the target transcript.

---

## 6. Analytical validity

| | |
|---|---|
| Network rows read | 64,516, all naming CollecTRI, none undirected |
| Signed edges retained after the frozen rules | **58,991** |
| — dropped sign-indeterminate | 2,590 |
| — dropped, cytochrome target (**guard 1**) | 811 |
| — duplicate rows collapsed | 2,004 |
| — dropped, duplicate resolving to opposite signs | 120 |
| Admissible probes refused because the array symbol is a cytochrome (**guard 2**) | 100 |
| Admissible probes refused because the *resolved* gene is a cytochrome (**guard 3**) | **4** |
| Scoring universe | 21,580 symbols; 21,549 with a usable probe in fold 0 |
| Cytochrome symbols held separately as predictor columns, never in a score | 51 |
| Donors analysed | 397 (CYP3A4M/T, CYP2A6, CYP2D6), 398, or 387 (CYP2C8) |
| Folds | 5 outer × 2 repeats × 3 inner, seed 20261109 |
| Fold fingerprint | `0ab8377a8db3959e718bb55cf937c761d3500285ca1819587d426b0b992d7bb1` |
| Activity settings computed per fold | 102 |
| Estimator | H2a's own ridge, imported by path, its constants asserted before the run |
| Individual-level data written | none |

**Determinism, verified rather than assumed.** The experiment was executed four
times: once before the contrast-guard fix of §7, and three times after it. Each
consecutive pair of post-fix artefacts is **identical in every field except the
timestamp, the timing fields and — between the last two — the recorded
checksum of the prespecification**, which changed because the change-log row
1.0.1 was added to that document in between. Nothing in the pipeline depends on
an unseeded source of randomness, and the diff was computed rather than
assumed: of the whole report, exactly one key differs, and it is a checksum of
a file, not a result.

That checksum deserves one sentence of its own, because read carelessly it
looks circular. The committed artefact carries the digest of prespecification
**v1.0.1**, whose change log describes defects found *after* the run. The
**frozen rules are v1.0.0 at commit `dc37476`**, made before any model was
fitted, and that commit is where the freeze is verifiable. v1.0.1 adds a
change-log row and changes no rule.

**A timing figure that would have misled.** One of these executions reported
`runtime_seconds` 6657.8. Folds one to nine took 213 s and fold ten took
6,436 s for identical work: the process was suspended, not computing. The
artefact now carries `fold_cumulative_seconds` per fold alongside the wall
clock, and the committed run's honest cost is **258.4 s total, 252.8 s in the
fold loop**. No result depends on it.

The per-fold timings are written as a mapping rather than a list, and the
reason is worth recording: this repository has a guard that refuses any numeric
array of ten or more entries in a committed artefact, as a crude proxy for
participant-level data. Ten fold timings tripped it. The guard is deliberately
crude and is worth more than the convenience of a list, so the artefact's shape
changed rather than the guard's threshold.

---

## 7. Two defects in this experiment's own execution, recorded

Both were in the reporting path, neither touched a fitted value, and both are
recorded because finding them after the freeze is exactly the situation the
freeze exists to make visible.

1. **The first execution reported the primary contrast as `NOT_EVALUABLE`.**
   Every arm had been fitted — the loss table was complete — but the guard in
   `contrast()` looked for the key `"real"` *inside* the dictionary of arms
   instead of the dictionary of settings, so the comparator was never found.
   The bug was silent: it produced a well-formed refusal rather than an error.
   Fixed, seven tests added (`TestContrastAssembly`), and the run repeated. No
   fitted value changed, because none of them depended on the guard, and this
   was checked rather than assumed: **all 47 arm losses across the ten
   endpoints agree to two decimal places before and after the fix.** The
   pre-fix artefact was overwritten by the repeat, so the comparison is at the
   precision at which those losses were read, not bit-for-bit.
2. **The freeze anticipates a contrast it never defined.** §4.5 states that
   CYP1A2 "is evaluable for `F4p` against `F0` alone and is reported there,
   separately", but §7's contrast list contains no `F4p`-versus-`F0` entry, so
   no interval was computed for it. The loss table gives `F4p` 465.26 against
   `F0` 492.23 on CYP1A2, a 5.5 % reduction, and that is **reported
   descriptively without an interval**: computing one now would be a post-hoc
   test presented as a planned one. The gap is recorded here rather than
   repaired in place.

---

## 8. What this does and does not establish

**Establishes, on 397 donors, one cohort, exploratory and not upgradeable:**

- Adding signed CollecTRI transcription-factor activities to a liver CYP3A4
  transcript **does not produce the gain the prespecification required**. The
  point estimate is a 1.36 % degradation and the data **exclude** an
  improvement of the prespecified 10 %. **Amended in 1.1.0:** they do *not*
  exclude a smaller benefit — the relative interval is [−1.21 %, +3.99 %], so a
  reduction of up to about 1.2 % remains compatible with this run. The earlier
  wording, "does not reduce the error", overstated that; what is established is
  the absence of a gain of the specified size, not the absence of any gain.
- Those activities are **significantly worse** than the target transcript when
  used in its place (+13.55 %, interval excluding zero), and on three further
  enzymes.
- Regulon-weighted aggregation is **not distinguishable** from putting the same
  regulators' transcripts in the model.
- **No network-specific advantage was detected** for the activity-only arm: a
  size-and-sign-matched random network scores 830.38 on average over 20
  replicates, range [815.68, 845.50], and the real network's 827.56 lies inside
  that range. **Amended in 1.1.0:** the earlier wording inferred from this that
  "the information those scores carry is not specific to the network's gene
  membership". Twenty replicates do not support that inference. They are an
  **exploratory** control whose finest one-sided rank statistic is 1/21 ≈ 0.048,
  so a network-specific effect smaller than that resolution cannot be excluded;
  and the `F5` arm points the other way, where the real network beats all twenty
  nulls (rank p ≤ 1/21, nominally significant, uncorrected, at the resolution
  floor). The defensible statement is "not detected at this resolution", not
  "absent".
- The scores are nonetheless **not noise**: destroying the expression-to-donor
  link costs 55 MAE.
- **Individual coefficient signs on these columns are not interpretable**:
  effective rank 7.16 of 49, condition number 1.66 × 10³³, at least one
  duplicated regulon.

**Does not establish:**

- That CollecTRI's signs are wrong. The coefficient-sign result is a property
  of a singular design; this experiment measured no sign independently.
- That a functional layer cannot help in general. It tested **one** resource,
  **one** scoring statistic, **one** tissue, **post-mortem**, against an **in
  vitro** endpoint, with an **mRNA** proxy for transcription-factor protein
  activity. Any of those five could be the binding constraint.
- Anything about clinical pharmacokinetics. The endpoint is
  `functional_in_vitro`. Closing that chain needs individual transcriptome
  **and** individual concentration-time pharmacokinetics in the same publicly
  downloadable subjects; the project's registry contains no such dataset, and
  that step remains `ABSTAIN` with the data requirement named rather than
  scheduled.

---

## 9. What follows for the architecture this came from

The proposal that prompted this work placed a signed-network layer between
omics and pharmacokinetic parameters, with CYP3A4 on the Human Liver Cohort as
its demonstrator. The demonstrator has now been run under a freeze, and the
layer did not earn its place at this step.

Three specific consequences, stated so they can be argued with:

1. **Do not build the remaining layers on this evidence.** SIGNOR, PROGENy,
   COSMOS and the scaling step were declared out of scope in the
   prespecification for reasons measured beforehand (SIGNOR has zero
   inflammation-to-CYP edges; PROGENy has no xenobiotic footprint). Nothing in
   this result weakens those reasons.
2. **The binding constraint is probably not the network.** The activity scores
   carry donor-specific information that a random network reproduces, which is
   the signature of a representation too coarse to be enzyme-specific, not of a
   wrong graph. A layer that could help would have to be enzyme-specific by
   construction — which the target transcript already is, and which is why `F1`
   wins.
3. **The real gap is the data, not the method.** The one thing that would change
   what is answerable here is a public cohort with a transcriptome and clinical
   pharmacokinetics in the same subjects. Until one exists, every functional
   layer tested on HLC is tested against an in vitro endpoint and cannot speak
   to exposure.

A negative result under a published freeze is a result. This one is reported
with its controls, its defects and its limits, and it is not upgradeable.

---

## 10. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-09 | First report of the frozen run (commit `d491152`). |
| 1.1.0 | 2026-11-09 | **Two claims in §8 narrowed after external review, both because they said more than the measurement supports. No number changed and no conclusion reversed.** (1) "Does not reduce the error" is withdrawn: the relative interval [−1.21 %, +3.99 %] excludes the prespecified 10 % gain but **not** a benefit of up to about 1.2 %. The established statement is the absence of a gain of the specified size. (2) The inference that the activity scores carry no network-specific information is withdrawn. Twenty rewired-null replicates are an exploratory control with a finest one-sided rank statistic of 1/21 ≈ 0.048; an effect below that resolution cannot be excluded, and the `F5` arm points the other way (real better than all twenty nulls, rank p ≤ 1/21, uncorrected). "Not detected at this resolution" replaces "not specific". Both corrections tighten the negative conclusion rather than weakening it: the primary contrast still fails, and §9's consequences are unaffected. |
