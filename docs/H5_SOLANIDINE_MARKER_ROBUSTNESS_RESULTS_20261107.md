# H5 results — robustness of the solanidine CYP2D6 markers

**Date:** 2026-11-07
**Cohort:** Müller 2026 standardized potato meal, n=20 (doi `10.3389/fphar.2026.1886427`)
**Script:** `scripts/run_h5_solanidine_marker_robustness.py`
**Source checksum:** `68d53071…1950`, verified at run time
**Direction check against:** ARIANA n=47 (doi `10.3389/fphar.2026.1812007`, checksum `7397b13e…fe37`)

**Revised 2026-11-07** after external review. Three claims in the first version were too strong and are corrected here; one of them is **reversed**. The change log is at the end.

## Prespecification, stated exactly

The first version of this document claimed more than is true. Precisely:

**Frozen in protocol v1.0.0 (commit `7632e9e`), before this dataset was known to exist** — the H5 *objectives* and the *criterion families*: "intra-individual coefficients of variation where data permit", "ranking stability", "robustness to matrix and time", "95 % intervals", "quantified versus semi-quantitative identification".

**Not frozen — chosen at analysis time on 2026-11-07** — Spearman rather than Pearson, the intraclass correlation, Kendall's W, and the reference times.

**One choice was made after seeing a result**, and is labelled accordingly: the 3 h reference time was fixed before any run, but **the 29 h reference time was added after the first run showed 29 h to be the informative timepoint**. Statistics anchored on 29 h are therefore post-hoc and exploratory. The association table itself is *not* selective — all nine timepoints are computed and all nine reported, none dropped.

So the accurate formulation is: **objectives and criterion families prespecified, statistical operationalisation added later.** Not "nothing here is tuned on the results".

**This is not pharmacokinetics.** No drug was administered. Nothing here predicts exposure, and nothing here replicates the ARIANA solanidine-to-trough association.

---

## Headline result

**The activity-score association is detected only in the late sample.** The association is significant at 29 h for both markers and in both censoring scenarios; no timepoint in the 0–5 h window reaches significance after correction.

**What this does not establish.** The early-window associations are **not established with the available precision** — that is the whole of what can be said. Two things must not be read into it:

- **The 80 %-power correlation is not an upper bound.** At n≈19 under Bonferroni correction the test would have 80 % power at ρ ≈ 0.76, but a non-significant result does not imply the true correlation lies below that. The first revision of this document made exactly that error.
- **The genuine compatibility bounds are the confidence intervals.** Across 0–5 h the widest upper limit is **+0.69** for OH-solanidine and **+0.79** for SSDA. The data are compatible with early associations up to those values — which **overlap the 29 h estimates** (+0.709 and +0.658). So this design does not distinguish the early window from the late sample on association strength at all.

Claiming otherwise would be the "difference in significance is not significance of difference" error.

Multiplicity is handled over the **primary family of 18 comparisons** (2 markers × 9 timepoints), threshold 0.002778. The three censoring scenarios are sensitivity analyses of the same comparison, not additional hypotheses; the most conservative count (54) is also reported.

| Marker | scenario | ρ at 29 h | p | survives ×18 | survives ×54 |
|---|---|--:|--:|:--:|:--:|
| OH-solanidine | all values | +0.753 | 0.00013 | ✓ | ✓ |
| OH-solanidine | below-LOQ excluded | +0.709 | 0.00068 | ✓ | ✓ |
| SSDA | all values | +0.711 | 0.00044 | ✓ | ✓ |
| SSDA | below-LOQ excluded | +0.658 | 0.00218 | ✓ | ✗ |

Below-LOQ excluded, across all nine timepoints:

| Time | OH-solanidine/solanidine ρ (p) | SSDA/solanidine ρ (p) |
|---|---|---|
| 0 h | −0.07 (0.80) | +0.23 (0.41) |
| 0.5 h | +0.16 (0.58) | +0.36 (0.18) |
| 1 h | +0.22 (0.43) | +0.37 (0.17) |
| 1.5 h | −0.00 (1.00) | +0.34 (0.19) |
| 2 h | +0.27 (0.33) | +0.35 (0.20) |
| 3 h | +0.28 (0.32) | +0.47 (0.074) |
| 4 h | +0.18 (0.48) | +0.46 (0.066) |
| 5 h | −0.14 (0.60) | +0.37 (0.14) |
| **29 h** | **+0.71 (0.0007)** ✓ | **+0.66 (0.0022)** ✓ |

Only the 29 h sample survives correction, for either marker. This is the central H5 finding and it is a **negative result about marker usability**: a sample drawn in the hours after a potato-containing meal is not interpretable as a CYP2D6 phenotype measure.

**Mechanistically coherent, and coherent with ARIANA.** Soon after the meal, absorbed parent solanidine dominates the denominator, so the ratio reflects absorption rather than metabolism. By 29 h the system has moved toward a metabolism-determined ratio. ARIANA sampled at trough, far from a meal, which is the condition resembling 29 h — consistent with its strong published correlations, and a reason not to generalise them to arbitrary sampling times.

### Reproducible but not valid — the subtle part

Early timepoints agree strongly **with each other** while agreeing poorly with the informative late sample. Spearman agreement of each time's subject ranking with the 3 h versus the 29 h ranking (below-LOQ excluded):

| Marker | with 3 h ranking | with 29 h ranking |
|---|---|---|
| OH-solanidine/solanidine | +0.46 to +0.90 | +0.00 to +0.55 |
| SSDA/solanidine | +0.75 to +0.96 | +0.43 to +0.75 |

A repeated early sample would look highly reproducible and would still rank patients differently from the one timepoint that tracks genotype. **Internal consistency is not validity**, and a reproducibility study confined to the early window would have produced a falsely reassuring answer.

---

## Reversed: the ICC drop is not cell-level imputation, it is subject inclusion at the floor

The first version of this document attributed the fall from ICC 0.942 to 0.400 to the below-LOQ *imputation*. That attribution was wrong, and the review was right to challenge it. Comparing "all values" with "below-LOQ excluded" confounds three things: the subject set changes, the number of observations changes, and the postprandial time trend sits inside the within-subject term.

Restricting both arms to a **common subject set** and removing the fixed time effect separates them:

| OH-solanidine/solanidine | all values | below-LOQ excluded | Δ |
|---|--:|--:|--:|
| headline, all 20 subjects, raw | 0.942 | 0.400 | **+0.542** |
| common 17 subjects, raw | 0.349 | 0.400 | −0.051 |
| common 17 subjects, time effect removed | 0.393 | 0.430 | **−0.037** |

| SSDA/solanidine | all values | below-LOQ excluded | Δ |
|---|--:|--:|--:|
| headline, all 20 subjects, raw | 0.886 | 0.754 | +0.132 |
| common 18 subjects, raw | 0.705 | 0.754 | −0.049 |
| common 18 subjects, time effect removed | 0.779 | 0.798 | **−0.020** |

**On a comparable population, censoring does not lower the ICC at all — it raises it slightly.** The part attributable to censoring alone is −0.037 for OH-solanidine and −0.020 for SSDA, i.e. in the opposite direction to the headline.

The real mechanism: **three subjects wholly pinned at the imputed floor**, participant 1 among them, inflate the between-subject variance when included. Excluding them from *both* arms leaves ICC ≈ 0.35–0.43 for OH-solanidine. So the finding is still about how below-LOQ subjects are handled — but by *inclusion of floor-pinned subjects*, not by corruption of retained subjects' values.

Removing the fixed time effect raises the ICC in every cell (0.349 → 0.393; 0.705 → 0.779), confirming that the postprandial trajectory was inflating the within-subject term, exactly as the review predicted.

**What survives, and is now better supported than before.** On a common subject set, with and without the time effect removed, OH-solanidine's ICC is **0.35–0.43** while SSDA's is **0.71–0.80**. The marker difference is not an artefact of differing subject sets or of the time trend. It holds in all four cells.

### The marker difference on the exact intersection, with uncertainty

The per-marker attribution above still uses each marker's *own* retained subjects — 17 for OH-solanidine, 18 for SSDA — so a difference between them was partly a difference of population. Restricting both to the **exact intersection** (n=17) and bootstrapping over subjects (2,000 resamples, subjects as the independent unit, timepoints not resampled):

| | OH-solanidine | SSDA | difference | bootstrap 95 % |
|---|--:|--:|--:|---|
| raw | 0.400 | 0.762 | **+0.362** | [+0.084, +0.596] — excludes zero |
| time effect removed | 0.430 | 0.811 | **+0.381** | [+0.135, +0.605] — excludes zero |

The difference is positive in favour of SSDA in both specifications, on identical subjects, with an interval excluding zero. **This is the strongest single result of the H5 analysis**: the semi-quantitative, counts-based analyte is the more reliable of the two, and that conclusion now survives population matching, removal of the postprandial trend, and explicit uncertainty quantification.

Three subjects are excluded from the intersection; they are the floor-pinned subjects, so this comparison is deliberately made in the population where neither marker is at the assay floor.

**What still stands about reading the font.** A value-only pipeline would have reported ICC 0.942 for OH-solanidine and called the marker highly reliable, when a comparable-population estimate is near 0.4. The headline figure is misleading; the reason is floor-pinned subjects rather than imputation per se.

### Descriptive figures, retained for the record

These are the raw all-subjects numbers. **Read them with the attribution above**: the apparent collapse is driven by excluding floor-pinned subjects, not by removing imputed cells from retained subjects.

The source imputes below-LOQ values and flags them by bold font only. Recomputing with those cells removed:

| | OH-solanidine/solanidine | SSDA/solanidine |
|---|---|---|
| ICC, all values incl. imputed | 0.942 | 0.886 |
| **ICC, below-LOQ excluded** | **0.400** | **0.754** |
| ICC, poor metaboliser excluded | 0.516 | 0.706 |
| between-subject SD (ln), all values | 1.480 | 1.638 |
| **between-subject SD (ln), below-LOQ excluded** | **0.284** | **0.945** |
| within-subject geometric CV, below-LOQ excluded | 35.9 % | 58.2 % |
| Kendall's W, all values → below-LOQ excluded | 0.682 → 0.663 | 0.825 → 0.830 |

For OH-solanidine the between-subject spread falls by a factor of ~5 on the ln SD and the ICC from 0.94 to 0.40. SSDA loses far less and its rank concordance is unchanged at 0.83.

**Correct reading:** these two columns are computed on different subject sets (20 versus 17), so the difference is not a censoring effect. The attribution analysis above shows the censoring-only effect is −0.037, in the opposite direction. The sentence in the first version of this document — "most of its apparent discriminative power came from imputed floor values creating artificial separation" — was wrong as a statement about imputation and right only as a statement about *including floor-pinned subjects*.

A pipeline reading values only, without inspecting the font, would still have reported ICC 0.94 for OH-solanidine and called the marker highly reliable, when a comparable-population estimate is near 0.4.

---

## The semi-quantitative analyte is the more robust one

SSDA is reported in **raw counts**, OH-solanidine in **ng/ml**. The expectation was that the quantified analyte would be the sounder marker. The opposite holds on every robustness measure: ICC after censoring (0.754 vs 0.400), Kendall's W (0.83 vs 0.66), rank agreement with the informative timepoint (+0.43…+0.75 vs +0.00…+0.55), and association sign consistency across all nine times (SSDA positive throughout; OH-solanidine crosses zero three times).

**This does not make SSDA the better marker outright.** At the informative timepoint OH-solanidine has the marginally *stronger* association, in both cohorts (Müller 29 h: 0.71 vs 0.66; ARIANA: 0.489 vs 0.405). The honest statement is that the two answer different questions: OH-solanidine is slightly more associated when sampling time is controlled, SSDA degrades less when it is not. Units are not interchangeable between cohorts, so only directions and orderings are compared.

---

## Food effect — inconclusive, not negative

Subject-level mean ln ratio, 3-day potato fasting (n=15) versus not (n=5):

| Marker | scenario | Δln (fold) | descriptive 95 % | Welch p |
|---|---|---|---|---|
| OH-solanidine/solanidine | all values | −0.556 (×0.57) | [−1.450, +0.337] | 0.24 |
| OH-solanidine/solanidine | below-LOQ excluded | −0.054 (×0.95) | [−0.403, +0.294] | 0.77 |
| SSDA/solanidine | all values | −0.931 (×0.39) | [−2.107, +0.244] | 0.14 |
| SSDA/solanidine | below-LOQ excluded | −0.724 (×0.48) | [−1.592, +0.144] | 0.14 |

Every interval covers zero. The direction is consistent — fasted subjects show lower ratios — but with 5 subjects in the comparison group this cannot be distinguished from noise. Reported as **`INCONCLUSIVE`**, explicitly not as evidence of no food effect. Note also that the fasting flag and the sampling-time effect are not separable here: both bear on recent dietary precursor supply.

---

## Direction consistency with ARIANA

| Contrast: ln ratio vs CYP2D6 activity score | ρ | 95 % |
|---|---|---|
| ARIANA n=47, OH-solanidine | +0.489 | [+0.234, +0.680] |
| ARIANA n=47, SSDA | +0.405 | [+0.134, +0.620] |
| Müller 29 h, OH-solanidine | +0.710 | [+0.38, +0.88] |
| Müller 29 h, SSDA | +0.660 | [+0.29, +0.86] |

Same sign, overlapping intervals, and the same ordering between the two markers. **Direction consistency holds — for the late sample only.** At 0–5 h the Müller associations are null, a regime ARIANA cannot test because it documents neither absolute sampling time nor repeated sampling.

**This is weaker evidence than it looks.** Sarömba, Stingl and Just author both studies, so the cohorts share a research group, laboratory and analytical platform. Shared method bias would produce exactly this agreement. The participants are different people, so no donor is double-counted, but this is a different population measured by the same group — not an independent laboratory.

---

## Status assignments, per the protocol vocabulary

| H5 item | Status | Basis |
|---|---|---|
| Sampling-time effect | **`PARTIALLY_SUPPORTED`** *(downgraded from `SUPPORTED`)* | Rank instability across timepoints is directly measured and robust. The association-strength comparison is **not** established: early-window intervals reach +0.69 and +0.79, overlapping the 29 h estimates, so this design does not separate the windows. |
| Intra-individual temporal stability | **`PARTIALLY_SUPPORTED`** | Marker-dependent, not a property of "the solanidine axis". On the exact intersection with the time effect removed, SSDA ICC 0.811 versus OH-solanidine 0.430. |
| Quantified versus semi-quantitative identification | **`SUPPORTED`**, direction opposite to expectation | On the exact intersection of subjects, ICC difference +0.362 raw and +0.381 time-adjusted, bootstrap 95 % [+0.084, +0.596] and [+0.135, +0.605], both excluding zero. The strongest single H5 result. |
| Sensitivity to censoring handling | **`SUPPORTED`** | The headline ICC moves from 0.942 to ≈0.40 depending on whether floor-pinned subjects are included. The attribution is subject inclusion, not cell-level imputation. |
| Food effect | **`INCONCLUSIVE`** | All intervals cover zero with n=5 in the comparison group. |
| Poor-metaboliser discrimination | **`ABSTAIN`** | The cohort's only poor metaboliser is below the limit of quantification at all 9 timepoints for both metabolites. |
| Direction consistency with ARIANA | **`PARTIALLY_SUPPORTED`** | Consistent sign and overlapping intervals at 29 h only, with shared laboratory and platform. |
| Between-day reproducibility | **`ABSTAIN`** | Nine samples follow a single meal. Within-occasion dispersion is not between-day reproducibility. |

## Change log

| Version | Date | Change |
|---|---|---|
| 1.0 | 2026-11-07 | First version. |
| 1.2 | 2026-11-07 | Two final statistical precisions after review. (a) The 80 %-power correlation ρ ≈ 0.76 was misused as an upper bound; it is not one. Replaced by the genuine compatibility bounds, the per-timepoint confidence intervals, whose widest early-window upper limits are +0.69 and +0.79 — overlapping the 29 h estimates, so the design does not separate the windows on association strength. (b) The marker ICC comparison used 17 versus 18 subjects. Recomputed on the **exact intersection** (n=17) with a 2,000-resample subject bootstrap: difference +0.362 raw and +0.381 time-adjusted, intervals [+0.084, +0.596] and [+0.135, +0.605], both excluding zero. The marker conclusion is strengthened, not weakened. |
| 1.1 | 2026-11-07 | Revised after external review. (a) The prespecification claim is corrected: objectives and criterion families were frozen in v1.0.0, the statistics were not, and the 29 h reference time was chosen after seeing a result. (b) **The ICC attribution is reversed**: on a common subject set with the fixed time effect removed, censoring *raises* the ICC slightly (Δ −0.037 and −0.020); the headline drop is caused by including three floor-pinned subjects. The marker-difference conclusion survives and is strengthened. (c) The sampling-time status is downgraded to `PARTIALLY_SUPPORTED` and a minimum-detectable-effect calculation added: the early-window nulls bound any association below ≈0.76 rather than excluding one. (d) Multiplicity restated over the correct primary family of 18. |

---

## What this means for the programme

The mission asks where mechanism-guided omics adds value and where it fails. This is a **where it fails** result, and a useful one for the manuscript's argument:

A biomarker can be mechanistically correct, strongly correlated with genotype under favourable conditions, and still behave unreliably when an external variable the mechanism itself predicts is left uncontrolled — here, recent dietary precursor intake. **Stated at the strength the data support:** which patient looks like a poor metaboliser changes with sampling time (rank agreement with the 29 h ranking as low as +0.00 for OH-solanidine), and the association with genotype is demonstrable only in the late sample. Whether the marker is genuinely less informative soon after a meal is *not* settled: the early-window intervals reach +0.69 and +0.79 and overlap the late-sample estimates, so n=20 cannot separate the two windows on association strength at all. The weaker claim — rank instability, demonstrated directly — is enough for the argument, and the stronger one is not made.

The source's own stool incubations support the same reading: solanidine is formed only with an active microbial community, so precursor supply is both diet- and microbiota-dependent, a variance source upstream of CYP2D6 that genotype cannot capture.

Two consequences for the remaining work:

1. **The ARIANA H1 analysis inherits an untested assumption.** ARIANA sampled at trough but documents neither absolute sampling time nor time since last meal. Its exposure-prediction result, favourable or not, will rest on the unverifiable assumption that its samples sit in the informative regime. This must be stated as a limitation of the H1 result rather than discovered afterwards.
2. **Censoring handling is not a detail, and neither is the subject set.** Had the imputation flags gone unread, OH-solanidine would have been reported at ICC 0.94 when a comparable-population estimate is ≈0.4. The lesson is sharper than "imputation matters": a reliability statistic computed over a population that includes subjects pinned at an assay floor is not a reliability statistic. Every later analysis using these analytes reports its below-LOQ handling *and* its subject set explicitly, and any ICC-style quantity is reported with the fixed time effect removed.

## Reproducing this

```bash
python scripts/run_h5_solanidine_marker_robustness.py \
    --muller "Solanidine standardized potato meal.xlsx" \
    --ariana "Solanidine_metoprolol_ratios_PGx.xlsx"
```

Both checksums are verified at run time and the script exits non-zero on mismatch. No individual measurement is printed and neither workbook is committed.
