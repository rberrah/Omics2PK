# H5 results — robustness of the solanidine CYP2D6 markers

**Date:** 2026-11-07
**Cohort:** Müller 2026 standardized potato meal, n=20 (doi `10.3389/fphar.2026.1886427`)
**Script:** `scripts/run_h5_solanidine_marker_robustness.py`
**Source checksum:** `68d53071…1950`, verified at run time
**Direction check against:** ARIANA n=47 (doi `10.3389/fphar.2026.1812007`, checksum `7397b13e…fe37`)

**Prespecification.** The H5 items executed here were frozen in protocol **v1.0.0, before this dataset was known to exist** — verifiable from commit `7632e9e`. The 3 h and 29 h reference times and Spearman as the primary association measure are protocol choices, not result-driven ones.

**This is not pharmacokinetics.** No drug was administered. Nothing here predicts exposure, and nothing here replicates the ARIANA solanidine-to-trough association.

---

## Headline result

**The markers are informative only in the late sample. In the 0–5 h window after a potato meal, neither ratio carries detectable CYP2D6 genotype information.**

With below-LOQ values excluded and Bonferroni correction across the 9 timepoints (threshold 0.00556):

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

## OH-solanidine's apparent reliability is largely an imputation artefact

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

For OH-solanidine the between-subject spread collapses by a factor of ~5 on the ln SD, roughly 27-fold in variance, and the ICC falls from 0.94 to 0.40. Most of its apparent discriminative power came from imputed floor values creating artificial separation. SSDA loses far less and its rank concordance is unchanged.

A pipeline reading values only — without inspecting the font — would have reported ICC 0.94 for OH-solanidine and called the marker highly reliable.

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
| Sampling-time effect | **`SUPPORTED`** | A strong time dependence exists; only the 29 h sample survives Bonferroni correction, for either marker. |
| Intra-individual temporal stability | **`PARTIALLY_SUPPORTED`** | SSDA retains ICC 0.754 after censoring; OH-solanidine falls to 0.400. Stability is marker-dependent, not a property of "the solanidine axis". |
| Quantified versus semi-quantitative identification | **`SUPPORTED`**, direction opposite to expectation | The counts-based analyte is more robust on every stability measure. |
| Food effect | **`INCONCLUSIVE`** | All intervals cover zero with n=5 in the comparison group. |
| Poor-metaboliser discrimination | **`ABSTAIN`** | The cohort's only poor metaboliser is below the limit of quantification at all 9 timepoints for both metabolites. |
| Direction consistency with ARIANA | **`PARTIALLY_SUPPORTED`** | Consistent sign and overlapping intervals at 29 h only, with shared laboratory and platform. |
| Between-day reproducibility | **`ABSTAIN`** | Nine samples follow a single meal. Within-occasion dispersion is not between-day reproducibility. |

---

## What this means for the programme

The mission asks where mechanism-guided omics adds value and where it fails. This is a **where it fails** result, and a useful one for the manuscript's argument:

A biomarker can be mechanistically correct, strongly correlated with genotype under favourable conditions, and still be unusable without controlling an external variable the mechanism itself predicts — here, recent dietary precursor intake. The solanidine ratio is not simply "a CYP2D6 phenotype measure"; it is one **conditional on sampling far from a potato-containing meal**. The source's own stool incubations support the same reading: solanidine is formed only with an active microbial community, so precursor supply is both diet- and microbiota-dependent, a variance source upstream of CYP2D6 that genotype cannot capture.

Two consequences for the remaining work:

1. **The ARIANA H1 analysis inherits an untested assumption.** ARIANA sampled at trough but documents neither absolute sampling time nor time since last meal. Its exposure-prediction result, favourable or not, will rest on the unverifiable assumption that its samples sit in the informative regime. This must be stated as a limitation of the H1 result rather than discovered afterwards.
2. **Censoring handling is not a detail.** Had the imputation flags gone unread, OH-solanidine would have been reported at ICC 0.94 instead of 0.40. Every later analysis using these analytes reports its below-LOQ handling explicitly.

## Reproducing this

```bash
python scripts/run_h5_solanidine_marker_robustness.py \
    --muller "Solanidine standardized potato meal.xlsx" \
    --ariana "Solanidine_metoprolol_ratios_PGx.xlsx"
```

Both checksums are verified at run time and the script exits non-zero on mismatch. No individual measurement is printed and neither workbook is committed.
