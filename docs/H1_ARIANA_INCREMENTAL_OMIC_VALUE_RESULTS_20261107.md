# H1 results — incremental value of omics in ARIANA

**Date:** 2026-11-07
**Cohort:** ARIANA n=47 (doi `10.3389/fphar.2026.1812007`), checksum `7397b13e…fe37` verified at run time
**Script:** `scripts/run_h1_ariana_incremental_omic_value.py`
**Endpoint:** observed pre-dose metoprolol concentration, ng/mL. **Never** called AUC, clearance or exposure.
**Status ceiling:** `exploratory`. This project has analysed ARIANA repeatedly, so no result here may be reported as confirmatory however favourable it is.

## Prespecification, stated exactly

**Frozen in protocol v1.0.0 (commit `7632e9e`):** the arm definitions, the primary contrast and its loss, the loss hierarchy with its reversal rule, the repeated subject-disjoint nested CV scheme, the mandatory influence reporting, the prohibition on `lnMR OH-metoprolol/metoprolol` as a predictor, the MAP-BE exclusion, and the ≥10 % exploratory success threshold.

**Operationalised at analysis time, not frozen:** fitting on the log of the trough with back-transformation, the penalty grids, the feature encodings, the bootstrap count.

Validation: 5 outer folds × 2 repeats, 3 inner folds, identical folds across every arm and model. Standardisation and hyperparameter selection inside the training fold only. Bootstrap over the 47 patients; fold repetition is never counted as additional subjects.

## The primary contrast

**Ridge M3 (clinical + PGx + omic) versus Ridge M1 (clinical + PGx), MAE, identical folds.** Ridge deliberately: it estimates its coefficients from the cohort and cannot borrow strength from literature-fixed constants, so it is the harder test.

| | value | 95 % descriptive interval |
|---|--:|---|
| **ΔMAE** | **−16.97 ng/mL (−32.1 %)** | **[−33.76, −4.69]** — excludes zero |
| Δ log MAE | −0.286 | [−0.415, −0.153] — excludes zero |
| Δ trimmed MAE (10 %) | −8.86 (24.48 vs 33.34) | — |

**Both mandatory robustness co-criteria agree with the primary loss.** No reversal, so the status is not capped by §2.1.

**The frozen exploratory success criterion — ≥10 % relative MAE improvement with a descriptive interval excluding zero — is met**, at 32.1 %.

The prespecified dose-weighted secondary target reproduces it: ΔMAE −10.80 (−31.5 %), 95 % [−20.15, −2.89].

**This is not circular.** The predictor is a solanidine-derived ratio, measured from a dietary precursor; the target is the metoprolol concentration. They share no analyte and no denominator. The prohibited `lnMR OH-metoprolol/metoprolol` column, which *does* share a denominator with the target, is absent from every arm and a test enforces that.

## The more interesting result: the omic substitutes for the genotype

Full arm table, all three model families on identical folds. MAE in ng/mL.

| model | M0 clinical | M1 clinical+PGx | M2 clinical+omic | M3 clinical+PGx+omic |
|---|--:|--:|--:|--:|
| Ridge | 55.64 | 52.81 | **34.32** | 35.84 |
| Elastic Net | 56.71 | 54.23 | 35.12 | **35.05** |
| XGBoost | 57.48 | 58.25 | **37.53** | 38.60 |

Two readings follow, and the second matters more than the first.

1. **Omics carry substantial incremental information over clinical + PGx.** The effect is large, consistent across three model families, and robust to the loss used.

2. **The genotype adds nothing once the omic is present — and slightly hurts.** M2 beats M3 for Ridge (34.32 vs 35.84) and for XGBoost (37.53 vs 38.60), and they are level for Elastic Net. Meanwhile M1 improves on M0 only marginally (55.64 → 52.81 for Ridge), and for XGBoost not at all (57.48 → 58.25).

So the measured metabolic ratio **substitutes for** the genotype-derived activity score rather than complementing it. This is biologically coherent: the ratio is a *functional* read-out that already integrates genotype, phenoconversion by comedication, and whatever else modulates enzyme activity, so once it is in the model the activity score is redundant and costs a degree of freedom.

It is also the same direction as the Jarmusch finding, where a prior functional probe beat the CYP3A genotype and adding genotype on top of the functional measure *worsened* MAE. Two public cohorts, two different enzymes, same ordering: **measured function > genotype, and genotype adds nothing on top of measured function.**

This is a secondary, exploratory observation within the H1 family. It is the hypothesis NetPK exists to test, and it is not established by these two cohorts.

## Influence analysis, mandatory under §2.2

| | value |
|---|--:|
| patients improved | 34 of 47 |
| patients worsened | 13 of 47 |
| net gain | +797.4 ng/mL summed |
| share of net gain from the five most-improved patients | **73.7 %** |
| mean gain by endpoint quartile Q1 / Q2 / Q3 / Q4 | +19.5 / +4.5 / +11.3 / +31.6 |

**The gain is heavily concentrated**: five patients out of 47 carry about three quarters of it. That must be reported next to the headline, not buried.

Two things nevertheless distinguish this from a pure outlier artefact: the gain is **positive in all four endpoint quartiles**, so it is not confined to one end of the concentration range; and a clear majority of patients (34 of 47) improve individually. The concentration reflects that absolute errors scale with concentration — the endpoint spans 7 to 572 ng/mL — which is exactly why the log and trimmed losses are mandatory co-criteria, and both of them agree.

## Limitations that bound this result

1. **Exploratory only.** Prior analysis exposure caps the claim; freezing the protocol stops further result-guided choices but cannot turn a re-analysis into independent validation.
2. **No replication is available.** No public cohort carries the CYP2D6/solanidine axis with drug exposure, so this cannot be upgraded in this phase.
3. **The H5 finding bounds it.** ARIANA sampled at trough but documents neither absolute sampling time nor time since last meal. H5 showed that which patient looks like a poor metaboliser changes with sampling time. This H1 result therefore rests on the unverifiable assumption that ARIANA's samples sit in the informative regime. Favourable as the result is, that assumption is not testable in this cohort.
4. **Gain concentration**, as above.
5. **n=47, single centre, single laboratory, single assay.**
6. **Multiplicity.** One primary contrast, reported unadjusted. Every other contrast here, including the M2-versus-M3 decomposition, is secondary within the H1 family and requires Holm adjustment before any claim.

## Status

| H1 item | Status |
|---|---|
| Incremental value of omics over clinical + PGx, ARIANA | **`SUPPORTED`** as an exploratory success: prespecified criterion met (32.1 % ≥ 10 %), interval excluding zero, both robustness co-criteria agreeing. |
| Independent confirmation of that gain | **`ABSTAIN`** — no public replication cohort exists for this axis. |
| Incremental value of genotype over clinical | **`NOT_SUPPORTED`** — marginal for Ridge and Elastic Net, absent for XGBoost. |
| Incremental value of genotype *on top of* the omic | **`NOT_SUPPORTED`** — M3 is worse than M2 for two of three model families, level for the third. |

`SUPPORTED` here means the prespecified criterion was actually met, on the primary contrast, with the robustness co-criteria agreeing — not that omics are established as useful for pharmacokinetics in general.

## What follows

The decisive question is now H2, and the result above sharpens it. The omic gain is real in this cohort, but **two variables** produced it. H2 must establish whether a mechanism-guided selection finds those variables for better reasons than a statistical selector does, at equal panel dimension, including against a deliberately wrong mechanistic panel — and whether NetPK's advantage, if any, comes from its selection rather than from pharmacological coefficients fixed from the literature. The primary contrast was placed on Ridge precisely so that this question cannot be answered by borrowed strength.

## Reproducing this

```bash
python scripts/run_h1_ariana_incremental_omic_value.py \
    --source "Solanidine_metoprolol_ratios_PGx.xlsx"
python scripts/run_h1_ariana_incremental_omic_value.py \
    --source "Solanidine_metoprolol_ratios_PGx.xlsx" \
    --target "c(metoprolol)weighted [ng/ml]"
```

The checksum is verified at run time and the script exits non-zero on mismatch. No individual prediction is printed and the workbook is never committed. `reports/omics2pk_results.Rmd` renders these JSON outputs into the manuscript tables and figures.
