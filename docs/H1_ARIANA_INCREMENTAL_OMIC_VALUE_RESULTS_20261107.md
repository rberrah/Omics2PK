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

## Corrected: no additional benefit from PGx was *detected* once the omics are in

**Revision 1.1, 2026-11-07.** The first version of this section claimed the omic "substitutes for the genotype". That is stronger than the data support, and the review was right to reject it. **Every contrast involving the PGx block has a 95 % interval covering zero**, including PGx added to clinical alone:

| model | M3 vs M2 (PGx on top of omic) | M1 vs M0 (PGx on clinical) |
|---|---|---|
| Ridge | +1.52 [−2.14, +5.04] | −2.83 [−6.46, +0.75] |
| Elastic Net | −0.08 [−3.19, +2.80] | −2.48 [−6.02, +0.92] |
| XGBoost | +1.07 [−0.30, +2.72] | +0.77 [−1.23, +3.15] |

So the admissible statement is: **no additional benefit from pharmacogenetics was detected once the omic ratios were included** — not that the omic replaces the genotype. "Not detected" is not "absent", and protocol §10 forbids reading a zero-crossing interval as evidence of equivalence.

To avoid stopping at "non-significant", each contrast now also reports whether its interval can exclude a benefit of the **prespecified size** (≥10 % relative MAE improvement). That distinguishes a genuine negative at the frozen margin from a design that simply cannot tell.

Full arm table, all three model families on identical folds. MAE in ng/mL.

| model | M0 clinical | M1 clinical+PGx | M2 clinical+omic | M3 clinical+PGx+omic |
|---|--:|--:|--:|--:|
| Ridge | 55.64 | 52.81 | **34.32** | 35.84 |
| Elastic Net | 56.71 | 54.23 | 35.12 | **35.05** |
| XGBoost | 57.48 | 58.25 | **37.53** | 38.60 |

Two readings follow, and the second matters more than the first.

1. **Omics carry substantial incremental information over clinical + PGx.** The effect is large, consistent across three model families, robust to the loss used, and its intervals exclude zero.

2. **No PGx benefit was detected, in either position.** The point estimates are small and inconsistent in sign across model families, and every interval covers zero. A plausible mechanism is that the ratio is a *functional* read-out already integrating genotype, phenoconversion by comedication and other modulation, so the activity score adds little while costing a degree of freedom — but that mechanism is a hypothesis consistent with the data, not something these data establish.

### The PGx block was badly specified, and separating it answers the attribution

The review caught a real design fault: `M1` and `M3` bundled the genotype-derived activity score with the comedication flag, so an effect could not be attributed to either. Separated (MAE, ng/mL):

| model | clinical | + genotype | + comedication | + both | + omic | + omic + genotype | + omic + both |
|---|--:|--:|--:|--:|--:|--:|--:|
| Ridge | 55.64 | 53.50 | 55.19 | 52.81 | **34.32** | 36.27 | 35.84 |
| Elastic Net | 56.71 | 53.62 | 56.52 | 54.23 | **35.12** | 36.25 | 35.05 |
| XGBoost | 57.48 | 59.41 | 57.07 | 58.25 | **37.53** | 38.37 | 38.60 |

The attribution is clear once the two are apart:

- **The comedication flag contributes essentially nothing.** Added to clinical: −0.45 [−1.87, +0.98] for Ridge, −0.18 [−1.58, +1.22] for Elastic Net, −0.42 [−1.80, +0.52] for XGBoost. A ≥10 % benefit is **excluded in all three model families**. Added on top of the omic, likewise excluded in all three. XGBoost in fact never used it at all — the contrast is exactly 0.00 [0.00, 0.00], because at depth 1–2 with the omic present the binary flag was never chosen as a split.
- **Whatever small signal the bundled block carried came from the activity score**, not the comedication flag: genotype alone gives −2.14 [−5.67, +1.51] for Ridge and −3.09 [−7.12, +0.97] for Elastic Net, close to the bundled −2.83 and −2.48.
- **The genotype's own contribution is `INCONCLUSIVE`, not negative**, for the linear models: not detected, but the intervals reach −10.2 % and −12.6 % relative, so a benefit of the prespecified size cannot be ruled out at n=47. For XGBoost it is +1.92 [−1.04, +6.01], where a ≥10 % benefit *is* excluded.
- **Adding the genotype on top of the omic is `NOT_SUPPORTED` at the frozen margin** in all three families: +1.95 [−1.07, +4.96] Ridge, +1.13 [−1.13, +3.25] Elastic Net, +0.84 [−0.25, +1.97] XGBoost. Not detected, *and* a ≥10 % benefit is incompatible with every interval. This is the one genuinely negative statement in the decomposition, and it is now attributed to the genotype specifically rather than to a bundled block.

**The Jarmusch parallel was overstated and is withdrawn.** The first version claimed "two public cohorts, two enzymes, same ordering". That does not hold: Jarmusch's "functional probe" is a **prior administered-probe PK measurement of the same drug**, which is close to a direct measurement of the outcome process, not an independent omic biomarker like a solanidine ratio. Its panels, sample sizes and predictor type all differ. Jarmusch provides a separate negative H1 test on untargeted metabolomics; it does **not** independently corroborate anything about genotype versus measured function.

This whole section is secondary and exploratory within the H1 family and requires Holm adjustment before any claim. It is the hypothesis H2 exists to test, and it is not established here.

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

Statuses are assigned on two distinct questions: was an effect **detected** (does the interval exclude zero), and can a benefit of the **prespecified size** be **excluded** (is the most favourable value in the interval short of −10 %). Only the second licenses `NOT_SUPPORTED`; failing both means `INCONCLUSIVE`.

| H1 item | Status |
|---|---|
| Incremental value of omics over clinical + PGx, ARIANA | **`SUPPORTED`** as an exploratory success: prespecified criterion met (32.1 % ≥ 10 %), interval excluding zero, both robustness co-criteria agreeing. |
| Independent confirmation of that gain | **`ABSTAIN`** — no public replication cohort exists for this axis. |
| Incremental value of PGx over clinical alone | **`INCONCLUSIVE`** *(revised from `NOT_SUPPORTED`)* — not detected in any model family, but the Ridge and Elastic Net intervals reach −11.6 % and −10.6 %, so a benefit of the prespecified size cannot be excluded. The design cannot tell. |
| — of which: the **genotype** alone over clinical | **`INCONCLUSIVE`** for Ridge and Elastic Net; **`NOT_SUPPORTED`** at the margin for XGBoost. |
| — of which: the **comedication flag** alone over clinical | **`NOT_SUPPORTED`** — a ≥10 % benefit is excluded in all three model families. XGBoost never used the variable. |
| Incremental value of the **genotype** on top of the omic | **`NOT_SUPPORTED`** — not detected, and a ≥10 % benefit is excluded in all three model families (+1.95, +1.13, +0.84). A negative at the frozen margin, attributed to the genotype specifically rather than to a bundled block. |
| Incremental value of the **comedication flag** on top of the omic | **`NOT_SUPPORTED`** — excluded in all three. |

`SUPPORTED` means the prespecified criterion was actually met on the primary contrast with the robustness co-criteria agreeing — not that omics are established as useful for pharmacokinetics in general.

## What follows — and why ARIANA cannot carry H2 on its own

The decisive question is now H2, and this result sharpens it while also showing where it cannot be answered.

**ARIANA has exactly two candidate omic variables.** It demonstrates that those two carry predictive information; it is **not a demanding test of selection**. Choosing 2 of 2 is not selection, and no selector — mechanistic, SHAP, Elastic Net or random — can distinguish itself on a panel that size. Running the H2 comparator set on ARIANA alone would produce a near-tie by construction and would be close to meaningless.

H2 therefore needs the wide panels, and they are available: HLC (whole-transcriptome, admissible since the v1.1.0 gate ruling), EPA, Wegler, and the public metabolomic sets. Those carry hundreds to tens of thousands of candidates, which is where a selection method can actually be tested — and where the false-positive and wrong-mechanism controls become informative.

**The endpoint classes must stay separate.** HLC, EPA and Wegler carry `functional_in_vitro` endpoints, so a selection result there is the *functional analogue* of H2 and is never reported in the same effect column as a clinical-PK result. ARIANA remains the only cohort with clinical exposure, PGx and named quantified biomarkers matched per individual, and Jarmusch the only other clinical-PK test — already negative on untargeted metabolomics.

Two distinct H2 experiments follow, and they answer different questions:

| Question | Design |
|---|---|
| Does NetPK **select** better? | Same model (Ridge, then XGBoost), **same number of variables**, selection by NetPK versus SHAP, Elastic Net, nested univariate, repeated random panel, a deliberately wrong mechanistic panel, and no selection |
| Does mechanistic **integration** then help? | Identical NetPK variables, statistical integration versus parameterisation as an enzyme-activity score feeding a PopPK structure |

Separating them is the point: the first isolates the value of the selection, the second the value of the mechanistic integration. Conflating them is how a method takes credit for borrowed pharmacological coefficients. The primary contrast was placed on Ridge for the same reason.

The mission's imperative check applies throughout: verify that the selection is **actually produced by the algorithm**, not written into the code after inspecting results.

## Revision log

| Version | Date | Change |
|---|---|---|
| 1.0 | 2026-11-07 | First version. |
| 1.1 | 2026-11-07 | Three corrections after review. (a) "The omic substitutes for the genotype" is **withdrawn**: every PGx contrast has an interval covering zero, so the admissible claim is that no additional PGx benefit was *detected*. Each contrast now also reports whether a benefit of the prespecified size can be excluded, so `NOT_SUPPORTED` is distinguished from `INCONCLUSIVE`. (b) The PGx block bundled the activity score with the comedication flag; **decomposed** into `M1a`/`M1b`/`M3a`/`M3b`, which attributes the small bundled signal to the activity score and shows the comedication flag contributes nothing in any model family. (c) The **Jarmusch parallel is withdrawn**: its "functional probe" is a prior administered-probe PK measurement of the same drug, not an independent omic biomarker, so it does not corroborate anything here. Also added: why ARIANA, with two candidate variables, cannot carry H2. The primary contrast is unchanged. |

## Reproducing this

```bash
python scripts/run_h1_ariana_incremental_omic_value.py \
    --source "Solanidine_metoprolol_ratios_PGx.xlsx"
python scripts/run_h1_ariana_incremental_omic_value.py \
    --source "Solanidine_metoprolol_ratios_PGx.xlsx" \
    --target "c(metoprolol)weighted [ng/ml]"
```

The checksum is verified at run time and the script exits non-zero on mismatch. No individual prediction is printed and the workbook is never committed. `reports/omics2pk_results.Rmd` renders these JSON outputs into the manuscript tables and figures.
