# Measured donor-matched hepatic RNA → CYP metabolic activity benchmark

**2026-10-16 | Stage E2: real in-vitro phenotype, NOT patient PK or clinical PopPK validation.**

## Original paired source
Liu J et al., Xenobiotica 2021, DOI https://doi.org/10.1080/00498254.2020.1867929 . Official US EPA deposited raw underlying workbook https://doi.org/10.23719/1519197 ; workbook SHA256 **cfc9dc334064ea61e99ca6fd39c2c8ef7c621d000a37c4da2c95f26ca85e36dd** (60,345 bytes).

The original EPA xlsx has three sheets: Figure 1 (expression overview), Figure 2 (in-vitro CYP metabolic activity with KUMC and Xenotech IDs) and Table 3 (RT-qPCR CYP mRNA with both identical source ID fields). The two latter worksheets contained **exactly 100 matched original biological donors, zero unmatched dual IDs and zero duplicated subject keys**, despite Figure 2's extra figure/layout rows. Original individual records are used transiently and never committed or included in CI artifacts/logs.

Eight prespecified matched transcript / microsomal enzyme assay pairs are: CYP1A2, CYP2B6, CYP2C9, CYP2C19, CYP2D6, CYP2E1, CYP3A4 and CYP4A11. This is **targeted qPCR** mRNA, not global RNA sequencing. CYP3A4 activity may include CYP3A5 effects; CYP4A11 may be CYP4A9/11 combined; interpret as functional proxy, not isoform-specific protein measurement.

### Frozen method
- Unit of independence = **human donor** (not eight CYP outcomes or technical repetitions).
- Leave-one-donor-out for each enzyme separately, with 94–100 complete paired donors depending on recorded missing assays.
- Two fixed ridge models, alpha=1: reference **age + sex**, experimental **age + sex + same CYP qPCR mRNA**.
- Train-only feature centering/scaling, same folds for reference/experimental, log1p enzyme activity target and log1p RNA feature, backtransform predictions into original activity units (pmol/min/mg microsomal protein).
- Secondary reference training log-mean and direct RNA alone; no fitted hyperparameter selection, no PGx, no network feature discovery.
- Paired bootstrap of donor absolute errors (1,500 resamples, seed 20261016): exploratory intervals, no eight-CYP multiplicity adjustment, no refitted nested bootstrap.

## First original-source donor-held-out results

**All MAEs and differences are pmol/min/mg protein.** Negative delta means adding RNA improves prediction versus age+sex.

| CYP measured probe activity | Complete paired donors | MAE age+sex | MAE age+sex+RNA | MAE delta | Relative reduction | Paired donor bootstrap 95% |
|---|---:|---:|---:|---:|---:|---|
| CYP1A2 | 98 | 31.249 | 25.660 | −5.589 | 17.9% | [−8.901,−2.341] |
| CYP2B6 | 100 | 72.297 | 60.540 | −11.757 | 16.3% | [−18.144,−5.511] |
| CYP2C19 | 94 | 72.397 | 70.867 | −1.530 | 2.1% | [−3.839,+0.735] |
| CYP2C9 | 100 | 648.412 | 628.740 | −19.673 | 3.0% | [−54.395,+18.420] |
| CYP2D6 | 100 | 144.330 | 126.327 | −18.003 | 12.5% | [−32.056,−2.980] |
| CYP2E1 | 100 | 953.020 | 957.697 | **+4.677 (worse)** | −0.5% | [−25.219,+37.950] |
| CYP3A4/5 activity | 99 | 2297.027 | 1959.749 | −337.278 | 14.7% | [−560.067,−131.913] |
| CYP4A9/11 activity | 98 | 611.976 | 595.807 | −16.170 | 2.6% | [−41.746,+8.862] |

**Seven of eight** prespecified endpoints showed lower held-donor-out MAE upon adding the corresponding transcript; four have *uncorrected* paired bootstrap 95% intervals excluding zero. CYP2E1 **worsened**, an important prespecified null/negative result. Do not average absolute MAE across different enzyme activity scales or treat eight outcomes as eight independent donor cohorts.

Code: src/omics2pk/epa_liver_cyp_validation.py and scripts/validate_epa_liver_rna_cyp.py, with original data SHA and strict dual-donor provenance. Source RNA and activity matrices are read ephemerally from official original EPA URL. Public CI runs the **actual 100-liver data** and prints only aggregate QC and model results; synthetic tests separately check subject matching, missingness, demographics and train-only validation.

## Independent numerical replication on the SAME measured donor data

The original-source scientific pipeline was **independently reproduced using scikit-learn 1.7.2 StandardScaler + Ridge(alpha=1)**, a separately implemented numerical estimator and normalization pipeline, under GitHub Actions run https://github.com/rberrah/Omics2PK/actions/runs/38930709329. All **eight CYP endpoints** reproduced the reference-age+sex and age+sex+RNA LO(D)O MAE values to the six published decimals; eight maximum absolute discrepancies **0.000000** on printed results. Actual source-original SHA-256 and donor crosswalk rechecked independently; no original donor values committed or printed. This is **independent estimator parity on the SAME cohort**, NOT biological external replication or an independent population.

## Limits and decision
- **E2 preclinical functional metabolism evidence:** original individual RNA does explain some additional variability in measured CYP activity beyond age/sex for this dataset.
- **NOT a proof of clinical PK:** no actual administered dose, patient concentration, systemic CL, AUC or PopPK model.
- **NOT a PGx-adjusted baseline:** CYP genotype does not appear in this joined original public workbook. Biological causality, nonlinearity, sample-source confounding and multiple endpoints require further testing.
- **Next:** independent estimator replication, multiplicity and outlier/assay sensitivity, multiomic independent source, then real molecular↔clinical PK joint data; BIOMARGIN alone contains no PK, any ABIS linkage remains hypothetical.
- **Website:** future nonclinical educational/research demonstration only once reviewed; not a dosing or MIPD deployment.

GitHub run: https://github.com/rberrah/Omics2PK/pull/19 .
