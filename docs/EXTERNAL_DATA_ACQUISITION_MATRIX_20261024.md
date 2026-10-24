# Omics2PK — external data acquisition matrix

**Frozen:** 24 October 2026

## Purpose

The current Omics2PK software and source-validation layer is sufficiently developed for the four remaining open scientific questions. New generic model features should not be added simply because source data are unavailable.

This document defines the exact next external acquisition for each remaining issue and the condition under which work should stop rather than fabricate a dataset.

| Priority | GitHub issue | Scientific question | Minimum next source | First executable gate | Stop / ABSTAIN condition |
|---|---|---|---|---|---|
| P0 | #31 tacrolimus | Does a molecular phenotype improve tacrolimus exposure prediction beyond clinical + CYP3A5 PGx? | **Zhu n137/693:** deidentified subject + biospecimen key, processed metabolomics, C0, daily dose, weight/postoperative time, interacting drugs; exact dose events/timestamps for prospective E4. **Phapale n29** remains the preferred predose true C(t)/AUC source if the original 29×1256 matrix + PK profiles can be obtained. | `irregular_repeated_tdm_metabolomics` for Zhu; relevant private source gate for Phapale. | No deterministic same-patient molecular↔PK rows, or only published model coefficients/figures. |
| P0 | #14 BIOMARGIN↔ABIS | Can BIOMARGIN molecular patients be deterministically linked to real tacrolimus/MPA TDM? | **Aggregate counts first only**: exact BIOMARGIN↔ABIS matches, timed dose/concentration records, pre-outcome omics, CYP3A5, repeated omics, holdout eligibility + governance attestations. | `biomargin_abis_feasibility` using the zero-filled JSON template. | Exact crosswalk count = 0, linkage not authorized, or no timed drug PK/TDM rows. |
| P0 | #63 Zhang/Ohtsuki | Does quantitative CYP protein add predictive information beyond RNA for measured CYP activity in the same liver donor? | Original Zhang n100 or Ohtsuki n17 donor table with stable donor IDs, RNA, protein, activity/CLint, units and BLQ/missing states. | Existing `multicyp_rna_protein_activity_source_gate` and frozen author/source manifest. | No legitimate donor-level rows; aggregate correlations, plots, row-order reconstruction or cross-cohort fusion are not substitutes. |
| P1 | #50 rosuvastatin | Does predose metabolomics improve measured rosuvastatin AUC/Cmax in the same 40 volunteers? | Author-provided five processed CSVs or equivalent deidentified 40-person processed omics + original AUC/Cmax + stable subject crosswalk and units. | `qualify_local_original_rosuvastatin_csvs.py`, then existing private E3 held-out protocol. | Google Drive remains private/inaccessible and no author-authorized equivalent is supplied. Raw mzML alone cannot reconstruct the missing patient PK vector. |

## Required comparator for tacrolimus

A future tacrolimus result is not positive merely because omics improves a weak clinical-only model.

The minimum reference is:

1. clinical/dose/TDM information available before the prediction landmark;
2. hematocrit and relevant organ-function/interacting-drug covariates;
3. **CYP3A5 PGx when genuinely measured**;
4. other CYP3A PGx only when justified and observed in the target population.

Yang 2022 / GSE196772 is frozen as a comparator-definition source. Its published R2 values are explanatory in-sample results, not held-out performance.

## HLC SWEET result still pending execution

The first remaining **non-clinical functional** experiment that can directly test the value of inductive SWEET is already implemented:

`benchmark_authorized_hlc_original --sweet-all`

It compares, across ten HLC CYP activity endpoints:

- age + sex;
- direct source-annotated CYP RNA;
- train-only selected RNA;
- train-only PCA4;
- inductive SWEET v1;
- topology-permuted SWEET negative control.

The source archives must remain private. Only aggregate results should be committed.

Do not rerun this benchmark without the exact qualified HLC archives. The absence of mounted source bytes is not a reason to replace them with another cohort or synthetic outcomes.

## Development stop rule

Until one of the external source gates above receives real rows, the following are **not** justified as priorities:

- additional SWEET variants;
- alternative neural encoders;
- larger unsupervised latent spaces;
- automated hyperparameter search;
- browser-side patient modelling;
- synthetic replacement of missing clinical endpoints;
- merging unrelated cohorts into a pseudo-multiomics subject.

Software work remains appropriate only when it is required to ingest a newly obtained legitimate source, prevent a demonstrated leakage/identity error, or execute a preregistered comparison.

## Decision hierarchy after data arrive

For each new source:

1. source identity / SHA / provenance;
2. deterministic biological identity and chronology;
3. privacy/governance gate;
4. fixed simple reference model;
5. direct molecular representation;
6. SWEET/mechanistic representation only when the modality supports it;
7. negative controls;
8. patient/donor-held-out evaluation;
9. independent cohort/site/time validation before transportable claims.

A negative incremental result is a valid Omics2PK result.
