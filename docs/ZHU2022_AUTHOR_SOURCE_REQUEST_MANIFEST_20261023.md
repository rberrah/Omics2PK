# Zhu 2022 tacrolimus repeated metabolomics — author/source acquisition manifest

**Date:** 23 October 2026  
**Study:** Zhu et al., *Frontiers in Pharmacology* 2022  
**DOI:** `10.3389/fphar.2022.1045843`  
**PMID:** `36386159`  
**PMCID:** `PMC9659571`

## Purpose

Obtain the minimum de-identified original source required to test whether a
measured molecular state adds held-patient predictive information for
tacrolimus exposure beyond clinical/TDM information.

The publication reports **137 liver-transplant recipients / 693 repeated TDM
samples**, with untargeted plasma metabolomics performed from residual material
from the clinical tacrolimus TDM biospecimen.

The public article reports:
- development set: **112 recipients / 570 samples**;
- validation set: **25 recipients / 123 samples**;
- approximately 3-6 repeated samples per recipient;
- C0 obtained within one hour before the morning dose after approximately three
  days at a fixed dose;
- daily tacrolimus dose, body weight, concomitant medicines and clinical
  biochemistry collected around C0;
- C0/D defined as C0 divided by **daily dose normalized to actual body weight**.

The public NLM/PMC supplementary catalogue contains four TIF figures and one
small DOCX (`DataSheet1`), not a source-qualified 693-sample analysis matrix.

## Preferred source format

The preferred acquisition is the **native source tables and processed
metabolomics matrix used by the authors**, with their original pseudonymous
keys and data dictionary.

Do not ask authors to manually recreate a matrix if the original analysis
files can be shared under appropriate institutional terms.

No direct identifiers, names, dates of birth, addresses or medical-record
numbers are needed.

## Level 1 — minimum package for contemporaneous E3

This is sufficient to test:

`clinical/TDM reference -> same-sample C0/D`

versus:

`clinical/TDM reference + same-biospecimen metabolomics -> same-sample C0/D`

Minimum sample-level fields:

1. stable pseudonymous `subject_id`;
2. stable `biospecimen_id` shared by the TDM result and metabolomics sample;
3. metabolomics `sample_id`;
4. sample/order grouping used for the published development vs validation split;
5. postoperative sampling day;
6. measured tacrolimus C0 and physical unit;
7. actual daily tacrolimus dose and unit;
8. actual body weight used for dose normalization;
9. formulation/brand if available;
10. interacting medicines at that sampling occasion, especially:
    voriconazole, fluconazole, Wuzhi capsule, omeprazole/esomeprazole,
    caspofungin and DHP calcium-channel blockers;
11. clinical chemistry used in the published analysis, preferably including
    hematocrit, creatinine/eGFR, albumin, bilirubin and liver enzymes;
12. processed metabolomics features for the exact same biospecimen;
13. feature annotation, ion mode, batch/run and QC metadata;
14. missing/BLQ/not-measured status rather than implicit zero imputation.

For this E3 use, exact dose-event timestamps are **preferred but not mandatory**
if the source table itself identifies the observation as the protocol-qualified
pre-dose TDM C0 described in the publication.

This level supports a contemporaneous exposure-phenotype analysis only. It must
not be described as prospective prediction or full PopPK.

## Level 2 — preferred package for prospective E4

To test whether an **earlier molecular state predicts a later tacrolimus
exposure**, additionally request:

1. timezone-aware collection date/time for every TDM/metabolomics biospecimen;
2. actual tacrolimus administration history around each observation;
3. explicit preceding and next dose-event IDs/times when available;
4. dose amount, route and formulation for those events;
5. treatment interruptions or dose changes;
6. the exact longitudinal order of all retained samples per participant;
7. previous C0 and previous dose history available at each prediction landmark;
8. any genuinely measured recipient/donor CYP3A5 or other PGx variables.

The Omics2PK strict E4 gate never infers dose timestamps from a nominal twice
daily regimen.

## Suggested exchange tables

If native files cannot be shared, three de-identified tables plus one feature
matrix are enough.

### A. Sample/TDM table

Use the header-only template:

`docs/templates/zhu2022_tdm_sample_template.csv`

### B. Metabolomics long table

Use:

`docs/templates/zhu2022_metabolomics_long_template.csv`

A native wide feature matrix is equally acceptable and preferable if that is
the original analysis format.

### C. Dose-event table — E4 preferred

Use:

`docs/templates/zhu2022_dose_event_template.csv`

## Exact primary comparisons after acquisition

### E3 — contemporaneous

Use whole-patient-disjoint evaluation.

Reference:
- prior/current clinically available tacrolimus information;
- daily dose and body weight;
- postoperative time;
- relevant interacting drugs;
- appropriate laboratory variables;
- PGx only if genuinely observed.

Then add molecular representations in a prespecified sequence:
1. published/targeted metabolite panel;
2. train-only low-dimensional metabolomic representation;
3. broader train-frozen representation if sample size supports it.

Primary outcome may reproduce source `log2(C0/D)`, where
`D = daily dose / actual body weight`.

### E4 — prospective

Preferred first landmark:

`molecular + clinical state at T_i -> C0/D at T_(i+1)`

Compare against the same clinical history available at `T_i`, especially
previous C0 and dose.

All observations from one recipient stay in one outer partition.

No feature selection, normalization, dimensionality reduction or tuning may use
the held-out patient or future observations.

## Privacy and source-integrity rules

Before fitting:
- hash every received original file;
- record transfer/source provenance and permitted research use;
- preserve original pseudonymous keys;
- never join by row order;
- never reconstruct patient rows from figures or published coefficients;
- never publish source individual rows, metabolomic matrices or predictions to
  GitHub;
- report only aggregate source counts, hashes, missingness and model metrics.

## Current public-source conclusion

Targeted searches of the article, PMC/NLM supplementary catalogue and major
metabolomics repositories have not identified a public source-qualified
693-sample patient-level matrix.

The NLM object named `DataSheet1` is the small DOCX supplementary file, not
the patient-level metabolomics/TDM table.

Therefore retain:

`ABSTAIN_ORIGINAL_693_SAMPLE_MATRIX_NOT_ACQUIRED`

until a legitimate source package is obtained.

## Public corresponding contacts

The article lists:
- Yizhun Zhu — `yzzhu@must.edu.mo`
- Weihong Ge — `geweihong.nju@gmail.com`

A request should ask for the **minimum de-identified research extract above**,
not identifiable clinical records.
