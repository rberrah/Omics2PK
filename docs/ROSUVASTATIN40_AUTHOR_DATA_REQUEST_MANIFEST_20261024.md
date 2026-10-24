# Rosuvastatin n40 — author data request manifest

**Frozen:** 24 October 2026
**Study:** Silveira et al. 2021, Frontiers in Pharmacology
**DOI:** `10.3389/fphar.2021.752960`
**Public omics accession:** `MTBLS3128`

## Purpose

The public raw LC-MS files are real, but the original processed five CSVs referenced by the author notebook are on a private/inaccessible Google Drive and the public MetaboLights deposit does not expose a verified 40-person AUC/Cmax vector.

Further automated attempts to retrieve the private Drive files are **not** a useful next step. The next legitimate path is an author-authorized deidentified research extract.

## Minimum dataset requested

One stable pseudonymous biological subject key must bridge all files.

### 1. Subject / PK table

Required columns:

`subject_id`
`rosuvastatin_dose`
`dose_unit`
`formulation`
`auc_value`
`auc_definition`
`auc_unit`
`cmax_value`
`cmax_unit`

Preferred additional columns when available:

`dose_time`
`pk_sampling_times_or_schedule_id`
`sex`
`age`
`weight`
`renal_function`
`relevant_concomitant_medication`

Requirements:
- original observed/NCA individual AUC and Cmax, not values digitized from figures;
- clear AUC definition such as AUC0-t or AUC0-inf;
- consistent units and BLQ/missingness information;
- exactly one biological subject identity per participant.

### 2. Predose metabolomics matrices

Preferred: the same five processed CSVs used by the published author notebook, or an equivalent reproducible processed export.

Required metadata:
- `subject_id` or an explicit deterministic crosswalk to subject_id;
- predose/baseline sample flag;
- ESI mode;
- technical replicate / pooled QC identity;
- feature identifiers;
- quantitative feature values;
- missing/nonquantifiable state;
- processing version or description.

Do not collapse technical replicates into new biological subjects.

### 3. Crosswalk

If PK and metabolomics use different pseudonyms, provide a two-column deterministic crosswalk:

`pk_subject_id, metabolomics_subject_id`

No row-order matching is acceptable.

## Optional but useful provenance

- source file names and checksums;
- ProgenesisQI export version/settings if preserved;
- internal-standard / pooled-QC normalization description;
- batch/run-order metadata;
- feature annotation level and confidence;
- exact relationship between ESI+ and ESI− acquisition rows.

## What is not required

- names, initials, dates of birth, addresses or other direct identifiers;
- raw medical-record identifiers;
- public redistribution permission if analysis can remain on approved/private infrastructure.

## First local execution after receipt

1. preserve the received source files unchanged;
2. compute SHA-256 checksums;
3. run `scripts/qualify_local_original_rosuvastatin_csvs.py`;
4. verify exactly 40 biological humans and a deterministic PK↔omics crosswalk;
5. confirm predose status and AUC/Cmax units;
6. only then run the existing private E3 held-out benchmark.

## Fixed scientific comparison

Primary question:

`clinical/reference -> clinical/reference + predose metabolomics`

on the **same 40 people**, with all preprocessing, feature selection and tuning learned inside training folds.

Published Elastic-Net R2 values are author results and must not be reused as independent validation.

## Stop condition

If the investigators cannot provide the original individual AUC/Cmax vector linked to the processed predose metabolomics, keep:

`ROSUVASTATIN_E3_ABSTAIN__ORIGINAL_PATIENT_PK_VECTOR_NOT_ACQUIRED`

Raw mzML alone is scientifically useful for MultiOmics method validation but cannot reconstruct the missing individual clinical PK target.
