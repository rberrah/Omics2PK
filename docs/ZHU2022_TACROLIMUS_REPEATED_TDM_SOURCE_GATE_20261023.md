# Zhu 2022 — irregular repeated tacrolimus TDM + metabolomics source gate

**Updated:** 23 October 2026  
**Study:** Zhu et al., *Frontiers in Pharmacology* 2022  
**DOI:** `10.3389/fphar.2022.1045843`  
**PMID:** `36386159`  
**PMCID:** `PMC9659571`

## Why this cohort is unusually relevant

The study reports **137 liver-transplant recipients and 693 repeated plasma
samples**. The metabolomics aliquots were derived from residual biospecimens
collected for routine tacrolimus therapeutic-drug monitoring.

For each TDM occasion the source study reports collection of clinical data
around the C0 measurement, including:
- tacrolimus daily dose;
- body weight;
- concomitant interacting medicines;
- clinical biochemistry;
- postoperative sampling time.

Tacrolimus C0 was measured from blood obtained within one hour before the
morning dose after a fixed dose for about three days.

The publication split the cohort into:
- development: **112 recipients / 570 samples**;
- validation: **25 recipients / 123 samples**.

This is a stronger Omics2PK candidate than a study where metabolomics and TDM
come from unrelated visits.

## Important limitation of the published analysis

The paper itself labels the design **repeated cross-sectional** and explicitly
states that prospective prediction from metabolite levels measured before a
future tacrolimus exposure still needs longitudinal verification.

Therefore:

- same-biospecimen metabolomics -> same C0 is **contemporaneous exposure
  phenotyping**, not prospective prediction;
- an actual E4 test requires an earlier metabolomics occasion to predict a
  later C0/C0-D for the same patient;
- all train/test splitting must be by whole patient.

## Why the He D7/M1/M3 source contract is not sufficient

Zhu has an irregular number and timing of TDM occasions (median around four,
reported range 3-6), rather than three fixed nominal visits.

Do not coerce the source into D7/M1/M3 labels.

The dedicated
`src/omics2pk/irregular_repeated_tdm_metabolomics.py` gate orders observations
using authenticated timestamps only.

## Minimum source extract

### Molecular table

```text
sample_id,biospecimen_id,subject_id,modality,collected_at
```

The metabolomics feature matrix is separately keyed by `sample_id`.

### Tacrolimus dose-event table

```text
dose_id,subject_id,drug,administered_at,amount,dose_unit,route
```

Actual administered events are required around every modeled C0. A nominal
twice-daily schedule is not an acceptable substitute.

### TDM observation table

```text
observation_id,biospecimen_id,subject_id,drug,
previous_dose_id,next_dose_id,sampled_at,
concentration,concentration_unit,
daily_dose,daily_dose_unit,weight_kg,postoperative_day
```

The `biospecimen_id` must be exactly shared with the metabolomics record.

This permits exact reconstruction of the study's dose/weight-normalized trough
target without inferring values from group summaries.

## Gate semantics

The gate requires:

1. exact same-biospecimen molecular/C0 identity;
2. same biological subject across omics, trough and bounding dose events;
3. pre-dose chronology:
   `previous dose < C0/metabolomics < next dose`;
4. one consistent concentration unit and one daily-dose unit;
5. observed daily dose and body weight;
6. real postoperative timing;
7. no duplicated trough time within one patient;
8. no fixed nominal visit schedule.

It returns only aggregate:
- linked subjects and sample count;
- distribution of samples per subject;
- number of patients with repeated qualified observations;
- number of possible earlier-omics -> later-C0 pairs;
- time-range/QC summaries.

No source patient IDs, concentrations, doses or feature values are emitted.

## Preregistered analyses after source acquisition

### E3-contemporaneous

For each same TDM biospecimen:

`clinical + dose + interacting drugs -> C0/D`

versus

`clinical + dose + interacting drugs + metabolomics -> C0/D`

Whole-patient grouped validation is mandatory.

### E4-prospective

Preferred landmark design:

`state at T_i -> C0/D at T_(i+1)`

with a clinical reference containing all information already available at
`T_i`, including previous C0 and dose history.

More distant horizons may be secondary analyses.

No later molecular measurement may enter an earlier landmark.

## Public-source status

The article/PMC public supplementary package exposes supplementary figures and
a small DOCX methodological supplement, but no source-qualified tabular
693-sample patient x metabolomics x C0/dose matrix has been established from
the public package.

Therefore the real numerical status remains:

`ABSTAIN_ORIGINAL_693_SAMPLE_MATRIX_NOT_ACQUIRED`

The preferred next acquisition is an author/institution-provided deidentified
sample-level extract retaining the exact biological and TDM keys above.
