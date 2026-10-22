# He 2022 — longitudinal tacrolimus metabolomics + pre-dose C0 source gate

**Updated:** 22 October 2026  
**Study:** He et al., *Molecules* 2022  
**DOI:** `10.3390/molecules27113517`  
**PMID:** `35684454`

## Why this cohort matters

The paper reports **109 kidney-transplant recipients** followed at **D7, M1 and
M3**, with plasma untargeted metabolomics and tacrolimus trough concentrations
(`C0`) obtained from fasting blood **before the morning dose**.

The source article reports:
- 2,695 initially detected LC-MS features;
- 1,514 retained after missingness/QC preprocessing;
- longitudinal plasma collection at D7, M1 and M3;
- a primary pharmacodynamic objective based on T-lymphocyte response;
- raw untargeted metabolomic MS data as Supplementary Materials;
- study data additionally available on request from the corresponding authors.

This makes the cohort a strong candidate for a repeated molecular-state /
tacrolimus-exposure study, but **not yet a qualified Omics-to-PK dataset**.
The public article does not establish that the downloadable supplement contains
a complete subject × visit × dose × C0 crosswalk.

## Critical temporal distinction

A tacrolimus trough is **pre-dose**.

For a valid source event, the chronology must be:

`previous tacrolimus dose < C0 / molecular blood draw < next morning dose`

The generic Omics2PK post-dose concentration linkage contract therefore cannot
be reused unchanged. The dedicated
`src/omics2pk/predose_trough_link_contract.py` requires explicit bounding dose
events and refuses a C0 linked to the morning dose as though it were a
post-dose concentration.

## Minimum authorized source package

Three linked tables are sufficient.

### 1. Molecular samples

```text
sample_id,subject_id,visit,modality,collected_at
```

Requirements:
- `visit` frozen to D7, M1 or M3;
- `modality=METABOLITE`;
- timezone-aware collection timestamp;
- one modeled molecular sample per subject × visit after technical replicate
  resolution;
- a separate feature matrix keyed by `sample_id`;
- original raw/processed provenance, ion mode, batch/QC and preprocessing
  version retained privately.

### 2. Tacrolimus dose events

```text
dose_id,subject_id,drug,administered_at,amount,dose_unit,route
```

Requirements:
- actual administered event times, not only a group mean or daily-dose summary;
- enough history to identify the dose immediately preceding each C0 and the
  next morning dose;
- formulation and complete relevant dosing history should be retained in the
  richer clinical source even though the minimal identity-time gate uses the
  fields above.

### 3. Trough observations

```text
observation_id,subject_id,visit,drug,previous_dose_id,next_dose_id,sampled_at,concentration,concentration_unit
```

Requirements:
- exact same subject key as molecular and dose tables;
- one modeled C0 per subject × visit;
- concentration > 0 with one consistent physical unit;
- sampled strictly between the explicitly identified preceding and next dose;
- no inferred dose IDs from row order.

No direct identifiers are required in the research extract.

## What the gate returns

The gate reports aggregate-only:
- number of linked independent subjects;
- linked D7/M1/M3 counts;
- number of subjects with at least two qualified visits;
- number of possible earlier-omics → later-C0 pairs;
- maximum molecular/C0 sampling-time separation.

It does **not** export patient IDs, timestamps, doses or concentrations.

## Two scientifically distinct analyses

### A. Same-visit molecular state → same-visit C0

This is a **concurrent exposure-sensing association**.

Because the molecular plasma sample and C0 come from essentially the same
pre-dose visit, a good prediction does not establish that the metabolome
prospectively determines future tacrolimus clearance. The metabolome may partly
reflect prior tacrolimus exposure, renal function, corticosteroids, MMF,
post-transplant recovery and other treatment effects.

It can still be useful as an E3-style exposure phenotype benchmark, but must be
labeled contemporaneous.

### B. Earlier molecular state → later-visit C0

This is the preferred longitudinal test.

Examples:
- D7 omics → M1 C0;
- D7 omics → M3 C0;
- M1 omics → M3 C0.

All train/test splitting must be by **whole patient**, never by visit, so no
other visit from a held-out patient enters training.

The primary clinical comparator should use the same patients and available
pre-outcome information, ideally including:
- previous C0;
- actual dose/formulation/history;
- age/sex/body size;
- hematocrit;
- renal/liver function;
- corticosteroid and interacting co-medications;
- CYP3A5/other PGx only if genuinely measured.

The incremental question is whether the earlier molecular state improves
held-patient prediction beyond that reference.

## What is not allowed

- assuming 109 × 3 = 327 complete visits;
- treating a group-average tacrolimus dose as a patient dose history;
- linking by row order;
- using a post-outcome molecular sample to predict an earlier C0;
- splitting visits from the same patient across train and test;
- treating same-visit metabolomics → C0 as prospective causality;
- calling C0-only prediction a full PopPK concentration-time or AUC model;
- inferring missing dose timestamps from a nominal twice-daily schedule without
  source confirmation.

## Current source status

The paper and PMC record confirm the study design and state that raw untargeted
MS data are available as Supplementary Materials. The publisher supplement
route is publicly referenced, but the current automated environment has not
materialized a verified original file containing the complete individual
metabolomics + dose + C0 identity/time crosswalk.

Therefore the clinical numerical gate remains:

`ABSTAIN_ORIGINAL_LONGITUDINAL_C0_CROSSWALK_NOT_ACQUIRED`

The next useful acquisition is not another model. It is the original
pseudonymized subject/visit crosswalk, the metabolomics feature matrix, the
individual C0 observations and actual bounding tacrolimus dose events.
