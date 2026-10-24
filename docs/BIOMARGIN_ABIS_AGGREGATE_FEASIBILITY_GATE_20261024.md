# BIOMARGIN ↔ ABIS aggregate feasibility gate

**Updated:** 24 October 2026

## What changed

The earlier Omics2PK issue treated an ABIS linkage as an unverified possibility.
A fresh public-source audit clarifies the situation:

1. **BIOMARGIN is a real multicentre kidney-transplant omics program**, with
   blood/urine/biopsy sampling in four European transplant centres including
   CHU Limoges.
2. **ABIS is a real clinical therapeutic-drug-monitoring database/service** used
   for Bayesian immunosuppressant dose adaptation across transplant centres.
3. The national **REXETRIS** project already performs authorized pseudonymized
   linkage of **ABIS + CRISTAL + SNDS**. CHU Limoges reports 47,842 kidney
   transplant recipients in that linked research database for transplants
   between 2005 and 2021.
4. None of these public sources proves that a BIOMARGIN subject identifier is
   currently linked to an ABIS patient record.

Therefore the BIOMARGIN ↔ ABIS route is no longer a vague data-source idea, but
it remains an **institutional feasibility question**, not a public dataset.

## Important distinction

BIOMARGIN itself does **not** contain a qualified tacrolimus or MPA PK dataset.

ABIS can contain clinical TDM/dose-adaptation records, but the existence of an
ABIS record for a BIOMARGIN participant must be established by an authorized
deterministic patient crosswalk.

Do not match by:
- age;
- sex;
- centre;
- transplant date;
- biopsy date;
- combinations of quasi-identifiers.

Only source-system identifiers or an approved pseudonymization/linkage service
may establish the match.

## First query to the data custodians

The first exchange should contain **aggregate counts only**.

No patient IDs or omic/PK values are needed to decide whether the project is
worth pursuing.

Frozen counts:

```text
biomargin_molecular_patients
exact_biomargin_abis_crosswalk_patients

tacrolimus_tdm_patients
tacrolimus_timed_dose_concentration_patients
tacrolimus_preoutcome_omics_patients
tacrolimus_cyp3a5_patients
tacrolimus_repeated_omics_patients

mpa_tdm_patients
mpa_timed_dose_concentration_patients
mpa_preoutcome_omics_patients

independent_site_or_time_holdout_patients
```

Required yes/no attestations:

```text
deterministic_same_person_crosswalk_confirmed
authorized_linkage_and_research_use_confirmed
omics_collection_time_available
no_direct_identifiers_in_aggregate_exchange
```

The machine-readable template is:
`docs/templates/biomargin_abis_aggregate_feasibility.json`.

## Interpretation

### Stop immediately

If:

`exact_biomargin_abis_crosswalk_patients = 0`

then BIOMARGIN cannot be the paired Omics2PK clinical proof cohort.

It remains useful for immunological / graft-injury omics and Omics2PD work, but
not for the proposed Omics→PK comparison.

### Crosswalk exists, but no timed PK

If a deterministic match exists but there are no original timed dose +
concentration records, do not call the cohort PK-ready.

A medication list, current daily dose or isolated historical trough without a
usable timing context is not sufficient for the intended PopPK question.

### E3 feasibility

For tacrolimus, an E3 feasibility set requires at minimum:
- exact same patient;
- molecular sample collection time;
- actual dose/formulation history appropriate for the PK endpoint;
- measured concentration with sampling time;
- relevant clinical covariates;
- CYP3A5 when genuinely observed for the comparator.

A pre-outcome molecular sample is preferable. A molecular sample obtained after
the exposure endpoint cannot be used as a causal predictor of that earlier
exposure.

### E4 feasibility

A prospective longitudinal candidate additionally requires repeated molecular
states and later PK/TDM outcomes, with all splits grouped by patient.

The relevant question becomes:

`molecular state at T_i -> PK/TDM at T_(i+1)`

not merely:

`molecular state measured together with C0 -> same C0`.

## Comparator

Any tacrolimus molecular model must be compared on the same people against the
clinical + PGx reference. The public Yang 2022 / GSE196772 audit is now frozen
as a useful comparator definition:
- clinical age + sex + HCT;
- CYP3A5*3 minimum PGx layer when measured;
- additional ancestry-relevant CYP3A variants only when justified.

The published Yang R2 values are not external performance targets; they simply
show that PGx carries substantial explanatory information and cannot be omitted
from the baseline.

## Governance boundary

The public evidence that REXETRIS links ABIS, CRISTAL and SNDS demonstrates that
authorized linkage infrastructure exists in the same institutional ecosystem.
It does **not** authorize or prove BIOMARGIN linkage.

The following still require confirmation from the data controller / ethics and
data-governance chain:
- whether BIOMARGIN identities can be linked to ABIS;
- the legal basis / consent compatibility;
- which centres are represented in ABIS;
- whether historical dose/concentration timestamps are retained;
- whether CYP3A5 is available;
- whether linked data can be analyzed on approved infrastructure.

## Current status

Public evidence supports:

`ABIS_IS_REAL_CLINICAL_TDM_DATABASE__AUTHORIZED_LINKAGE_INFRASTRUCTURE_EXISTS`

Public evidence does **not** support:

`BIOMARGIN_ABIS_CROSSWALK_CONFIRMED`

Therefore Omics2PK remains:

`ABSTAIN_BIOMARGIN_ABIS_PATIENT_LINKAGE_NOT_YET_CONFIRMED`
