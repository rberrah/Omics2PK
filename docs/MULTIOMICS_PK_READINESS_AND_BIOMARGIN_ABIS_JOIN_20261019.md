# MultiOmics-to-PK readiness: original person / specimen / clock / PGx coverage

**Date: 2026-10-19. Research tooling only. The output is source availability, NOT evidence of improved PK.**

## Why this check is essential

The central research comparison in Omics2PK is whether genuine individual RNA, proteins, metabolites, their combinations, or biology-informed scores improve **measured drug PK** over clinical variables and authentic pharmacogenetics **on the exact same biological people**.

A six-view omics study may have 131 overall patients yet only 35 complete six-view specimens; gene expression data in GEO can contain 384 assays without proven exact person keys linked to the 131-person subset, and **neither observation proves linked tacrolimus PK**. Counting assay rows as people, assigning the same chronological row order across sources, or fabricating time-zero ISO timestamps would make an apparent positive biomarker result invalid.

`src/omics2pk/pk_readiness.py` therefore reports *per-person original source intersections and missing provenance*, before any model fitting.

## Command-line check with no private clinical sidecar

Install the existing dependency-free Python package and run:

```bash
omics2pk audit-pk-readiness /path/to/private_or_public_multiomics_export/exchange
```

The command reads the original `DOWNSTREAM_EXCHANGE_V1` `manifest.json` and SHA-256-verified `samples.csv` / source features and matrices, and reports:

- Distinct independent biological `subject_id` values and total assay rows **separately**.
- Number of unique people with exactly ONE source specimen/assay for each omic.
- Number of people with multiple assays in the same modality, excluded until visits or replicates are explicitly frozen.
- Intersections for each original modality, each pair, and all three modalities.
- An explicit **unknown** value for real observed PK and true PGx availability when no trusted clinical sidecar has been supplied, instead of inventing target labels or qualifying eligibility.
- No subject/sample IDs, original molecular intensities, patient records, PK outcomes or actual source values printed.

**Important:** sample intersections are **metadata coverage only**, not an assurance that transcriptomes represent hepatic CYP function, that matched proteins are present, or that assays were collected before a PK prediction landmark.

## Private clinical availability audit (optional Python API)

In an authorized original-data environment, the researcher/data steward may supply a metadata-only in-memory list with keys:

```python
{
    "subject_id": "<real original matched pseudonym; never commit>",
    "has_clinical": True,
    "has_measured_pk": True,
    "has_pgx": True,
    "has_dose_timing": True,
    "assay_ids_by_omic": {"transcriptomics": "<exact source assay ID>"},
    "sample_ids_by_omic": {"transcriptomics": "<exact source specimen ID>"},
    "collected_at_by_omic": {"transcriptomics": "<genuine ISO time WITH offset>"},
    "index_dose_at": "<original dose / prediction landmark ISO time WITH offset>",
}
```

This is an **example schema**, not a released original observation. Exact source original specimen timestamps must ALSO be present in the MultiOmics `samples.csv` `collection_time` column. Only exact chronological matches that precede the original dose/landmark within 30 days (configurable) enter temporal qualification. A label such as `T0` never qualifies as an absolute timestamp.

```python
from omics2pk.pk_readiness import assess_pk_readiness

aggregate = assess_pk_readiness(
    "/protected/multiomics/exchange",
    original_authorized_private_availability_records,
    lookback_days=30,
)
# Only aggregate per-modality and paired intersections; no original IDs or values.
```

The method reports **three different levels** without conflating them:

1. Original source RNA/protein/metabolite sample availability and independent biological people.
2. Attested clinical + measured PK + original PGx reference availability for an eligible person (authenticity of PK/PGx values is **not independently certified**; none are loaded).
3. More restrictive exact assay + **specimen** source join and absolute source-verified molecular-before-dose clocks. The final counts measure *potential eligibility to test future PK hypothesis*, not a biomarker effect.

For any single modality or combination, `n_predose_qualified_by_modalities` is the intersection of **the same eligible people**, not the number of molecule peaks or independent timepoints. Counts must be interpreted alongside source missingness, stratification and sample-use permissions. Subgroup reporting by centers could reveal participants if strata are tiny and therefore is intentionally not emitted by this API.

## Strengthened actual native PK adapter

In `run_multimodal_exchange_pk`, an original patient `subject_id` now also requires the authentic `sample_ids_by_omic` in addition to `assay_ids_by_omic`, matching the same individual's source assay and specimen. An original-source-declared study must provide a source ISO `collection_time`; a sidecar may **not** invent an absolute collection clock from a generic `T0`. Both sides must agree and precede the PK prediction landmark. Existing tests use only 12–14 **fictional** people; no patient source uploads.

## Direct consequence for BIOMARGIN ↔ ABIS

Before requesting files with thousands of patient rows, request only an authorized **aggregate availability table** with original biological-person intersection counts stratified into:

- full clinical+verified CYP3A4/CYP3A5 PGx+actual tacrolimus dose/PK source;
- those same people with each original independent RNA/protein/metabolite layer;
- those same people with exact **two/three modalities** and valid pre-landmark sampling clocks;
- cases with matching original sample ID, assay ID and longitudinal visit keys;
- counts and missingness by clinical post-transplant period and center, subject to minimal-cell privacy.

The public BIOMARGIN/ROCKET transcriptomics sources do **not** automatically supply actual PK, genotypes or original `GEO ↔ BIOS ↔ ABIS` source identity. No original patient rows, files or private alias keys belong in GitHub, Actions, PR discussions or third-party data logs.

## Not yet supported

The audit does NOT fit MOFA2, SWEET, XGBoost or PopPK. It is **not** a PK assay, a clinical dose recommendation or a medical device. It does not independently confirm that provided original clinical PK, genotypes and collection clocks are true; this remains source-governance review. For a valid incremental PK comparison, still require prospectively frozen patient-heldout folds, original measured drug concentration/time or qualified individual AUC, actual dosing records, and independent external validation.
