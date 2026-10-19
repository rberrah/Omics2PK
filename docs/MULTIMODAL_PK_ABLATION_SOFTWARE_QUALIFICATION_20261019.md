# Fair single- vs multi-omics → observed PK: implementation qualification

**19 October 2026 — exploratory research methods, NOT a tacrolimus or other real-drug validation**

## What is newly implemented

`src/omics2pk/multimodal_pk_ablation.py` has a public Python API `compare_modalities_pk(...)` for **1 to 3 molecular assay modalities** measured in the SAME genuine biological people with an independently measured single scalar PK endpoint. `src/omics2pk/multimodal_exchange_pk.py` adds a private-source adapter `run_multimodal_exchange_pk(...)` from original versioned MultiOmics `DOWNSTREAM_EXCHANGE_V1` assay matrices, with a separately authorized **in-memory** person/Pgx/PK sidecar.

This is a software-comparison protocol only. There is **NO original tacrolimus omics plus authentic PGx plus measured drug PK linked cohort inside GitHub or a new positive result**. It does not run an NLME fit on a concentration–time likelihood, infer any causal enzyme activity or prove the original SWEET algorithm.

## Scientific estimand

Compare out-of-patient predictive MAE for:
`clinical + authentic PGx + omic representation` against
`clinical + authentic PGx` **on the same participants and same PK endpoint**.
If true source PGx is unavailable, label the baseline `clinical`; never claim genotype superiority or impute a fabricated genotype.

The PK target may be drug AUC, observed CL/F, a qualified numerical summary or another **actual originally measured in-vivo** PK endpoint. Its physical unit must be declared, but the code cannot independently verify that a caller's unit or biomarker baseline classification is correct. Regressing a pre-estimated PK parameter is NOT the same as fitting measured individual drug concentrations with original dosing schedules; that requires a separate NLME/PopPK engine.

### Prespecified benchmark arms

| Name | Representation | Modeling rule |
|---|---|---|
| `clinical` | Original clinical covariates | Ridge alpha=1 on log(observed PK) |
| `clinical_pgx` | Clinical + genuinely observed variant dosages 0/1/2 | Same Ridge. PRIMARY when present |
| `<omic>_pca2` | Each modality independently, 2 train-fitted PC scores | Same Ridge + clinical/PGx |
| `<omic>_direct_ridge` | Train-standardized original numeric features | Same Ridge + clinical/PGx |
| `<omic>_fixed_score` | **Externally prespecified**, signed fixed weights on train-scaled original feature IDs | Same Ridge; no score weights estimated on PK |
| `early_fusion_pca2` | Concatenated separate train-only PCs for 2–3 modalities | Same Ridge + clinical/PGx |
| `late_fusion_fixed_average` | **Fixed unweighted average** of held-out single-omic PCA+PGx predictions | Not tuned/stacked on test folds |
| `reference_xgboost` and `omic_pca_xgboost` | Baseline vs multimodal PCA factors | Optional deterministic low-capacity XGBoost, **blocked for n<80**, no hyperparameter grid search |

No biological winner is assumed. XGBoost n80 is an explicit **minimum software guard, not a mathematical power guarantee**; a much larger cohort and genuine external validation may be needed.

### Per-participant source requirements

- **Independent biological subject** one true specimen/assay per selected omic for a specified visit. No assay-row or repeated-timepoint pseudoreplication.
- Original **source feature IDs**, exact source `assay_id`, `subject_id`, and `sample_id`; no position-based source alignment.
- `collected_at_by_omic` with original absolute timezone-aware timestamp (one per selected modality), `index_dose_at` for the actual PK prediction administration/landmark, original tissue/compartment, actual collected molecular state **before** the index and within prespecified 30-day lookback by default. A categorical `T0` label is NOT an ISO clock and cannot be accepted.
- Original per-person `measured_pk`, named target and unit, original numeric clinical covariates, optional actual genotype dosage vector. Absent PGx => no PGx arm. Partial PGx **fails closed**.
- MultiOmics source V1 `manifest.json` with checksummed `samples.csv`, `features.csv` and `matrices/<omic>.csv`. True clinical sidecar is NOT bundled, logged, committed or an Actions artifact. For `source_kind='original_source_declared'`, upstream V1 status must be `qualified_for_declared_scope`; this **is an upstream declared status, not independent ethical/source validation**.
- If some persons lack one selected modality, defaults to **ABSTAIN**. Explicit `allow_complete_case=True` allows intersection only, reports the number excluded, and reruns **all compared models on the same complete-case people**. This is prone to selection bias, which needs sensitivity analysis beyond this engine.

**Measurement scale qualification:** source normalized RNA/proteins accepted as externally preprocessed values; authentic raw RNA counts REJECTED. For MS1/peak-area only, the code can apply per-person TIC fraction ×1e6 followed by log1p, with any inter-person standardisation recomputed only on each training fold. This is **not guaranteed to be an appropriate analytical normalization** for every LC-MS experiment. Upstream global batch correction or high-dimensional feature alignment may already have been transductive, and the program deliberately does **not** claim it prevented upstream leakage. Matrix source IDs do not imply molecules are positively identified.

## Cross-validation and reporting

- Canonical original participant ID sort, with identical **5-fold held-person partitions** (customizable 2–10).
- Entire modality StandardScaler, PCA2, Ridge/XGBoost fits computed on each outer training partition. No heldout PK outcome enters feature projection. Fixed score weights require prior definition; cannot prove preregistration by declaration alone.
- Per-model: **aggregate-only** MAE, RMSE, signed bias, P20 and incremental MAE; paired 999 bootstrap-of-heldout-loss intervals **WITHOUT model refitting** (descriptive, not uncertainty on fitted methodology and not independent external validation).
- No IDs, raw sample values, original clinical records, genotype calls or individual heldout predictions in return objects, GitHub logs or tests.
- Original measured PK outcome unit, timing and source provenance are caller-attested and must be independently qualified. Scores are research covariates, not bedside recommendation.

## Offline local/private use

Install Python 3.11+ and optional research benchmark dependencies:

```bash
python -m pip install -e '.[benchmark]'
python -m unittest discover -s tests -p 'test_multimodal_pk_ablation.py' -v
python -m unittest discover -s tests -p 'test_multimodal_exchange_pk.py' -v
```

The primary pure-library function is for **original same-subject patient arrays kept in an approved in-memory research process**:

```python
from omics2pk.multimodal_pk_ablation import compare_modalities_pk

# These keys are only an API sketch: do not create fake evidence by
# replacing participant measurements with values invented for this call.
report = compare_modalities_pk(
    subject_ids=authorized_original_subject_ids,
    modalities={
        "transcriptomics": {
            "feature_ids": original_gene_ids,
            "values": source_normalized_person_by_gene_expression,
            "collected_at": original_omic_iso_times,
            "tissue": "source_verified_tissue",
            "value_type": "normalized_log2_expression",
        },
        # Optionally proteomics and/or metabolomics, same genuine people.
    },
    clinical=original_patient_clinical_covariates,
    genotype_dosages=original_actual_variant_dosages,  # or None
    measured_pk=original_measured_pk_outcomes,
    index_dose_at=original_iso_dose_times,
    pk_target_name="observed_drug_AUC",
    pk_target_unit="original_source_verified_units",
    source_kind="original_source_declared",
    source_sha256=original_source_verified_sha256,
)
# report contains aggregate statistical error summaries only.
```

For the native upstream handoff, replace the abstract above with:

```python
from omics2pk.multimodal_exchange_pk import run_multimodal_exchange_pk

report = run_multimodal_exchange_pk(
    "/private/authorized_source/exchange",
    authorized_original_in_memory_pk_sidecar,
    selected_omics=("transcriptomics", "metabolomics"),
    pk_target_name="observed_drug_AUC",
    pk_target_unit="original_source_verified_units",
    source_kind="original_source_declared",
    source_sha256=original_source_verified_sha256,
)
```

Never upload protected clinical/omics source ZIPs to GitHub, CI artifacts or user-visible issue threads. Execute locally or on approved institutional infrastructure.

## Next step not blocked by unavailable tacrolimus source

- Finish qualified synthetic source/temporal/reordered-identity tests and monitor the dedicated CI on every change.
- Use approved public source matched RNA/proteomics/metabolomics to **qualify the representation software**, without inventing an observed PK label.
- Once a real **paired same-patient** tacrolimus or other drug original source is obtained, preregister a small number of score definitions and exact PGx reference. Compare at original patient/site-heldout scale and obtain an independent external cohort for biological conclusions.
- After positive externally qualified scalar-PK evidence, integrate actual dose–concentration–time PopPK/NLME with appropriate IIV, uncertainty and identifiability diagnostics, as a distinct experimental phase.
