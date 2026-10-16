"""Optional future V1 individual molecular signature export contract.

Results with group effects, p-values or fold changes are NOT a patient-level
covariate. An inferred MOFA/SWEET/pathway score intended for a held-out PK
individual must have been produced out-of-sample with a frozen reference.
No signature model or PK estimation is performed by this reader.
"""
from __future__ import annotations

from collections import Counter
from math import isfinite
from typing import Mapping, Sequence


class OmicSignatureError(ValueError):
    pass


def qualify_individual_signatures(rows: Sequence[Mapping[str, object]]) -> dict:
    """Review a research score export without outputting individual vectors."""
    if not rows:
        raise OmicSignatureError("Individual scores are not exported")
    ids = set()
    subjects = set()
    samples = set()
    modality_counts = Counter()
    signatures = set()
    scopes = set()
    references = set()
    for row in rows:
        required = ("subject_id", "sample_id", "assay_id", "omic",
                    "signature_id", "score", "score_scale", "fit_scope",
                    "reference_model_sha256", "feature_definition_sha256",
                    "source_exchange_analysis_id")
        if any(not str(row.get(k, "")).strip() for k in required):
            raise OmicSignatureError("Missing score provenance, assay identity, value or model SHA")
        modality=str(row["omic"])
        if modality not in ("transcriptomics","proteomics","metabolomics","integrated"):
            raise OmicSignatureError("Unknown individual score modality")
        scope=str(row["fit_scope"])
        if scope not in ("training_frozen_test_projection","out_of_fold_training"):
            raise OmicSignatureError("In-sample, transductive or group-only score cannot be used in PK held-out comparison")
        for key in ("reference_model_sha256","feature_definition_sha256"):
            sh=str(row[key])
            if len(sh)!=64 or any(c not in "0123456789abcdef" for c in sh):
                raise OmicSignatureError("Model and feature set need precise source SHA-256")
        x=row["score"]
        if isinstance(x,bool):
            raise OmicSignatureError("Boolean score invalid")
        try: val=float(x)
        except (TypeError,ValueError,OverflowError):
            raise OmicSignatureError("Numeric observed individual score required") from None
        if not isfinite(val):
            raise OmicSignatureError("No undefined individual model scores")
        key=(str(row["assay_id"]),str(row["signature_id"]))
        if key in ids:
            raise OmicSignatureError("Duplicate per-assay signature")
        ids.add(key)
        subjects.add(str(row["subject_id"]))
        samples.add(str(row["sample_id"]))
        signatures.add(str(row["signature_id"]))
        scopes.add(scope)
        modality_counts[modality]+=1
        references.add(str(row["reference_model_sha256"]))
    return {
        "status":"QUALIFIED_INDIVIDUAL_SIGNATURE_SCHEMA_NOT_PREDICTIVE_EFFECT",
        "n_unique_independent_biological_subjects_at_most":len(subjects),
        "n_original_unique_specimens":len(samples),
        "n_assay_signature_rows":len(rows),
        "n_named_signatures":len(signatures),
        "n_frozen_reference_models":len(references),
        "n_by_omic":dict(modality_counts),
        "eligible_fit_scopes":sorted(scopes),
        "no_patient_or_score_values_exported":True,
        "note":"Real PK gain requires frozen train/test patient grouping, biomarker collection times, age/PGx reference and measured dose-concentration/NCA outcome sources",
    }
