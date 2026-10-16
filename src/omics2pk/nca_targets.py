"""Clinical-PK independent NCA-target contract for a MultiOmics research exchange.

Receives only REAL previously computed NCA endpoints with complete units,
exposure interval, original original-PK source SHA and subject/occasion IDs.
No PK outcomes are generated or fabricated by MultiOmics, and no source values
are exported by the QC functions. This E3 contract is only an input gate.
"""
from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict
from math import isfinite
from pathlib import Path
from typing import Mapping, Sequence

from .exchange import ExchangeError, _csv, _file, inspect_exchange
from .multiomics_pk_readiness import _temporal_hour


class NCATargetError(ValueError):
    """Unverifiable PK NCA endpoint, subject join or prospective state."""


UNITS = {
    "Cmax": "mg/L",
    "Tmax": "h",
    "AUC0last": "mg*h/L",
}
REQUIRED = (
    "subject_id", "occasion_id", "dose_id", "drug", "specimen", "analyte",
    "dose_time_h", "window_end_h", "endpoint", "value", "unit",
    "n_observed_concentrations", "pk_source_sha256",
)


def _number(v, field, *, strictly_positive=False):
    if isinstance(v, bool):
        raise NCATargetError(f"{field} may not be boolean")
    try:
        n = float(v)
    except (ValueError, TypeError, OverflowError):
        raise NCATargetError(f"{field} must be numeric") from None
    if not isfinite(n) or n < 0 or (strictly_positive and n == 0):
        raise NCATargetError(f"{field} must be finite and positive/zero per protocol")
    return n


def qualify_nca_targets(rows: Sequence[Mapping[str, object]]) -> dict:
    if not rows:
        raise NCATargetError("No original measured NCA targets")
    seen = set()
    subjects = set()
    occasions = set()
    endpoints = Counter()
    intervals = defaultdict(set)
    sources = set()
    for r in rows:
        if not all(str(r.get(k, "")).strip() for k in REQUIRED):
            raise NCATargetError("Missing required source NCA field or original measurement")
        endpoint = str(r["endpoint"])
        if endpoint not in UNITS:
            raise NCATargetError("Only Cmax/Tmax/AUC0last explicit targets accepted")
        if str(r["unit"]) != UNITS[endpoint]:
            raise NCATargetError("Original PK endpoint physical unit mismatch")
        sha = str(r["pk_source_sha256"])
        if len(sha) != 64 or any(c not in "0123456789abcdef" for c in sha):
            raise NCATargetError("Original PK dose/concentration source SHA256 required")
        sources.add(sha)
        sid = str(r["subject_id"])
        occ = str(r["occasion_id"])
        dose = str(r["dose_id"])
        drug = str(r["drug"])
        key = (sid, occ, dose, drug, str(r["analyte"]), endpoint)
        if key in seen:
            raise NCATargetError("Duplicated donor/occasion/endpoint is not a new independent PK measurement")
        seen.add(key)
        start = _number(r["dose_time_h"], "dose_time_h")
        stop = _number(r["window_end_h"], "window_end_h", strictly_positive=True)
        if stop <= start:
            raise NCATargetError("NCA time window cannot end before the dose")
        samples = _number(r["n_observed_concentrations"], "n_observed_concentrations", strictly_positive=True)
        if samples != int(samples) or samples < 3:
            raise NCATargetError("NCA target must be based on at least three original measured PK concentrations")
        value = _number(r["value"], "measured_endpoint_value")
        if endpoint == "Tmax" and value > stop - start:
            raise NCATargetError("Tmax exceeds original observed NCA window")
        if endpoint == "AUC0last" and value <= 0:
            raise NCATargetError("Nonpositive AUC0-last cannot be a measured exposure")
        if endpoint == "AUC0last" and str(r.get("integration_method", "")) not in ("linear_trapezoid", "linear_up_log_down"):
            raise NCATargetError("AUC requires an explicit predeclared integration method")
        period = (start, stop)
        intervals[(sid, occ, dose, drug, str(r["analyte"]))].add(period)
        subjects.add(sid)
        occasions.add((sid, occ))
        endpoints[endpoint] += 1
    if any(len(x) != 1 for x in intervals.values()):
        raise NCATargetError("Cmax, Tmax and AUC labels use incompatible PK time windows")
    return {
        "status": "MEASURED_NCA_TARGET_SCHEMA_ONLY_NOT_POPPK_PROOF",
        "n_unique_original_subjects":len(subjects),
        "n_independent_subject_occasions":len(occasions),
        "n_endpoint_rows":len(rows),
        "n_original_pk_sources":len(sources),
        "n_by_parameter":dict(endpoints),
        "no_PK_or_omics_measurements_exported":True,
        "requires_original_observed_dose_Ct_timecourse":True,
    }


def qualify_baseline_omic_nca_join(
    exchange_dir: str | Path,
    targets: Sequence[Mapping[str, object]],
    *,
    cohort_time_axis_confirmed: bool,
    modality: str = "transcriptomics",
) -> dict:
    """One original pre-dose molecular specimen and matched NCA, by subject.

    This does NOT create feature vectors or estimate any parametric effect.
    A result must NOT be called matched when the subject ID namespaces or
    enrollment anchors do not come from the same documented clinical study.
    """
    if not cohort_time_axis_confirmed:
        raise NCATargetError("Explicit clinical protocol crosswalk/shared time origin required")
    if modality not in ("transcriptomics","proteomics","metabolomics"):
        raise NCATargetError("Unknown molecular modality")
    inspect_exchange(exchange_dir)
    source = Path(exchange_dir).resolve()
    manifest = json.loads((source/"manifest.json").read_text(encoding="utf-8"))
    _, samples = _csv(_file(source, manifest["files"]["samples"]),
                      ["assay_id","sample_id","subject_id","omic","collection_time","collection_time_unit"])
    report = qualify_nca_targets(targets)
    per_subject = defaultdict(list)
    for sample in samples:
        if sample["omic"] != modality:
            continue
        hour = _temporal_hour(sample["collection_time"],sample["collection_time_unit"])
        if hour is None:
            raise NCATargetError("Unquantified or categorical MultiOmics time cannot join original PK chronology")
        per_subject[sample["subject_id"]].append((hour,sample["sample_id"]))
    n_linked = 0
    unmatched_subjects = set()
    n_ambiguous = 0
    has_same_time = False
    for target in targets:
        subject = str(target["subject_id"])
        start = float(target["dose_time_h"])
        # strictly before administered dose; same-time biomarker is not
        # deemed pre-dose without an original specimen ordering protocol.
        matches = [(t,sid) for t,sid in per_subject.get(subject,()) if t < start]
        if not matches:
            unmatched_subjects.add(subject)
            if any(t==start for t,_ in per_subject.get(subject,())):
                has_same_time = True
            continue
        recent = max(t for t,_ in matches)
        sample_ids = {sample for t,sample in matches if t==recent}
        if len(sample_ids) != 1:
            n_ambiguous += 1
        else:
            n_linked += 1
    return {
        "status": "BASELINE_MOLECULAR_NCA_LINK_CHECK_ONLY_NOT_AN_EFFECT_ESTIMATE",
        "n_nca_endpoint_rows":report["n_endpoint_rows"],
        "n_original_independent_subjects_with_nca":report["n_unique_original_subjects"],
        "n_nca_endpoints_linked_to_pre_dose_omic":n_linked,
        "n_subjects_without_pre_dose_omic":len(unmatched_subjects),
        "n_ambiguous_omics_sample_timestamp":n_ambiguous,
        "same_time_omic_without_order_abstained":has_same_time,
        "clinical_shared_time_origin_confirmed_by_caller":True,
        "group_effect_or_qvalue_is_not_an_individual_PK_covariate":True,
        "requires_subject_grouped_outer_validation":True,
        "no_individual_source_records_exported":True,
    }
