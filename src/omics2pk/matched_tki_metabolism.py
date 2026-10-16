"""Local-only, source-gated same-donor two-drug CYP→metabolism benchmark.

Tucatinib DOI 10.1016/j.dmd.2025.100061 and imatinib
DOI 10.1016/j.dmd.2025.100196 use the SAME original 21 HLM lots.
No published donor-level data are bundled or uploaded to GitHub.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
from statistics import mean

from omics2pk.tucatinib_local_validation import (
    qualify as qualify_tucatinib,
    verify_published_table_correlations, _spearman, _predict_ridge,
)


class MatchedMetabolismError(ValueError):
    pass


IMATINIB_DOI = "10.1016/j.dmd.2025.100196"
EXPECTED_IMATINIB_CORRELATIONS = {
    "cyp2c8_protein": 0.909,
    "cyp3a_protein": 0.902,
    "cyp2c8_activity": 0.664,
    "cyp3a_activity": 0.863,
}
FIELDS = (
    "donor_id", "row_kind", "source_doi",
    "imatinib_cl_app_ul_min_mg", "ndesmethyl_imatinib_peak_area_ratio",
    "cyp2c8_protein_pmol_mg", "cyp3a_protein_pmol_mg",
    "cyp2c8_activity_pmol_min_mg", "cyp3a_activity_pmol_min_mg",
)


def _strict_positive(value, label):
    try:
        v = float(value)
    except (ValueError, TypeError):
        raise MatchedMetabolismError(f"{label} must be a numerical observation") from None
    if not math.isfinite(v) or v <= 0:
        raise MatchedMetabolismError(f"{label} must be positive and finite")
    return v


def qualify_pair(tucatinib_rows, imatinib_rows):
    """Return only validated donor-indexed tuples; never match rows by order."""
    tuc = qualify_tucatinib(tucatinib_rows)
    verify_published_table_correlations(tuc)
    if len(imatinib_rows) != 22:
        raise MatchedMetabolismError("Expected 21 original HLM individual lots plus 150-donor pool")
    imatinib = {}
    pools = 0
    for row in imatinib_rows:
        if not set(FIELDS).issubset(row):
            raise MatchedMetabolismError("Missing required original imatinib source columns")
        if str(row["source_doi"]).strip() != IMATINIB_DOI:
            raise MatchedMetabolismError("Imatinib provenance mismatch")
        kind = str(row["row_kind"]).strip().lower()
        if kind == "pool":
            pools += 1
            continue
        if kind != "individual":
            raise MatchedMetabolismError("Explicit individual/pool label required")
        donor = str(row["donor_id"]).strip()
        if not donor or donor in imatinib:
            raise MatchedMetabolismError("Duplicate or missing original imatinib donor ID")
        predictors = (
            _strict_positive(row["cyp2c8_protein_pmol_mg"], "CYP2C8 protein"),
            _strict_positive(row["cyp3a_protein_pmol_mg"], "CYP3A protein"),
            _strict_positive(row["cyp2c8_activity_pmol_min_mg"], "CYP2C8 activity"),
            _strict_positive(row["cyp3a_activity_pmol_min_mg"], "CYP3A activity"),
        )
        cl_raw = str(row["imatinib_cl_app_ul_min_mg"]).strip()
        cl = None if cl_raw in ("", "-") else _strict_positive(cl_raw, "Imatinib CL")
        metabolite = _strict_positive(row["ndesmethyl_imatinib_peak_area_ratio"], "N-desmethyl peak-area ratio")
        imatinib[donor] = (cl, metabolite, predictors)
    if pools != 1 or len(imatinib) != 21:
        raise MatchedMetabolismError("Imatinib has 21 genuine donors and one excluded pooled lot")
    if set(imatinib) != {r[0] for r in tuc}:
        raise MatchedMetabolismError("Donor IDs differ across the two source publications: no fuzzy join")
    if sum(v[0] is None for v in imatinib.values()) != 4:
        raise MatchedMetabolismError("Four original imatinib clearances must remain missing, not zero")
    paired = []
    for donor, tuc_cl, tuc_prot, tuc_activity in tuc:
        im_cl, formation, other = imatinib[donor]
        original = (*tuc_prot, *tuc_activity)
        if any(abs(x - y) > 0.15 for x, y in zip(original, other)):
            raise MatchedMetabolismError("Cross-study donor CYP phenotypes disagree: no join")
        paired.append((donor, tuc_cl, im_cl, formation, tuc_prot, tuc_activity))
    matched = [r for r in paired if r[2] is not None]
    target = [r[2] for r in matched]
    checks = {
        "cyp2c8_protein": _spearman([r[4][0] for r in matched], target),
        "cyp3a_protein": _spearman([r[4][1] for r in matched], target),
        "cyp2c8_activity": _spearman([r[5][0] for r in matched], target),
        "cyp3a_activity": _spearman([r[5][1] for r in matched], target),
    }
    if any(abs(checks[k] - published) > 0.014
           for k, published in EXPECTED_IMATINIB_CORRELATIONS.items()):
        raise MatchedMetabolismError("Independent published imatinib Spearman validation failed")
    return paired, {k: round(v, 6) for k, v in checks.items()}


def _lodo(records, outcome_index):
    """Fixed, pre-specified same-donor LODO without hidden imputation."""
    selected = [r for r in records if r[outcome_index] is not None]
    if len(selected) < 12:
        raise MatchedMetabolismError("Not enough original measured donors")
    errors = {"baseline": [], "protein": [], "activity": []}
    for i, held in enumerate(selected):
        training = [x for j, x in enumerate(selected) if i != j]
        y = held[outcome_index]
        predictions = {
            "baseline": mean(x[outcome_index] for x in training),
            "protein": _predict_ridge([(x[4], x[outcome_index]) for x in training], held[4]),
            "activity": _predict_ridge([(x[5], x[outcome_index]) for x in training], held[5]),
        }
        for name, pred in predictions.items():
            errors[name].append(abs(pred - y))
    return {
        "n_independent_biological_donors": len(selected),
        "heldout_donor_MAE": {model: round(mean(err), 6) for model, err in errors.items()},
    }


def validate_pair(tucatinib_rows, imatinib_rows):
    paired, spearman = qualify_pair(tucatinib_rows, imatinib_rows)
    return {
        "status": "SOURCE_VERIFIED_TWO_DRUG_SHARED_DONORS_IN_VITRO_NOT_PATIENT_PK",
        "original_donors": 21,
        "imatinib_missing_clearance_not_imputed": 4,
        "pool_excluded": True,
        "source_dois": ["10.1016/j.dmd.2025.100061", IMATINIB_DOI],
        "imatinib_published_spearman_replicated": spearman,
        "tucatinib_cl_ul_min_mg": _lodo(paired, 1),
        "imatinib_observed_cl_ul_min_mg": _lodo(paired, 2),
        "ndesmethyl_imatinib_peak_area_ratio": _lodo(paired, 3),
        "models": "train-only mean vs fixed alpha1 two-CYP protein ridge vs enzyme-activity ridge",
        "no_cross_study_independence": True,
        "no_clinical_PopPK_or_PK_concentrations": True,
        "no_original_donor_measurements_exported": True,
    }


def from_local_csv(tucatinib_path, imatinib_path):
    def load(path):
        raw = Path(path).read_bytes()
        if len(raw) > 500_000:
            raise MatchedMetabolismError("Unexpected original donor table size")
        return list(csv.DictReader(raw.decode("utf-8-sig").splitlines())), hashlib.sha256(raw).hexdigest()
    t, th = load(tucatinib_path)
    i, ih = load(imatinib_path)
    return {**validate_pair(t, i), "original_local_csv_sha256": {
        "tucatinib": th, "imatinib": ih}}
