"""Privacy-preserving qualification of locally supplied original rosuvastatin data.

Scientific source: Silveira et al, Front Pharmacol 2021,
DOI 10.3389/fphar.2021.752960, MTBLS3128, 40 biological human volunteers.

No original person/peak/PK rows bundled, uploaded to Actions or committed.
No row-order joins, automatic guessing of person-ID columns, hidden imputation,
PK-target-as-feature leakage, or retrospective outcome reconstructions.

A positive result establishes ONLY a *structural* 40-person sample↔PK join.
Original authorship, specimen predose time, units/assay quality and manuscript
source must be independently verified before any clinical E3 model fitting.
"""
from __future__ import annotations

import csv
import hashlib
import io
import math
from pathlib import Path
import re

ORIGINAL_DOI = "10.3389/fphar.2021.752960"
ORIGINAL_MTBLS = "MTBLS3128"
EXPECTED_HUMANS = 40
MAX_BYTES = 64_000_000
MAX_COLUMNS = 40_000
MAX_ROWS = 100_000
PK_OUTCOME_TOKEN = re.compile(
    r"(^|[\W_])(?:"
    r"auc|auc0(?:to|t|\d)*|area\s*under\s*(?:the\s*)?curve|"
    r"cmax|c[\W_]*max|pk\s*(?:target|response|outcome|label)|"
    r"trough|cmin|tmax|clearance|"
    r"rosuvastatin[\W_]*(?:plasma|conc|concentration|auc|cmax)"
    r")($|[\W_])",
    re.IGNORECASE,
)
EMPTY_TOKENS = frozenset({"", "NA", "N/A", "NAN", "NULL", "-"})
SAMPLE_METADATA_TOKENS = (
    "subject", "participant", "volunteer", "patient", "sample", "donor", 
    "batch", "run", "retention", "time", "file", "label", "tissue", 
    "group", "gender", "sex", "date", "ionization", "polarity",
)


class SourceGateError(ValueError):
    """Only source-contract descriptions; never includes source cell values."""


def _source_error(message: str):
    raise SourceGateError(message)


def original_csv(path: str | Path) -> tuple[list[str], list[dict[str, str]], str]:
    source = Path(path)
    if not source.is_file():
        _source_error("Original source CSV not found")
    if source.stat().st_size > MAX_BYTES:
        _source_error("Source exceeds private local CSV size cap; do not import mzML as CSV")
    payload = source.read_bytes()
    if b"\x00" in payload[:1000]:
        _source_error("Binary content, not a source CSV")
    try:
        text = payload.decode("utf-8-sig")
    except UnicodeDecodeError:
        _source_error("Original CSV must be UTF-8; recode privately, preserve source copy")
    if text.lstrip().lower().startswith(("<!doctype", "<html")):
        _source_error("Google Drive returned HTML rather than an original data CSV")
    try:
        d = csv.Sniffer().sniff(text[:16000], delimiters=",;\t")
    except csv.Error:
        d = csv.excel
    reader = csv.DictReader(io.StringIO(text), dialect=d)
    headers = reader.fieldnames or []
    if not headers or len(headers) > MAX_COLUMNS:
        _source_error("Invalid or excessively wide original CSV header")
    if len(headers) != len(set(headers)) or any(not h.strip() for h in headers):
        _source_error("Duplicate or blank original CSV column headers")
    records = []
    for row in reader:
        if len(records) >= MAX_ROWS:
            _source_error("Too many source rows")
        if None in row:
            _source_error("Malformed original CSV row width: no silent column truncation")
        records.append(row)
    if not records:
        _source_error("No original data rows")
    return headers, records, hashlib.sha256(payload).hexdigest()


def _ids(rows: list[dict], name: str, kind: str) -> set[str]:
    values = [str(row[name]).strip() for row in rows]
    if not all(values):
        _source_error(f"Absent original {kind} ID in some original source rows")
    if len(values) != len(set(values)):
        _source_error(f"Duplicate original {kind} ID: independent humans must not be pseudoreplicated")
    return set(values)


def _all_numeric(values, reason: str) -> int:
    n = 0
    for v in values:
        s = str(v).strip()
        if s.upper() in EMPTY_TOKENS:
            _source_error(f"Missing original {reason}; no original data imputation")
        try:
            a = float(s)
        except ValueError:
            _source_error(f"Non-numeric original {reason}; no synthetic target")
        if not math.isfinite(a) or a <= 0.0:
            _source_error(f"Original {reason} must be finite and positive")
        n += 1
    return n


def summarize_local_sources(paths: list[str | Path]) -> dict:
    """Summarize real uploaded CSV dimensions; never print names, headers or data."""
    files = []
    for i, p in enumerate(paths, 1):
        headers, rows, digest = original_csv(p)
        normalized = [" ".join(h.casefold().replace("_", " ").split()) for h in headers]
        files.append({
            "file_position": i,
            "source_sha256": digest,
            "n_data_rows": len(rows),
            "n_columns": len(headers),
            "headers_suggest_auc": any(PK_OUTCOME_TOKEN.search(x) and "auc" in x for x in normalized),
            "headers_suggest_cmax": any(PK_OUTCOME_TOKEN.search(x) and ("cmax" in x or "c max" in x) for x in normalized),
            "status": "LOCAL_STRUCTURE_ONLY_NO_PERSON_JOIN_CLAIM",
            "individual_values_or_original_header_exported": False,
        })
    return {
        "study": ORIGINAL_MTBLS,
        "article": ORIGINAL_DOI,
        "n_original_source_files": len(files),
        "files": files,
        "qualified_original_PK_omic_person_join": False,
        "clinical_model_evaluated": False,
    }


def qualify_exact_subject_link(
    omics_path: str | Path,
    pk_path: str | Path,
    *,
    omics_id_col: str,
    pk_id_col: str,
    auc_col: str,
    cmax_col: str,
    expected_subjects: int = EXPECTED_HUMANS,
    author_crosswalk_path: str | Path | None = None,
    crosswalk_sample_col: str = "omics_sample_id",
    crosswalk_person_col: str = "pk_subject_id",
) -> dict:
    """Fail-closed strictly matched original data *contract* — never model fit.

    Original study source rights, pre-dose timing and PK numerical assay must
    be verified independently; this function can confirm only join shape.
    """
    if not 5 <= expected_subjects <= 10_000:
        _source_error("Expected biological person count must be an explicit sensible number")
    omics_cols, omics, ohash = original_csv(omics_path)
    pk_cols, pk, phash = original_csv(pk_path)
    if not (omics_id_col and pk_id_col and auc_col and cmax_col):
        _source_error("Explicit original sample/person/PK column names required")
    if omics_id_col not in omics_cols or pk_id_col not in pk_cols:
        _source_error("Exact original ID column missing; never infer ID from order")
    if not {auc_col, cmax_col}.issubset(pk_cols):
        _source_error("Quantitative source AUC and Cmax columns both required")
    # Hard-block model leakage: PK outcomes, rosuvastatin concentration
    # and other postdose descriptors must NEVER enter predose omics features.
    for col in omics_cols:
        if col != omics_id_col and PK_OUTCOME_TOKEN.search(col):
            _source_error("PK outcome or postdose concentration embedded among omics features")
    omics_ids = _ids(omics, omics_id_col, "omics sample")
    pk_ids = _ids(pk, pk_id_col, "PK human")
    _all_numeric((row[auc_col] for row in pk), "source quantitative AUC")
    _all_numeric((row[cmax_col] for row in pk), "source quantitative Cmax")
    if len(pk_ids) != expected_subjects:
        _source_error("Published original biological-person count does not match PK labels")

    source_crosswalk_hash = None
    if author_crosswalk_path is not None:
        map_cols, pairs, mhash = original_csv(author_crosswalk_path)
        if not {crosswalk_sample_col, crosswalk_person_col}.issubset(map_cols):
            _source_error("Author sample↔person crosswalk columns missing")
        mapping = {}
        for row in pairs:
            original = str(row[crosswalk_sample_col]).strip()
            person = str(row[crosswalk_person_col]).strip()
            if not original or not person or original in mapping:
                _source_error("Crosswalk contains blank or repeated original sample key")
            mapping[original] = person
        if set(mapping) != omics_ids or len(set(mapping.values())) != len(mapping):
            _source_error("Crosswalk must be complete and 1:1: no inferred row-order join")
        source_crosswalk_hash = mhash
        mapped = set(mapping.values())
    else:
        mapped = omics_ids

    if len(mapped) != expected_subjects or mapped != pk_ids:
        _source_error("Original biological donor IDs are not an exact strict 40-person match")
    # Require measured original numeric features rather than empty MAF annotations.
    exclude = {omics_id_col}
    candidates = [
        col for col in omics_cols
        if col not in exclude and
        not any(tok in col.casefold() for tok in SAMPLE_METADATA_TOKENS)
    ]
    usable = []
    for feature in candidates:
        n_numeric = 0
        n_filled = 0
        for row in omics:
            value = str(row[feature]).strip()
            if value.upper() in EMPTY_TOKENS:
                continue
            n_filled += 1
            try:
                z = float(value)
            except ValueError:
                continue
            if math.isfinite(z):
                n_numeric += 1
        if n_numeric >= max(5, int(expected_subjects * 0.8)) and n_numeric == n_filled:
            usable.append(feature)
    if len(usable) < 2:
        _source_error("Source does not contain enough original quantitative omics features")
    if Path(omics_path).resolve() == Path(pk_path).resolve():
        _source_error("Independent original PK and omics tables required")
    return {
        "study": ORIGINAL_MTBLS,
        "article_doi": ORIGINAL_DOI,
        "source_status": "STRICT_ORIGINAL_40_PERSON_STRUCTURE_JOIN_ONLY",
        "n_biological_subjects": expected_subjects,
        "n_original_omics_sample_rows": len(omics),
        "n_independent_PK_subjects": len(pk),
        "n_usable_quantitative_omic_features": len(usable),
        "n_original_PK_auc_and_cmax_both_numeric": len(pk),
        "source_SHA256_omics": ohash,
        "source_SHA256_PK": phash,
        "source_SHA256_author_crosswalk_if_supplied": source_crosswalk_hash,
        "original_source_provenance_independently_authenticated": False,
        "predose_timing_independently_verified": False,
        "quantitative_PK_units_and_assay_independently_verified": False,
        "model_fit_performed": False,
        "heldout_performance_claim": False,
        "source_patient_rows_exported": False,
        "next_step": "Human-verify original author source, pre-dose time, PK units and IDs before nested biological-subject CV",
    }


def qualify_combined_source(
    combined_path: str | Path,
    *,
    person_id_col: str,
    auc_col: str,
    cmax_col: str,
    expected_subjects: int = EXPECTED_HUMANS,
) -> dict:
    """For original authors' *one combined CSV* with features and true AUC/Cmax.

    Explicit exclusion of the two PK targets and every PK-related column from
    molecular predictors is required before any nested model validation.
    Proves only dimensional integrity; neither study authenticity, timing nor
    numerical target provenance can be established solely from CSV contents.
    """
    if not 5 <= expected_subjects <= 10_000:
        _source_error("Expected biological person count must be an explicit sensible number")
    headers, rows, digest = original_csv(combined_path)
    if not person_id_col or person_id_col not in headers:
        _source_error("Explicit original one-human-per-row ID required for combined source")
    if not auc_col or not cmax_col or auc_col == cmax_col:
        _source_error("Require two distinct actual quantitative PK AUC/Cmax columns")
    if not {auc_col,cmax_col}.issubset(headers):
        _source_error("Original combined table does not contain explicit PK targets")
    humans = _ids(rows, person_id_col, "combined original human")
    if len(humans) != expected_subjects:
        _source_error("Original combined table is not one row per expected human")
    _all_numeric((r[auc_col] for r in rows), "original quantitative rosuvastatin AUC")
    _all_numeric((r[cmax_col] for r in rows), "original quantitative rosuvastatin Cmax")
    excluded_pk_columns = [
        c for c in headers if c not in (auc_col,cmax_col,person_id_col)
        and PK_OUTCOME_TOKEN.search(c)
    ]
    if excluded_pk_columns:
        _source_error("Combined dataset has extra potential postdose PK columns; manually define a predose-only molecular view before fitting")
    usable_features = 0
    for col in headers:
        if col in (auc_col,cmax_col,person_id_col):
            continue
        if any(marker in col.lower() for marker in SAMPLE_METADATA_TOKENS):
            continue
        vals = [str(r[col]).strip() for r in rows]
        numeric = 0
        nonmissing = 0
        for value in vals:
            if value.upper() in EMPTY_TOKENS:
                continue
            nonmissing += 1
            try:
                if math.isfinite(float(value)):
                    numeric += 1
            except ValueError:
                pass
        if numeric >= max(5,int(expected_subjects*.8)) and numeric==nonmissing:
            usable_features += 1
    if usable_features < 2:
        _source_error("Combined source contains no substantive original quantitative molecular feature matrix")
    return {
        "study": ORIGINAL_MTBLS,
        "article": ORIGINAL_DOI,
        "source_status": "COMBINED_ONE_ROW_PER_HUMAN_AND_TWO_ORIGINAL_PK_TARGETS_STRUCTURE_ONLY",
        "n_unique_original_biological_humans": len(humans),
        "n_numeric_auc": len(rows),
        "n_numeric_cmax": len(rows),
        "n_candidate_numeric_features_excluding_explicit_PK_targets": usable_features,
        "original_combined_sha256": digest,
        "original_combined_source_authorship_verified": False,
        "raw_omics_feature_specimen_predose_verified": False,
        "original_AUC_Cmax_PK_target_assay_independently_verified": False,
        "cross_modal_ISA_sample_ID_link_verified": False,
        "target_columns_explicitly_excluded_from_feature_candidates": True,
        "source_patients_or_PK_values_exported": False,
        "model_fit_or_predictive_performance_claim": False,
        "next_step": "Independently verify published source provenance, predose metabolomics, rosuvastatin quantified AUC/Cmax and true person/sample ID crosswalk",
    }
