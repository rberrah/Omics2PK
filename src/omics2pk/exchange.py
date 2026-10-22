"""Strict offline reader for MultiOmics DOWNSTREAM_EXCHANGE_V1. Research only."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

OMICS = frozenset({"transcriptomics", "proteomics", "metabolomics"})
SCHEMA_VERSION = "1.0.0"


class ExchangeError(ValueError):
    """Invalid MultiOmics exchange; never repair errors silently."""


def _file(root: Path, record: dict[str, str]) -> Path:
    rel = record.get("path", "")
    if not isinstance(rel, str) or not rel or Path(rel).is_absolute():
        raise ExchangeError("File path must be nonempty and relative")
    path = (root / rel).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ExchangeError(f"Unsafe or missing file: {rel}")
    expected = record.get("sha256", "")
    if not isinstance(expected, str) or len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
        raise ExchangeError(f"Invalid SHA-256 declaration: {rel}")
    with path.open("rb") as stream:
        actual = hashlib.file_digest(stream, "sha256").hexdigest()
    if actual != expected:
        raise ExchangeError(f"Checksum mismatch: {rel}")
    return path


def _csv(path: Path, required: list[str]) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
            raise ExchangeError(f"Invalid/duplicate CSV headers: {path.name}")
        if any(c not in reader.fieldnames for c in required):
            raise ExchangeError(f"Missing CSV columns: {path.name}")
        rows = list(reader)
        if not rows or any(None in row for row in rows):
            raise ExchangeError(f"Empty or ragged CSV: {path.name}")
    return list(reader.fieldnames), rows


def inspect_exchange(directory: str | Path) -> dict[str, Any]:
    """Verify SHA, assay-subject linkage, namespace and numeric matrix values."""
    root = Path(directory).resolve()
    mf = root / "manifest.json"
    if not mf.is_file():
        raise ExchangeError("Expected exchange/manifest.json")
    manifest = json.loads(mf.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ExchangeError("Unsupported schema version")
    study = manifest.get("study", {})
    if study.get("data_access") not in {"public", "private", "controlled"}:
        raise ExchangeError("Missing data access classification")
    if study.get("biological_unit") not in {"participant", "animal", "culture", "organoid", "other"}:
        raise ExchangeError("Missing biological unit")
    if not manifest.get("software", {}).get("commit"):
        raise ExchangeError("Missing producer commit")
    status = manifest.get("validation", {}).get("status")
    if status not in {"exploratory", "reference_computed_unqualified", "qualified_for_declared_scope"}:
        raise ExchangeError("Invalid upstream validation status")

    files = manifest["files"]
    _, samples = _csv(_file(root, files["samples"]), ["assay_id", "sample_id", "subject_id", "omic"])
    _, features = _csv(_file(root, files["features"]), ["omic", "feature_id", "id_namespace"])
    assays: dict[tuple[str, str], dict[str, str]] = {}
    for row in samples:
        omic, assay = row["omic"], row["assay_id"]
        if omic not in OMICS or not assay or not row["sample_id"] or not row["subject_id"]:
            raise ExchangeError("Invalid assay/sample/subject/omic")
        key = (omic, assay)
        if key in assays:
            raise ExchangeError(f"Duplicate assay {key}")
        assays[key] = row
    feature_ids: set[tuple[str, str]] = set()
    for row in features:
        key = (row["omic"], row["feature_id"])
        if row["omic"] not in OMICS or not row["feature_id"] or not row["id_namespace"] or key in feature_ids:
            raise ExchangeError(f"Invalid/duplicate feature {key}")
        feature_ids.add(key)
    matrices: dict[str, dict] = {}
    entries = manifest.get("matrices", [])
    if not isinstance(entries, list) or not entries:
        raise ExchangeError("No matrices")
    for entry in entries:
        omic = entry.get("omic")
        if omic not in OMICS or omic in matrices or entry.get("orientation") != "feature_rows_sample_columns":
            raise ExchangeError("Invalid/duplicate omic matrix or orientation")
        if not entry.get("value_type"):
            raise ExchangeError("Missing value type")
        header, rows = _csv(_file(root, entry["file"]), ["feature_id"])
        if header[0] != "feature_id" or len(header) < 2:
            raise ExchangeError("Matrix requires feature_id and assays")
        assay_ids = header[1:]
        if set(assay_ids) != {a for o, a in assays if o == omic}:
            raise ExchangeError(f"Sample/assay mismatch in {omic}")
        seen = set()
        for row in rows:
            f = row["feature_id"]
            if not f or f in seen or (omic, f) not in feature_ids:
                raise ExchangeError(f"Unknown/duplicate feature {omic}/{f}")
            seen.add(f)
            for assay in assay_ids:
                val = row[assay]
                if val == "NA":
                    continue
                try:
                    num = float(val)
                except (ValueError, TypeError):
                    raise ExchangeError(f"Non-numeric value {omic}/{f}/{assay}") from None
                if not math.isfinite(num):
                    raise ExchangeError("Infinite/NaN matrix value")
        if seen != {f for o, f in feature_ids if o == omic}:
            raise ExchangeError(f"Feature/matrix mismatch: {omic}")
        matrices[omic] = {"n_features": len(rows), "n_assays": len(assay_ids), "value_type": entry["value_type"]}
    if {r["omic"] for r in samples} != set(matrices) or {r["omic"] for r in features} != set(matrices):
        raise ExchangeError("Inconsistent omic coverage")
    for lock in manifest.get("software", {}).get("dependency_lockfiles", []):
        _file(root, lock)
    return {"schema_version": SCHEMA_VERSION, "analysis_id": manifest.get("analysis_id"),
            "data_access": study["data_access"], "validation": status,
            "subjects": len({r["subject_id"] for r in samples}), "assays": len(samples),
            "matrices": matrices, "clinical_use": False, "network_inferred": False}


def _sweet_expression_preflight(
    rows: list[dict[str, str]], assays: list[str]
) -> dict[str, Any]:
    """Reject expression shapes that pinned upstream SWEET handles unsafely.

    The original Step 2 drops a gene when its *sum* is exactly zero rather
    than when all values are zero, while NumPy correlation returns non-finite
    values for constant gene/sample profiles.  Omics2PK fails before writing
    an upstream input instead of silently accepting these cases.
    """
    if not rows or not assays:
        raise ExchangeError("SWEET expression preflight requires genes and assays")
    matrix: list[list[float]] = []
    gene_ids: list[str] = []
    for row in rows:
        gene_id=row["feature_id"]
        values=[]
        for assay in assays:
            raw=row[assay]
            if raw=="NA" or raw.strip()=="":
                raise ExchangeError("SWEET input cannot contain missing values; handle missingness upstream")
            try:
                value=float(raw)
            except (TypeError,ValueError):
                raise ExchangeError(f"SWEET input contains non-numeric expression: {gene_id}/{assay}") from None
            if not math.isfinite(value):
                raise ExchangeError("SWEET input cannot contain NaN or infinite expression")
            values.append(value)
        spread=max(values)-min(values)
        if spread<=1e-12:
            raise ExchangeError(
                f"SWEET upstream correlation is undefined for constant gene: {gene_id}")
        magnitude=max(1.0,math.fsum(abs(x) for x in values))
        if abs(math.fsum(values))<=1e-12*magnitude:
            raise ExchangeError(
                f"SWEET upstream Step 2 deletes zero-sum gene even when nonconstant: {gene_id}")
        matrix.append(values)
        gene_ids.append(gene_id)
    for col,assay in enumerate(assays):
        values=[row[col] for row in matrix]
        if max(values)-min(values)<=1e-12:
            raise ExchangeError(
                f"SWEET sample-weight correlation is undefined for constant sample: {assay}")
    return {
        "finite_expression":True,
        "constant_genes":0,
        "constant_samples":0,
        "zero_sum_nonconstant_genes":0,
        "n_genes_checked":len(gene_ids),
        "n_samples_checked":len(assays),
        "upstream_missing_value_coercion_used":False,
    }


def prepare_sweet(directory: str | Path, output_dir: str | Path) -> dict[str, Any]:
    """Create official SWEET input format; does NOT infer networks.

    One assay per subject and tissue; callers must separately validate an
    inductive SWEET method before prospective individual prediction.
    """
    report = inspect_exchange(directory)
    if "transcriptomics" not in report["matrices"]:
        raise ExchangeError("SWEET needs transcriptomics")
    root = Path(directory).resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    entry = next(e for e in manifest["matrices"] if e["omic"] == "transcriptomics")
    vt = entry["value_type"].strip().lower()
    if not any(s in vt for s in ("normalized", "log2", "log-transformed", "vst", "rlog", "tpm", "fpkm")) or ("raw" in vt and "count" in vt):
        raise ExchangeError(f"SWEET requires normalized/transformed expression, not {entry['value_type']}")
    _, all_samples = _csv(_file(root, manifest["files"]["samples"]), ["assay_id", "sample_id", "subject_id", "omic"])
    samples = [s for s in all_samples if s["omic"] == "transcriptomics"]
    if len({s["subject_id"] for s in samples}) != len(samples):
        raise ExchangeError("Repeated subject: define explicit design or aggregation before SWEET")
    if len({s.get("tissue", "") for s in samples}) > 1:
        raise ExchangeError("Mixed tissues need explicit strata")
    _, features = _csv(_file(root, manifest["files"]["features"]), ["omic", "feature_id", "id_namespace"])
    genes = [f for f in features if f["omic"] == "transcriptomics"]
    if len({f["id_namespace"] for f in genes}) != 1:
        raise ExchangeError("Mixed gene ID namespaces require explicit mapping")
    _, rows = _csv(_file(root, entry["file"]), ["feature_id"])
    assays = [s["assay_id"] for s in samples]
    preflight=_sweet_expression_preflight(rows,assays)
    out = Path(output_dir).resolve()
    if out == root or root in out.parents:
        raise ExchangeError("Output must be outside source exchange")
    out.mkdir(parents=True, exist_ok=True)
    with (out / "expression.txt").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene"] + assays)
        for row in rows:
            writer.writerow([row["feature_id"]] + [row[a] for a in assays])
    (out / "patient.txt").write_text("\n".join(assays) + "\n", encoding="utf-8")
    (out / "gene.txt").write_text("\n".join(r["feature_id"] for r in rows) + "\n", encoding="utf-8")
    return {"expression_file": str(out / "expression.txt"), "sample_count": len(assays),
            "gene_count": len(rows), "sweet_input_preflight":preflight,
            "requires_independent_sweet_qualification": True,
            "message": "Inputs only. SWEET not run; upstream edge-case and out-of-sample qualification remain separate."}
