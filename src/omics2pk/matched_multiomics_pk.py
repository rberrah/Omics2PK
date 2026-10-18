"""Private in-memory MultiOmics exchange V1 -> original clinical RNA/PK ablation.

Never print patient/assay identifiers or per-subject predictions. All original
source clinical sidecars stay caller-provided and local/authorized. A user-provided
identity/time crosswalk is NOT independently source-confirmed by software.
"""
from __future__ import annotations

import json
from pathlib import Path

from .exchange import ExchangeError, _csv, _file, inspect_exchange
from .predose_representation_benchmark import compare_predose_rna_pk


class MatchedExchangePKError(ValueError):
    pass


def run_matched_rna_exchange_pk(
    exchange_dir, original_clinical_sidecar, *, target_name, target_unit,
    source_kind, source_sha256, prespecified_gene_weights=None
):
    """Run a single-assay-per-person RNA/PK study without saving person values.

    original_clinical_sidecar is a user-authorized, in-memory sequence of:
      subject_id, assay_id, collected_at (absolute), first_dose_at (absolute),
      measured_pk (actual in-vivo endpoint), clinical_values (numeric list),
      and optional genotype_dosages (genuine 0/1/2 variants).
    MultiOmics RNA expression is expected to be normalized and gene-wise.
    """
    root = Path(exchange_dir).resolve()
    report = inspect_exchange(root)  # SHA-256, file containment, subject/assay links
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest["study"]["biological_unit"] != "participant":
        raise MatchedExchangePKError("Patient PK comparison needs actual human participants")
    if source_kind == "original_source_declared" and report["validation"] != "qualified_for_declared_scope":
        raise MatchedExchangePKError("Original-data RNA export not qualified upstream")
    entries = [m for m in manifest["matrices"] if m["omic"] == "transcriptomics"]
    if len(entries) != 1:
        raise MatchedExchangePKError("Exactly one transcriptomics assay matrix required")
    entry = entries[0]
    vt = str(entry.get("value_type", "")).lower()
    if not (any(x in vt for x in ("normalized", "log2", "vst", "rlog", "tpm", "fpkm", "log-transformed"))
            and not ("raw" in vt and "count" in vt)):
        raise MatchedExchangePKError("Original raw counts/unspecified expression not fit-ready")

    _, sample_rows = _csv(_file(root, manifest["files"]["samples"]),
                          ["assay_id", "sample_id", "subject_id", "omic"])
    rna = [x for x in sample_rows if x["omic"] == "transcriptomics"]
    if not rna or len({r["subject_id"] for r in rna}) != len(rna):
        raise MatchedExchangePKError("Repeated RNA person; freeze visit/replicate selection first")
    if any(not r.get("tissue", "").strip() for r in rna):
        raise MatchedExchangePKError("Original RNA specimen tissue required")
    tissue = {r["tissue"] for r in rna}
    if len(tissue) != 1:
        raise MatchedExchangePKError("Mixed source tissues cannot be pooled")
    assays = {r["assay_id"]: r for r in rna}
    people = {r["subject_id"]: r for r in rna}
    if not isinstance(original_clinical_sidecar, (tuple, list)) or not original_clinical_sidecar:
        raise MatchedExchangePKError("Authorized original clinical sidecar required")
    clinical_by_person = {}
    clinical_assays = set()
    for original in original_clinical_sidecar:
        if not isinstance(original, dict):
            raise MatchedExchangePKError("Private source clinical record is not a mapping")
        person, assay = original.get("subject_id"), original.get("assay_id")
        if not isinstance(person, str) or not isinstance(assay, str):
            raise MatchedExchangePKError("Missing exact original patient and molecular assay join")
        if person in clinical_by_person or assay in clinical_assays:
            raise MatchedExchangePKError("Duplicate original biological person or assay")
        clinical_by_person[person] = original
        clinical_assays.add(assay)
    if set(clinical_by_person) != set(people) or clinical_assays != set(assays):
        raise MatchedExchangePKError("Clinical and molecular subjects or assay IDs do not match exactly")
    if any(people[p]["assay_id"] != c["assay_id"] for p, c in clinical_by_person.items()):
        raise MatchedExchangePKError("Crossed original sample-person identity")
    _, features = _csv(_file(root, manifest["files"]["features"]),
                       ["omic", "feature_id", "id_namespace"])
    rna_features = [f for f in features if f["omic"] == "transcriptomics"]
    if len({f["id_namespace"] for f in rna_features}) != 1:
        raise MatchedExchangePKError("Mixed RNA gene namespaces")
    header, rows = _csv(_file(root, entry["file"]), ["feature_id"])
    if header[0] != "feature_id" or set(header[1:]) != set(assays):
        raise MatchedExchangePKError("RNA matrix/sample assay columns differ")
    matrix_genes = [r["feature_id"] for r in rows]
    if set(matrix_genes) != {f["feature_id"] for f in rna_features}:
        raise MatchedExchangePKError("Transcript identity lost across source matrices")
    if any(v == "NA" for row in rows for a in assays for v in [row[a]]):
        raise MatchedExchangePKError("Missing omics values must be handled in outer train folds")
    # Canonical deterministic order is for reproducibility only, NOT a row-order
    # molecular-to-clinical join. All matches above are exact original source keys.
    ids = sorted(people)
    assay_order = [people[p]["assay_id"] for p in ids]
    matched = [clinical_by_person[p] for p in ids]
    clinical = [r.get("clinical_values") for r in matched]
    if any(not isinstance(x, (tuple, list)) for x in clinical):
        raise MatchedExchangePKError("Each person must have observed clinical baseline vector")
    genotype = [r.get("genotype_dosages") for r in matched]
    if all(x is None for x in genotype):
        genotype = None
    elif any(x is None or not isinstance(x, (tuple, list)) for x in genotype):
        raise MatchedExchangePKError("Do not compare incomplete or fabricated PGx")
    return compare_predose_rna_pk(
        subject_ids=ids,
        gene_ids=matrix_genes,
        expression=[[row[a] for a in assay_order] for row in rows],
        clinical=clinical,
        genotype_dosages=genotype,
        measured_pk=[r.get("measured_pk") for r in matched],
        collected_at=[r.get("collected_at") for r in matched],
        first_dose_at=[r.get("first_dose_at") for r in matched],
        tissue=next(iter(tissue)), target_name=target_name, target_unit=target_unit,
        source_kind=source_kind, source_sha256=source_sha256,
        prespecified_gene_weights=prespecified_gene_weights,
    )
