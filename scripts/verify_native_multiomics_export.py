#!/usr/bin/env python3
"""Validate a real JS-generated MultiOmics v1 exchange, synthetic demo only."""
import csv
import json
import sys
from pathlib import Path

from omics2pk.exchange import ExchangeError, inspect_exchange, prepare_sweet
from omics2pk.multiomics_pk_readiness import inspect_pk_readiness


def main(root: Path) -> dict:
    result = inspect_exchange(root)
    assert result["subjects"] == 8, result
    assert result["assays"] == 49, result
    assert set(result["matrices"]) == {
        "transcriptomics", "proteomics", "metabolomics"
    }, result
    assert result["matrices"]["transcriptomics"]["n_assays"] == 17
    assert result["matrices"]["proteomics"]["n_assays"] == 16
    assert result["matrices"]["metabolomics"]["n_assays"] == 16
    # The source exporter intentionally marks exchange as private by default;
    # this is valid even when the test fixture comes from publicly shipped demo.
    assert result["data_access"] == "private"
    manifest = json.loads((root / "manifest.json").read_text())
    assert manifest["software"]["commit"] == "unknown", (
        "If producer commit capture is implemented, update this explicit provenance gate"
    )
    with (root / "samples.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    rna = [r for r in rows if r["omic"] == "transcriptomics"]
    assert len(rna) == 17
    assert len({r["subject_id"] for r in rna}) == 8, (
        "Longitudinal public demo has only eight distinct subjects"
    )
    counts = {}
    for row in rna:
        counts[row["subject_id"]] = counts.get(row["subject_id"], 0) + 1
    assert sorted(counts.values()) == [2, 2, 2, 2, 2, 2, 2, 3], counts
    assert len({r["sample_id"] for r in rna}) == 16, (
        "One RNA sample must have a technical repeat assay"
    )
    # Never relabel original counts as SWEET-compatible normalized expression.
    # In the unchanged real producer package, prepare-sweet must refuse.
    blocked = False
    try:
        prepare_sweet(root, root.parent / "should-not-write-sweet")
    except ExchangeError as exc:
        blocked = "normalized/transformed" in str(exc)
    assert blocked, "Raw-count clinical/research input must be rejected for SWEET"
    assert not (root.parent / "should-not-write-sweet").exists()
    pk = inspect_pk_readiness(root)
    # Real current MultiOmics original demo labels samples "T0"/"T12"
    # in input 'timepoint', but its V1 producer writes only collection_time,
    # collection_time_unit; both are empty on the NATIVE exchange CSV.
    assert pk["n_independent_biological_subjects"] == 8
    assert pk["n_omic_assays_missing_collection_time"] == 49
    assert pk["n_subjects_with_two_or_more_distinct_numeric_collection_times"] == 0
    assert not pk["temporal_numeric_metadata_ready"]
    assert pk["n_original_unique_specimens"] == 16
    assert pk["n_extra_technical_repeat_assays_not_independent_biology"] == 1
    assert not pk["has_original_PK_dose_concentration_history"]
    assert not pk["has_structured_person_level_NCA_targets"]
    assert "PRESERVE_CATEGORICAL_TIMEPOINT_FROM_ORIGINAL_MULTIOMICS_SOURCE" in pk["missing_for_E3_and_E4"]
    report = {
        "result": "PASS_native_multiomics_three_omic_demo_contract",
        "source_data": "bundled synthetic MultiOmics samples only",
        "subjects": result["subjects"],
        "assays": result["assays"],
        "omics": sorted(result["matrices"]),
        "source_data_access": result["data_access"],
        "producer_commit": manifest["software"]["commit"],
        "longitudinal_rna_subjects": 8,
        "rna_technical_replicate_assays": 1,
        "sweet_raw_counts_guard": "PASS_BLOCKED",
        "real_native_exchange_PK_time_loss": "CONFIRMED_TIMEPOINT_FIELD_NOT_PROPAGATED",
        "real_native_original_omic_assays_with_missing_quantitative_time": pk["n_omic_assays_missing_collection_time"],
        "native_source_PK_E3_E4_status": pk["status"],
        "clinical_validation": False,
        "scientific_analysis": False,
        "no_patient_export": True,
    }
    print(json.dumps(report, sort_keys=True, indent=2))
    return report


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: verify_native_multiomics_export.py <real-export-dir>")
    main(Path(sys.argv[1]))
