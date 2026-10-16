#!/usr/bin/env python3
"""Validate a real JS-generated MultiOmics v1 exchange, synthetic demo only."""
import csv
import json
import sys
from pathlib import Path

from omics2pk.exchange import ExchangeError, inspect_exchange, prepare_sweet


def main(root: Path) -> dict:
    result = inspect_exchange(root)
    assert result["subjects"] == 16, result
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
    assert len({r["subject_id"] for r in rna}) == 16, (
        "Public demo includes a repeat RNA assay for one subject"
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
    report = {
        "result": "PASS_native_multiomics_three_omic_demo_contract",
        "source_data": "bundled synthetic MultiOmics samples only",
        "subjects": result["subjects"],
        "assays": result["assays"],
        "omics": sorted(result["matrices"]),
        "source_data_access": result["data_access"],
        "producer_commit": manifest["software"]["commit"],
        "rna_repeated_subject": True,
        "sweet_raw_counts_guard": "PASS_BLOCKED",
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
