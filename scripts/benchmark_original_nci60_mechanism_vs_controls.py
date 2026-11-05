#!/usr/bin/env python3
"""Original 47 NCI60 source-concordant RNA/doxorubicin GI50, research only.

Reuses original mixOmics + IntLIM 47-line source-parity R export. No source
rows or heldout values leave ephemeral runner. RNA evidence is not function.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from scripts.benchmark_public_nci60_twoomics import load_dataset
from omics2pk.mechanism_matched_control import benchmark_prior_against_controls


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--controls", type=int, default=64)
    a = p.parse_args()
    ids, classes, y, rna_names, _, rna, _ = load_dataset(a.data)
    if len(ids) != 47:
        raise ValueError("Original 47 concordant NCI-60 source lines required")
    result = benchmark_prior_against_controls(
        ids, classes, y, rna_names, rna, repeats=3, n_controls=a.controls)
    result["n_original_source_transcriptomic_features"] = len(rna_names)
    result["original_cell_line_source_identity_provenance"] = (
        "mixOmics doxorubicin GI50 plus IntLIM original exact line keys "
        "with MDA-MB-435 source histology mismatch excluded")
    result["source_rna_and_drug_assay_already_not_same_specimen"] = True
    a.out.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
