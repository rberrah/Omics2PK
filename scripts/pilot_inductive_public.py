#!/usr/bin/env python3
"""Public SWEET expression-only leakage and network-stability diagnostic.

No PK endpoints; no patient identifiers or individual expression exported.
Only small predeclared first-gene panel: gene selection is NOT evidence-based.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from omics2pk.inductive_sweet import (
    build_outer_split_features, fit_frozen_reference, predict_network,
)


def pilot(upstream: Path) -> dict:
    path = upstream / "example" / "expression.txt"
    if not path.is_file():
        raise ValueError(f"Missing official SWEET expression fixture: {path}")
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        if header[0] != "gene":
            raise ValueError("Expected upstream first column gene")
        samples = tuple(header[1:])
        dataset = list(reader)
    if len(samples) < 80 or len(dataset) < 25:
        raise ValueError("Unexpected upstream fixture size")
    # This fixed rule precedes all calculations. No held-out data or endpoint
    # is used for feature selection, though it is not biologically informed.
    all_genes = tuple(r[0] for r in dataset)
    all_values = np.array([[float(v) for v in r[1:]] for r in dataset])
    chosen = dataset[:25]
    genes = tuple(r[0] for r in chosen)
    values = all_values[:25, :]
    if not np.isfinite(values).all() or values.shape != (25, len(samples)):
        raise ValueError("Invalid published expression")
    if np.any(np.std(values[:, :60], axis=1) < 1e-12):
        raise ValueError("Frozen 25-gene panel contains constant training gene")
    if len(set(samples)) != len(samples):
        raise ValueError("Repeated sample in example")
    train, test = samples[:60], samples[60:]
    split = build_outer_split_features(
        values, genes=genes, subject_ids=samples,
        train_ids=train, test_ids=test, tissue="TCGA-archived-unspecified", folds=0,
        weight_expression=all_values, weight_genes=all_genes)
    ref = fit_frozen_reference(
        values[:, :60], genes=genes, subject_ids=train,
        tissue="TCGA-archived-unspecified",
        weight_expression=all_values[:, :60], weight_genes=all_genes)
    first_graph, a = predict_network(
        ref, values[:, 60], subject_id=test[0], genes=genes,
        tissue="TCGA-archived-unspecified",
        weight_query=all_values[:, 60], weight_genes=all_genes)
    # Mutating a DIFFERENT test subject must not change this subject's result.
    replaced = values.copy()
    replaced[:, -1] += np.linspace(-3, 4, 25)
    other_graph, b = predict_network(
        ref, replaced[:, 60], subject_id=test[0], genes=genes,
        tissue="TCGA-archived-unspecified",
        weight_query=all_values[:, 60], weight_genes=all_genes)
    if not np.array_equal(first_graph, other_graph):
        raise ValueError("Held-out/test contamination detected")
    if a["reference_fingerprint"] != b["reference_fingerprint"]:
        raise ValueError("Training reference drifted across held-out queries")
    if len(split["training"]) != 60 or len(split["test"]) != 20:
        raise ValueError("Unexpected cohort partition")
    # Official SWEET transductive contrast for the same patient, with ALL
    # 80 subjects in the reference. The difference is expected (not a failure).
    subject_corr = np.corrcoef(all_values.T)
    means = (subject_corr.sum(axis=1) - 1) / (len(samples) - 1)
    all_weights = ((means - means.min() + 0.01) /
                   (means.max() - means.min() + 0.01)) * 0.1 * len(samples)
    full_ref_corr = np.corrcoef(values)
    full_duplicate_corr = np.corrcoef(np.column_stack([values, values[:, 60]]))
    transductive = full_ref_corr + all_weights[60] * (
        full_duplicate_corr - full_ref_corr)
    np.fill_diagonal(transductive, 0.0)
    rel_difference = float(np.linalg.norm(first_graph - transductive) /
                           max(np.linalg.norm(transductive), 1e-12))
    clipped = sum(x["weight_clamped_outside_reference"] for x in split["test"])
    abs_strength = np.array([list(r["features"].values()) for r in split["test"]])
    return {
        "upstream": "https://github.com/SysMednet/SWEET",
        "dataset": "official 80-sample expression example; outcome-free",
        "panel_genes": len(genes),
        "sample_weight_genes": len(all_genes),
        "sample_weight_genes_include_panel": set(genes).issubset(set(all_genes)),
        "relative_frobenius_difference_vs_transductive_full_cohort": rel_difference,
        "training_subjects": 60,
        "test_subjects": 20,
        "train_reference_exclusion": "leave-one-out; each training sample excluded",
        "training_network_reference_n": 59,
        "test_network_reference_n": 60,
        "all_heldout_predictions_share_frozen_fingerprint": (
            len({r["reference_fingerprint"] for r in split["test"]}) == 1),
        "heldout_patient_invariant_to_other_heldout_profile": True,
        "clipped_test_weights": clipped,
        "test_abs_strength_min": float(np.min(abs_strength)),
        "test_abs_strength_max": float(np.max(abs_strength)),
        "test_abs_strength_mean": float(np.mean(abs_strength)),
        "result": "PASS_public_expression_only_no_PK_target",
        "not_validated": [
            "published SWEET equivalence",
            "genome-wide scale (100 official example genes, not full RNA-seq)",
            "cross-fold preprocessing learned on training only",
            "tissue-specific drug activity",
            "ML benefit", "PK/PD predictive benefit",
        ],
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--upstream", type=Path, required=True)
    p.add_argument("--report", type=Path, required=True)
    args = p.parse_args()
    result = pilot(args.upstream)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
