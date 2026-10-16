"""Strict research-only reader of raw SWEET edge lists.

This module DOES NOT infer SWEET networks and does not interpret co-expression
as a molecular activity or causal biological mechanism.
"""
from __future__ import annotations

import csv
import math
import re
import statistics
from pathlib import Path

MAX_EDGES = 200_000  # POC memory/CPU safety bound, not genome-wide qualification
SAFE_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


class SWEETOutputError(ValueError):
    """Fail closed on corrupt or incomplete network data."""


def _read_edges(path: Path):
    if not path.is_file():
        raise SWEETOutputError(f"Missing edge file: {path.name}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene1", "gene2", "raw_edge_score"]:
            raise SWEETOutputError(f"Unexpected SWEET edge schema: {path.name}")
        seen = set()
        nodes = set()
        totals = {}
        signed = {}
        weights = []
        for record in reader:
            a = record.get("gene1", "")
            b = record.get("gene2", "")
            if not a or not b or a == b:
                raise SWEETOutputError(f"Invalid gene pair in {path.name}")
            key = tuple(sorted((a, b)))
            if key in seen:
                raise SWEETOutputError(f"Repeated edge in {path.name}: {a}/{b}")
            seen.add(key)
            if len(seen) > MAX_EDGES:
                raise SWEETOutputError(
                    f"POC edge safety cap {MAX_EDGES} exceeded: select a fixed subset of genes")
            try:
                weight = float(record["raw_edge_score"])
            except (TypeError, ValueError) as err:
                raise SWEETOutputError(f"Non-numeric edge in {path.name}") from err
            if not math.isfinite(weight):
                raise SWEETOutputError(f"Non-finite SWEET edge in {path.name}")
            nodes.update((a, b))
            for g in (a, b):
                totals[g] = totals.get(g, 0.0) + abs(weight)
                signed[g] = signed.get(g, 0.0) + weight
            weights.append(weight)
    n = len(nodes)
    if n < 2 or len(seen) != n * (n - 1) // 2:
        raise SWEETOutputError(f"Incomplete non-sparse network: {path.name}")
    return {
        "genes": sorted(nodes),
        "n_edges": len(seen),
        "edge_mean": statistics.fmean(weights),
        "absolute_edge_mean": statistics.fmean(abs(w) for w in weights),
        "negative_edge_fraction": sum(w < 0 for w in weights) / len(weights),
        "edge_sd": statistics.pstdev(weights),
        "abs_strength": {g: totals[g] / (n - 1) for g in sorted(nodes)},
        "signed_strength": {g: signed[g] / (n - 1) for g in sorted(nodes)},
    }


def summarize_networks(edges_dir: str | Path, patient_list: str | Path,
                       output_dir: str | Path) -> dict:
    """Create tabular network features for M1bis, strictly local.

    The original SWEET edge files must already have been produced. This does
    not validate whether their reference cohort included the held-out subject.
    """
    folder = Path(edges_dir).resolve()
    patients_path = Path(patient_list).resolve()
    out = Path(output_dir).resolve()
    if out == folder or folder in out.parents:
        raise SWEETOutputError("Output must not overwrite the SWEET edge directory")
    if not patients_path.is_file():
        raise SWEETOutputError("Missing SWEET patient list")
    ids = [x.strip() for x in patients_path.read_text(encoding="utf-8").splitlines()]
    if not ids or any(not SAFE_ID.fullmatch(x) or x in {".", ".."} for x in ids):
        raise SWEETOutputError("Invalid or unsafe SWEET patient identifier")
    if len(set(ids)) != len(ids):
        raise SWEETOutputError("Duplicate patient identifier")
    reports = []
    genes = None
    for patient_id in ids:
        record = _read_edges(folder / f"{patient_id}.txt")
        if genes is None:
            genes = record["genes"]
        elif genes != record["genes"]:
            raise SWEETOutputError("Different gene universes across patient networks")
        reports.append((patient_id, record))
    assert genes is not None
    out.mkdir(parents=True, exist_ok=True)
    table = [
        ["patient_id", "n_genes", "n_edges", "edge_mean", "absolute_edge_mean",
         "negative_edge_fraction", "edge_sd"]
    ]
    strength_rows = [["patient_id"] + [f"abs_strength::{g}" for g in genes]]
    signed_rows = [["patient_id"] + [f"signed_strength::{g}" for g in genes]]
    for sample, rec in reports:
        table.append([sample, len(genes), rec["n_edges"], rec["edge_mean"],
                      rec["absolute_edge_mean"], rec["negative_edge_fraction"], rec["edge_sd"]])
        strength_rows.append([sample] + [rec["abs_strength"][g] for g in genes])
        signed_rows.append([sample] + [rec["signed_strength"][g] for g in genes])
    for filename, rows in [
        ("network_summary.csv", table),
        ("gene_absolute_strength.csv", strength_rows),
        ("gene_signed_strength.csv", signed_rows),
    ]:
        with (out / filename).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerows(rows)
    return {
        "n_patients": len(ids), "n_genes": len(genes),
        "n_edges_per_patient": reports[0][1]["n_edges"],
        "files": [str(out / name) for name in
                  ("network_summary.csv", "gene_absolute_strength.csv", "gene_signed_strength.csv")],
        "research_only": True,
        "sweet_inference_executed": False,
        "out_of_sample_sweet_validated": False,
        "notes": "Network structure statistics only; not CYP activity or measured PK.",
    }
