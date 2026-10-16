#!/usr/bin/env python3
"""Reproduce the *official* SWEET example with unmodified pinned upstream scripts.

Requires a separate checkout of SysMednet/SWEET@cf064eb5... and numpy. This
checks mathematical/software reproduction, NOT clinical validity, transductive
leakage, or inference of an unseen patient.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
import sys
import tempfile
from pathlib import Path

UPSTREAM_SHA = "cf064eb5f5b996c40bba3cce507384f6c0fdecb5"
ATOL = 1e-8
RTOL = 1e-8
SCRIPTS = [
    "1.SWEET_sample_weight_calculating.py",
    "2.SWEET_edge_score_calculating.py",
    "3.SWEET_calculating_mean_std_zscore.py",
]


class ReproductionError(RuntimeError):
    pass


def run_script(cmd: list[str], cwd: Path) -> None:
    result = subprocess.run(cmd, cwd=cwd, text=True, capture_output=True,
                            timeout=180, check=False)
    if result.returncode:
        raise ReproductionError(
            f"Original SWEET failed ({result.returncode}): {cmd[1]}\n"
            f"{result.stdout[-800:]}\n{result.stderr[-1500:]}")


def rows(path: Path, expected: list[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != expected:
            raise ReproductionError(f"Unexpected schema: {path}")
        data = [r for r in reader if r.get(expected[0])]
    if not data:
        raise ReproductionError(f"Empty official fixture: {path}")
    return data


def compare_csv(actual: Path, reference: Path, fields: list[str],
                numeric_field: str, allow_discrepancy: bool = False) -> dict:
    observed = rows(actual, fields)
    expected = rows(reference, fields)
    keys = [k for k in fields if k != numeric_field]
    def extract(data):
        out = {}
        for row in data:
            key = tuple(row[k] for k in keys)
            if key in out:
                raise ReproductionError(f"Duplicate key {key}")
            try:
                out[key] = float(row[numeric_field])
            except ValueError as exc:
                raise ReproductionError(f"Non-numeric value: {key}") from exc
        return out
    obs, ref = extract(observed), extract(expected)
    if set(obs) != set(ref):
        raise ReproductionError(
            f"Official keys do not match for {numeric_field}: "
            f"observed={len(obs)} reference={len(ref)}")
    worst = 0.0
    mismatches = 0
    for key, target in ref.items():
        value = obs[key]
        if not (math.isfinite(target) and math.isfinite(value)):
            raise ReproductionError(f"Non-finite value: {key}")
        diff = abs(value - target)
        worst = max(diff, worst)
        if diff > ATOL + RTOL * abs(target):
            mismatches += 1
            if not allow_discrepancy:
                raise ReproductionError(
                    f"Numerical mismatch in {numeric_field}, {key}: {value} != {target}")
    return {"compared": len(ref), "mismatched": mismatches,
            "maximum_absolute_difference": worst, "matches_archive": mismatches == 0}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=None)
    parser.add_argument("--strict-archived-zscores", action="store_true",
                        help="Treat upstream archived z-score discrepancy as a blocking error")
    args = parser.parse_args()
    upstream = args.upstream.resolve()
    if not (upstream / ".git").exists():
        raise ReproductionError("A separate git checkout of SWEET is required")
    head = subprocess.check_output(
        ["git", "-C", str(upstream), "rev-parse", "HEAD"], text=True).strip()
    if head != UPSTREAM_SHA:
        raise ReproductionError(f"Unpinned SWEET revision: {head}")
    for name in SCRIPTS:
        if not (upstream / name).is_file():
            raise ReproductionError(f"Missing upstream source: {name}")
    example = upstream / "example"
    patient_ids = [p.strip() for p in (example / "patient.txt").read_text().splitlines()
                   if p.strip()]
    if len(patient_ids) != len(set(patient_ids)) or not patient_ids:
        raise ReproductionError("Duplicated/empty published patient list")
    with (example / "expression.txt").open("r", encoding="utf-8") as stream:
        header = stream.readline().strip().split("\t")
    if not header or header[0] != "gene" or len(set(header[1:])) != len(header[1:]):
        raise ReproductionError("Invalid published expression header")
    cohort_sample_count = len(header) - 1
    if not set(patient_ids).issubset(set(header[1:])):
        raise ReproductionError("Target patient(s) missing from expression cohort")

    with tempfile.TemporaryDirectory(prefix="omics2pk-official-sweet-") as d:
        temp = Path(d)
        edge_dir = temp / "edges"
        edge_dir.mkdir()
        weight = temp / "weight.txt"
        stats = temp / "mean_std.txt"

        run_script([sys.executable, str(upstream / SCRIPTS[0]),
                    "-f", str(example / "expression.txt"),
                    "-s", str(weight)], temp)
        weights = compare_csv(
            weight, example / "weight.txt",
            ["patient", "sample_weight"], "sample_weight")
        if weights["compared"] != cohort_sample_count:
            raise ReproductionError("Official weights do not cover complete expression reference cohort")

        run_script([sys.executable, str(upstream / SCRIPTS[1]),
                    "-f", str(example / "expression.txt"),
                    "-w", str(weight),
                    "-p", str(example / "patient.txt"),
                    "-g", str(example / "gene.txt"),
                    "-s", str(edge_dir)], temp)
        if any(not (edge_dir / f"{patient}.txt").is_file()
               for patient in patient_ids):
            raise ReproductionError("Missing generated individual network(s)")
        first = patient_ids[0]
        edges = compare_csv(
            edge_dir / f"{first}.txt",
            example / f"{first}.txt",
            ["gene1", "gene2", "raw_edge_score"], "raw_edge_score")

        # Upstream argparse type=bool: deliberately pass True (not 'False').
        run_script([sys.executable, str(upstream / SCRIPTS[2]),
                    "-p", str(example / "patient.txt"),
                    "-l", str(edge_dir),
                    "-s", str(stats),
                    "-z", "True"], temp)
        def mean_std(file):
            result = {}
            for line in file.read_text().splitlines():
                if not line.strip():
                    continue
                k, v = line.split("\t")
                result[k] = float(v)
            return result
        expected_stats = mean_std(example / "mean_std.txt")
        observed_stats = mean_std(stats)
        for key in ("mean", "std"):
            if key not in expected_stats or key not in observed_stats:
                raise ReproductionError(f"Missing official {key}")
            if abs(expected_stats[key] - observed_stats[key]) > ATOL + RTOL * abs(expected_stats[key]):
                raise ReproductionError(
                    f"Official {key} mismatch: {observed_stats[key]} != {expected_stats[key]}")
        zs = compare_csv(
            edge_dir / f"{first}_zscore.txt",
            example / f"{first}_zscore.txt",
            ["gene1", "gene2", "z_score"], "z_score",
            allow_discrepancy=True)
        raw_rows = rows(edge_dir / f"{first}.txt", ["gene1", "gene2", "raw_edge_score"])
        computed_rows = rows(edge_dir / f"{first}_zscore.txt",
                             ["gene1", "gene2", "z_score"])
        raw_map = {(r["gene1"], r["gene2"]): float(r["raw_edge_score"]) for r in raw_rows}
        recomputed_z_checks = 0
        for record in computed_rows:
            key = (record["gene1"], record["gene2"])
            if key not in raw_map:
                raise ReproductionError("Generated z-score key missing raw edge")
            z_expected = (raw_map[key] - observed_stats["mean"]) / observed_stats["std"]
            z_calculated = float(record["z_score"])
            if abs(z_calculated - z_expected) > ATOL + RTOL * abs(z_expected):
                raise ReproductionError("Original script internally inconsistent with its mean/std")
            recomputed_z_checks += 1
        if recomputed_z_checks != len(raw_map):
            raise ReproductionError("Incomplete internally computed z-scores")
        if args.strict_archived_zscores and zs["mismatched"]:
            raise ReproductionError(
                f"Archived SWEET z-scores differ from algorithm output: "
                f"{zs['mismatched']} / {zs['compared']} mismatches")

    report = {
        "upstream_repository": "https://github.com/SysMednet/SWEET",
        "upstream_commit": UPSTREAM_SHA,
        "reference_cohort_sample_count": cohort_sample_count,
        "requested_network_count": len(patient_ids),
        "sample_weights": weights,
        "first_patient_edge_scores": edges,
        "first_patient_archived_z_scores": zs,
        "recomputed_zscore_formula_verified_edges": recomputed_z_checks,
        "archived_zscore_discrepancy": zs["mismatched"] > 0,
        "published_global_mean": expected_stats["mean"],
        "published_global_std": expected_stats["std"],
        "numerical_absolute_tolerance": ATOL,
        "numerical_relative_tolerance": RTOL,
        "result": ("CORE_REPRODUCED_BUT_ARCHIVED_ZSCORES_INCONSISTENT" if zs["mismatched"]
                   else "FULL_OFFICIAL_EXAMPLE_REPRODUCED"),
        "does_not_establish": [
            "out-of-sample network inference", "lack of transductive leakage",
            "complete archived z-score reproduction when discrepancy flag is true",
            "biological significance", "PK or PD prediction", "clinical validation",
        ],
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
