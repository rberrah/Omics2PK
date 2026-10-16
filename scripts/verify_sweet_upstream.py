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
                numeric_field: str) -> dict:
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
    for key, target in ref.items():
        value = obs[key]
        if not (math.isfinite(target) and math.isfinite(value)):
            raise ReproductionError(f"Non-finite value: {key}")
        diff = abs(value - target)
        worst = max(diff, worst)
        if diff > ATOL + RTOL * abs(target):
            raise ReproductionError(
                f"Numerical mismatch in {numeric_field}, {key}: {value} != {target}")
    return {"compared": len(ref), "maximum_absolute_difference": worst}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--report", type=Path, default=None)
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
        if weights["compared"] != len(patient_ids):
            raise ReproductionError("Official weights do not cover all samples")

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
            ["gene1", "gene2", "z_score"], "z_score")

    report = {
        "upstream_repository": "https://github.com/SysMednet/SWEET",
        "upstream_commit": UPSTREAM_SHA,
        "sample_count": len(patient_ids),
        "sample_weights": weights,
        "first_patient_edge_scores": edges,
        "first_patient_z_scores": zs,
        "published_global_mean": expected_stats["mean"],
        "published_global_std": expected_stats["std"],
        "numerical_absolute_tolerance": ATOL,
        "numerical_relative_tolerance": RTOL,
        "result": "PASS_official_example_only",
        "does_not_establish": [
            "out-of-sample network inference", "lack of transductive leakage",
            "biological significance", "PK or PD prediction", "clinical validation",
        ],
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
