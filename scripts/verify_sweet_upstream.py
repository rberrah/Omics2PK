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




def compare_archived_legacy_zscore(raw_path: Path, archived_z_path: Path,
                                   mean_value: float, std_value: float) -> dict:
    """Identify archived values created by the historical precedence bug.

    Older SWEET source computed raw - mean/std before upstream corrected the
    line to (raw - mean)/std. This identifies provenance; it never endorses
    the legacy values as valid z-scores.
    """
    raw_rows = rows(raw_path, ["gene1", "gene2", "raw_edge_score"])
    archived_rows = rows(archived_z_path, ["gene1", "gene2", "z_score"])
    raw_map = {(r["gene1"], r["gene2"]): float(r["raw_edge_score"]) for r in raw_rows}
    archived_map = {(r["gene1"], r["gene2"]): float(r["z_score"]) for r in archived_rows}
    if set(raw_map) != set(archived_map):
        raise ReproductionError("Archived z-score keys do not match raw edge keys")
    legacy_offset = mean_value / std_value
    mismatches = 0
    worst = 0.0
    for key, raw_value in raw_map.items():
        expected = raw_value - legacy_offset
        observed = archived_map[key]
        diff = abs(observed - expected)
        worst = max(worst, diff)
        if diff > ATOL + RTOL * abs(expected):
            mismatches += 1
    return {
        "compared": len(raw_map),
        "mismatched_legacy_formula": mismatches,
        "matches_legacy_raw_minus_mean_over_std": mismatches == 0,
        "legacy_offset_mean_over_std": legacy_offset,
        "maximum_absolute_difference_from_legacy_formula": worst,
    }


def _run_capture(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        cmd, cwd=cwd, text=True, capture_output=True,
        timeout=60, check=False,
    )


def _write_expression(path: Path, samples: list[str], rows_: list[tuple[str, list[object]]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer=csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene",*samples])
        for gene,values in rows_:
            writer.writerow([gene,*values])


def _write_simple_weight(path: Path, samples: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer=csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(["patient","sample_weight"])
        for sample in samples:
            writer.writerow([sample,"1.0"])


def audit_upstream_edge_cases(upstream: Path) -> dict:
    """Execute synthetic edge cases on the unmodified pinned upstream scripts.

    These are characterization tests, not proposed fixes to upstream SWEET.
    No biological or patient data are used.
    """
    py=sys.executable
    step1=str(upstream/SCRIPTS[0])
    step2=str(upstream/SCRIPTS[1])
    step3=str(upstream/SCRIPTS[2])
    samples=["S1","S2","S3","S4"]
    with tempfile.TemporaryDirectory(prefix="omics2pk-sweet-edge-cases-") as d:
        root=Path(d)

        # Step 1 does not run the missing-token coercion used by Step 2.
        na=root/"na"; na.mkdir()
        _write_expression(na/"expression.txt",samples,[
            ("g1",[1,"NA",3,4]),("g2",[2,3,5,8]),("g3",[4,6,7,9]),
        ])
        na_run=_run_capture(
            [py,step1,"-f",str(na/"expression.txt"),"-s",str(na/"weight.txt")],na)

        # The literal token "nan" is accepted by float conversion, after which
        # NumPy correlations/weights become non-finite without a non-zero exit.
        nan_case=root/"nan"; nan_case.mkdir()
        _write_expression(nan_case/"expression.txt",samples,[
            ("g1",[1,"nan",3,4]),("g2",[2,3,5,8]),("g3",[4,6,7,9]),
        ])
        nan_run=_run_capture(
            [py,step1,"-f",str(nan_case/"expression.txt"),
             "-s",str(nan_case/"weight.txt")],nan_case)
        nan_weights=[]
        if (nan_case/"weight.txt").is_file():
            for row in rows(nan_case/"weight.txt",["patient","sample_weight"]):
                nan_weights.append(float(row["sample_weight"]))
        step1_nan_nonfinite=(
            nan_run.returncode==0 and bool(nan_weights)
            and any(not math.isfinite(x) for x in nan_weights)
        )

        # Shared Step 2 fixtures.
        def run_step2(case: Path, expression_rows, genes=("g1","g2","g3")):
            case.mkdir()
            _write_expression(case/"expression.txt",samples,expression_rows)
            _write_simple_weight(case/"weight.txt",samples)
            (case/"patient.txt").write_text("S1\n",encoding="utf-8")
            (case/"gene.txt").write_text("\n".join(genes)+"\n",encoding="utf-8")
            edges=case/"edges"; edges.mkdir()
            proc=_run_capture([
                py,step2,"-f",str(case/"expression.txt"),
                "-w",str(case/"weight.txt"),"-p",str(case/"patient.txt"),
                "-g",str(case/"gene.txt"),"-s",str(edges),
            ],case)
            output=(edges/"S1.txt").read_text(encoding="utf-8") if (edges/"S1.txt").is_file() else ""
            return proc,output

        missing2_run,missing2_edges=run_step2(root/"missing2",[
            ("g1",["NA",2,3,4]),("g2",[1,3,5,8]),("g3",[2,5,4,9]),
        ])
        step2_missing_to_zero=(
            missing2_run.returncode==0
            and "assigned to 0" in missing2_run.stdout
            and bool(missing2_edges)
        )

        zerosum_run,zerosum_edges=run_step2(root/"zerosum",[
            ("g1",[-1,0,1,0]),("g2",[1,3,5,8]),("g3",[2,5,4,9]),
        ])
        zero_sum_gene_deleted=(
            zerosum_run.returncode==0
            and "g1" not in zerosum_edges
            and "g2" in zerosum_edges and "g3" in zerosum_edges
        )

        constant_run,constant_edges=run_step2(root/"constant",[
            ("g1",[5,5,5,5]),("g2",[1,3,5,8]),("g3",[2,5,4,9]),
        ])
        constant_gene_nonfinite=(
            constant_run.returncode==0 and "nan" in constant_edges.casefold()
        )

        # argparse type=bool treats the non-empty string "False" as True.
        boolcase=root/"bool"; boolcase.mkdir()
        (boolcase/"patient.txt").write_text("S1\n",encoding="utf-8")
        edges=boolcase/"edges"; edges.mkdir()
        (edges/"S1.txt").write_text(
            "gene1\tgene2\traw_edge_score\n"
            "g1\tg2\t0.1\n"
            "g1\tg3\t0.2\n"
            "g2\tg3\t0.3\n",encoding="utf-8")
        bool_run=_run_capture([
            py,step3,"-p",str(boolcase/"patient.txt"),"-l",str(edges),
            "-s",str(boolcase/"mean_std.txt"),"-z","False",
        ],boolcase)
        false_string_enables_zscore=(
            bool_run.returncode==0 and (edges/"S1_zscore.txt").is_file()
        )

        # Global edge SD=0 makes z-scoring undefined/non-finite.
        zerostd=root/"zerostd"; zerostd.mkdir()
        (zerostd/"patient.txt").write_text("S1\n",encoding="utf-8")
        zedges=zerostd/"edges"; zedges.mkdir()
        (zedges/"S1.txt").write_text(
            "gene1\tgene2\traw_edge_score\n"
            "g1\tg2\t0.5\n"
            "g1\tg3\t0.5\n"
            "g2\tg3\t0.5\n",encoding="utf-8")
        zerostd_run=_run_capture([
            py,step3,"-p",str(zerostd/"patient.txt"),"-l",str(zedges),
            "-s",str(zerostd/"mean_std.txt"),"-z","True",
        ],zerostd)
        stats={}
        if (zerostd/"mean_std.txt").is_file():
            for line in (zerostd/"mean_std.txt").read_text().splitlines():
                if line.strip():
                    key,value=line.split("\t")
                    stats[key]=float(value)
        z_nonfinite=False
        if (zedges/"S1_zscore.txt").is_file():
            for row in rows(zedges/"S1_zscore.txt",["gene1","gene2","z_score"]):
                if not math.isfinite(float(row["z_score"])):
                    z_nonfinite=True
                    break
        zero_std_undefined=(
            zerostd_run.returncode==0 and stats.get("std")==0.0 and z_nonfinite
        )

    observed={
        "step1_NA_token_returns_nonzero":na_run.returncode!=0,
        "step1_literal_nan_can_emit_nonfinite_weights":step1_nan_nonfinite,
        "step2_missing_token_is_coerced_to_zero":step2_missing_to_zero,
        "step2_nonconstant_zero_sum_gene_is_deleted":zero_sum_gene_deleted,
        "step2_constant_nonzero_gene_can_emit_nan_edges":constant_gene_nonfinite,
        "step3_dash_z_False_string_still_generates_zscores":false_string_enables_zscore,
        "step3_zero_global_edge_sd_produces_nonfinite_zscores":zero_std_undefined,
    }
    if not all(observed.values()):
        raise ReproductionError(
            "Pinned SWEET edge-case behavior changed; review qualification before proceeding: "
            +json.dumps(observed,sort_keys=True))
    return {
        **observed,
        "synthetic_only":True,
        "qualification_interpretation":{
            "missingness":"Step 1 and Step 2 are inconsistent; Omics2PK forbids missing expression before SWEET.",
            "zero_sum":"Step 2 tests row sum==0, not all-zero/variance; Omics2PK rejects nonconstant zero-sum genes.",
            "constant_profiles":"Undefined Pearson correlations are possible; Omics2PK rejects constant genes/samples.",
            "z_flag":"Do not pass '-z False'; argparse bool treats non-empty strings as true.",
            "zero_edge_sd":"Z-score standardization is undefined; raw edge summaries must fail closed on degenerate networks.",
        },
    }


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
        legacy_archive_checks = {}
        for patient in patient_ids:
            legacy_archive_checks[patient] = compare_archived_legacy_zscore(
                edge_dir / f"{patient}.txt",
                example / f"{patient}_zscore.txt",
                expected_stats["mean"], expected_stats["std"])
        archived_files_match_legacy_bug = all(
            item["matches_legacy_raw_minus_mean_over_std"]
            for item in legacy_archive_checks.values())
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

    edge_case_audit=audit_upstream_edge_cases(upstream)

    report = {
        "upstream_repository": "https://github.com/SysMednet/SWEET",
        "upstream_commit": UPSTREAM_SHA,
        "reference_cohort_sample_count": cohort_sample_count,
        "synthetic_edge_case_audit":edge_case_audit,
        "requested_network_count": len(patient_ids),
        "sample_weights": weights,
        "first_patient_edge_scores": edges,
        "first_patient_archived_z_scores": zs,
        "recomputed_zscore_formula_verified_edges": recomputed_z_checks,
        "archived_zscore_discrepancy": zs["mismatched"] > 0,
        "archived_zscore_legacy_precedence_provenance": legacy_archive_checks,
        "archived_zscore_files_match_legacy_bug": archived_files_match_legacy_bug,
        "legacy_bug_formula": "raw_edge_score - mean/std",
        "correct_current_formula": "(raw_edge_score - mean)/std",
        "upstream_history_note": (
            "Archived example z-scores match the historical operator-precedence bug; "
            "current upstream source has corrected parentheses, but archived z-score files were not regenerated."
            if archived_files_match_legacy_bug else
            "Archived z-score provenance not fully explained by the known historical precedence bug."
        ),
        "published_global_mean": expected_stats["mean"],
        "published_global_std": expected_stats["std"],
        "numerical_absolute_tolerance": ATOL,
        "numerical_relative_tolerance": RTOL,
        "result": (
            "CORE_REPRODUCED_ARCHIVED_ZSCORES_IDENTIFIED_AS_LEGACY_PRECEDENCE_BUG"
            if zs["mismatched"] and archived_files_match_legacy_bug
            else ("CORE_REPRODUCED_BUT_ARCHIVED_ZSCORES_UNEXPLAINED" if zs["mismatched"]
                  else "FULL_OFFICIAL_EXAMPLE_REPRODUCED")
        ),
        "does_not_establish": [
            "out-of-sample network inference", "lack of transductive leakage",
            "validity of the legacy archived z-score values as z-scores",
            "biological significance", "PK or PD prediction", "clinical validation",
        ],
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
