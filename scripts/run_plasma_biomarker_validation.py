#!/usr/bin/env python
"""Execute the frozen Lee / Jackson validation: replication, and one new contrast.

Protocol:
``docs/OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_PRESPECIFICATION_20261110.md``

Input: the per-donor table written OUTSIDE this repository by
``scripts/audit_lee_jackson_cyp3a_panel.py``, whose own run verified the
extraction against three correlations the source article publishes.

No donor identifier, value or prediction is written to the artefact.

Usage:
    python -I scripts/run_plasma_biomarker_validation.py \
        --per-donor <outside-repo>/lee_jackson_cyp3a_panel_per_donor.csv \
        --out reports/data/plasma_biomarker_validation.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8")

from omics2pk import plasma_biomarker_validation as PBV  # noqa: E402
from omics2pk import protein_clint_pilot as PCP  # noqa: E402

PRESPEC = (ROOT / "docs"
           / "OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_PRESPECIFICATION_20261110.md")
PROTOCOL_VERSION = "1.0.0"

#: What the Wegler pilot found for the same relationship, to be read beside
#: this cohort's answer. Both are hepatocyte systems and midazolam is the probe
#: in both; the endpoints differ, so this is concordance of DIRECTION and
#: never a pooled estimate.
WEGLER_COMPARATOR = {
    "source": "docs/OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_RESULTS_20261110.md",
    "cohort": "Wegler, 15 liver donors, Uppsala",
    "system": "hepatocytes",
    "endpoint": "midazolam CLint by equation 7",
    "finding": ("the target CYP's own abundance predicted its substrate's "
                "intrinsic clearance out of sample for CYP3A4 and CYP2C19 "
                "after Bonferroni; midazolam's geometric fold error fell from "
                "1.55 to 1.18"),
    "why_this_is_not_a_pooled_estimate": (
        "different donors, different laboratory and a different endpoint -- a "
        "fitted intrinsic clearance there, a single-concentration formation "
        "rate here. Two cohorts agreeing in direction is a replication of the "
        "direction, not one estimate measured twice."),
}


class RunError(RuntimeError):
    pass


def check_imported_machinery():
    """Refuse to run if what is imported is not what was frozen."""
    if PBV.FROZEN["ridge_alphas"] != PCP.FROZEN["ridge_alphas"]:
        raise RunError("the penalty grid no longer matches the frozen pilot's")
    if PBV.FROZEN["seed"] != PCP.FROZEN["seed"]:
        raise RunError("the seed no longer matches the frozen pilot's")
    for name, spec in PBV.CONTRASTS.items():
        for key in ("arm", "comparator"):
            if spec[key] not in PBV.ARMS:
                raise RunError(f"{name} names a missing arm {spec[key]!r}")
        added = set(PBV.ARMS[spec["arm"]]) - set(PBV.ARMS[spec["comparator"]])
        if set(spec["permuted_block"]) != added:
            raise RunError(
                f"{name} permutes {spec['permuted_block']} but the arm adds "
                f"{sorted(added)}; the null would not test the contrast")
    # The guard must refuse every predictor as an endpoint. Exercising it
    # against itself is the only way to know it is still a guard.
    for forbidden in PBV.FORBIDDEN_ENDPOINTS:
        try:
            PBV.assert_endpoint_is_admissible(forbidden)
        except PBV.PlasmaBiomarkerValidationError:
            continue
        raise RunError(
            f"the circularity guard accepted {forbidden!r} as an endpoint, "
            "which is the one thing it exists to refuse")
    expected = (len(PBV.CONTRASTS)
                * len({PBV.FROZEN["primary_endpoint"],
                       PBV.FROZEN["secondary_endpoint"]}))
    if PBV.FROZEN["declared_family_size"] != expected:
        raise RunError(
            f"the declared family size is {PBV.FROZEN['declared_family_size']} "
            f"but the design has {expected} primary contrasts")


def read_per_donor(path):
    rows, order = {}, []
    with open(path, newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle):
            donor = record.pop("donor")
            parsed = {}
            for key, value in record.items():
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    parsed[key] = value
            parsed["cyp3a_total_abundance"] = (parsed["cyp3a4_abundance"]
                                               + parsed["cyp3a5_abundance"])
            rows[donor] = parsed
            order.append(donor)
    if len(order) != 15:
        raise RunError(f"expected fifteen donors, read {len(order)}")
    return rows, order


def run_endpoint(rows, donors, endpoint, family_size):
    out = {}
    for name in PBV.CONTRASTS:
        block, _, _ = PBV.contrast(rows, donors, endpoint, name,
                                   family_size=family_size)
        null = PBV.permutation_null(rows, donors, endpoint, name)
        observed = block["stored_error_bootstrap"]["mean_difference_log10"]
        block["permutation_control"] = PCP.permutation_rank(observed, null)
        block["permutation_control"]["null_deciles"] = PBV.deciles(null)
        block["permutation_control"].pop("null_min", None)
        block["permutation_control"].pop("null_max", None)
        survives = (
            block["refitting_out_of_bag_bootstrap"]["interval_excludes_zero"]
            and block["refitting_out_of_bag_bootstrap"][
                "mean_difference_log10"] < 0.0
            and block["permutation_control"]["rank_p_value"] <= 0.05
        )
        block["verdict"] = "SUPPORTED" if survives else "NOT_SUPPORTED"
        block["verdict_rule"] = (
            "SUPPORTED requires all three: the REFITTING interval excludes "
            "zero, its point estimate is an improvement, and the permutation "
            "rank p-value is at most 0.05. The stored-error interval alone is "
            "never sufficient.")
        out[name] = block
    return out


def run_sensitivity(rows, donors, endpoint):
    """CYP3A4 + CYP3A5 in place of CYP3A4, declared and not Bonferroni-pooled."""
    results = {}
    for arm, columns in PBV.SENSITIVITY_ARMS.items():
        score, errors = PBV.score_arm(rows, donors, endpoint, columns)
        results[arm] = {"score": score}
        results[arm]["score"].pop("penalty_by_fold", None)
    base, base_errors = PBV.score_arm(rows, donors, endpoint,
                                      PBV.ARMS["V0_mean"])
    total_score, total_errors = PBV.score_arm(
        rows, donors, endpoint, PBV.SENSITIVITY_ARMS["V1s_cyp3a_total"])
    plus_score, plus_errors = PBV.score_arm(
        rows, donors, endpoint,
        PBV.SENSITIVITY_ARMS["V3s_cyp3a_total_plus_plasma"])
    results["replication_with_cyp3a_total_vs_mean"] = PCP.paired_bootstrap(
        total_errors, base_errors, family_size=1)
    results["incremental_value_with_cyp3a_total"] = PCP.paired_bootstrap(
        plus_errors, total_errors, family_size=1)
    results["how_to_read"] = (
        "A declared sensitivity on the choice of the protein predictor. It is "
        "not a second primary analysis and carries no multiplicity "
        "adjustment of its own, so it may corroborate the primary reading and "
        "may not replace it.")
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-donor", required=True)
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "plasma_biomarker_validation.json"))
    args = parser.parse_args(argv)

    if not PRESPEC.is_file():
        raise RunError(f"the frozen protocol is missing: {PRESPEC}")
    check_imported_machinery()

    per_donor = pathlib.Path(args.per_donor).resolve()
    if ROOT == per_donor.parent or ROOT in per_donor.parents:
        raise RunError(
            f"the per-donor table is inside this repository ({per_donor}); "
            "this source is CC BY-NC-ND and individual rows are not committed")
    rows, donors = read_per_donor(per_donor)

    family = PBV.FROZEN["declared_family_size"]
    endpoints = {
        "primary": PBV.FROZEN["primary_endpoint"],
        "secondary": PBV.FROZEN["secondary_endpoint"],
    }
    results = {role: run_endpoint(rows, donors, endpoint, family)
               for role, endpoint in endpoints.items()}
    sensitivity = run_sensitivity(rows, donors, endpoints["primary"])

    supported = [
        f"{role}::{name}"
        for role, block in results.items()
        for name, cell in block.items()
        if cell["verdict"] == "SUPPORTED"
    ]

    payload = {
        "run": "axis_b2_plasma_biomarker_validation",
        "generated": "2026-11-10",
        "protocol": PRESPEC.name,
        "protocol_version": PROTOCOL_VERSION,
        "protocol_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
        "per_donor_table_sha256": hashlib.sha256(
            per_donor.read_bytes()).hexdigest(),
        "per_donor_table_is_outside_this_repository": str(per_donor),
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in PBV.FROZEN.items()},
        "arms": {k: list(v) for k, v in PBV.ARMS.items()},
        "endpoints": endpoints,
        "n_donors": len(donors),
        "declared_family_size": family,
        "results": results,
        "sensitivity_cyp3a_total": sensitivity,
        "contrasts_supported": supported,
        "n_contrasts_supported": len(supported),
        "wegler_comparator": WEGLER_COMPARATOR,
        "what_this_cannot_claim": PBV.what_this_cannot_claim(),
        "status_ceiling": PBV.FROZEN["status_ceiling"],
        "no_donor_identifiers_values_or_predictions_exported": True,
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    for role, block in results.items():
        print(f"\n=== {role}: {endpoints[role]} ===")
        for name, cell in block.items():
            stored = cell["stored_error_bootstrap"]
            refit = cell["refitting_out_of_bag_bootstrap"]
            null = cell["permutation_control"]
            print(f"  {name}")
            print(f"    arm {cell['arm']:<24} fold error "
                  f"{cell['arm_score']['geometric_fold_error']:.3f}  "
                  f"| comparator {cell['comparator']:<22} fold error "
                  f"{cell['comparator_score']['geometric_fold_error']:.3f}")
            print(f"    stored-error  d={stored['mean_difference_log10']:+.4f} "
                  f"[{stored['interval'][0]:+.4f}, {stored['interval'][1]:+.4f}]"
                  f"  excludes zero {stored['interval_excludes_zero']}")
            print(f"    REFITTING     d={refit['mean_difference_log10']:+.4f} "
                  f"[{refit['interval'][0]:+.4f}, {refit['interval'][1]:+.4f}]"
                  f"  excludes zero {refit['interval_excludes_zero']}"
                  f"  ({refit['n_draws_used']} draws)")
            print(f"    null median {null['null_median']:+.4f}  "
                  f"rank p {null['rank_p_value']:.4f}  -> {cell['verdict']}")
    print(f"\nsupported contrasts: {len(supported)} of "
          f"{family} -> {supported}")
    print(f"artefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
