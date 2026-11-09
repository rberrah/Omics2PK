#!/usr/bin/env python
"""Merge endpoint batch files into one result, refusing any inconsistency.

The endpoints were run in parallel batches to fit the wall clock. That is a
scheduling decision and must have no effect on the result, so this script
verifies it had none: **every top-level field except the endpoint results and
the run date must be identical in every batch**, and no endpoint may appear
twice. A mismatch is an error, not something to reconcile.

The all-fields rule replaced a named list of fields in v1.4.5, for two reasons:
it cannot fall out of date when a report gains a field -- the qualified pool is
recorded inside `universe`, and a merge across two different pools would have
passed the old list unnoticed -- and it makes the gate independent of which
analysis produced the batches, so an H2d batch set is checked as strictly as an
H2a one.

    python scripts/merge_h2a_parts.py --parts reports/data/h2a_part_*.json \
        --out reports/data/h2a_hlc_selection.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

EXPECTED_ENDPOINTS = (
    "CYP1A2_activity", "CYP2A6_activity", "CYP2B6_activity", "CYP2C8_activity",
    "CYP2C9_activity", "CYP2C19_activity", "CYP2D6_activity", "CYP2E1_activity",
    "CYP3A4M_activity", "CYP3A4T_activity",
)

# Everything outside this set must be identical in every batch. A batch's own
# endpoint results obviously differ, and a run that crosses midnight carries a
# different date without meaning anything.
PER_BATCH = ("endpoints", "run_on")

# Fields no batch may omit, whatever the analysis produced it.
REQUIRED = ("endpoint_class", "frozen", "universe")


class MergeError(RuntimeError):
    pass


def merge(paths, *, require_all=True, expected=EXPECTED_ENDPOINTS):
    parts = []
    for path in paths:
        text = Path(path).read_text(encoding="utf-8")
        parts.append((str(path), json.loads(text)))
    if not parts:
        raise MergeError("no batch files given")

    reference_name, reference = parts[0]
    for name, part in parts:
        absent = [field for field in REQUIRED if field not in part]
        if absent:
            raise MergeError(f"{name} declares no {absent}")
    for name, part in parts[1:]:
        for field in sorted((set(reference) | set(part)) - set(PER_BATCH)):
            if part.get(field) != reference.get(field):
                raise MergeError(
                    f"{name} disagrees with {reference_name} on '{field}'; "
                    "the batches were not run under the same specification")

    merged = {key: value for key, value in reference.items() if key != "endpoints"}
    merged["endpoints"] = {}
    merged["merged_from"] = []
    for name, part in parts:
        merged["merged_from"].append(
            {"file": Path(name).name,
             "endpoints": sorted(part["endpoints"]),
             "run_on": part.get("run_on")})
        for endpoint, result in part["endpoints"].items():
            if endpoint in merged["endpoints"]:
                raise MergeError(f"{endpoint} appears in more than one batch")
            merged["endpoints"][endpoint] = result

    present = set(merged["endpoints"])
    unexpected = sorted(present - set(expected))
    if unexpected:
        raise MergeError(f"unexpected endpoints: {unexpected}")
    missing = sorted(set(expected) - present)
    if missing and require_all:
        raise MergeError(f"incomplete: {len(missing)} endpoints absent: {missing}")
    merged["endpoints_missing"] = missing

    designs = {result["fold_design"] for result in merged["endpoints"].values()}
    if len(designs) > 1:
        raise MergeError(f"batches used different fold designs: {sorted(designs)}")
    repeats = {result["repeats"] for result in merged["endpoints"].values()}
    if len(repeats) > 1:
        raise MergeError(f"batches used different repeat counts: {sorted(repeats)}")

    merged["merge_note"] = (
        "Run in parallel batches for wall-clock reasons only. Every top-level "
        "field except the endpoint results and the run date was checked to be "
        "identical across batches -- the frozen parameters, the network, the "
        "candidate universe including its qualified pool, the fold design and "
        "the repeat count among them -- and no endpoint appears twice.")
    merged["endpoints_merged"] = len(present)
    if any("first_fold_selection" in result
           for result in merged["endpoints"].values()):
        merged["wrong_mechanism_control_audit"] = audit_wrong_mechanism_control(
            merged["endpoints"])
    else:
        merged["wrong_mechanism_control_audit"] = {
            "not_applicable": (
                "These batches carry no per-fold selection table, so the "
                "wrong-mechanism rotation cannot be audited from them. The "
                "audit is absent rather than reported as passed.")}
    return merged


def audit_wrong_mechanism_control(endpoints):
    """Check that the deliberately wrong seed really is a different gene.

    The rotation was fixed before running and maps each endpoint to the next
    one's enzyme. It does not guarantee a *different* enzyme: where a cohort
    measures one enzyme with two substrates, the next label can resolve to the
    same gene, and the control then compares an arm against itself. That is a
    void control, not a result, and it has to be declared rather than noticed
    by eye in a table of near-identical numbers.
    """
    audit = {"valid": [], "void_same_gene": [], "undetermined": []}
    for endpoint, result in sorted(endpoints.items()):
        selections = result.get("first_fold_selection") or {}
        target = (selections.get("netpk_target_seeded") or {}).get("1") or []
        wrong = (selections.get("netpk_wrong_mechanism") or {}).get("1") or []
        if not target or not wrong:
            audit["undetermined"].append(endpoint)
        elif target[0] == wrong[0]:
            audit["void_same_gene"].append(
                {"endpoint": endpoint, "gene": target[0],
                 "wrong_seed_label": result.get("wrong_seed_label")})
        else:
            audit["valid"].append(
                {"endpoint": endpoint, "target_gene": target[0],
                 "wrong_gene": wrong[0]})
    audit["how_to_read"] = (
        "A void control means the wrong-mechanism arm was seeded with the same "
        "gene as the target arm, so for that endpoint the control says nothing "
        "and its numbers must not be read as evidence that a wrong mechanism "
        "costs nothing. The frozen rotation is not re-run to repair this: the "
        "defect is reported, and the endpoints with a valid control carry the "
        "finding.")
    audit["valid_count"] = len(audit["valid"])
    audit["void_count"] = len(audit["void_same_gene"])
    return audit


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parts", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--allow-incomplete", action="store_true",
                        help="Merge what exists; the result records what is absent.")
    parser.add_argument("--expected", nargs="+", default=list(EXPECTED_ENDPOINTS),
                        help="The endpoints this merge must end up with; any "
                             "other one is refused. Defaults to the ten "
                             "cohort endpoints.")
    args = parser.parse_args(argv)
    try:
        merged = merge(args.parts, require_all=not args.allow_incomplete,
                       expected=tuple(args.expected))
    except MergeError as exc:
        print(f"MERGE REFUSED: {exc}", file=sys.stderr)
        return 2
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(merged, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    print(f"wrote {out} with {merged['endpoints_merged']} endpoints")
    if merged["endpoints_missing"]:
        print(f"  absent: {merged['endpoints_missing']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
