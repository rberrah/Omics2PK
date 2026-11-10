"""Does an independent cohort reproduce the protein-explains-activity figure?

The project's strongest funding number so far rested on ONE published source:
Grangeon 2025, 16 microsomal donors, median protein-to-activity Spearman
**0.944**, against a median transcript-to-activity Spearman of **0.508** in the
Human Liver Cohort. Three explanations were recorded as inseparable from the
public record: biology, assay precision, and a shared microsomal preparation.

Wegler 2021 Table S3 now supplies the missing test. It publishes, per donor,
the matching CYP concentration **and** the intrinsic clearance of its probe
substrate in the same 15 donors, in microsomes and in hepatocytes. Unlike
Grangeon's, its correlations can be recomputed, and they were: all twenty
anchors of its own Table S4 reproduce
(`scripts/audit_wegler_table_s3_per_donor_kinetics.py`).

This script computes, for each of the five matched enzyme-substrate pairs, the
protein-to-clearance Spearman correlation with a donor bootstrap, and asks two
questions with a yes-or-no answer:

* does **Grangeon's** published value for the same enzyme fall inside the
  interval this cohort supports?
* does the **HLC transcript** correlation for the same enzyme fall inside it?

No model, no selection, no tuning. Aggregates only; the per-donor table stays
in the local data directory.

Run with ``python -I``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RNA_ARTEFACT = REPO / "reports" / "data" / "rna_versus_protein_activity_correlation.json"

#: SHA-256 of the canonical per-donor table, as emitted by the audit script.
EXPECTED_TABLE_SHA256 = (
    "750a43095e223bff"  # prefix only; the audit report carries the full value
)

#: Wegler substrate -> (its CYP, the HLC endpoint key for the same enzyme).
PAIRS = {
    "Midazolam": ("CYP3A4", "CYP3A4M_activity"),
    "Omeprazole": ("CYP2C19", "CYP2C19_activity"),
    "Diclofenac": ("CYP2C9", "CYP2C9_activity"),
    "Bupropion": ("CYP2B6", "CYP2B6_activity"),
    "Bufuralol": ("CYP2D6", "CYP2D6_activity"),
}

SYSTEMS = {
    "HLM": ("cyp_pmol_per_mg_HLM", "clint_mic_eq7"),
    "HH": ("cyp_pmol_per_mg_HH", "clint_hep_eq7"),
}

N_BOOTSTRAP = 2000
SEED = 20261110


class ComparisonError(RuntimeError):
    pass


def _ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    index = 0
    while index < len(order):
        stop = index
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[index]]:
            stop += 1
        average = (index + stop) / 2.0 + 1.0
        for position in range(index, stop + 1):
            ranks[order[position]] = average
        index = stop + 1
    return ranks


def _pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return None
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    return sxy / math.sqrt(sxx * syy)


def spearman(xs, ys):
    return _pearson(_ranks(list(xs)), _ranks(list(ys)))


def bootstrap_interval(xs, ys, n_boot=N_BOOTSTRAP, seed=SEED):
    import numpy as np

    rng = np.random.default_rng(seed)
    n = len(xs)
    draws = []
    for _ in range(n_boot):
        picked = rng.integers(0, n, size=n)
        resample_x = [xs[i] for i in picked]
        resample_y = [ys[i] for i in picked]
        value = spearman(resample_x, resample_y)
        if value is not None:
            draws.append(value)
    if len(draws) < n_boot // 2:
        return None
    draws.sort()
    low = draws[int(0.025 * len(draws))]
    high = draws[int(0.975 * len(draws)) - 1]
    return [float(low), float(high)]


def load_per_donor(path):
    rows = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["substrate"], row["donor"])
            rows[key] = row
    if not rows:
        raise ComparisonError(f"no rows read from {path}")
    return rows


def run(per_donor_csv, audit_json):
    audit = json.loads(Path(audit_json).read_text(encoding="utf-8"))
    if not audit["verification"]["all_anchors_agree"]:
        raise ComparisonError(
            "the per-donor table has not been verified against the article; "
            "refusing to build a comparison on an unverified transcription"
        )
    if not audit["per_donor_table_sha256"].startswith(EXPECTED_TABLE_SHA256):
        raise ComparisonError(
            "the audit report's per-donor table checksum does not match the one "
            "this comparison was written against"
        )
    rows = load_per_donor(per_donor_csv)
    rna = json.loads(RNA_ARTEFACT.read_text(encoding="utf-8"))

    results = {}
    for substrate, (cyp, hlc_key) in PAIRS.items():
        entry = {"cyp": cyp, "hlc_endpoint": hlc_key, "systems": {}}
        for system, (concentration, clearance) in SYSTEMS.items():
            xs, ys = [], []
            not_determined = 0
            for donor_index in range(1, 16):
                row = rows.get((substrate, f"UU{donor_index}"))
                if row is None:
                    raise ComparisonError(f"{substrate} UU{donor_index} absent")
                x, y = row[concentration], row[clearance]
                if "ND" in (x, y):
                    not_determined += 1
                    continue
                xs.append(float(x))
                ys.append(float(y))
            if len(xs) < 3:
                entry["systems"][system] = {
                    "status": "NOT_EVALUABLE",
                    "n_usable": len(xs),
                    "n_not_determined": not_determined,
                }
                continue
            rho = spearman(xs, ys)
            interval = bootstrap_interval(xs, ys)
            entry["systems"][system] = {
                "n_donors": len(xs),
                "n_not_determined": not_determined,
                "spearman_rho": round(float(rho), 4),
                "bootstrap_95_interval": (
                    [round(v, 4) for v in interval] if interval else None
                ),
            }

        hlc = rna["per_endpoint"].get(hlc_key)
        if hlc is None:
            raise ComparisonError(
                f"the committed mRNA-versus-protein artefact has no endpoint "
                f"{hlc_key!r}, so {substrate} cannot be compared"
            )
        grangeon = hlc.get("published_protein", {}).get("hlm_r")
        transcript_rho = hlc.get("hlc_transcript_activity_spearman")
        transcript_interval = (
            [hlc["hlc_ci_low"], hlc["hlc_ci_high"]]
            if "hlc_ci_low" in hlc and "hlc_ci_high" in hlc
            else None
        )

        microsomal = entry["systems"].get("HLM", {})
        interval = microsomal.get("bootstrap_95_interval")
        entry["comparison_in_microsomes"] = {
            "grangeon_published_hlm_r": grangeon,
            "grangeon_inside_this_cohorts_interval": (
                None
                if grangeon is None or interval is None
                else bool(interval[0] <= grangeon <= interval[1])
            ),
            "hlc_transcript_rho": transcript_rho,
            "hlc_transcript_interval_in_its_own_cohort": transcript_interval,
            "hlc_transcript_n_donors": hlc.get("hlc_n_donors"),
            "hlc_transcript_inside_this_cohorts_interval": (
                None
                if transcript_rho is None or interval is None
                else bool(interval[0] <= transcript_rho <= interval[1])
            ),
            "intervals_overlap_with_the_transcript": (
                None
                if transcript_interval is None or interval is None
                else bool(
                    interval[0] <= transcript_interval[1]
                    and transcript_interval[0] <= interval[1]
                )
            ),
            "why_microsomes": (
                "Grangeon's figure is microsomal, so the microsomal system is "
                "the like-for-like comparison. The hepatocyte system is "
                "reported beside it and is a different preparation."
            ),
        }
        results[substrate] = entry

    def median(values):
        ordered = sorted(values)
        if not ordered:
            return None
        middle = len(ordered) // 2
        return (
            ordered[middle]
            if len(ordered) % 2
            else (ordered[middle - 1] + ordered[middle]) / 2.0
        )

    microsomal_rhos = [
        results[s]["systems"]["HLM"]["spearman_rho"]
        for s in PAIRS
        if "spearman_rho" in results[s]["systems"].get("HLM", {})
    ]
    hepatocyte_rhos = [
        results[s]["systems"]["HH"]["spearman_rho"]
        for s in PAIRS
        if "spearman_rho" in results[s]["systems"].get("HH", {})
    ]
    inside_grangeon = [
        s
        for s in PAIRS
        if results[s]["comparison_in_microsomes"]["grangeon_inside_this_cohorts_interval"]
        is True
    ]
    outside_grangeon = [
        s
        for s in PAIRS
        if results[s]["comparison_in_microsomes"]["grangeon_inside_this_cohorts_interval"]
        is False
    ]
    inside_transcript = [
        s
        for s in PAIRS
        if results[s]["comparison_in_microsomes"][
            "hlc_transcript_inside_this_cohorts_interval"
        ]
        is True
    ]

    report = {
        "label": "INDEPENDENT_REPLICATION_ATTEMPT",
        "question": (
            "Does an independent cohort with recomputable per-donor data "
            "reproduce the published protein-explains-activity figure?"
        ),
        "this_cohort": (
            "Wegler et al., Mol Pharm 2021, 18(4):1792-1805, Table S3: 15 "
            "donors, five probe substrates, CYP concentration and intrinsic "
            "clearance in the same donors, in microsomes and hepatocytes."
        ),
        "verification_of_this_cohorts_transcription": {
            "anchors_agreeing": audit["verification"]["n_anchors_agreeing"],
            "anchors_total": audit["verification"]["n_anchors"],
            "median_check": audit["verification"]["median_midazolam_clint_mic_eq7"],
        },
        "endpoint_difference": (
            "Grangeon and the HLC compare protein or transcript against a probe "
            "ACTIVITY; this cohort compares protein against an intrinsic "
            "CLEARANCE measured by substrate depletion. Those are related but "
            "not the same quantity, and the comparison is descriptive."
        ),
        "statistic": "Spearman rho, average ranks for ties, no model",
        "interval": f"{N_BOOTSTRAP}-draw donor bootstrap, seed {SEED}",
        "per_pair": results,
        "summary": {
            "n_pairs": len(PAIRS),
            "microsomal_median_spearman": median(microsomal_rhos),
            "microsomal_min": min(microsomal_rhos) if microsomal_rhos else None,
            "microsomal_max": max(microsomal_rhos) if microsomal_rhos else None,
            "hepatocyte_median_spearman": median(hepatocyte_rhos),
            "hepatocyte_min": min(hepatocyte_rhos) if hepatocyte_rhos else None,
            "hepatocyte_max": max(hepatocyte_rhos) if hepatocyte_rhos else None,
            "grangeon_published_microsomal_median": rna["summary"][
                "published_protein_rho_median"
            ],
            "hlc_transcript_median": rna["summary"]["hlc_transcript_rho_median"],
            "pairs_where_grangeon_falls_inside_this_cohorts_interval": inside_grangeon,
            "pairs_where_grangeon_falls_outside": outside_grangeon,
            "pairs_where_the_hlc_transcript_falls_inside_this_cohorts_interval":
                inside_transcript,
        },
        "how_to_read": (
            "n = 15 gives wide intervals, so this cohort rarely EXCLUDES a "
            "value. A published value falling outside the interval is "
            "therefore a strong statement; a transcript value falling inside "
            "it means this cohort cannot tell protein from transcript as an "
            "explanation of function for that enzyme."
        ),
        "what_this_is_not": [
            "Not a model, not a prediction, not a validation.",
            "Not a like-for-like replication: the endpoints differ (activity "
            "against intrinsic clearance), the probes differ for some enzymes, "
            "and the preparations differ.",
            "Not a statement about clinical pharmacokinetics.",
            "Not a correction of Grangeon's arithmetic: Grangeon publishes no "
            "per-donor value, so its correlations cannot be recomputed at all.",
        ],
        "no_person_identifiers_values_or_predictions_exported": True,
    }
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--per-donor-csv", required=True,
                        help="the LOCAL csv emitted by the audit script")
    parser.add_argument("--audit-json", required=True,
                        help="reports/data/wegler_table_s3_per_donor_kinetics_audit.json")
    parser.add_argument("--out", default=None)
    args = parser.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        report = run(args.per_donor_csv, args.audit_json)
    except ComparisonError as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    text = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
        print(f"written: {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
