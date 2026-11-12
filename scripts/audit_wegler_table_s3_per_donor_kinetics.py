"""Recover Wegler 2021 Table S3 — per-donor kinetics — and verify it against the article.

The commissioning note asked whether the individual CLint values of the 15
Wegler donors are *really* accessible and matchable. They are: Table S3 of the
supporting information, "Kinetics and proteomics data for the 15 donors",
publishes for each of the five probe substrates and each of the 15 donors the
matching CYP concentration in microsomes, hepatocytes and homogenate, the CYP
amount in the incubation, Kpuu, and the intrinsic clearance under three
equations in microsomes and two in hepatocytes.

**They cannot be read with a text dump.** ``pdftotext -layout`` breaks the row
correspondence of this table — the clearance cells are vertically offset from
the batch labels, so a text-order reading silently assigns one donor's
clearance to another. This script therefore reconstructs the table from word
coordinates and then **verifies the reconstruction against values the article
itself published**: the twenty correlation coefficients of Table S4 that pair
each substrate with its own CYP, plus the median midazolam CLint,mic quoted in
the Results. A reconstruction that does not reproduce them is refused.

Two outputs, deliberately separated:

* ``--local-out`` receives the per-donor table, as a CSV, and must be a
  directory **outside** this repository. Donor-level rows from a third party's
  article do not belong in a public repository, whatever its licence.
* ``--out`` receives a JSON of **aggregates and verification results only**,
  plus the SHA-256 of the canonical per-donor table, which is enough to prove
  a later run recovered the same values without publishing them.

Run with ``python -I``; the PDF is a downloaded file and is passed by path.
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

#: SHA-256 of ``mp1c00053_si_001.pdf`` as recovered from the Europe PMC
#: supplementary archive for PMC8041379 on 2026-11-10.
EXPECTED_PDF_SHA256 = (
    "5beaab1117a32bb281734504dc19237bd264b5db0857837dfd794a9319892bab"
)

SUBSTRATES = ("Midazolam", "Omeprazole", "Diclofenac", "Bupropion", "Bufuralol")

#: Which CYP each substrate's block reports, from the block headers.
SUBSTRATE_CYP = {
    "Midazolam": "CYP3A4",
    "Omeprazole": "CYP2C19",
    "Diclofenac": "CYP2C9",
    "Bupropion": "CYP2B6",
    "Bufuralol": "CYP2D6",
}

#: Column names in table order, after the batch label.
COLUMNS = (
    "cyp_pmol_per_mg_HLM",
    "cyp_pmol_per_mg_HH",
    "cyp_pmol_per_mg_HL",
    "cyp_pmol_in_incubation_HLM",
    "cyp_pmol_in_incubation_HH",
    "kpuu_HH",
    "clint_mic_eq7",
    "clint_mic_eq8",
    "clint_mic_eq9",
    "clint_hep_eq7",
    "clint_hep_eq8",
)

#: Left-edge bands for the eleven value columns, read off the page geometry.
#: Numbers are right-aligned, so a wider number starts further left; the bands
#: are wide enough to absorb that and narrow enough never to merge two columns.
COLUMN_BANDS = (
    (100.0, 150.0),
    (150.0, 180.0),
    (180.0, 212.0),
    (212.0, 240.0),
    (240.0, 272.0),
    (272.0, 305.0),
    (305.0, 345.0),
    (345.0, 385.0),
    (385.0, 425.0),
    (425.0, 465.0),
    (465.0, 520.0),
)

DONORS = tuple(f"UU{index}" for index in range(1, 16))

#: Table S4, verbatim: Pearson on log-transformed values and Spearman, for the
#: pairing of each substrate with its OWN CYP. HH pairs the hepatocyte
#: concentration with the hepatocyte clearance, HLM the microsomal with the
#: microsomal. Twenty published anchors.
PUBLISHED_CORRELATIONS = {
    ("Midazolam", "HH", "r"): 0.94,
    ("Midazolam", "HH", "rs"): 0.89,
    ("Midazolam", "HLM", "r"): 0.84,
    ("Midazolam", "HLM", "rs"): 0.76,
    ("Omeprazole", "HH", "r"): 0.88,
    ("Omeprazole", "HH", "rs"): 0.90,
    ("Omeprazole", "HLM", "r"): 0.32,
    ("Omeprazole", "HLM", "rs"): 0.47,
    ("Diclofenac", "HH", "r"): 0.73,
    ("Diclofenac", "HH", "rs"): 0.75,
    ("Diclofenac", "HLM", "r"): 0.18,
    ("Diclofenac", "HLM", "rs"): 0.11,
    ("Bupropion", "HH", "r"): 0.40,
    ("Bupropion", "HH", "rs"): 0.37,
    ("Bupropion", "HLM", "r"): 0.33,
    ("Bupropion", "HLM", "rs"): 0.15,
    ("Bufuralol", "HH", "r"): 0.13,
    ("Bufuralol", "HH", "rs"): 0.36,
    ("Bufuralol", "HLM", "r"): 0.14,
    ("Bufuralol", "HLM", "rs"): 0.41,
}

#: Quoted in the Results: "The median midazolam CLint,mic was 375 mL/min/kg
#: bodyweight (bw)".
PUBLISHED_MEDIAN_MIDAZOLAM_CLINT_MIC = 375.0

#: How close a recomputed correlation must be to the published one. The
#: published values carry two decimals, so rounding alone can move a value by
#: 0.005; 0.05 leaves room for the authors' exact log base and tie handling
#: while still failing a mis-assigned row, which moves a correlation by far
#: more than that.
CORRELATION_TOLERANCE = 0.05
MEDIAN_TOLERANCE = 1.0

#: Minimum number of the twenty anchors that must agree. Set to all twenty:
#: a single disagreement means the reconstruction is not trustworthy.
REQUIRED_ANCHORS = len(PUBLISHED_CORRELATIONS)


class AuditError(RuntimeError):
    pass


def sha256_of_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def cluster_rows(page, tolerance=2.5):
    """Group a page's words into rows by the vertical centre of each word."""
    clusters = []
    for x0, y0, x1, y1, text, *_ in page.get_text("words"):
        if not text.strip():
            continue
        centre = (y0 + y1) / 2.0
        for cluster in clusters:
            if abs(cluster["centre"] - centre) <= tolerance:
                cluster["words"].append((x0, text))
                count = len(cluster["words"])
                cluster["centre"] = (cluster["centre"] * (count - 1) + centre) / count
                break
        else:
            clusters.append({"centre": centre, "words": [(x0, text)]})
    clusters.sort(key=lambda c: c["centre"])
    for cluster in clusters:
        cluster["words"].sort(key=lambda w: w[0])
    return clusters


def _value(text):
    if text.strip().upper() == "ND":
        return None
    try:
        return float(text)
    except ValueError:
        raise AuditError(f"unparseable cell {text!r} in Table S3") from None


def parse_table_s3(pdf_path):
    """Reconstruct Table S3 from word coordinates across its three pages."""
    import fitz

    rows = {}
    pages_used = []
    with fitz.open(pdf_path) as document:
        for number in range(document.page_count):
            page = document[number]
            if "Table S3" not in page.get_text():
                continue
            pages_used.append(number + 1)
            current = None
            for cluster in cluster_rows(page):
                words = cluster["words"]
                labels = [text for _, text in words]
                if len(labels) == 1 and labels[0] in SUBSTRATES:
                    current = labels[0]
                    continue
                if not labels or labels[0] not in DONORS:
                    continue
                if current is None:
                    raise AuditError(
                        f"donor row {labels[0]} on page {number + 1} precedes any "
                        "substrate heading, so its block is unknown"
                    )
                donor = labels[0]
                cells = [None] * len(COLUMNS)
                for x0, text in words[1:]:
                    for index, (low, high) in enumerate(COLUMN_BANDS):
                        if low <= x0 < high:
                            if cells[index] is not None:
                                raise AuditError(
                                    f"two values fell in column {COLUMNS[index]} "
                                    f"for {current} {donor}: the column bands do "
                                    "not match this page's geometry"
                                )
                            cells[index] = _value(text)
                            break
                    else:
                        raise AuditError(
                            f"value {text!r} at x={x0:.0f} for {current} {donor} "
                            "falls outside every column band"
                        )
                rows[(current, donor)] = dict(zip(COLUMNS, cells))
                current_row = rows[(current, donor)]
                missing = [
                    name for name in COLUMNS if name not in current_row
                ]
                if missing:
                    raise AuditError(f"{current} {donor}: columns absent {missing}")
    if sorted(pages_used) != pages_used or not pages_used:
        raise AuditError("Table S3 pages were not found in order")
    expected = {(substrate, donor) for substrate in SUBSTRATES for donor in DONORS}
    found = set(rows)
    if found != expected:
        raise AuditError(
            f"Table S3 is incomplete: {len(found)} of {len(expected)} donor rows; "
            f"missing {sorted(expected - found)[:5]}"
        )
    return rows, pages_used


# --- correlations, written out so the verification owes nothing to a library -
def _pearson(xs, ys):
    n = len(xs)
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / math.sqrt(sxx * syy)


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


def _spearman(xs, ys):
    return _pearson(_ranks(xs), _ranks(ys))


def verify_against_the_article(rows):
    """Recompute the twenty Table S4 anchors and the quoted median."""
    checks = []
    agreed = 0
    for substrate in SUBSTRATES:
        for system, concentration, clearance in (
            ("HH", "cyp_pmol_per_mg_HH", "clint_hep_eq7"),
            ("HLM", "cyp_pmol_per_mg_HLM", "clint_mic_eq7"),
        ):
            pairs = []
            for donor in DONORS:
                row = rows[(substrate, donor)]
                x = row[concentration]
                y = row[clearance]
                if x is None or y is None or x <= 0 or y <= 0:
                    continue
                pairs.append((x, y))
            if len(pairs) < 3:
                checks.append({
                    "substrate": substrate, "system": system,
                    "status": "NOT_EVALUABLE",
                    "reason": f"only {len(pairs)} usable pairs",
                })
                continue
            xs = [math.log(x) for x, _ in pairs]
            ys = [math.log(y) for _, y in pairs]
            recomputed = {
                "r": _pearson(xs, ys),
                "rs": _spearman([x for x, _ in pairs], [y for _, y in pairs]),
            }
            for kind, value in recomputed.items():
                published = PUBLISHED_CORRELATIONS[(substrate, system, kind)]
                difference = None if value is None else abs(value - published)
                ok = difference is not None and difference <= CORRELATION_TOLERANCE
                agreed += int(ok)
                checks.append({
                    "substrate": substrate,
                    "system": system,
                    "coefficient": kind,
                    "published_in_table_s4": published,
                    "recomputed": None if value is None else round(value, 4),
                    "absolute_difference": None if difference is None else round(difference, 4),
                    "n_pairs": len(pairs),
                    "agrees": bool(ok),
                })

    midazolam = sorted(
        rows[("Midazolam", donor)]["clint_mic_eq7"]
        for donor in DONORS
        if rows[("Midazolam", donor)]["clint_mic_eq7"] is not None
    )
    median = midazolam[len(midazolam) // 2] if len(midazolam) % 2 else (
        (midazolam[len(midazolam) // 2 - 1] + midazolam[len(midazolam) // 2]) / 2.0
    )
    median_ok = abs(median - PUBLISHED_MEDIAN_MIDAZOLAM_CLINT_MIC) <= MEDIAN_TOLERANCE

    return {
        "n_anchors": len(PUBLISHED_CORRELATIONS),
        "n_anchors_agreeing": agreed,
        "required": REQUIRED_ANCHORS,
        "tolerance": CORRELATION_TOLERANCE,
        "all_anchors_agree": agreed >= REQUIRED_ANCHORS,
        "checks": checks,
        "median_midazolam_clint_mic_eq7": {
            "recomputed": round(median, 3),
            "published_in_results": PUBLISHED_MEDIAN_MIDAZOLAM_CLINT_MIC,
            "agrees": bool(median_ok),
            "n_donors": len(midazolam),
        },
        "why_this_verification": (
            "pdftotext -layout breaks the row correspondence of Table S3, so a "
            "text-order reading assigns one donor's clearance to another. "
            "Reproducing twenty correlation coefficients the article itself "
            "published, each of which would move far more than the tolerance "
            "under a single mis-assigned row, is what makes the coordinate "
            "reconstruction trustworthy."
        ),
    }


def canonical_text(rows):
    """A stable serialisation of the per-donor table, for checksumming."""
    lines = ["substrate,donor," + ",".join(COLUMNS)]
    for substrate in SUBSTRATES:
        for donor in DONORS:
            row = rows[(substrate, donor)]
            cells = [
                "ND" if row[name] is None else f"{row[name]:.10g}" for name in COLUMNS
            ]
            lines.append(f"{substrate},{donor}," + ",".join(cells))
    return "\n".join(lines) + "\n"


def aggregate(rows):
    out = {}
    for substrate in SUBSTRATES:
        entry = {"cyp": SUBSTRATE_CYP[substrate], "columns": {}}
        for name in COLUMNS:
            values = [
                rows[(substrate, donor)][name]
                for donor in DONORS
                if rows[(substrate, donor)][name] is not None
            ]
            missing = len(DONORS) - len(values)
            if not values:
                entry["columns"][name] = {"n": 0, "n_not_determined": missing}
                continue
            ordered = sorted(values)
            middle = len(ordered) // 2
            median = (
                ordered[middle]
                if len(ordered) % 2
                else (ordered[middle - 1] + ordered[middle]) / 2.0
            )
            entry["columns"][name] = {
                "n": len(values),
                "n_not_determined": missing,
                "median": round(median, 4),
                "min": round(ordered[0], 4),
                "max": round(ordered[-1], 4),
                "fold_range": (
                    round(ordered[-1] / ordered[0], 2) if ordered[0] > 0 else None
                ),
            }
        out[substrate] = entry
    return out


def write_local_csv(rows, local_out):
    directory = Path(local_out).resolve()
    if REPO == directory or REPO in directory.parents:
        raise AuditError(
            "the per-donor table must not be written inside this repository: "
            f"{directory} is under {REPO}. Donor-level rows from a third "
            "party's article stay in the local data directory."
        )
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "wegler2021_table_s3_per_donor_kinetics.csv"
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["substrate", "cyp", "donor", *COLUMNS])
        for substrate in SUBSTRATES:
            for donor in DONORS:
                row = rows[(substrate, donor)]
                writer.writerow([
                    substrate, SUBSTRATE_CYP[substrate], donor,
                    *["ND" if row[name] is None else row[name] for name in COLUMNS],
                ])
    return target


def audit(pdf_path, local_out):
    checksum = sha256_of_file(pdf_path)
    if checksum != EXPECTED_PDF_SHA256:
        raise AuditError(
            "supporting-information PDF checksum mismatch\n"
            f"  expected {EXPECTED_PDF_SHA256}\n  actual   {checksum}"
        )
    rows, pages = parse_table_s3(pdf_path)
    verification = verify_against_the_article(rows)
    if not verification["all_anchors_agree"]:
        disagreeing = [
            check for check in verification["checks"] if check.get("agrees") is False
        ]
        raise AuditError(
            "the coordinate reconstruction does not reproduce the article's own "
            f"published correlations ({verification['n_anchors_agreeing']} of "
            f"{verification['n_anchors']} agree). Refusing to emit a table that "
            f"cannot be verified. First disagreements: {disagreeing[:3]}"
        )
    if not verification["median_midazolam_clint_mic_eq7"]["agrees"]:
        raise AuditError(
            "the recomputed median midazolam CLint,mic does not match the value "
            "quoted in the Results"
        )
    text = canonical_text(rows)
    written = write_local_csv(rows, local_out) if local_out else None
    return {
        "source": "Wegler et al., Mol Pharm 2021, 18(4):1792-1805, Table S3",
        "doi": "10.1021/acs.molpharmaceut.1c00053",
        "pmcid": "PMC8041379",
        "licence_recorded_by_europe_pmc": "cc by",
        "licence_caveat": (
            "Europe PMC records the licence as CC BY. The supporting-information "
            "PDF carries no licence statement of its own, so this is the "
            "repository's assertion and not a term read from the document. "
            "Silence about terms is not permission: the per-donor rows are kept "
            "out of this repository regardless."
        ),
        "supporting_information_pdf_sha256": checksum,
        "table_s3_pages": pages,
        "n_substrates": len(SUBSTRATES),
        "n_donors": len(DONORS),
        "substrate_to_cyp": dict(SUBSTRATE_CYP),
        "columns": list(COLUMNS),
        "n_values_recovered": len(rows) * len(COLUMNS),
        "per_donor_table_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "per_donor_table_written_to": str(written) if written else None,
        "per_donor_rows_are_not_in_this_repository": True,
        "verification": verification,
        "aggregates_only": aggregate(rows),
        "what_table_s1_does_not_give": (
            "Table S1 reports donor characteristics as aggregates only: age mean "
            "64 (39-79), BMI mean 25.7 (20.1-32.9), 10 male and 5 female, 13 "
            "colorectal cancer, 1 breast cancer, 1 GIST. There is no per-donor "
            "age, sex, BMI or diagnosis, so no per-donor demographic baseline "
            "block can be built for these 15 donors."
        ),
        "why_a_text_dump_will_not_do": (
            "pdftotext -layout misaligns the clearance cells against the batch "
            "labels in this table. Any transcription that does not verify "
            "against the article's own published statistics should be assumed "
            "to have mis-assigned donors."
        ),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pdf", required=True, help="mp1c00053_si_001.pdf")
    parser.add_argument(
        "--local-out",
        default=None,
        help="directory OUTSIDE this repository for the per-donor CSV",
    )
    parser.add_argument("--out", default=None, help="aggregate-only JSON")
    args = parser.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        report = audit(args.pdf, args.local_out)
    except AuditError as exc:
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
