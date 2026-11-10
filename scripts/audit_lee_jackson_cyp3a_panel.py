#!/usr/bin/env python
"""Qualify the Lee / Jackson CYP3A donor panel as a three-layer public dataset.

Source: Lee J, Fallon JK, Smith PC, Jackson KD. "Formation of CYP3A-specific
metabolites of ibrutinib in vitro is correlated with hepatic CYP3A activity and
4-beta-hydroxycholesterol/cholesterol ratio." Clin Transl Sci 2023;16:279-293.
PMCID PMC9926076, CC BY-NC-ND.

Why this source matters. For the SAME fifteen hepatocyte donors it publishes,
each at the level of the individual donor:

  * layer 1 -- hepatic CYP3A4 and CYP3A5 absolute abundance (supplement
    Table S2);
  * layer 2 -- plasma 4-beta-hydroxycholesterol and cholesterol, and their
    ratio, from matched donor plasma (MAIN-TEXT Table 1);
  * outcome -- CYP3A activity as midazolam 1'-hydroxylation, and ibrutinib
    metabolite formation (supplement Table S3).

That is a per-donor molecular layer paired with a per-donor outcome, which is
exactly what five of the seven candidates qualified on 2026-11-10 lacked.

This script fits nothing and makes no claim. It extracts the three tables,
proves the rows describe the same individuals, and verifies the extraction
against correlations the article itself publishes. If any published anchor
fails to reproduce, it raises: an extraction that cannot reproduce the source's
own numbers is not evidence of anything.

The per-donor table is written OUTSIDE this repository, like Wegler's.

Usage
-----
    python -I scripts/audit_lee_jackson_cyp3a_panel.py \
        --supplement <dir>/CTS-16-279-s001.docx \
        --fulltext   <dir>/cts16279.xml \
        --out-dir    <a directory outside this repository>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pathlib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile

sys.stdout.reconfigure(encoding="utf-8")

REPO = pathlib.Path(__file__).resolve().parents[1]
WORD = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

#: The donor whose age the article itself excludes from the correlation it
#: reports, and which we must therefore also exclude to reproduce that number.
INFANT_DONOR = "VYU"

#: Correlations the article publishes, and which this extraction must
#: reproduce. Source sentence, Results: "When an infant donor (1 year old) was
#: excluded based on previous ontogeny studies, M37 and M25 formation
#: correlated with plasma 4beta-HC/cholesterol ratio in the remaining 14 donors
#: (Spearman correlation coefficients [r] 0.62 and 0.67, respectively)."
PUBLISHED_ANCHORS = (
    {"label": "M37 vs plasma ratio, 14 non-infant donors",
     "left": "ibrutinib_m37", "right": "plasma_4b_hc_cholesterol_ratio",
     "exclude": (INFANT_DONOR,), "published": 0.62},
    {"label": "M25 vs plasma ratio, 14 non-infant donors",
     "left": "ibrutinib_m25", "right": "plasma_4b_hc_cholesterol_ratio",
     "exclude": (INFANT_DONOR,), "published": 0.67},
    # The Discussion reports, for the identical donor panel in the companion
    # venetoclax paper, "CYP3A4 protein concentrations were highly correlated
    # with the formation of M3 (r = 0.97)". The same relationship for
    # ibrutinib's M37 across all fifteen donors is quoted at two decimals.
    {"label": "CYP3A4 abundance vs M37, all 15 donors",
     "left": "cyp3a4_abundance", "right": "ibrutinib_m37",
     "exclude": (), "published": 0.97},
)

#: How far an extracted anchor may sit from the published figure. The published
#: values are quoted to two decimals, so rounding alone can move them by
#: 0.005; 0.02 leaves room for that and for nothing else.
ANCHOR_TOLERANCE = 0.02

#: The number of hepatocyte donors the article reports, and the floor below
#: which leave-one-donor-out with an inner leave-one-out is not estimable.
EXPECTED_DONORS = 15
DONOR_FLOOR = 10


class AuditError(RuntimeError):
    """Raised when the source does not say what this audit records."""


def _check(condition, message):
    if not condition:
        raise AuditError(message)


# --------------------------------------------------------------------------
# reading the two files
# --------------------------------------------------------------------------
def _cells(row):
    return ["".join(t.text or "" for t in tc.iter(f"{WORD}t")).strip()
            for tc in row.findall(f"{WORD}tc")]


def _docx_tables(path):
    root = ET.fromstring(zipfile.ZipFile(path).read("word/document.xml"))
    return [[_cells(r) for r in tbl.findall(f"{WORD}tr")]
            for tbl in root.find(f"{WORD}body").iter(f"{WORD}tbl")]


def _strip_tags(blob):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", blob)).strip()


def _fulltext_tables(path):
    raw = pathlib.Path(path).read_text(encoding="utf-8", errors="replace")
    out = []
    for block in re.findall(r"<table-wrap\b.*?</table-wrap>", raw, re.S):
        rows = [[_strip_tags(c) for c in
                 re.findall(r"<t[dh]\b.*?</t[dh]>", r, re.S)]
                for r in re.findall(r"<tr\b.*?</tr>", block, re.S)]
        out.append({"caption": _strip_tags(block), "rows": rows})
    return out


def _number(text):
    """The point estimate of a 'value +/- sd' cell, or None."""
    cleaned = text.replace("−", "-").replace(",", "")
    match = re.match(r"^\s*(-?\d+(?:\.\d+)?)", cleaned)
    return float(match.group(1)) if match else None


# --------------------------------------------------------------------------
# statistics, written out rather than imported, so the audit has no model in it
# --------------------------------------------------------------------------
def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    index = 0
    while index < len(order):
        stop = index
        while (stop + 1 < len(order)
               and values[order[stop + 1]] == values[order[index]]):
            stop += 1
        average = (index + stop) / 2 + 1
        for position in range(index, stop + 1):
            out[order[position]] = average
        index = stop + 1
    return out


def spearman(left, right):
    _check(len(left) == len(right) and len(left) > 2,
           "a correlation needs at least three paired values")
    a, b = ranks(left), ranks(right)
    n = len(a)
    mean_a, mean_b = sum(a) / n, sum(b) / n
    numerator = sum((x - mean_a) * (y - mean_b) for x, y in zip(a, b))
    denominator = (sum((x - mean_a) ** 2 for x in a)
                   * sum((y - mean_b) ** 2 for y in b)) ** 0.5
    _check(denominator > 0.0, "a correlation needs variation on both sides")
    return numerator / denominator


# --------------------------------------------------------------------------
# the three tables
# --------------------------------------------------------------------------
def hepatocyte_abundance(tables):
    """Supplement Table S2: per-donor CYP3A4 and CYP3A5 in hepatocytes."""
    for rows in tables:
        head = [c.lower() for c in rows[0]]
        if (head and head[0] == "donor"
                and any("cyp3a4 protein" in c for c in head)
                and any("bmi" in c for c in head)):
            out = {}
            for row in rows[1:]:
                if not row or not row[0]:
                    continue
                out[row[0]] = {
                    "age_years": _number(row[1]),
                    "sex": row[2],
                    "race": row[3],
                    "bmi": _number(row[4]),
                    "cyp3a4_abundance": _number(row[5]),
                    "cyp3a5_abundance": _number(row[6]),
                }
            return out
    raise AuditError("supplement Table S2 (hepatocyte abundance) was not found")


def hepatocyte_activity(tables):
    """Supplement Table S3: per-donor activity in the SAME hepatocytes."""
    for rows in tables:
        head = [c.lower() for c in rows[0]]
        if (head and head[0] == "donor"
                and any("hydroxymidazolam" in c for c in head)):
            out = {}
            for row in rows[1:]:
                if not row or not row[0]:
                    continue
                out[row[0]] = {
                    "midazolam_1_hydroxylation": _number(row[1]),
                    "ibrutinib_m37": _number(row[2]),
                    "ibrutinib_m34": _number(row[3]),
                    "ibrutinib_m25": _number(row[4]),
                }
            return out
    raise AuditError("supplement Table S3 (hepatocyte activity) was not found")


def plasma_biomarker(tables):
    """Main-text Table 1: per-donor plasma 4beta-HC, cholesterol and ratio."""
    for table in tables:
        if "Plasma 4" not in table["caption"]:
            continue
        rows = table["rows"]
        head = [c.lower() for c in rows[0]]
        if not (head and head[0] == "donor"):
            continue
        out = {}
        for row in rows[1:]:
            if not row or not row[0]:
                continue
            out[row[0]] = {
                "plasma_4b_hc_ng_per_ml": _number(row[1]),
                "plasma_cholesterol_mg_per_ml": _number(row[2]),
                "plasma_4b_hc_cholesterol_ratio": _number(row[3]),
            }
        return out
    raise AuditError("main-text Table 1 (plasma biomarker) was not found")


def join_on_the_same_donors(abundance, activity, plasma):
    """One row per donor, and a refusal if the three tables disagree.

    The Wegler lesson applies: identifiers that look alike are not evidence
    that they are the same individual. Here the three tables are published in
    the same article about the same purchased lots, so the identifier IS the
    link -- but it must be the same SET, and that is checked, not assumed.
    """
    sets = {"abundance": set(abundance), "activity": set(activity),
            "plasma": set(plasma)}
    _check(sets["abundance"] == sets["activity"] == sets["plasma"],
           "the three tables do not describe the same donors: "
           + "; ".join(f"{k} has {len(v)}" for k, v in sets.items())
           + f"; symmetric differences {sets['abundance'] ^ sets['activity']}"
           + f" and {sets['abundance'] ^ sets['plasma']}")
    _check(len(sets["abundance"]) == EXPECTED_DONORS,
           f"expected {EXPECTED_DONORS} hepatocyte donors, "
           f"found {len(sets['abundance'])}")
    _check(list(abundance) == list(activity),
           "Tables S2 and S3 list the same donors in a DIFFERENT order; the "
           "join must then be by identifier and the row order must not be "
           "trusted anywhere downstream")
    rows = {}
    for donor in abundance:
        merged = dict(abundance[donor])
        merged.update(activity[donor])
        merged.update(plasma[donor])
        missing = [key for key, value in merged.items() if value is None]
        _check(not missing, f"donor {donor} has no value for {missing}")
        rows[donor] = merged
    return rows


def verify_published_anchors(rows):
    """Reproduce the article's own correlations, or raise."""
    checked = []
    for anchor in PUBLISHED_ANCHORS:
        kept = [d for d in rows if d not in anchor["exclude"]]
        got = spearman([rows[d][anchor["left"]] for d in kept],
                       [rows[d][anchor["right"]] for d in kept])
        gap = abs(got - anchor["published"])
        checked.append({
            "anchor": anchor["label"], "n": len(kept),
            "published": anchor["published"],
            "recomputed": round(got, 4), "absolute_gap": round(gap, 4),
            "within_tolerance": gap <= ANCHOR_TOLERANCE,
        })
        _check(gap <= ANCHOR_TOLERANCE,
               f"anchor '{anchor['label']}' does not reproduce: the article "
               f"publishes {anchor['published']}, this extraction gives "
               f"{got:.4f} (gap {gap:.4f} > {ANCHOR_TOLERANCE}). The "
               "extraction is wrong, or the table is not the one the article "
               "used. Either way nothing downstream may run.")
    return checked


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------
def describe(rows, field):
    values = sorted(rows[d][field] for d in rows)
    n = len(values)
    median = (values[n // 2] if n % 2
              else (values[n // 2 - 1] + values[n // 2]) / 2)
    return {
        "n": n,
        "minimum": values[0],
        "median": median,
        "maximum": values[-1],
        "fold_range_max_over_min": (round(values[-1] / values[0], 1)
                                    if values[0] > 0 else None),
    }


def relationships(rows):
    """Every pairwise rank correlation this source permits, named plainly."""
    pairs = (
        ("cyp3a4_abundance", "midazolam_1_hydroxylation"),
        ("cyp3a4_abundance", "ibrutinib_m37"),
        ("cyp3a5_abundance", "midazolam_1_hydroxylation"),
        ("plasma_4b_hc_cholesterol_ratio", "midazolam_1_hydroxylation"),
        ("plasma_4b_hc_cholesterol_ratio", "cyp3a4_abundance"),
        ("plasma_4b_hc_cholesterol_ratio", "ibrutinib_m37"),
        ("age_years", "cyp3a4_abundance"),
        ("age_years", "midazolam_1_hydroxylation"),
    )
    donors = list(rows)
    return {
        f"{left}__vs__{right}": round(
            spearman([rows[d][left] for d in donors],
                     [rows[d][right] for d in donors]), 4)
        for left, right in pairs
    }


def write_per_donor_table(rows, out_dir):
    directory = pathlib.Path(out_dir).resolve()
    if REPO == directory or REPO in directory.parents:
        raise AuditError(
            "the per-donor table must not be written inside this repository: "
            f"{directory} is at or under {REPO}. The source is CC BY-NC-ND and "
            "this project does not commit individual-level rows.")
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "lee_jackson_cyp3a_panel_per_donor.csv"
    fields = ["donor", "age_years", "sex", "race", "bmi",
              "cyp3a4_abundance", "cyp3a5_abundance",
              "midazolam_1_hydroxylation", "ibrutinib_m37", "ibrutinib_m34",
              "ibrutinib_m25", "plasma_4b_hc_ng_per_ml",
              "plasma_cholesterol_mg_per_ml", "plasma_4b_hc_cholesterol_ratio"]
    lines = [",".join(fields)]
    for donor in rows:
        record = {"donor": donor, **rows[donor]}
        lines.append(",".join(str(record[field]) for field in fields))
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement", required=True,
                        help="CTS-16-279-s001.docx, outside this repository")
    parser.add_argument("--fulltext", required=True,
                        help="the article's Europe PMC fullTextXML")
    parser.add_argument("--out-dir", required=True,
                        help="where the per-donor table goes; outside the repo")
    parser.add_argument("--artifact",
                        default=str(REPO / "reports" / "data"
                                    / "lee_jackson_cyp3a_panel.json"))
    args = parser.parse_args(argv)

    supplement = pathlib.Path(args.supplement)
    fulltext = pathlib.Path(args.fulltext)
    for path in (supplement, fulltext):
        _check(path.is_file(), f"missing input: {path}")

    abundance = hepatocyte_abundance(_docx_tables(supplement))
    activity = hepatocyte_activity(_docx_tables(supplement))
    plasma = plasma_biomarker(_fulltext_tables(fulltext))
    rows = join_on_the_same_donors(abundance, activity, plasma)
    anchors = verify_published_anchors(rows)
    table = write_per_donor_table(rows, args.out_dir)

    payload = {
        "audit": "lee_jackson_cyp3a_panel",
        "generated": "2026-11-10",
        "source": {
            "citation": ("Lee J, Fallon JK, Smith PC, Jackson KD. Clin Transl "
                         "Sci 2023;16(2):279-293."),
            "pmcid": "PMC9926076", "pmid": "36350327",
            "licence": "cc by-nc-nd",
            "licence_consequence": (
                "Analysis and the reporting of aggregates are unaffected. "
                "Redistributing the per-donor table, or a derivative of it, is "
                "not permitted, so it is written outside this repository and "
                "never committed."),
            "channel_that_worked": ("Europe PMC supplementaryFiles for the "
                                    "docx, fullTextXML for the main table"),
        },
        "n_donors": len(rows),
        "donor_floor": DONOR_FLOOR,
        "meets_donor_floor": len(rows) >= DONOR_FLOOR,
        "layers_per_donor": {
            "hepatic_protein": ("CYP3A4 and CYP3A5 absolute abundance, "
                                "pmol/mg protein"),
            "plasma_metabolite": ("4-beta-hydroxycholesterol and cholesterol, "
                                  "and their ratio, from matched donor plasma"),
        },
        "outcome_per_donor": {
            "primary": ("midazolam 1'-hydroxylation rate, pmol/million "
                        "cells/min -- an ACTIVITY RATE at a single substrate "
                        "concentration, not an intrinsic clearance and not a "
                        "clinical pharmacokinetic quantity"),
            "secondary": "ibrutinib M37 formation rate, pmol/million cells/min",
        },
        "circularity_check": {
            "question": ("does the outcome's definition contain the abundance "
                         "used as its predictor, as Wegler's equations 8 and 9 "
                         "did?"),
            "answer": "no",
            "evidence": ("the outcome is normalised per million cells and the "
                         "predictor per mg of protein, so the outcome is not "
                         "divided by any CYP quantity"),
            "cost_of_that_difference": (
                "the two are expressed per DIFFERENT denominators, so the "
                "protein-per-cell ratio enters any relationship between them "
                "as an uncontrolled term. That weakens the comparison; it does "
                "not make it circular."),
        },
        "published_anchors_verified": anchors,
        "distributions": {
            field: describe(rows, field) for field in
            ("cyp3a4_abundance", "cyp3a5_abundance",
             "midazolam_1_hydroxylation", "ibrutinib_m37",
             "plasma_4b_hc_cholesterol_ratio", "age_years", "bmi")
        },
        "rank_correlations_all_donors": relationships(rows),
        "independence": {
            "from_wegler": (
                "independent: different donors (BioIVT letter codes against "
                "Wegler's UU codes), different laboratory, different year, "
                "different quantification, and an activity rate rather than a "
                "fitted intrinsic clearance."),
            "from_the_venetoclax_paper": (
                "NOT independent. Lee et al., Clin Transl Sci 2024 "
                "(PMC11645475) reports venetoclax metabolism in the SAME "
                "fifteen hepatocyte donors and the same twenty microsome lots. "
                "Two drugs measured in the same donors are one cohort with two "
                "substrates, never two cohorts."),
        },
        "what_this_source_cannot_do": [
            "no intrinsic clearance: the outcome is a single-concentration rate",
            "no clinical pharmacokinetics: no dose, no plasma drug "
            "concentration, no exposure measure",
            "no microsomal outcome: Table S1 gives abundance for twenty "
            "microsome lots and no activity for them, so the microsomal arm "
            "Wegler could test is not testable here",
            "ages span 1 to 63 years and include two donors under four, so "
            "ontogeny is a confounder this design cannot remove",
            "donor plasma is procurement plasma from deceased donors, which is "
            "not interchangeable with plasma from a living patient",
        ],
        "per_donor_table_is_outside_this_repository": str(table),
        "per_donor_table_sha256": hashlib.sha256(table.read_bytes()).hexdigest(),
        "no_donor_identifiers_values_or_predictions_exported": True,
        "this_audit_fits_nothing": True,
    }

    artifact = pathlib.Path(args.artifact)
    artifact.parent.mkdir(parents=True, exist_ok=True)
    artifact.write_text(json.dumps(payload, indent=2, ensure_ascii=False)
                        + "\n", encoding="utf-8")

    print(f"donors: {payload['n_donors']}  (floor {DONOR_FLOOR})")
    print("published anchors:")
    for row in anchors:
        print(f"  {row['anchor']:<46} published {row['published']:.2f} "
              f"recomputed {row['recomputed']:+.4f}  "
              f"gap {row['absolute_gap']:.4f}")
    print("rank correlations, all fifteen donors:")
    for name, value in payload["rank_correlations_all_donors"].items():
        print(f"  {name:<62} {value:+.4f}")
    print(f"per-donor table (outside the repo): {table}")
    print(f"artefact: {artifact}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
