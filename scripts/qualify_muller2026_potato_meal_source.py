"""Acquire and qualify the Mueller 2026 standardized-potato-meal source.

The source publication releases its data openly:
  Mueller JP et al. Solanidine metabolic ratios for CYP2D6 phenotyping after a
  standardized potato meal. Front Pharmacol. 2026 Oct 2;17:1886427.
  doi 10.3389/fphar.2026.1886427
  https://github.com/JulianPeterMueller/Solanidine-standardized-potato-meal
  "Solanidine standardized potato meal.xlsx"
  sha256 68d53071db79f961a63e433137a8a3e7dbd40daf8c8ad0f9f7f1b268bb3a8e07

This script verifies that checksum and qualifies the workbook for the public
data evidence registry. It is a **qualification** tool, not an analysis tool:
it reports structure, completeness and analytical quality, and deliberately
computes no H5 result, because the H5 criteria were frozen in protocol v1.0.0
before this dataset was known to exist and must stay that way.

The one thing it does that a value-only reader cannot: the source imputed
values below the limit of quantification and flagged them **only with bold
font**. A pandas read sees imputed floors as ordinary measurements. This
script extracts the font flags with openpyxl and counts them per analyte and
per participant, so the registry can record how much of the signal is floor.

No individual measurement is printed: only counts, structure and aggregate
quality statistics.

Usage
-----
    python scripts/qualify_muller2026_potato_meal_source.py --source PATH
    python scripts/qualify_muller2026_potato_meal_source.py --download DIR
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

DATA_URL = (
    "https://raw.githubusercontent.com/JulianPeterMueller/"
    "Solanidine-standardized-potato-meal/HEAD/"
    "Solanidine%20standardized%20potato%20meal.xlsx"
)
EXPECTED_SHA256 = "68d53071db79f961a63e433137a8a3e7dbd40daf8c8ad0f9f7f1b268bb3a8e07"
EXPECTED_BYTES = 20965

SHEET_PARTICIPANTS = "Participant data"
SHEET_PLASMA = "Plasma concentration data"
SHEET_GLYCOALKALOID = "Glycoalkaloid concentrations"
SHEET_MICROBIAL = "Microbial metabolism"
EXPECTED_SHEETS = (
    SHEET_PARTICIPANTS,
    SHEET_PLASMA,
    SHEET_GLYCOALKALOID,
    SHEET_MICROBIAL,
)

#: First data row (1-indexed, openpyxl) of each stacked analyte block, and the
#: nine sampling times in hours that head every block.
PLASMA_BLOCKS = {
    "solanidine_ng_ml": 3,
    "oh4_solanidine_ng_ml": 26,
    "ssda_counts": 49,
}
TIMEPOINTS_H = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 29.0)
N_PARTICIPANTS = 20

COL_PARTICIPANT = "Participant-No."
COL_GENOTYPE = "CYP2D6 genotype"
COL_PHENOTYPE = "CYP2D6 genotype-predicted phenotype"
COL_ACTIVITY_SCORE = "CYP2D6 activity score"
COL_FASTED = "Fasted on potato products for 3 days (Yes=1, No=0)"

#: SSDA is reported in raw counts, so it is semi-quantitative. Recorded here
#: because it is the same analytical weakness already logged for ARIANA.
SEMI_QUANTITATIVE_ANALYTES = ("ssda_counts",)

#: The source states its own imputation rule; pinned so a change is reviewable.
SOURCE_LOQ_RULE = (
    "Values below limit of quantification were imputed using half of the lower "
    "limit of quantification for solanidine (0.01 ng/ml) and 4-OH-solanidine "
    "(0.02 ng/ml) or the lowest measured values for SSDA."
)

#: This cohort scores *41 as 0.25 (post-2021 CPIC). The Rudesheim 2020 mapping
#: used by the mechanistic module is built on a grid that predates that change,
#: so activity scores are not comparable across cohorts without alignment.
RUDESHEIM_AS_GRID = (0.0, 0.5, 1.25, 1.5, 2.0, 3.0)


class QualificationError(RuntimeError):
    pass


def sha256_of(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_source(path):
    """Return the checksum, raising if it does not match the published file."""
    path = Path(path)
    actual = sha256_of(path)
    if actual != EXPECTED_SHA256:
        raise QualificationError(
            f"checksum mismatch for {path}\n"
            f"  expected {EXPECTED_SHA256}\n  actual   {actual}"
        )
    return actual


def download(target_dir):
    """Fetch the public workbook into *target_dir* and verify its checksum."""
    from urllib.request import urlopen

    target_dir = Path(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / "Solanidine standardized potato meal.xlsx"
    with urlopen(DATA_URL, timeout=90) as response:  # noqa: S310 - fixed https URL
        target.write_bytes(response.read())
    verify_source(target)
    return target


def _bold(cell):
    return bool(cell.font is not None and cell.font.bold)


def qualify_plasma_blocks(worksheet):
    """Count present, missing and LOQ-imputed cells per analyte block.

    The imputed cells are identifiable **only** by bold font, which is why
    this reads the worksheet through openpyxl rather than pandas.
    """
    out = {}
    for analyte, first_row in PLASMA_BLOCKS.items():
        ids, present, imputed, missing = [], 0, 0, 0
        per_participant = []
        for offset in range(N_PARTICIPANTS):
            row = first_row + offset
            participant = worksheet.cell(row=row, column=1).value
            ids.append(participant)
            imputed_here = 0
            for index in range(len(TIMEPOINTS_H)):
                cell = worksheet.cell(row=row, column=2 + index)
                if cell.value is None:
                    missing += 1
                    continue
                present += 1
                if _bold(cell):
                    imputed += 1
                    imputed_here += 1
            per_participant.append(
                {"participant": participant, "n_imputed": imputed_here}
            )
        expected = N_PARTICIPANTS * len(TIMEPOINTS_H)
        fully_imputed = [
            record["participant"]
            for record in per_participant
            if record["n_imputed"] == len(TIMEPOINTS_H)
        ]
        out[analyte] = {
            "participant_ids": ids,
            "participant_ids_are_1_to_n": ids == list(range(1, N_PARTICIPANTS + 1)),
            "cells_expected": expected,
            "cells_present": present,
            "cells_missing": missing,
            "complete": missing == 0 and present == expected,
            "loq_imputed_cells": imputed,
            "loq_imputed_percent": round(100.0 * imputed / max(present, 1), 1),
            "participants_fully_below_loq": fully_imputed,
            "semi_quantitative": analyte in SEMI_QUANTITATIVE_ANALYTES,
            "most_imputed_participants": sorted(
                per_participant, key=lambda r: -r["n_imputed"]
            )[:3],
        }
    return out


def qualify(source):
    """Return a JSON-serialisable qualification report."""
    import pandas as pd
    from openpyxl import load_workbook

    checksum = verify_source(source)
    size = Path(source).stat().st_size

    excel = pd.ExcelFile(source)
    missing_sheets = [s for s in EXPECTED_SHEETS if s not in excel.sheet_names]
    if missing_sheets:
        raise QualificationError(f"source workbook is missing sheets: {missing_sheets}")

    participants = excel.parse(SHEET_PARTICIPANTS)
    required = [COL_PARTICIPANT, COL_GENOTYPE, COL_PHENOTYPE,
                COL_ACTIVITY_SCORE, COL_FASTED]
    absent = [c for c in required if c not in participants.columns]
    if absent:
        raise QualificationError(f"participant sheet is missing columns: {absent}")

    workbook = load_workbook(source, data_only=True)
    plasma = qualify_plasma_blocks(workbook[SHEET_PLASMA])

    activity_scores = sorted(float(v) for v in participants[COL_ACTIVITY_SCORE].unique())
    fasted = participants[COL_FASTED].value_counts().to_dict()

    report = {
        "source_sha256": checksum,
        "source_bytes": size,
        "source_bytes_match": size == EXPECTED_BYTES,
        "doi": "10.3389/fphar.2026.1886427",
        "sheets": list(excel.sheet_names),
        "endpoint_class": "biomarker_only",
        "no_drug_administered": True,
        "cannot_support": [
            "H1",
            "any pharmacokinetic claim",
            "replication of the ARIANA solanidine-to-exposure association",
        ],
        "objectives_supported": ["H5"],
        "participants": {
            "n": int(len(participants)),
            "n_genotyped": int(participants[COL_GENOTYPE].notna().sum()),
            "phenotype_counts": {
                str(k): int(v)
                for k, v in participants[COL_PHENOTYPE].value_counts().items()
            },
            "activity_score_values_present": activity_scores,
            "fasted_on_potato_products_counts": {str(k): int(v) for k, v in fasted.items()},
        },
        "sampling": {
            "timepoints_h": list(TIMEPOINTS_H),
            "observations_per_subject": len(TIMEPOINTS_H),
            "absolute_times_documented": True,
            "repeated_sampling": True,
        },
        "plasma_blocks": plasma,
        "loq_imputation": {
            "present": True,
            "source_rule": SOURCE_LOQ_RULE,
            "flagged_only_by_bold_font": True,
            "why_this_matters": (
                "Imputed cells carry no data flag. A value-only reader treats the "
                "imputed floor as a measurement, which would overstate the "
                "marker's resolution at the low end."
            ),
        },
        "activity_score_convention": {
            "this_cohort_scores_star41_as": 0.25,
            "rudesheim_2020_grid": list(RUDESHEIM_AS_GRID),
            "comparable_with_ariana_without_alignment": False,
            "note": (
                "Post-2021 CPIC scores *41 as 0.25; the Rudesheim mapping predates "
                "that change. Aligning the conventions is a reported decision, never "
                "a silent transformation."
            ),
        },
        "independence_limitation": (
            "Saromba, Stingl and Just author both this study and ARIANA, so this is a "
            "different population measured by the same group, laboratory and platform "
            "- not an independent laboratory. Participants are different people, so "
            "donors are not double-counted."
        ),
    }

    report["qualified"] = (
        report["source_bytes_match"]
        and not missing_sheets
        and all(block["complete"] for block in plasma.values())
        and all(block["participant_ids_are_1_to_n"] for block in plasma.values())
        and report["participants"]["n"] == N_PARTICIPANTS
        and report["participants"]["n_genotyped"] == N_PARTICIPANTS
    )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--source", help="path to the public workbook")
    group.add_argument("--download", metavar="DIR", help="download it into DIR first")
    args = parser.parse_args(argv)

    try:
        source = download(args.download) if args.download else args.source
        report = qualify(source)
    except QualificationError as exc:
        parser.exit(2, f"INVALID: {exc}\n")

    print(json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True))
    return 0 if report["qualified"] else 1


if __name__ == "__main__":
    sys.exit(main())
