#!/usr/bin/env python
"""Qualify and execute the venetoclax clinical-PK contrast.

Protocol: ``docs/OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_PRESPECIFICATION_20261110.md``

Reads the per-patient dose and venetoclax concentration from the article's
main-text table (Europe PMC ``fullTextXML``) and the per-patient plasma
4-beta-hydroxycholesterol, cholesterol and ratio from supplement Table S3
(retrieved through BioStudies and the EBI FTP mirror). Verifies the extraction
against three exposure figures the article publishes, then runs the frozen
contrast.

No patient identifier, value or prediction is written to the artefact.

Usage:
    python -I scripts/run_venetoclax_clinical_pk_contrast.py \
        --fulltext <dir>/PMC11645475.xml \
        --supplement <dir>/CTS-17-e70106-s001.docx \
        --out-dir <a directory outside this repository>
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

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8")

from omics2pk import venetoclax_clinical_pk_contrast as VCP  # noqa: E402

WORD = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
PRESPEC = (ROOT / "docs"
           / "OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_PRESPECIFICATION_20261110.md")

#: Values the article marks as not determined. They are recorded and excluded
#: by a rule declared before the run, never dropped silently.
NOT_DETERMINED = {"ND", "NA", "Completed treatment", ""}


class RunError(RuntimeError):
    pass


def _clean(blob):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", blob)).strip()


def _number(text):
    cleaned = (text or "").replace("−", "-").replace(",", "")
    match = re.match(r"^\s*(-?\d+(?:\.\d+)?)", cleaned)
    return float(match.group(1)) if match else None


def read_exposure(fulltext):
    """Per-patient dose and venetoclax concentration from the main-text table.

    The table carries the patient number only on its first row per patient, so
    the reader must carry it forward -- the same shape of trap that mis-assigned
    Wegler's donors.
    """
    raw = pathlib.Path(fulltext).read_text(encoding="utf-8", errors="replace")
    blocks = [b for b in re.findall(r"<table-wrap\b.*?</table-wrap>", raw, re.S)
              if "daily dose" in _clean(b)]
    if not blocks:
        raise RunError("the per-patient dose and concentration table "
                       "was not found in the full text")
    rows, current, excluded = {}, None, []
    for tr in re.findall(r"<tr\b.*?</tr>", blocks[0], re.S):
        cells = [_clean(c) for c in
                 re.findall(r"<t[dh]\b.*?</t[dh]>", tr, re.S)]
        if len(cells) < 4 or cells[0].lower().startswith("patient"):
            continue
        if cells[0] and cells[0].isdigit() and cells[1] in {"1", "2"}:
            current, visit = cells[0], int(cells[1])
            dose, concentration = cells[2], cells[3]
        elif cells[0] in {"1", "2"} and current is not None:
            visit = int(cells[0])
            dose, concentration = cells[1], cells[2]
        else:
            continue
        if dose in NOT_DETERMINED or concentration in NOT_DETERMINED:
            excluded.append({"visit": visit,
                             "reason": f"dose={dose!r} concentration="
                                       f"{concentration!r}"})
            continue
        dose_value, concentration_value = _number(dose), _number(concentration)
        if dose_value is None or concentration_value is None:
            excluded.append({"visit": visit,
                             "reason": f"unparsable dose={dose!r} "
                                       f"concentration={concentration!r}"})
            continue
        rows[(current, visit)] = {"dose_mg_per_day": dose_value,
                                  "concentration_nm": concentration_value}
    return rows, excluded


def read_biomarker(supplement):
    """Per-patient plasma 4beta-HC, cholesterol and ratio from Table S3."""
    root = ET.fromstring(
        zipfile.ZipFile(supplement).read("word/document.xml"))
    for table in root.find(f"{WORD}body").iter(f"{WORD}tbl"):
        rows = [["".join(t.text or "" for t in tc.iter(f"{WORD}t")).strip()
                 for tc in tr.findall(f"{WORD}tc")]
                for tr in table.findall(f"{WORD}tr")]
        head = [c.lower() for c in rows[0]] if rows else []
        if not (head and head[0] == "patient" and "visit" in head[1]):
            continue
        out, current = {}, None
        for cells in rows[1:]:
            if len(cells) < 5:
                continue
            if cells[0]:
                current = cells[0]
            if current is None:
                continue
            visit = int(cells[1])
            ratio = _number(cells[4])
            if ratio is None:
                continue
            out[(current, visit)] = {
                "plasma_4b_hc_ng_per_ml": _number(cells[2]),
                "plasma_cholesterol_mg_per_ml": _number(cells[3]),
                "plasma_4b_hc_cholesterol_ratio": ratio,
            }
        return out
    raise RunError("supplement Table S3 (per-patient plasma biomarker) "
                   "was not found")


def join(exposure, biomarker):
    keys = sorted(set(exposure) & set(biomarker))
    if not keys:
        raise RunError("the exposure and biomarker tables share no "
                       "patient-visit key")
    rows = {k: {**exposure[k], **biomarker[k]} for k in keys}
    for key, record in rows.items():
        missing = [f for f, v in record.items() if v is None]
        if missing:
            raise RunError(f"patient-visit {key} has no value for {missing}")
    return rows, keys


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fulltext", required=True)
    parser.add_argument("--supplement", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "venetoclax_clinical_pk_contrast.json"))
    args = parser.parse_args(argv)

    if not PRESPEC.is_file():
        raise RunError(f"the frozen protocol is missing: {PRESPEC}")

    exposure, excluded = read_exposure(args.fulltext)
    biomarker = read_biomarker(args.supplement)
    rows, keys = join(exposure, biomarker)
    anchors = VCP.verify_published_anchors(rows)

    visit_one = [k for k in keys if k[1] == 1]
    if len(visit_one) < 10:
        raise RunError(f"the primary analysis has {len(visit_one)} patients, "
                       "below the declared floor of ten")

    primary = VCP.contrast(rows, visit_one,
                           family_size=VCP.FROZEN["declared_family_size"])
    sensitivity = VCP.contrast(rows, keys,
                               family_size=VCP.FROZEN["declared_family_size"])

    out_dir = pathlib.Path(args.out_dir).resolve()
    if ROOT == out_dir or ROOT in out_dir.parents:
        raise RunError("the per-patient table must not be written inside "
                       "this repository")
    out_dir.mkdir(parents=True, exist_ok=True)
    table = out_dir / "venetoclax_clinical_pk_per_patient.csv"
    fields = ["patient", "visit", "dose_mg_per_day", "concentration_nm",
              "plasma_4b_hc_ng_per_ml", "plasma_cholesterol_mg_per_ml",
              "plasma_4b_hc_cholesterol_ratio"]
    lines = [",".join(fields)]
    for key in keys:
        record = {"patient": key[0], "visit": key[1], **rows[key]}
        lines.append(",".join(str(record[f]) for f in fields))
    table.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "run": "axis_c2_venetoclax_clinical_pk_contrast",
        "generated": "2026-11-10",
        "protocol": PRESPEC.name,
        "protocol_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
        "source": {
            "citation": "Lee J, et al. Clin Transl Sci 2024;17:e70106.",
            "pmcid": "PMC11645475", "licence": "cc by-nc",
            "channel_that_worked": (
                "main text via Europe PMC fullTextXML; supplement Table S3 via "
                "BioStudies S-EPMC11645475 redirecting to ftp.ebi.ac.uk. The "
                "Europe PMC supplementaryFiles endpoint returns only images "
                "for this article, which is how an earlier note wrongly "
                "concluded the per-patient biomarker was figure-only."),
        },
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in VCP.FROZEN.items()},
        "arms": {k: list(v) for k, v in VCP.ARMS.items()},
        "n_patients_primary": len(visit_one),
        "n_observations_all_visits": len(keys),
        "observations_excluded_by_the_declared_rule": excluded,
        "published_anchors_verified": anchors,
        "primary_visit_1_only": primary,
        "sensitivity_both_visits": sensitivity,
        "what_this_cannot_claim": VCP.what_this_cannot_claim(),
        "status_ceiling": VCP.FROZEN["status_ceiling"],
        "per_patient_table_is_outside_this_repository": str(table),
        "per_patient_table_sha256": hashlib.sha256(
            table.read_bytes()).hexdigest(),
        "no_patient_identifiers_values_or_predictions_exported": True,
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    print("published anchors:")
    for key, value in anchors.items():
        if key != "what_this_proves":
            print(f"  {key:<42} {value}")
    for label, block in (("PRIMARY (visit 1 only)", primary),
                         ("SENSITIVITY (both visits)", sensitivity)):
        stored = block["stored_error_bootstrap"]
        refit = block["refitting_patient_level_bootstrap"]
        null = block["permutation_control"]
        print(f"\n{label}: {block['arm_score']['n_patients']} patients, "
              f"{block['arm_score']['n_observations']} observations")
        print(f"  biomarker arm fold error "
              f"{block['arm_score']['geometric_fold_error']:.3f}  |  mean "
              f"{block['comparator_score']['geometric_fold_error']:.3f}")
        print(f"  stored-error d={stored['mean_difference_log10']:+.4f} "
              f"[{stored['interval'][0]:+.4f}, {stored['interval'][1]:+.4f}]")
        print(f"  REFITTING    d={refit['mean_difference_log10']:+.4f} "
              f"[{refit['interval'][0]:+.4f}, {refit['interval'][1]:+.4f}]  "
              f"excludes zero {refit['interval_excludes_zero']}")
        print(f"  null median {null['null_median']:+.4f}  rank p "
              f"{null['rank_p_value']:.4f}  ->  {block['verdict']}")
    print(f"\nper-patient table (outside the repo): {table}")
    print(f"artefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
