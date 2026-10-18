"""Full-sheet structural search for genuine individual rosuvastatin AUC/Cmax.

Official article PMC8633954 Europe PMC original supplementary package. The
previous source qualification checked only first-row Excel labels; this audit
checks ALL worksheet cells without publishing any cell value, header, row,
sample ID or PK number. It does not infer patient joins from aggregate tables.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile, BadZipFile

PMCID = "PMC8633954"
DOI = "10.3389/fphar.2021.752960"
URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
MAX_PACKAGE = 12_000_000
MAX_MEMBER = 8_000_000
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
AUC = re.compile(r"(?<![a-z])(?:auc(?:0?[-– ]\d+)?|area\s+under\s+(?:the\s+)?curve)(?![a-z])", re.I)
CMAX = re.compile(r"(?<![a-z])c[\s_\-]*max(?![a-z])|maximum\s+(?:plasma\s+)?concentration", re.I)
PERSON = re.compile(r"\b(?:patient|subject|volunteer|participant|individual|code|id|number)\b", re.I)


def _value(cell, strings: list[str]) -> str:
    t = cell.get("t", "")
    v = cell.find(NS + "v")
    if t == "s":
        try:
            k = int(v.text) if v is not None and v.text else -1
            return strings[k] if 0 <= k < len(strings) else ""
        except (TypeError, ValueError):
            return ""
    if t == "inlineStr":
        return "".join(cell.itertext())
    return v.text if v is not None and v.text else ""


def _is_numeric(s: str) -> bool:
    try:
        x = float(s)
        return x == x and abs(x) < float("inf")
    except (TypeError, ValueError, OverflowError):
        return False


def inspect_original_xlsx(payload: bytes) -> dict:
    """Parse all original worksheet cells transiently, emitting metadata only."""
    if not payload.startswith(b"PK"):
        raise ValueError("Original supplement not a workbook")
    result = []
    with ZipFile(BytesIO(payload)) as z:
        names = z.namelist()
        strings = []
        if "xl/sharedStrings.xml" in names:
            if z.getinfo("xl/sharedStrings.xml").file_size > MAX_MEMBER:
                raise ValueError("Unsafe giant sharedStrings")
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            strings = ["".join(si.itertext()) for si in root.iter(NS + "si")]
        for name in sorted(x for x in names if re.fullmatch(r"xl/worksheets/sheet\d+\.xml", x)):
            if z.getinfo(name).file_size > MAX_MEMBER:
                raise ValueError("Unsafe large worksheet")
            root = ET.fromstring(z.read(name))
            rows = []
            for row in root.iter(NS + "row"):
                contents = []
                for cell in row.iter(NS + "c"):
                    ref = cell.get("r", "")
                    col = re.match(r"[A-Z]+", ref)
                    if not col:
                        continue
                    contents.append((col.group(0), _value(cell, strings)))
                rows.append(contents)
            n_auc = 0
            n_cmax = 0
            n_subject_labels = 0
            candidate_pairs = []
            all_cells = { (ri,c):s for ri,row in enumerate(rows) for c,s in row }
            for ri,row in enumerate(rows):
                auc_columns = [c for c,s in row if AUC.search(s)]
                cmax_columns = [c for c,s in row if CMAX.search(s)]
                n_auc += len(auc_columns)
                n_cmax += len(cmax_columns)
                n_subject_labels += sum(bool(PERSON.search(s)) for _,s in row)
                if auc_columns and cmax_columns:
                    for a in auc_columns:
                        for b in cmax_columns:
                            if a != b:
                                # Same header line with PK outcome label columns.
                                # Count only *numeric* cells below each header.
                                n_following_auc = sum(
                                    _is_numeric(all_cells.get((k,a),"")) for k in range(ri+1, len(rows))
                                )
                                n_following_cmax = sum(
                                    _is_numeric(all_cells.get((k,b),"")) for k in range(ri+1, len(rows))
                                )
                                candidate_pairs.append((n_following_auc, n_following_cmax))
            max_pair = max((min(a,b) for a,b in candidate_pairs), default=0)
            result.append({
                "n_rows": len(rows),
                "n_cells": sum(len(row) for row in rows),
                "n_auc_keyword_cells_anywhere": n_auc,
                "n_cmax_keyword_cells_anywhere": n_cmax,
                "has_possible_person_label_anywhere": bool(n_subject_labels),
                "n_rows_with_auc_cmax_as_different_headers": len(candidate_pairs),
                "max_both_PK_outcome_numeric_rows_below_same_header": max_pair,
                "candidate_40_original_subject_PK_outcome_table": max_pair >= 40,
                "individual_original_values_or_headers_exported": False,
            })
    return {"n_sheets": len(result), "worksheet_aggregate_structure": result}


def inspect_supplement_package(payload: bytes) -> dict:
    if len(payload) > MAX_PACKAGE or not payload.startswith(b"PK"):
        raise ValueError("Unsafe original Europe PMC package")
    reports = []
    with ZipFile(BytesIO(payload)) as outer:
        for item in outer.infolist():
            if item.is_dir() or not item.filename.lower().endswith(".xlsx"):
                continue
            if item.file_size > MAX_MEMBER or ".." in item.filename.split("/"):
                raise ValueError("Unsafe supplement workbook")
            data = outer.read(item)
            reports.append({
                "source_book_position": len(reports) + 1,
                "original_source_SHA256": sha256(data).hexdigest(),
                **inspect_original_xlsx(data),
            })
    return {
        "original_article": DOI,
        "original_supplement_SHA256": sha256(payload).hexdigest(),
        "n_original_excel_files": len(reports),
        "source_workbooks": reports,
        "no_patient_table_or_omics_PK_link_claim": True,
        "source_values_never_exported": True,
    }


def main():
    try:
        with urlopen(Request(URL, headers={"User-Agent": "Omics2PK-MTBLS3128-publisher-full-sheet-metadata-audit/1.0"}),timeout=60) as response:
            b = response.read(MAX_PACKAGE+1)
        print(json.dumps(inspect_supplement_package(b), indent=2, sort_keys=True))
    except (OSError, ValueError, BadZipFile, ET.ParseError) as err:
        print(json.dumps({"status": "ORIGINAL_SOURCE_UNAVAILABLE_ABSTAIN", "error_type": type(err).__name__}))


if __name__ == "__main__":
    main()
