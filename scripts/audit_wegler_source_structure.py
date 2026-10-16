#!/usr/bin/env python3
"""Public 2021 Wegler supplement *structure* audit, not donor CLint extraction.

Downloads only original SHA-pinned PMC AWS source documents into transient
memory. Prints only document structure and aggregate named-table indicators.
No donor IDs, numerical clearance values, or proteome rows in outputs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import io
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

from scripts.audit_public_pk_sources import (
    SOURCES, EXPECTED_PMC_SUPPLEMENTS, bounded_http,
    official_pmc_media_https, pmc_cloud_metadata,
)

SHEET_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
OFFICE_REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PACKAGE_REL_NS = "{http://schemas.openxmlformats.org/package/2006/relationships}"
EXPECTED_HASHES = {
    "wegler_si_methods_pdf": "5beaab1117a32bb281734504dc19237bd264b5db0857837dfd794a9319892bab",
    "wegler_s1_proteomics": "8ce5df4f3e167915bbb963c8706b747360cf8ea9dd9538e7184d2a454d9c6ea6",
    "wegler_s2_bioinformatics": "623b6f4813193c348aa4a213ea9c4f5d380230a8932145cab472420059b3b213",
}


def pdf_structure(binary: bytes) -> dict:
    from pypdf import PdfReader  # only needed in the dedicated research CI
    doc = PdfReader(io.BytesIO(binary))
    term_counts = Counter()
    populated_pages = 0
    possible_tables = set()
    possible_figures = set()
    for page in doc.pages:
        # PDF native text only: NO OCR, NO raw assay values exported.
        content = page.extract_text() or ""
        if len(content.strip()) > 40:
            populated_pages += 1
        for term in ("clearance", "clint", "kinetics", "donor", "CYP3A4",
                     "midazolam", "bupropion", "omeprazole"):
            if term.casefold() in content.casefold():
                term_counts[term] += 1
        for match in re.finditer(r"\b(?:Table|TABLE)\s*S?([1-9]\d?)\b", content):
            possible_tables.add(match.group(1))
        for match in re.finditer(r"\b(?:Figure|FIGURE|Fig\.?)\s*S?([1-9]\d?)\b", content):
            possible_figures.add(match.group(1))
    return {
        "source": "supporting PDF native text, without extracting donor values",
        "pdf_pages": len(doc.pages),
        "pdf_pages_with_native_text": populated_pages,
        "pages_containing_research_terms": dict(sorted(term_counts.items())),
        "n_distinct_table_number_mentions": len(possible_tables),
        "n_distinct_figure_number_mentions": len(possible_figures),
        "donor_level_clint_table_verified": False,
        "verification_note": "Keyword and caption counts cannot establish donor-level numerical endpoint availability.",
    }


def workbook_structure_only(binary: bytes, chosen: tuple[str, ...]) -> dict:
    """List only sheet XML row/cell COUNTS, never individual cell values."""
    report = {}
    with ZipFile(io.BytesIO(binary)) as zf:
        wb = ET.fromstring(zf.read("xl/workbook.xml"))
        rel_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rels = {x.attrib.get("Id"): x.attrib.get("Target")
                for x in rel_root.findall(f"{PACKAGE_REL_NS}Relationship")}
        for sheet in wb.iter(f"{SHEET_NS}sheet"):
            name = sheet.attrib.get("name")
            if name not in chosen:
                continue
            rid = sheet.attrib.get(f"{OFFICE_REL_NS}id")
            target = rels.get(rid)
            if not target:
                raise ValueError("Missing worksheet relation")
            target = target.lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            if target not in zf.namelist():
                raise ValueError("Missing worksheet XML")
            row_count = 0
            nonempty_cells = 0
            first_row_counts = []
            # stream parsing: no workbook values or human-readable headers kept
            with zf.open(target) as fh:
                for event, el in ET.iterparse(fh, events=("end",)):
                    if el.tag == f"{SHEET_NS}row":
                        n = len(el.findall(f"{SHEET_NS}c"))
                        row_count += 1
                        nonempty_cells += n
                        if len(first_row_counts) < 4:
                            first_row_counts.append(n)
                        el.clear()
            report[name] = {
                "xml_rows": row_count,
                "xml_cells_with_elements": nonempty_cells,
                "first_four_row_cell_counts": first_row_counts,
                "note": "Cells may be formulae, headers, technical replicates or missing; not donor-level records.",
            }
    if not report:
        raise ValueError("No required compartments found")
    return report


def get_verified_supplements(fetch=bounded_http) -> dict:
    import hashlib
    cloud = json.loads(fetch(SOURCES["wegler_cloud_metadata"][0]).decode("utf-8"))
    meta = pmc_cloud_metadata(json.dumps(cloud).encode("utf-8"))
    if not meta["is_pmc_openaccess"]:
        raise ValueError("PMC source is not flagged open access")
    urls = cloud.get("media_urls") or []
    if not isinstance(urls, list):
        raise ValueError("Unexpected media listing")
    found = {}
    for name, (filename, _) in EXPECTED_PMC_SUPPLEMENTS.items():
        matches = [u for u in urls if isinstance(u, str)
                   and u.split("?")[0].endswith("/" + filename)]
        if len(matches) != 1:
            raise ValueError(f"Ambiguous/missing official media listing: {name}")
        blob = fetch(official_pmc_media_https(matches[0]))
        if hashlib.sha256(blob).hexdigest() != EXPECTED_HASHES[name]:
            raise ValueError(f"Primary source content changed: {name}")
        found[name] = blob
    return found


def audit_source_structures(fetch=bounded_http) -> dict:
    raw = get_verified_supplements(fetch=fetch)
    return {
        "purpose": "PUBLIC_IN_VITRO_SOURCE_STRUCTURE_ONLY",
        "population_pk_or_clinical_validation": False,
        "measured_donor_specific_clint_values_extracted": False,
        "pdf": pdf_structure(raw["wegler_si_methods_pdf"]),
        "protein_workbook": workbook_structure_only(
            raw["wegler_s1_proteomics"], ("Hepatocytes", "Microsomes", "Description")),
        "bioinformatics_workbook": workbook_structure_only(
            raw["wegler_s2_bioinformatics"],
            ("HigherInMic_Hep", "HigherInHep_Mic")),
        "source_sha256": EXPECTED_HASHES,
        "blocking_next_step": (
            "Independent manual confirmation of tabular donor-specific original CLint, "
            "assay units/status and donor-matched proteomics columns before regression."
        ),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    options = ap.parse_args()
    report = audit_source_structures()
    options.out.parent.mkdir(parents=True, exist_ok=True)
    options.out.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
