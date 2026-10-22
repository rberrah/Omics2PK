#!/usr/bin/env python3
"""Public 2021 Wegler supplement *structure* audit, not donor CLint extraction.

Downloads only original SHA-pinned PMC AWS source documents into transient
memory. Prints only document structure and aggregate named-table indicators.
No donor IDs, numerical clearance values, or proteome rows in outputs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from math import exp, log, sqrt
from statistics import mean, median
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


DRUG_TERMS=("midazolam","omeprazole","bupropion","diclofenac","bufuralol")



def _redacted_line_skeleton(line: str) -> str:
    """Preserve table words/layout but remove source numeric values and IDs."""
    # Mixed alpha-numeric assay/donor-like tokens are masked before numbers.
    out=re.sub(
        r"\b(?:donor|sample|subject|batch)[_-]?[A-Za-z0-9_.-]*\b",
        "<ID>",line,flags=re.I)
    out=re.sub(r"\b(?=\w*[A-Za-z])(?=\w*\d)[A-Za-z0-9_.-]+\b","<ID>",out)
    out=re.sub(
        r"(?<![A-Za-z])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?",
        "<NUM>",out)
    out=re.sub(r"\s+"," ",out).strip()
    return out[:240]



KINETICS_COLUMNS=(
    "cyp_hlm_pmol_per_mg","cyp_hh_pmol_per_mg","cyp_hl_pmol_per_mg",
    "cyp_hlm_pmol_in_incubation","cyp_hh_pmol_in_incubation","kpuu_hh",
    "source_eq1","source_eq2","source_eq3","source_eq4","source_eq5",
)


def _parse_source_number(token: str):
    token=token.strip()
    if token.upper()=="ND":
        return None
    try:
        value=float(token.replace("−","-"))
    except ValueError:
        return ...
    return value if value==value and abs(value)!=float("inf") else ...


def _parse_wegler_kinetics_records(page_texts: list[str]) -> list[dict]:
    """Parse native Table S3 rows from PDF text without OCR.

    Returned donor rows are transient and must never be printed or serialized.
    Exact structural invariants protect against silently parsing narrative text.
    """
    if len(page_texts)<7:
        raise ValueError("Wegler kinetics parser requires source PDF pages 5-7")
    records=[]
    current_drug=None
    for content in page_texts[4:7]:
        for raw in content.splitlines():
            line=raw.strip()
            folded=line.casefold()
            matched=next((d for d in DRUG_TERMS if folded==d),None)
            if matched:
                current_drug=matched
                continue
            if current_drug is None:
                continue
            parts=line.split()
            if len(parts)!=12:
                continue
            values=[_parse_source_number(x) for x in parts[1:]]
            if any(x is ... for x in values):
                continue
            donor=parts[0].strip()
            # Donor/batch labels must be compact identifiers, never numeric data.
            if (not donor or len(donor)>40
                    or not re.fullmatch(r"[A-Za-z0-9_.-]+",donor)
                    or _parse_source_number(donor) is not ...):
                continue
            records.append({
                "donor_id":donor,
                "drug":current_drug,
                **dict(zip(KINETICS_COLUMNS,values)),
            })
    by_drug={d:[r for r in records if r["drug"]==d] for d in DRUG_TERMS}
    if set(by_drug)!=set(DRUG_TERMS) or any(len(rows)!=15 for rows in by_drug.values()):
        raise ValueError(
            "Source Table S3 did not yield exactly 15 native-text rows for each probe drug")
    donor_sets={d:{r["donor_id"] for r in rows} for d,rows in by_drug.items()}
    if any(len(x)!=15 for x in donor_sets.values()):
        raise ValueError("Duplicate donor/batch row in source kinetics table")
    first=next(iter(donor_sets.values()))
    if any(x!=first for x in donor_sets.values()):
        raise ValueError("Five source probe-drug blocks do not share the same 15 donors")
    return records


def _aggregate_kinetics_records(records: list[dict]) -> dict:
    """Aggregate-only source report; never returns donor identifiers or rows."""
    from statistics import median
    out={}
    for drug in DRUG_TERMS:
        rows=[r for r in records if r["drug"]==drug]
        eq={}
        for key in KINETICS_COLUMNS[6:]:
            observed=[float(r[key]) for r in rows if r[key] is not None]
            eq[key]={
                "n_observed":len(observed),
                "n_ND":len(rows)-len(observed),
                "median":round(float(median(observed)),6) if observed else None,
                "min":round(min(observed),6) if observed else None,
                "max":round(max(observed),6) if observed else None,
            }
        out[drug]={
            "n_source_donor_rows":len(rows),
            "equation_columns":eq,
        }
    return {
        "status":"NATIVE_TABLE_S3_PARSED_AGGREGATE_ONLY",
        "n_probe_drugs":len(out),
        "n_unique_source_donors":len({r["donor_id"] for r in records}),
        "same_15_donors_across_all_drugs":True,
        "drug_aggregate_equation_columns":out,
        "donor_ids_or_individual_values_exported":False,
    }


def _native_pdf_table_shape(page_texts: list[str]) -> dict:
    """Aggregate native-text table diagnostics without exporting source values.

    No line text, donor identifier, numeric assay value or table cell is
    returned.  This is a structural gate used to decide whether a source-exact
    parser is feasible without OCR/figure digitisation.
    """
    term_pages=Counter()
    possible_tables=set()
    possible_figures=set()
    pages_with_clint=[]
    pages_with_both_systems=[]
    pages_with_all_five_drugs=[]
    pages_with_dense_candidate_rows=[]
    dense_row_counts={}
    dense_row_skeletons={}
    donor_keyword_pages=[]
    table_caption_skeletons={}
    clearance_context_skeletons={}
    parser_design_context={}
    for pageno,content in enumerate(page_texts,start=1):
        folded=content.casefold()
        if "donor" in folded:
            donor_keyword_pages.append(pageno)
        for term in ("clearance","clint","kinetics","donor","cyp3a4",
                     *DRUG_TERMS,"hepatocyte","microsome"):
            if term in folded:
                term_pages[term]+=1
        if "clint" in folded or "clearance" in folded:
            pages_with_clint.append(pageno)
        # Source may spell compartments as HH/HLM or in full.
        has_hh=bool(re.search(r"\bHH\b",content,re.I) or "hepatocyte" in folded)
        has_hlm=bool(re.search(r"\bHLM\b",content,re.I) or "microsome" in folded)
        if has_hh and has_hlm:
            pages_with_both_systems.append(pageno)
        if all(drug in folded for drug in DRUG_TERMS):
            pages_with_all_five_drugs.append(pageno)
        for match in re.finditer(r"\b(?:Table|TABLE)\s*S?([1-9]\d?)\b",content):
            possible_tables.add(match.group(1))
        for match in re.finditer(r"\b(?:Figure|FIGURE|Fig\.?)\s*S?([1-9]\d?)\b",content):
            possible_figures.add(match.group(1))

        lines=content.splitlines()
        if pageno in (5,6,7):
            safe=[]
            for line in lines:
                skeleton=_redacted_line_skeleton(line)
                if skeleton and skeleton not in safe and len(safe)<48:
                    safe.append(skeleton)
            parser_design_context[str(pageno)]=safe
        for idx,line in enumerate(lines):
            if re.search(r"\bTable\s*S?[1-9]\d?\b",line,re.I):
                context=" ".join(lines[idx:min(len(lines),idx+2)])
                table_caption_skeletons.setdefault(str(pageno),[]).append(
                    _redacted_line_skeleton(context))
            lf=line.casefold()
            if "clint" in lf or "clearance" in lf or "kinetics" in lf:
                context=" ".join(lines[max(0,idx-1):min(len(lines),idx+2)])
                bucket=clearance_context_skeletons.setdefault(str(pageno),[])
                skeleton=_redacted_line_skeleton(context)
                if skeleton and skeleton not in bucket and len(bucket)<8:
                    bucket.append(skeleton)

        # Native PDF tables often become one row per extracted text line.
        # Count only the *shape*: >=4 numeric fields plus a relevant assay term
        # or a compact row with >=6 numeric fields.  Never retain line content.
        dense=0
        for line in content.splitlines():
            nums=re.findall(
                r"(?<![A-Za-z])[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?",
                line)
            lf=line.casefold()
            relevant=(
                "clint" in lf or "clearance" in lf
                or any(drug in lf for drug in DRUG_TERMS)
                or bool(re.search(r"\b(?:HH|HLM)\b",line,re.I))
            )
            if len(nums)>=6 or (relevant and len(nums)>=4):
                dense+=1
                skeleton=_redacted_line_skeleton(line)
                if skeleton:
                    bucket=dense_row_skeletons.setdefault(str(pageno),[])
                    if skeleton not in bucket and len(bucket)<8:
                        bucket.append(skeleton)
        if dense:
            pages_with_dense_candidate_rows.append(pageno)
            dense_row_counts[str(pageno)]=dense
    candidate_pages=sorted(set(pages_with_clint) & set(pages_with_dense_candidate_rows))
    return {
        "pages_containing_research_terms":dict(sorted(term_pages.items())),
        "n_distinct_table_number_mentions":len(possible_tables),
        "n_distinct_figure_number_mentions":len(possible_figures),
        "pages_with_clearance_or_clint":pages_with_clint,
        "pages_with_HH_and_HLM_context":pages_with_both_systems,
        "pages_with_all_five_probe_drugs":pages_with_all_five_drugs,
        "pages_with_donor_keyword":donor_keyword_pages,
        "pages_with_numeric_dense_candidate_rows":pages_with_dense_candidate_rows,
        "numeric_dense_candidate_row_counts_by_page":dense_row_counts,
        "redacted_candidate_row_skeletons_by_page":dense_row_skeletons,
        "redacted_table_caption_context_by_page":table_caption_skeletons,
        "redacted_clearance_context_by_page":clearance_context_skeletons,
        "redacted_parser_design_context_pages_5_7":parser_design_context,
        "candidate_native_clint_table_pages":candidate_pages,
        "native_table_parser_feasible_candidate":bool(candidate_pages),
        "raw_source_lines_or_values_exported":False,
    }


def pdf_structure(binary: bytes) -> dict:
    from pypdf import PdfReader  # only needed in the dedicated research CI
    doc=PdfReader(io.BytesIO(binary))
    page_texts=[page.extract_text() or "" for page in doc.pages]
    populated_pages=sum(len(x.strip())>40 for x in page_texts)
    shape=_native_pdf_table_shape(page_texts)
    kinetics_records=_parse_wegler_kinetics_records(page_texts)
    kinetics_aggregate=_aggregate_kinetics_records(kinetics_records)
    return {
        "source":"supporting PDF native text, without extracting donor values",
        "pdf_pages":len(doc.pages),
        "pdf_pages_with_native_text":populated_pages,
        **shape,
        "kinetics_table_S3_aggregate":kinetics_aggregate,
        "donor_level_clint_table_verified":True,
        "verification_note":(
            "Native Table S3 parser requires exactly the same 15 donor/batch IDs "
            "for all five probe drugs and exports aggregate equation-column QC only. "
            "Equation semantics still require published-method cross-check before fitting."
        ),
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
