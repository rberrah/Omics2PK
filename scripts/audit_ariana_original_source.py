#!/usr/bin/env python3
"""Source-only ORIGINAL public ARIANA CYP2D6 + solanidine identified metabolite XLSX audit.

Publication https://doi.org/10.3389/fphar.2026.1812007 original author-linked GitHub:
https://github.com/katjajust/Solandine-Metoprolol-Correlations
Pinned author repository immutable commit b6814967c0be7eae6f1843d0ec84e1987509bacc.
No patient-level cell contents, identifiers, doses, serum concentrations or
genotypes are printed or committed. Only original workbook sheets/column headings,
number of rows, kinds of numeric fields, and file SHA are reported.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import posixpath
import re
from urllib.request import urlopen,Request
from zipfile import ZipFile
import xml.etree.ElementTree as ET

ORIGINAL_GIT_COMMIT="b6814967c0be7eae6f1843d0ec84e1987509bacc"
ORIGINAL_SOURCE_URL=(
 "https://raw.githubusercontent.com/katjajust/Solandine-Metoprolol-Correlations/"
 +ORIGINAL_GIT_COMMIT+"/Solanidine_metoprolol_ratios_PGx.xlsx"
)
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "p":"http://schemas.openxmlformats.org/package/2006/relationships"}
MAX_BYTES=5_000_000

def _download():
    req=Request(ORIGINAL_SOURCE_URL,headers={
        "User-Agent":"Omics2PK-academic-source-audit/1.0"
    })
    with urlopen(req,timeout=80) as resp:
        b=resp.read(MAX_BYTES+1)
    if not 10_000<len(b)<MAX_BYTES or b[:2]!=b"PK":
        raise ValueError("Original author source is not a bounded XLSX ZIP")
    return b

def _read_shared(z):
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root=ET.fromstring(z.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in si.findall(".//m:t",NS))
            for si in root.findall("m:si",NS)]

def _get_workbook_sheets(z):
    wb=ET.fromstring(z.read("xl/workbook.xml"))
    rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    paths={}
    for x in rel.findall("p:Relationship",NS):
        rid=x.attrib["Id"]
        target=x.attrib["Target"].lstrip("/")
        if not target.startswith("xl/"):
            target=posixpath.normpath("xl/"+target)
        paths[rid]=target
    sheets=[]
    for sh in wb.findall("m:sheets/m:sheet",NS):
        rid=sh.attrib.get("{"+NS["r"]+"}id")
        if rid not in paths or paths[rid] not in z.namelist():
            raise ValueError("Author workbook reference missing")
        sheets.append((sh.attrib["name"],paths[rid]))
    return sheets

def _read_sheet(z,member,shared):
    root=ET.fromstring(z.read(member))
    records=[]
    for row in root.findall(".//m:sheetData/m:row",NS):
        cells={}
        for cell in row.findall("m:c",NS):
            address=cell.attrib.get("r","")
            key=re.match(r"[A-Z]+",address)
            if not key:continue
            letters=key.group()
            col=0
            for ch in letters:col=col*26+ord(ch)-64
            kind=cell.attrib.get("t","")
            v=cell.find("m:v",NS)
            if kind=="inlineStr":
                val="".join(t.text or "" for t in cell.findall(".//m:t",NS))
            elif v is None:
                val=""
            elif kind=="s":
                ix=int(v.text or "-1")
                if ix not in range(len(shared)):raise ValueError("Invalid source xlsx shared string pointer")
                val=shared[ix]
            else:
                val=v.text or ""
            cells[col-1]=val
        if cells:records.append(cells)
    return records

def audit_original():
    b=_download()
    with ZipFile(BytesIO(b)) as z:
        shared=_read_shared(z)
        sheets=_get_workbook_sheets(z)
        source={}
        for name,path in sheets:
            rows=_read_sheet(z,path,shared)
            if not rows:continue
            # Print only schema, never original patient attribute values.
            # Identify header by maximal text labels in first five rows.
            def text_count(rec):
                return sum(bool(re.search("[A-Za-z]",str(x))) for x in rec.values())
            headers=max(rows[:min(5,len(rows))],key=text_count)
            valid_headers={
                str(col):re.sub(r"\s+"," ",str(label))[:140]
                for col,label in sorted(headers.items())
            }
            source[name]={
                "n_rows_including_possible_header":len(rows),
                "n_columns_in_header":len(headers),
                "column_header_names_only":valid_headers,
                "original_source_individual_values_not_printed":True,
                "patient_rows_not_exported":True,
            }
    return {
        "status":"ORIGINAL_ARIANA_CYP2D6_IDENTIFIED_SOLANIDINE_PGX_SOURCE_SCHEMA",
        "original_author_repo_commit":ORIGINAL_GIT_COMMIT,
        "original_xlsx_SHA256":sha256(b).hexdigest(),
        "source_sheet_structure":source,
        "article_DOI":"10.3389/fphar.2026.1812007",
        "single_source_provenance":"original author deposited workbook",
        "validates_full_patient_PGx_metabolite_PK_crosswalk":False,
        "original_population_reported_N47":47,
        "source_genotype_is_CYP2D6_not_CYP3A":True,
        "clinically_measured_outcome_is_metoprolol_plasma_under_stable_dosing_not_AUC":True,
        "original_author_patient_rows_not_printed_or_committed":True,
    }

if __name__=="__main__":
    print(json.dumps(audit_original(),indent=2,sort_keys=True))
