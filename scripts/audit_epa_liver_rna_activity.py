#!/usr/bin/env python3
"""Bounded official U.S. EPA original liver mRNA/activity workbook availability audit.

DO NOT infer matched individual CYP transcript/activity from the paper's n=100
or from figure-only aggregates. Output sheet metadata only. No measurements,
identifiers, cells, or per-sample values are published by this audit.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
from urllib.error import HTTPError, URLError
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile

SOURCE="https://pasteur.epa.gov/uploads/10.23719/1519197/Liu%20et%20al.%20P450%20Expression%20in%20Human%20Liver%20ScienceHub.xlsx"
SOURCE_DOI="10.23719/1519197"
PAPER_DOI="10.1080/00498254.2020.1867929"
MS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RELS="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
RELPKG="{http://schemas.openxmlformats.org/package/2006/relationships}"
CAP=28_000_000


def _get(url=SOURCE):
    if url!=SOURCE:
        raise ValueError("Allowlisted original official U.S. EPA source URL only")
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-source-structure-audit/1.0"}),timeout=45) as r:
        blob=r.read(CAP+1)
    if len(blob)>CAP: raise ValueError("EPA original document exceeds safe size cap")
    if not blob.startswith(b"PK"): raise ValueError("EPA response not an XLSX archive")
    return blob


def workbook_structure(blob:bytes):
    result=[]
    with ZipFile(io.BytesIO(blob)) as z:
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        relationships=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rels={x.attrib["Id"]:x.attrib.get("Target","") for x in relationships.findall(RELPKG+"Relationship")}
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            r=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(n.itertext()) for n in r.findall(MS+"si")]
        for sh in wb.iter(MS+"sheet"):
            name=sh.attrib.get("name","")
            target=rels.get(sh.attrib.get(RELS+"id"),"")
            target=target.lstrip("/")
            if not target.startswith("xl/"):target="xl/"+target
            if target not in z.namelist():raise ValueError("Broken original sheet reference")
            root=ET.fromstring(z.read(target))
            rows=root.findall(".//"+MS+"sheetData/"+MS+"row")
            n_rows=len(rows)
            nonempty=0
            preview=[]
            for row in rows[:1]:
                cellvals=[]
                for c in row.findall(MS+"c"):
                    val=c.find(MS+"v")
                    if val is None:
                        inline=c.find(MS+"is")
                        display="".join(inline.itertext()) if inline is not None else ""
                    elif c.attrib.get("t")=="s":
                        try:display=shared[int(val.text)]
                        except (ValueError,IndexError,TypeError):
                            raise ValueError("Unresolvable XLSX shared string index") from None
                    elif c.attrib.get("t") in ("str","inlineStr"):
                        display=val.text or ""
                    else:
                        display=""
                    if display.strip():
                        cellvals.append(display.strip()[:70])
                preview.append(cellvals[:10])
            for row in rows:
                nonempty+=len(row.findall(MS+"c"))
            flattened=" ".join(" ".join(r) for r in preview)
            result.append({
                "sheet_name":name[:100],
                "xml_row_count":n_rows,
                "xml_cell_count":nonempty,
                "header_text_tokens_only":preview,
                "mentions_cyp_in_first_five_rows":bool(re.search("CYP",flattened,re.I)),
                "mentions_activity_in_first_five_rows":bool(re.search("activit|metaboli|microsom",flattened,re.I)),
                "donor_specific_transcript_activity_crosswalk_verified":False,
            })
    if not result:raise ValueError("No workbook sheets")
    return result


def audit():
    try:
        blob=_get()
    except (HTTPError,URLError,TimeoutError,OSError,ValueError) as exc:
        return {
            "status":"REMOTE_SOURCE_NOT_CONFIRMED",
            "official_source_doi":SOURCE_DOI,
            "article_doi":PAPER_DOI,
            "error_type":type(exc).__name__,
            "individual_donor_omics_to_activity_verified":False,
            "note":"Metadata-only endpoint probe cannot establish phenotype availability or donor-level matching",
        }
    return {
        "status":"XLSX_STRUCTURE_VERIFIED_NO_DONOR_LINK_CLAIM",
        "official_source_doi":SOURCE_DOI,
        "article_doi":PAPER_DOI,
        "bytes":len(blob),
        "sha256":hashlib.sha256(blob).hexdigest(),
        "source_sheet_structure":workbook_structure(blob),
        "individual_donor_omics_to_activity_verified":False,
        "no_original_values_or_subjects_exported":True,
    }

if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
