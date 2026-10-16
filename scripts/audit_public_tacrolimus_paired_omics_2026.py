"""Original 2026 PLOS tacrolimus human PK×gut omics supplement inventory.

Only official DOI-scoped PLOS open supplementary spreadsheet. No patient rows
or measurements printed/persisted. A detected column label is NOT evidence
that a matching microbiome or metabolite matrix exists for those patients.
"""
from __future__ import annotations

import hashlib
import io
import json
import re
from urllib.request import Request, urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile

DOI="10.1371/journal.pone.0343817"
SUPPLEMENT="https://journals.plos.org/plosone/article/file?id="+DOI+".s001&type=supplementary"
MAX_BYTES=400_000
M="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
R="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
P="{http://schemas.openxmlformats.org/package/2006/relationships}"

def fetch_official():
    req=Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-Original-PLOS-Tac-2026-Source-QC/1.0"})
    with urlopen(req,timeout=40) as r: b=r.read(MAX_BYTES+1)
    if len(b)>MAX_BYTES or b[:2]!=b"PK":
        raise ValueError("Official PLOS supplement must be a small native XLSX")
    return b


def audit_spreadsheet(blob):
    sheets=[]
    with ZipFile(io.BytesIO(blob)) as zip:
        wb=ET.fromstring(zip.read("xl/workbook.xml"))
        rel=ET.fromstring(zip.read("xl/_rels/workbook.xml.rels"))
        targets={x.attrib["Id"]:x.attrib.get("Target","") for x in rel.findall(P+"Relationship")}
        strs=[]
        if "xl/sharedStrings.xml" in zip.namelist():
            roots=ET.fromstring(zip.read("xl/sharedStrings.xml"))
            strs=["".join(s.itertext()) for s in roots.findall(M+"si")]
        for sh in wb.findall(".//"+M+"sheet"):
            target=targets.get(sh.attrib[R+"id"],"").lstrip("/")
            if not target.startswith("xl/"):target="xl/"+target
            if target not in zip.namelist():raise ValueError("Original xlsx sheet target missing")
            doc=ET.fromstring(zip.read(target))
            rows=doc.findall(".//"+M+"sheetData/"+M+"row")
            text=[]
            numeric_cells=0
            for ix,row in enumerate(rows):
                for c in row.findall(M+"c"):
                    v=c.find(M+"v")
                    if v is None:
                        t=c.find(M+"is")
                        val="".join(t.itertext()) if t is not None else ""
                    elif c.attrib.get("t")=="s":
                        val=strs[int(v.text)]
                    elif c.attrib.get("t") in ("inlineStr","str"):
                        val=v.text or ""
                    else:
                        numeric_cells+=1
                        val=""
                    if ix<3 and val.strip():
                        text.append(val.strip()[:100])
            # Only FIRST row is the original XLSX scientific column header.
            header_text=[]
            if rows:
                for cell in rows[0].findall(M+"c"):
                    v=cell.find(M+"v")
                    if v is not None and cell.attrib.get("t")=="s":
                        header_text.append(strs[int(v.text)][:100])
                    elif v is not None and cell.attrib.get("t")=="str":
                        header_text.append((v.text or "")[:100])
                    else:
                        inline=cell.find(M+"is")
                        if inline is not None:header_text.append("".join(inline.itertext())[:100])
            # Only source scientific HEADER strings, never original IDs/values.
            permitted=[t for t in text if re.search(
                r"sample|group|patient|subject|person|tacro|concen|c/d|dose|sex|age|time|transplant|metabol|microbiome|clinical|conc",t,re.I)]
            sheets.append({
                "sheet_name":sh.attrib.get("name","")[:70],
                "xml_row_count":len(rows),
                "numeric_cell_count_without_original_values":numeric_cells,
                "only_first_original_row_header_labels":header_text[:25],
                "source_original_patient_row_count":max(0,len(rows)-1),
                "individual_subject_pk_omics_crosswalk_verified":False,
            })
    return sheets


def run():
    b=fetch_official()
    return {
        "status":"PLOS_ORIGINAL_SUPPLEMENT_FORMAT_AUDIT_NOT_YET_PAIRED_OMICS_PK",
        "article_doi":DOI,
        "official_source_kind":"CC_BY_PLOS_article_small_supplementary_xlsx",
        "size_bytes":len(b),
        "source_sha256":hashlib.sha256(b).hexdigest(),
        "worksheets":audit_spreadsheet(b),
        "article_reported_n_liver_recipients":36,
        "article_study_design":"human_cross_sectional_not_longitudinal",
        "original_SRA_metagenomic_bioproject":"PRJNA1416570",
        "does_not_prove_subject_join_to_raw_metabolomics":True,
        "does_not_prove_individual_concentration_time_PK":True,
        "patient_row_or_measurement_exported":False,
    }

if __name__=="__main__":
    print(json.dumps(run(),sort_keys=True,indent=2))
