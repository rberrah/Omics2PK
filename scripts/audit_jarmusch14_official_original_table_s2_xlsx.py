"""Source-only audit of actual OFFICIAL PMC Jarmusch Table S2 original XLSX.

DOI 10.1111/cts.12785, PMC7485946 official supplementary file:
https://pmc.ncbi.nlm.nih.gov/articles/instance/7485946/bin/CTS-13-972-s002.xlsx

This reports ONLY public original attachment byte SHA and aggregate workbook
structure and whether candidate PK assay HEADER concepts are visible.
It never emits participant IDs, original numeric values, named original
sample rows, personal source metadata or attached file bytes into GitHub.
Source xlsx attachment can have scientific PK stats without individual
concentration-time labels; n~14 alone is never proof of a real same-donor
omics/PK join.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile,BadZipFile

URL="https://pmc.ncbi.nlm.nih.gov/articles/instance/7485946/bin/CTS-13-972-s002.xlsx"
S="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def check_official_jarmusch_table_s2():
    result={"source":"PMC7485946 original supplementary Table S2 (original publisher XLSX)",
            "doi":"10.1111/cts.12785",
            "url":URL,
            "original_patient_PK_omics_crosswalk_verified":False,
            "original_drug_plasma_Ct_AUC_per_person_retrieved":False,
            "original_human_identifiers_or_numeric_values_exported":False}
    try:
        with urlopen(Request(URL,headers={
            "User-Agent":"Omics2PK-genuine-source-research-supplement-structure/1.0"
        }),timeout=35) as r:
            original=r.read(2_000_001)
            result["HTTP_response_content_type"]=r.headers.get("Content-Type","").split(";")[0]
        if len(original)>2_000_000:
            result["status"]="ORIGINAL_SUPPLEMENT_TOO_LARGE_ABSTAIN"
            return result
    except HTTPError as err:
        result["status"]="ORIGINAL_SUPPLEMENT_HTTP_UNAVAILABLE_ABSTAIN"
        result["HTTP_status"]=err.code
        return result
    except (URLError,OSError,TimeoutError) as err:
        result["status"]="ORIGINAL_SUPPLEMENT_TRANSPORT_ABSTAIN"
        result["error_class"]=type(err).__name__
        return result
    result["retrieved_original_attachment_bytes"]=len(original)
    result["retrieved_original_attachment_SHA256"]=sha256(original).hexdigest()
    if not original.startswith(b"PK\x03\x04"):
        result["status"]="OFFICIAL_SUPPLEMENT_URL_RETURNED_NON_XLSX_ABSTAIN"
        return result
    sheets=[]
    try:
        with ZipFile(BytesIO(original)) as z:
            strings=[]
            if "xl/sharedStrings.xml" in z.namelist():
                parsed=ET.fromstring(z.read("xl/sharedStrings.xml"))
                for si in parsed.findall(S+"si"):
                    strings.append("".join(t.text or "" for t in si.iter(S+"t")))
            for name in sorted(z.namelist()):
                if not re.fullmatch(r"xl/worksheets/sheet\d+\.xml",name):continue
                element=ET.fromstring(z.read(name))
                rownodes=element.findall(".//"+S+"sheetData/"+S+"row")
                header_cells=[]
                widths=[]
                n_numeric=0
                n_rows_numeric=0
                for i,row in enumerate(rownodes):
                    cells=row.findall(S+"c")
                    widths.append(len(cells))
                    k=0
                    for cell in cells:
                        v=cell.find(S+"v")
                        if v is None:continue
                        txt=v.text or ""
                        if cell.get("t")=="s":
                            try:val=strings[int(txt)]
                            except (IndexError,ValueError):val=""
                        else:val=txt
                        if i<4 and len(header_cells)<200:
                            header_cells.append(val[:100])
                        try:
                            x=float(txt) if cell.get("t")!="s" else float("nan")
                            if x==x and abs(x)<float("inf"):
                                k+=1
                        except ValueError:pass
                    n_numeric+=k
                    n_rows_numeric+=k>0
                header=" ".join(header_cells).lower()
                flags={concept:bool(pattern.search(header)) for concept,pattern in {
                    "drug_analyte":re.compile(r"midazolam|omeprazole|caffeine"),
                    "PK_AUC":re.compile(r"\bauc\b|area under|a\.u\.c"),
                    "PK_clearance":re.compile(r"clearance|\bcl/f\b"),
                    "PK_Cmax":re.compile(r"cmax|maximum concentration"),
                    "PK_concentration_time":re.compile(r"concentration.*time|time.*concentration"),
                    "subject_key_column":re.compile(r"subject|participant|patient|volunteer|sample.?id"),
                    "untargeted_omic":re.compile(r"metabolomic|microbiome|16s|gnps"),
                }.items()}
                sheets.append({
                    "n_table_rows":len(rownodes),
                    "max_table_columns":max(widths,default=0),
                    "n_table_cells_numeric":n_numeric,
                    "n_rows_with_at_least_one_number":n_rows_numeric,
                    "candidate_biological_phenotype_header_role_flags":flags,
                    "individual_values_or_ids_not_exported":True})
    except (BadZipFile,ET.ParseError,KeyError,ValueError) as err:
        result["status"]="ORIGINAL_XLSX_UNPARSEABLE_ABSTAIN"
        result["error_class"]=type(err).__name__
        return result
    result["status"]="OFFICIAL_ORIGINAL_PK_S2_XLSX_WORKBOOK_STRUCTURE_CHECKED_NOT_PATIENT_LINKED"
    result["n_original_sheets"]=len(sheets)
    result["source_original_sheets_structural_metadata_only"]=sheets
    # NO code path here can mark real PK-to-omic source numerical data ready.
    return result


if __name__=="__main__":
    print(json.dumps(check_official_jarmusch_table_s2(),indent=2,sort_keys=True))
