"""Original E3 clinical source audit: Jarmusch 2020 published targeted PK and
untargeted metabolomics on same individuals, but not necessarily shared public
machine-readable subject records.

ONLY official PMC OA S3 + Wiley publisher HTTPS source; never log subject
identifiers, raw concentrations, raw features or crosswalk pseudonyms.
All conclusions require numeric PK observations AND raw omic features keyed to
the SAME biologically identified subject. Source inventory != qualified join.
"""
from __future__ import annotations

from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from xml.etree import ElementTree as ET
from zipfile import ZipFile, BadZipFile

PMCID="PMC7485946"
DOI="10.1111/cts.12785"
OA_META="https://pmc-oa-opendata.s3.amazonaws.com/metadata/PMC7485946.1.json"
LIMIT=14_000_000
OOXML="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
SS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DOCX="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
PUBLISHER={
    "Supplementary methods DOCX":"https://ascpt.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fcts.12785&file=cts12785-sup-0001-Supinfo.docx",
    "Table S2 XLSX":"https://ascpt.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fcts.12785&file=cts12785-sup-0002-TableS2.xlsx",
}


def _read(url):
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK original-source qualification/1.0",
       "Accept":"application/json,application/octet-stream,*/*;q=.7"}),timeout=32) as conn:
        data=conn.read(LIMIT+1)
        if len(data)>LIMIT:
            raise ValueError("Source too large for restricted metadata-only audit")
        return data


def inspect_meta(data):
    record=json.loads(data)
    if record.get("pmcid")!=PMCID:
        raise ValueError("Incorrect scientific source PMCID")
    docs=[]
    for original in record.get("media_urls",[]):
        if not isinstance(original,str) or not original.startswith("s3://pmc-oa-opendata/"):
            continue
        filename=original.rsplit("/",1)[-1]
        if filename.lower().endswith((".xlsx",".docx")):
            docs.append((filename,original.replace("s3://pmc-oa-opendata/",
                       "https://pmc-oa-opendata.s3.amazonaws.com/",1)))
    return {"license_code":record.get("license_code"),
            "public_oa_version":record.get("version"),
            "n_official_original_office_supplements":len(docs)},docs


def inspect_xlsx_binary(data):
    if len(data)>LIMIT or not data.startswith(b"PK"):
        raise ValueError("Original supplement is not an Excel ZIP package")
    with ZipFile(BytesIO(data)) as z:
        names=z.namelist()
        if "xl/workbook.xml" not in names:
            raise ValueError("No actual spreadsheet inside source")
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=[n.attrib.get("name","") for n in wb.iter(SS+"sheet")]
        result=[]
        for fname in sorted(n for n in names if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",n)):
            root=ET.fromstring(z.read(fname))
            rows=list(root.iter(SS+"row"))
            header=[]
            if rows:
                for cell in rows[0].iter(SS+"c"):
                    value=cell.find(SS+"v")
                    inline=cell.find(SS+"is")
                    # header is not personally identifiable data, but do not leak
                    # any shared-string index as donor identifiers.
                    if inline is not None:
                        v=" ".join(inline.itertext())
                    elif value is not None and cell.attrib.get("t")!="s":
                        v=value.text or ""
                    else:
                        v="shared-string-header"
                    header.append(v[:55])
            n_cells=sum(1 for _ in root.iter(SS+"c"))
            n_num=sum(1 for c in root.iter(SS+"c")
                if c.attrib.get("t") not in ("s","str","inlineStr","b")
                and c.find(SS+"v") is not None)
            result.append({"sheet_xml":fname,"n_rows_including_header":len(rows),
                           "header_labels_only":header[:12],
                           "n_numeric_cells":n_num,
                           "total_cells":n_cells})
    return {"xlsx_sheet_names":sheets,
            "xlsx_sheets_metadata":result,
            "n_sheets":len(sheets),
            "subject_join_is_proven":False,
            "calibrated_drug_concentration_matrix_is_proven":False}


def inspect_docx_binary(data):
    if len(data)>LIMIT or not data.startswith(b"PK"):
        raise ValueError("Original supplement is not a Word OOXML package")
    with ZipFile(BytesIO(data)) as z:
        if "word/document.xml" not in z.namelist():
            raise ValueError("No Word document inside source")
        root=ET.fromstring(z.read("word/document.xml"))
    result=[]
    for table in root.iter(DOCX+"tbl"):
        rows=table.findall(DOCX+"tr")
        # first-row headings only (usually a table header), never participant data.
        header=[]
        if rows:
            for cell in rows[0].findall(DOCX+"tc"):
                val=" ".join(t.text or "" for t in cell.iter(DOCX+"t"))
                # Do not log even inadvertent original human values
                header.append(re.sub(r"\d{2,}","[digits]",val)[:45])
        result.append({"n_rows":len(rows),"n_columns_first_row":len(header),
                       "possible_column_names_only":header[:12]})
    return {"n_docx_tables":len(result),"docx_tables":result,
            "subject_join_is_proven":False}


def inspect_source(data,filename):
    return {
        "source_sha256":sha256(data).hexdigest(),
        "source_n_bytes":len(data),
        "filename":filename,
        "structural_evidence":(inspect_xlsx_binary(data) if filename.lower().endswith(".xlsx")
                              else inspect_docx_binary(data)),
    }


def audit():
    report={
        "primary_doi":DOI,"pmcid":PMCID,
        "design":"n14 humans: targeted LC-MS/MS plasma C(t), untargeted MS, fecal 16S and CYP genotype",
        "original_massive":"MSV000082493",
        "source_codes":"https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism",
        "same_human_protocol_evidence":True,
        "patient_level_targeted_PK_and_untargeted_MS_crosswalk_verified":False,
        "qualifies_for_clinical_E3_benchmark":False,
        "numeric_targeted_pk_results_extracted":False,
    }
    try:
        metadata,docs=inspect_meta(_read(OA_META))
        report["official_pmc_open_access_metadata"]=metadata
        report["official_oa_supplements"]=[]
        for filename,url in docs:
            try:
                report["official_oa_supplements"].append(
                    inspect_source(_read(url),filename))
            except (ValueError,HTTPError,URLError,BadZipFile,ET.ParseError) as err:
                report["official_oa_supplements"].append(
                    {"filename":filename,"status":"SOURCE_UNAVAILABLE",
                     "reason_type":type(err).__name__})
    except (ValueError,HTTPError,URLError,KeyError) as err:
        report["original_pmc_oa_status"]="SOURCE_UNAVAILABLE"
        report["original_pmc_oa_reason_type"]=type(err).__name__
    report["original_publisher_approach"]=[]
    for label,url in PUBLISHER.items():
        try:
            report["original_publisher_approach"].append(
                {"label":label,**inspect_source(_read(url),url.rsplit("file=",1)[-1])})
        except (ValueError,HTTPError,URLError,BadZipFile,ET.ParseError) as err:
            report["original_publisher_approach"].append(
                {"label":label,"status":"SOURCE_UNAVAILABLE","reason_type":type(err).__name__})
    return report


if __name__=="__main__":
    print(json.dumps(audit(),indent=2,sort_keys=True,ensure_ascii=False))
