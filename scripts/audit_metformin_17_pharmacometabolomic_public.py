"""Read-only qualification of Tee et al 2022 metformin n17 paired PK/metabolomics.

Original article doi:10.3390/pharmaceutics14061268, PMC9231303.
17 healthy humans × 15 plasma PK times = 255 assay samples.
Pre-dose + personalized Cmax *untargeted metabolomics* was measured.
Original supplement Table S5 claims participant PK profiles, but the raw
person×feature MS intensity matrix must NEVER be inferred from PCA, pathway
hits or manuscript presence alone.

Only official PMC OA metadata and small permitted supplement inventory.
No individual C(t) or MS features are exported to CI logs.
"""
from __future__ import annotations
import hashlib
import io
import json
from urllib.request import Request, urlopen
from urllib.error import HTTPError,URLError
from zipfile import ZipFile
from xml.etree import ElementTree as ET

PMCID="PMC9231303"
DOI="10.3390/pharmaceutics14061268"
SOURCE=f"https://pmc-oa-opendata.s3.amazonaws.com/metadata/{PMCID}.1.json"
MAX_METADATA=250_000
MAX_SUPP=6_000_000
DOCX_NS="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def source_metadata(data):
    if len(data)>MAX_METADATA:
        raise ValueError("Unexpected oversized official PMC article metadata")
    record=json.loads(data)
    if record.get("pmcid")!=PMCID:
        raise ValueError("Incorrect original metformin source PMID/PMCID")
    urls=record.get("media_urls") or []
    if not isinstance(urls,list):
        raise ValueError("Malformed public original supplementary listing")
    supp=[]
    for src in urls:
        if not isinstance(src,str) or not src.startswith("s3://pmc-oa-opendata/"):
            continue
        name=src.split("/")[-1]
        if name.lower().endswith((".docx",".xlsx",".csv",".zip",".pdf")):
            supp.append({"source_url":src.replace("s3://pmc-oa-opendata/",
                "https://pmc-oa-opendata.s3.amazonaws.com/",1),
                "filename":name})
    return {
        "pmcid":PMCID,
        "record_version":str(record.get("version","")),
        "license_code":str(record.get("license_code","")),
        "n_declared_supplementary_media":len(supp),
        "original_supplement_sources":supp,
    }


def supplement_structure(raw,name):
    if len(raw)>MAX_SUPP:
        raise ValueError("Unexpected oversized original supplementary object")
    result={"name":name,
            "original_source_sha256":hashlib.sha256(raw).hexdigest(),
            "original_bytes":len(raw),
            "has_person_level_omics_matrix_verified":False,
            "has_original_numeric_patient_PK_verified":False}
    if raw.startswith(b"PK") and name.lower().endswith(".docx"):
        with ZipFile(io.BytesIO(raw)) as z:
            if "word/document.xml" not in z.namelist():
                raise ValueError("DOCX missing original document")
            root=ET.fromstring(z.read("word/document.xml"))
        tables=[]
        for t in root.iter(DOCX_NS+"tbl"):
            rows=list(t.findall(DOCX_NS+"tr"))
            heads=[]
            if rows:
                for c in rows[0].findall(DOCX_NS+"tc"):
                    heads.append(" ".join("".join(c.itertext()).split())[:100])
            tables.append({"n_rows":len(rows),"n_columns_in_header":len(heads),
                           "head_labels_only":heads[:10]})
        result["docx_table_structure"]=tables
    elif raw.startswith(b"PK") and name.lower().endswith(".xlsx"):
        with ZipFile(io.BytesIO(raw)) as z:
            if "xl/workbook.xml" in z.namelist():
                root=ET.fromstring(z.read("xl/workbook.xml"))
                sheet_tag="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}sheet"
                result["xlsx_sheet_names_only"]=[e.attrib.get("name","") for e in root.iter(sheet_tag)]
    else:
        result["structure_status"]="NOT_A_MACHINE_READABLE_OMICS_MATRIX"
    return result


def _get(url,cap):
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-public-source-audit/1.0",
        "Accept":"application/json,application/zip,application/octet-stream"}),timeout=35) as response:
        payload=response.read(cap+1)
        if len(payload)>cap:
            raise ValueError("Source exceeds safety cap")
        return payload


def audit():
    report={
        "study":"Tee et al., 2022 metformin pharmacometabolomics, human n17",
        "doi":DOI,"pmcid":PMCID,
        "same_human_original_omics_and_PK_experiment":True,
        "original_PK_measurement_design":"17 humans, oral 1000 mg, 15 C(t) samples each, 255 samples",
        "molecular_design":"untargeted LC-MS pre-dose plasma and selected individual Tmax-related samples",
        "same_human_matched_predictive_matrix_PUBLIC_verified":False,
        "clinical_omics_to_PK_numeric_benchmark_qualified":False,
        "claim":"Source design verified, no donor-level molecular feature/dose/PK crosswalk yet",
    }
    try:
        original=_get(SOURCE,MAX_METADATA)
        source=source_metadata(original)
        report["source_metadata"]={
            k:v for k,v in source.items() if k!="original_supplement_sources"
        }
        reports=[]
        for item in source["original_supplement_sources"]:
            name=item["filename"]
            # Only Office structured original files; never log donor data.
            if not name.lower().endswith((".docx",".xlsx")):
                continue
            try:
                reports.append(supplement_structure(
                    _get(item["source_url"],MAX_SUPP),name))
            except (HTTPError,URLError,OSError,ValueError,ET.ParseError) as exc:
                reports.append({"name":name,"status":"ORIGINAL_SUPPLEMENT_UNAVAILABLE",
                                "reason_type":type(exc).__name__})
        report["official_supplement_inventory"]=reports
        report["access_status"]="OFFICIAL_PMC_METADATA_INSPECTED_NO_PAIRED_MATRIX_CLAIM"
    except (HTTPError,URLError,OSError,ValueError,KeyError) as exc:
        report["access_status"]="SOURCE_ABSTAIN_NO_ORIGINAL_OMICS_MATRIX"
        report["access_error_type"]=type(exc).__name__
    return report


if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,indent=2,sort_keys=True))
