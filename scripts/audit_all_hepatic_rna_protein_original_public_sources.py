"""Official-source, no-subject-leak qualification of SEVEN liver omics studies.

This script uses primary Europe PMC / NCBI/GEO endpoints and transiently
inspects only bounded released original supplement bytes. A literature claim
that original laboratory assayed n donors is *NOT* release of an aligned
patient-level n-row numeric table. Report source format/rows, source checksum,
and modality metadata flags ONLY, never per-donor values, sample labels,
original subject keys or measured source-column content.

A positive original-file structural screen is NOT donor crosswalk verification
or a qualification to run a predictive model. Only an independently checked
DONOR × RNA × PROTEIN × measured substrate turnover table can pass E2 fusion.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
from io import BytesIO
import gzip
import json
import re
from urllib.parse import quote
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from xml.etree import ElementTree as ET
from zipfile import ZipFile,BadZipFile

STUDIES={
    "Zhang_2016_10CYP_n100":{
        "doi":"10.1124/jpet.116.233635","published_n":100,
        "experiment":"human liver 10 CYP targeted mRNA+protein+probe substrate turnover",
    },
    "Tanner_2017_CYP2A6_n360":{
        "doi":"10.1124/jpet.116.237594","pmcid":"PMC5193072","published_n":360,
        "experiment":"CYP2A6 mRNA+protein, cotinine and coumarin metabolite formation",
    },
    "Bhatt_2018_UGT2B17_n455":{
        "doi":"10.1124/dmd.118.080952","pmcid":"PMC5938891","published_n":455,
        "published_n_RNA":230,"published_n_activity":346,
        "experiment":"UGT2B17 RNA+protein and androgen glucuronidation; different subsets",
    },
    "Ohtsuki_2012_CYP_n17":{
        "doi":"10.1124/dmd.111.042259","published_n":17,
        "experiment":"multiple CYP targeted RNA/protein and probe activities",
    },
    "Jamwal_2018_CYP3A4_n74":{
        "doi":"10.1021/acs.molpharmaceut.8b00159","pmcid":"PMC7065849",
        "published_n":74,
        "experiment":"CYP3A4 RNA/protein, midazolam 1-hydroxylation",
    },
    "Wenzel_2024_jejunum_liver_n11":{
        "doi":"10.1002/cpt.3055","published_n":11,
        "experiment":"same-person jejunum/liver targeted RNA+protein; NO measured drug clearance",
        "drug_turnover_assayed":False,
    },
}
GEO="GSE9588"
MAX_ARCHIVE=20_000_000
MAX_MEMBER=9_000_000
MAX_XML=2_000_000
W="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
S="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
# Structural header role checks only. Neither patient IDs nor actual headers
# from source will be printed. Unknown units/person pairing => ABSTAIN.
ROLES={
    "RNA":re.compile(r"\b(?:mRNA|RNA[- ]?seq|transcript|FPKM|qPCR)\b",re.I),
    "protein":re.compile(r"\b(?:protein|proteom|abundance|pmol\/mg)\b",re.I),
    "metabolic_activity":re.compile(
        r"\b(?:activity|metabolite|glucuronidation|hydroxylation|clearance|CLint|CYP)\b",re.I),
    "donor_id":re.compile(r"\b(?:donor|subject|liver|sample|patient)[\s_-]*(?:ID|no\.?|number|code)\b",re.I),
}


def get_bytes(url,cap=MAX_XML):
    try:
        with urlopen(Request(url,headers={"User-Agent":"Omics2PK source-validation research/1.0",
                                          "Accept":"application/octet-stream, application/json, */*"}),timeout=25) as r:
            body=r.read(cap+1)
            fmt=r.headers.get("Content-Type","").split(";")[0]
        return {"ok":len(body)<=cap,"bytes":body if len(body)<=cap else b"",
                "size":len(body),"type":fmt,"sha256":sha256(body).hexdigest()
                if len(body)<=cap else None}
    except (URLError,OSError,TimeoutError) as exc:
        return {"ok":False,"error_type":type(exc).__name__}


def pmc_metadata(doi,expected_pmcid=None):
    url=("https://www.ebi.ac.uk/europepmc/webservices/rest/search?"
         f"query={quote('DOI:'+doi)}&format=json&pageSize=3")
    r=get_bytes(url,500_000)
    if not r["ok"]:
        return {"doi_verified":False,"access_error":r.get("error_type","SIZE")}
    try:
        results=json.loads(r["bytes"])["resultList"]["result"]
    except (ValueError,KeyError):
        return {"doi_verified":False,"access_error":"NOT_JSON"}
    valid=[x for x in results if str(x.get("doi","")).casefold()==doi.casefold()]
    if not valid:
        return {"doi_verified":False,"access_error":"NO_MATCHING_DOI"}
    x=valid[0]
    pmcid=x.get("pmcid")
    conflict=bool(expected_pmcid and pmcid and expected_pmcid!=pmcid)
    return {"doi_verified":not conflict,"pmcid":pmcid or expected_pmcid,
            "pmcid_conflict":conflict,"is_open_access":x.get("isOpenAccess"),
            "has_fulltext":x.get("inEPMC"),"source_article_year":x.get("pubYear"),
            "original_study_pmid":x.get("id")}


def _row_metadata(rows):
    header=" | ".join(" ".join(x) for x in rows[:min(5,len(rows))])[:30000]
    return {
        "n_rows":len(rows),
        "max_cells_in_one_row":max((len(r) for r in rows),default=0),
        "n_source_data_rows_with_three_numeric_measurements":sum(
            sum(_numeric(v) for v in row)>=3 for row in rows[1:]),
        "n_source_header_roles":sum(bool(p.search(header)) for p in ROLES.values()),
        "source_header_role_flags":{name:bool(p.search(header))
                                    for name,p in ROLES.items()},
        "original_human_identifiers_or_measurements_emitted":False,
    }


def _numeric(value):
    try:
        v=float(value.replace(",","").strip())
        return v==v and abs(v)<float("inf")
    except (ValueError,TypeError,OverflowError):return False


def inspect_docx(payload):
    with ZipFile(BytesIO(payload)) as z:
        if "word/document.xml" not in z.namelist():
            return {"type":".docx","n_tables":0,"word_xml_missing":True}
        if z.getinfo("word/document.xml").file_size>MAX_MEMBER:
            return {"type":".docx","source_table_scan_abstain":"OVERSIZE"}
        root=ET.fromstring(z.read("word/document.xml"))
        summary=[]
        for t in root.iter(W+"tbl"):
            rows=[]
            for row in t.findall(W+"tr"):
                cells=[" ".join((c.text or "") for c in td.iter(W+"t"))
                       for td in row.findall(W+"tc")]
                rows.append(cells)
            summary.append(_row_metadata(rows))
        return {"type":".docx","n_tables":len(summary),
                "max_source_table_rows":max((d["n_rows"] for d in summary),default=0),
                "n_rough_original_paired_table_candidates":sum(
                    d["n_rows"]>=12 and d["n_source_header_roles"]==4
                    and d["n_source_data_rows_with_three_numeric_measurements"]>=10
                    for d in summary),
                "n_original_tables_with_50plus_numeric_rows":sum(
                    d["n_source_data_rows_with_three_numeric_measurements"]>=50 for d in summary),
                "original_table_values_emitted":False}


def inspect_xlsx(payload):
    rowsizes=[]
    with ZipFile(BytesIO(payload)) as z:
        infos=[i for i in z.infolist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",i.filename)]
        if any(i.file_size>MAX_MEMBER for i in infos):
            return {"type":".xlsx","source_table_scan_abstain":"OVERSIZE"}
        for info in infos:
            xml=ET.fromstring(z.read(info))
            rows=xml.findall(".//"+S+"row")
            rowsizes.append({"n_rows":len(rows),
                             "n_rows_with_3_or_more_numeric_cell_values":sum(
                                 sum(c.find(S+"v") is not None and
                                     _numeric(c.find(S+"v").text or "")
                                     for c in row.findall(S+"c"))>=3 for row in rows)})
    return {"type":".xlsx","n_worksheets":len(rowsizes),
            "max_worksheet_rows":max((x["n_rows"] for x in rowsizes),default=0),
            "max_worksheet_numeric_rows":max(
                (x["n_rows_with_3_or_more_numeric_cell_values"] for x in rowsizes),default=0),
            "original_table_values_emitted":False}


def inspect_package(body):
    """Official source only, never follow article HTML arbitrary external URLs."""
    if len(body)>MAX_ARCHIVE or not body.startswith(b"PK"):
        return {"archive_usable":False,
                "reason":"NOT_ZIP_OR_OFFICIAL_XML_ERROR_OR_TOO_LARGE"}
    files=[]
    extension_counts=Counter()
    try:
        with ZipFile(BytesIO(body)) as z:
            if "word/document.xml" in z.namelist():
                return {"archive_usable":True,"delivery":"direct_docx",
                        "documents":[inspect_docx(body)]}
            if "xl/workbook.xml" in z.namelist():
                return {"archive_usable":True,"delivery":"direct_xlsx",
                        "documents":[inspect_xlsx(body)]}
            for x in z.infolist():
                if x.is_dir():continue
                if ".." in x.filename.split("/"):
                    return {"archive_usable":False,"reason":"UNSAFE_SOURCE_MEMBER"}
                ext="."+x.filename.rsplit(".",1)[-1].lower() if "." in x.filename else "unknown"
                extension_counts[ext]+=1
                if ext not in (".docx",".xlsx",".csv",".tsv"):continue
                if x.file_size>MAX_MEMBER:
                    files.append({"type":ext,"source_table_scan_abstain":"OVERSIZE"})
                    continue
                payload=z.read(x)
                if ext==".docx":files.append(inspect_docx(payload))
                elif ext==".xlsx":files.append(inspect_xlsx(payload))
                else:
                    txt=payload.decode("utf-8-sig",errors="replace")
                    lines=txt.splitlines()
                    delim="\t" if ext==".tsv" else ","
                    sample=[x.split(delim) for x in lines[:min(len(lines),200)]]
                    head=" ".join(sample[0]) if sample else ""
                    files.append({"type":ext,"n_original_lines":len(lines),
                                 "n_column_header_roles":sum(bool(p.search(head)) for p in ROLES.values()),
                                 "original_subject_rows_emitted":False})
        return {"archive_usable":True,"delivery":"official_supplement_archive",
                "archive_member_extension_counts":dict(extension_counts),
                "documents":files}
    except (BadZipFile,ET.ParseError,ValueError) as e:
        return {"archive_usable":False,"reason":type(e).__name__}


def audit_study(name,study):
    metadata=pmc_metadata(study["doi"],study.get("pmcid"))
    result={
        "study":name,"doi":study["doi"],
        "published_total_human_donors":study["published_n"],
        "original_study_sample_count_not_proven_complete_cases":True,
        "original_donor_matched_RNA_protein_measured_metabolism_matrix_qualified":False,
        "original_RNA_vs_RNA_plus_protein_model_allowed":False,
        "why":"Source-original same-donor numerical rows and assay units NOT verified",
        "metadata":metadata
    }
    if "published_n_RNA" in study:result["reported_RNA_n"]=study["published_n_RNA"]
    if "published_n_activity" in study:result["reported_activity_n"]=study["published_n_activity"]
    if study.get("drug_turnover_assayed") is False:
        result["why"]="No actual metabolism or drug CL measurement in the source experiment"
    if not metadata.get("doi_verified") or not metadata.get("pmcid"):
        result["original_supplement"]={"status":"NO_DOI_VERIFIED_PMC_PUBLISHER_SUPPLEMENT_ENDPOINT"}
        return result
    resource=(f"https://www.ebi.ac.uk/europepmc/webservices/rest/"
              f"{metadata['pmcid']}/supplementaryFiles")
    src=get_bytes(resource,MAX_ARCHIVE)
    if not src["ok"]:
        result["original_supplement"]={"status":"NOT_RETRIEVABLE",
                                         "error_type":src.get("error_type","SIZE")}
        return result
    inspect=inspect_package(src["bytes"])
    result["original_supplement"]={
        "status":"SOURCE_ARCHIVE_INSPECTED" if inspect["archive_usable"] else
                 "OFFICIAL_ENDPOINT_NOT_NUMERICAL_ARCHIVE",
        "response_content_type":src["type"],
        "retrieved_body_size_bytes":src["size"],
        "download_SHA256":src["sha256"],
        **inspect,
    }
    # Do not silently equate DOCX table structural candidate with original
    # donor-level source join; this is solely an inspection step.
    return result


def audit_gse9588():
    url=f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO}&targ=self&form=text&view=quick"
    r=get_bytes(url,4_000_000)
    ret={"study":"Human_Liver_Cohort_GSE9588","accession":GEO,
         "GEO_RNA_source_released":False,"same_donor_CYP_activity_source_crosswalk_qualified":False,
         "same_donor_RNA_protein_activity_source_qualified":False}
    if not r["ok"]:
        ret["source_state"]="GEO_TEXT_UNAVAILABLE"
        ret["error_type"]=r.get("error_type","SIZE")
        return ret
    raw=r["bytes"].decode("utf-8",errors="replace")
    gsm=set(re.findall(r"\bGSM\d{6,9}\b",raw))
    ret["GEO_RNA_source_released"]=bool(gsm or "!Series_title" in raw)
    ret["original_GEO_METADATA_SHA256"]=r["sha256"]
    ret["n_GSM_accession_labels_in_GEO_summary"]=len(gsm)
    ret["GEO_summary_has_series_metadata"]="!Series_title" in raw
    ret["source_state"]="GEO_SERIES_METADATA_VERIFIED" if ret["GEO_RNA_source_released"] else "GEO_UNEXPECTED_TEXT"
    ret["note"]="Exact n427 RNA and n398 CYP activity donor ID join NOT obtained; GSM count is metadata coverage, NOT complete patients."
    return ret


def main():
    reports={name:audit_study(name,p) for name,p in STUDIES.items()}
    reports["GSE9588"]=audit_gse9588()
    return {
        "status":"ORIGINAL_SEVEN_LIVER_STUDY_SOURCE_LEVEL_AUDIT_NO_FALSE_MULTIOMIC_FIT",
        "source_studies_audited":len(STUDIES)+1,
        "checks":reports,
        "n_new_source_qualified_same_donor_RNA_protein_metabolic_function_matrices":sum(
            bool(d.get("original_donor_matched_RNA_protein_measured_metabolism_matrix_qualified"))
            for d in reports.values()),
        "no_patient_raw_rows_or_identifiers_exported":True,
    }


if __name__=="__main__":
    print(json.dumps(main(),sort_keys=True,indent=2))
