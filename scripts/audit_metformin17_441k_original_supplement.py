"""Tee et al. 2022 metformin original 441 KiB source supplement audit.

The public PMCID record reports one 441.4KiB ZIP, although the PMC OA metadata
API lists zero DOCX/XLSX media. This tests ONLY official source URLs and safely
reports member names and table dimensions, never subject rows or PK values.
Original 17 human donors × 15 calibrated concentrations (255 samples) and
pre-dose untargeted LCMS features are NOT a matched dataset until joined on
identical author-defined IDs with known units and collection timing.
"""
from __future__ import annotations
import hashlib
import io
import json
import re
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from zipfile import ZipFile,BadZipFile
from xml.etree import ElementTree as ET

PMCID="PMC9231303"
DOI="10.3390/pharmaceutics14061268"
SOURCE_URLS={
 "PMC canonical supplement":
   "https://pmc.ncbi.nlm.nih.gov/articles/PMC9231303/bin/pharmaceutics-14-01268-s001.zip",
 "PMC instance original supplement":
   "https://pmc.ncbi.nlm.nih.gov/articles/instance/9231303/bin/pharmaceutics-14-01268-s001.zip",
 "PMC Open Access S3 original media":
   "https://pmc-oa-opendata.s3.amazonaws.com/PMC9231303.1/media/pharmaceutics-14-01268-s001.zip",
 "Europe PMC official supplementaryFiles":
   "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9231303/supplementaryFiles",
 "MDPI article supporting information":
   "https://www.mdpi.com/article/10.3390/pharmaceutics14061268/s1",
}
MAX_SIZE=5_000_000
SS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DOC="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"

def get_original_source(url):
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-academic-original-supplement-audit/1.0",
       "Accept":"application/zip,application/octet-stream,*/*;q=.5"}),timeout=24) as rsp:
        body=rsp.read(MAX_SIZE+1)
        if len(body)>MAX_SIZE:
            raise ValueError("Original source exceeds declared small supplement limit")
        if body[:2]!=b"PK":
            raise ValueError("Original supplement not returned as a ZIP/OOXML binary")
        return body

def members_no_rows(data, *, depth=0):
    if not data.startswith(b"PK") or len(data)>MAX_SIZE or depth>2:
        raise ValueError("Not one allowed original supplement ZIP at approved nesting depth")
    out=[]
    with ZipFile(io.BytesIO(data)) as z:
        for info in z.infolist():
            if info.is_dir():
                continue
            if len(out)>=70:
                raise ValueError("Unexpectedly many ZIP members")
            path=info.filename
            if ".." in path.split("/") or path.startswith("/") or info.file_size>MAX_SIZE:
                raise ValueError("Unsafe original supplement archive")
            label=path.rsplit("/",1)[-1]
            # Never extract/print donor rows; only neutral file names and counts.
            summary={"filename":label[:150],"uncompressed_bytes":info.file_size}
            ext=label.lower()
            if ext.endswith(".xlsx"):
                buf=z.read(info)
                with ZipFile(io.BytesIO(buf)) as wb:
                    paths=wb.namelist()
                    if "xl/workbook.xml" in paths:
                        root=ET.fromstring(wb.read("xl/workbook.xml"))
                        summary["sheet_names"]=[s.attrib.get("name","")[:80]
                            for s in root.iter(SS+"sheet")]
                        summary["per_sheet_row_counts"]=[
                            sum(1 for _ in ET.fromstring(wb.read(p)).iter(SS+"row"))
                            for p in paths if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",p)]
            elif ext.endswith(".docx"):
                with ZipFile(io.BytesIO(z.read(info))) as doc:
                    if "word/document.xml" in doc.namelist():
                        root=ET.fromstring(doc.read("word/document.xml"))
                        summary["n_word_tables"]=sum(1 for _ in root.iter(DOC+"tbl"))
            elif ext.endswith(".csv"):
                # Count rows only. No values, raw feature names or IDs to logs.
                rows=z.read(info).count(b"\n")
                summary["csv_line_count_approx"]=rows
            elif ext.endswith(".zip"):
                # Europe PMC bundles article figures plus the actual MDPI
                # supplemental ZIP: audit that inner archive's member structure.
                summary["inner_supplement_members"]=members_no_rows(
                    z.read(info),depth=depth+1)
            out.append(summary)
    return out

def audit():
    result={"doi":DOI,"pmcid":PMCID,"n_original_human_subjects":17,
            "n_original_metformin_targeted_plasma_samples_reported":255,
            "article_lists_441k_supplementary_zip":True,
            "same_subject_omics_features_and_PK_join_verified":False,
            "qualified_for_individual_omic_to_PK_fit":False,
            "attempts":[]}
    for name,url in SOURCE_URLS.items():
        try:
            data=get_original_source(url)
            listing=members_no_rows(data)
            result["attempts"].append({
              "source":name,"status":"ORIGINAL_ZIP_AVAILABLE_METADATA_INSPECTED",
              "source_sha256":hashlib.sha256(data).hexdigest(),
              "n_archive_members":len(listing),"members":listing})
            # Once one legitimate source works, stop; do not use alternative
            # host or duplicate publisher content unnecessarily.
            break
        except (HTTPError,URLError,TimeoutError,ValueError,BadZipFile,
                ET.ParseError,OSError) as err:
            result["attempts"].append({"source":name,
              "status":"UNAVAILABLE_OR_NO_SOURCE_ZIP","error_type":type(err).__name__})
    # Even finding source numeric PK is only one half of an independent join:
    # values must be actively qualified in a separate private investigation.
    result["status"]="SUPPLEMENT_ARCHIVE_CHECK_ONLY_NO_E3_PREDICTIVE_CLAIM"
    return result

if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,indent=2,sort_keys=True))
