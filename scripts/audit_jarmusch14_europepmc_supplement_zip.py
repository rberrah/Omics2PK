"""Official Europe PMC exact-source Jarmusch 2020 original supplemental archive audit.

Unlike PMC OA S3 which listed 0 Office files, the Europe PMC supplementaryFiles
response can be a nested study ZIP including original publisher Office attachments.
Report exact source hashes, file names, dimensional structure, never observed
human PK, omic intensities, subject IDs, or model performance.
"""
from __future__ import annotations
from hashlib import sha256
from io import BytesIO
import json
from urllib.error import HTTPError,URLError
from urllib.request import Request,urlopen
from zipfile import ZipFile,BadZipFile
from xml.etree import ElementTree as ET

from scripts.audit_jarmusch14_original_omics_targeted_pk import (
    DOI,PMCID, inspect_source,
)

URL=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
MAX_BYTES=18_000_000
MAX_MEMBERS=120


def inspect_archive(blob:bytes,*,depth=0):
    if len(blob)>MAX_BYTES or not blob.startswith(b"PK") or depth>2:
        raise ValueError("Invalid or oversized official nested original archive")
    files=[]
    with ZipFile(BytesIO(blob)) as z:
        for f in z.infolist():
            if f.is_dir():
                continue
            if len(files)>=MAX_MEMBERS:
                raise ValueError("Unexpectedly large original archive")
            path=f.filename
            if ".." in path.split("/") or path.startswith("/") or f.file_size>MAX_BYTES:
                raise ValueError("Unsafe source ZIP member")
            name=path.rsplit("/",1)[-1]
            info={"file":name[:160],"uncompressed_bytes":f.file_size}
            if name.lower().endswith((".xlsx",".docx")):
                try:
                    source=inspect_source(z.read(f),name)
                    info["sha256"]=source["source_sha256"]
                    info["source_structure"]=source["structural_evidence"]
                except (BadZipFile,ValueError,ET.ParseError) as error:
                    info["status"]="UNREADABLE_ORIGINAL_OFFICE_DOCUMENT"
                    info["error_type"]=type(error).__name__
            elif name.lower().endswith(".zip"):
                info["nested_items"]=inspect_archive(z.read(f),depth=depth+1)
            files.append(info)
    return files


def audit():
    report={"doi":DOI,"pmcid":PMCID,
           "source_type":"Europe PMC original supplementaryFiles ZIP only",
           "human_targeted_PK_original_value_matrix_qualified":False,
           "human_omics_PK_exact_sample_crosswalk_qualified":False,
           "do_not_treat_archive_list_as_join":True}
    try:
        with urlopen(Request(URL,headers={"User-Agent":"Omics2PK-original-scientific-source-audit/1.1"}),timeout=35) as conn:
            data=conn.read(MAX_BYTES+1)
        if len(data)>MAX_BYTES:
            raise ValueError("Exceeds original supplement limit")
        report["source_zip_sha256"]=sha256(data).hexdigest()
        report["source_zip_bytes"]=len(data)
        report["structure"]=inspect_archive(data)
        report["access"]="OFFICIAL_ORIGINAL_ARCHIVE_RETRIEVED_STRUCTURE_ONLY"
    except (HTTPError,URLError,TimeoutError,OSError,ValueError,BadZipFile,ET.ParseError) as err:
        report["access"]="SOURCE_UNAVAILABLE_ABSTAIN"
        report["error_type"]=type(err).__name__
    return report


if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,sort_keys=True,indent=2))
