"""Stage 2: qualify ORIGINAL MetaboLights MTBLS3128 ISA sample↔assay structure.

Uses official unauthenticated MetaboLights metadata file inventory and file
endpoints, not random web mirrors or screenshots. No volunteer IDs, sample
names, molecular intensities or individual AUC/Cmax leave transient memory.

Even identical ISA Sample Names do NOT prove that individual PK outcome values
are present: this audit inspects whether a subject-specific, quantitative PK
table and unique participant key are provided.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
import csv
from io import StringIO
import json
from urllib.parse import quote
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

STUDY="MTBLS3128"
API=f"https://www.ebi.ac.uk/metabolights/ws/studies/{STUDY}"
MAX=6_000_000


def _get(url):
    with urlopen(Request(url,headers={"User-Agent":"Omics2PK-public-ISA-study-qualification/1.0"}),timeout=35) as conn:
        data=conn.read(MAX+1)
        if len(data)>MAX:raise ValueError("Original metadata too large for safe individual-study audit")
        return data


def find_metadata_filenames(payload):
    info=json.loads(payload)
    if not isinstance(info,dict):raise ValueError("Invalid official metadata file inventory")
    files=info.get("study")
    if not isinstance(files,list):
        raise ValueError("No official study metadata listing")
    output=[]
    for item in files:
        if not isinstance(item,dict):continue
        f=item.get("file")
        if not isinstance(f,str):continue
        typ=item.get("type","unknown")
        if f.startswith(("s_", "a_", "m_", "i_")) and f.lower().endswith((".txt",".tsv",".csv")):
            if "/" in f or ".." in f:raise ValueError("Unexpected metadata filename path")
            output.append((f,str(typ)))
    return output


def summarize_isa(data,category):
    raw=data.decode("utf-8-sig")
    if not raw.strip():raise ValueError("Empty original ISA file")
    table=csv.DictReader(StringIO(raw),delimiter="\t")
    headers=table.fieldnames or []
    if len(headers)>250:raise ValueError("Unexpectedly wide original ISA metadata")
    hmap=[x.casefold().strip() for x in headers]
    rx=lambda marker:any(marker in key for key in hmap)
    unique_sample=set()
    unique_source=set()
    n=0
    has_pk_sample=False
    for row in table:
        n+=1
        if n>6000:raise ValueError("Unexpectedly many rows in metadata-only ISA")
        sid=row.get("Sample Name","").strip()
        if sid:unique_sample.add(sid)
        source=row.get("Source Name","").strip()
        if source:unique_source.add(source)
        if sid and ("auc" in sid.casefold() or "cmax" in sid.casefold()):
            has_pk_sample=True
    return {
       "kind":category,"n_original_metadata_rows":n,
       "n_unique_isa_sample_names":len(unique_sample),
       "n_unique_source_names":len(unique_source),
       "has_sample_name_column":"Sample Name" in headers,
       "has_subject_id_column":rx("subject") or rx("volunteer") or rx("participant"),
       "has_pk_auc_column":rx("auc"),
       "has_pk_cmax_column":rx("cmax") or rx("c max"),
       "has_drug_concentration_column":rx("concentration"),
       "has_normalized_peak_area_column":rx("peak area") or rx("intensity"),
       "has_predose_factor":rx("time") or rx("predose") or rx("pre-dose"),
       "unusual_sample_labels_found":has_pk_sample,
       "original_file_sha256":sha256(data).hexdigest(),
       "sample_values_exported":False,
    },unique_sample


def audit():
    out={
       "study":STUDY,
       "article":"10.3389/fphar.2021.752960",
       "n_original_volunteers_from_article":40,
       "status":"ISA_METADATA_ONLY_NOT_CLINICAL_PREDICTION",
       "source_listing_available":False,
       "n_exact_matched_target_auc_pk_records_verified":0,
       "subject_level_omics_pk_predictive_benchmark_ready":False,
    }
    try:
        listing=find_metadata_filenames(_get(API+"/files?include_raw_data=false"))
        out["source_listing_available"]=True
        out["n_metabolights_ISA_source_files"]=len(listing)
    except (ValueError,HTTPError,URLError,TimeoutError,OSError) as err:
        out["source_error_type"]=type(err).__name__
        return out
    inspected=[]
    samples=[]
    assays=[]
    for filename,kind in listing:
        if filename[0] not in "sami":continue
        success=False
        for url in (API+"/download?file="+quote(filename),API+"/"+quote(filename)):
            try:
                b=_get(url)
                if b"\x00" in b[:100]:raise ValueError("Not a plaintext original metadata file")
                result,ids=summarize_isa(b,kind)
                result["is_original_sample_table"]=filename.startswith("s_")
                result["is_original_assay_table"]=filename.startswith("a_")
                inspected.append(result)
                if filename.startswith("s_"):samples.append(ids)
                if filename.startswith("a_"):assays.append(ids)
                success=True
                break
            except (ValueError,HTTPError,URLError,TimeoutError,OSError,UnicodeError,csv.Error):
                continue
        if not success:
            inspected.append({"kind":kind,"source_unavailable":True,
                              "is_original_assay_table":filename.startswith("a_")})
    # Sample-level reuse in original author assays; no IDs emitted.
    if len(samples)==1:
        out["n_unique_ISA_source_samples"]=len(samples[0])
        out["assay_sample_overlap_counts"]=[len(samples[0]&x) for x in assays]
    out["source_file_results"]=inspected
    # A quantified AUC/Cmax result needs an actual original numeric PK file
    # and same biological donor crosswalk; data availability alone isn't proof.
    out["source_IDs_or_PK_values_exported"]=False
    return out

if __name__=="__main__":
    print(json.dumps(audit(),indent=2,ensure_ascii=False,sort_keys=True))
