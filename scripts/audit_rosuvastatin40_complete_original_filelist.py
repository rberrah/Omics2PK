"""Complete public original MetaboLights MTBLS3128 file inventory audit.

Official Study files endpoint used by the public editor's Files view.
The previous audit used include_raw_data=false (six original ISA files).
This one tests include_raw_data=true and explicitly groups ALL source file
types, including whether ORIGINAL processed omic matrices or individual
clinical PK tables are already deposited among raw mzML acquisitions.

No filenames containing biological sample IDs, individual labels, metabolic
intensities or concentration values are exported to logs/reports.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
import json
from pathlib import PurePosixPath
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

STUDY="MTBLS3128"
BASE=f"https://www.ebi.ac.uk/metabolights/ws/studies/{STUDY}"
MAX=30_000_000
CANDIDATE_EXT={"csv","tsv","xlsx","xls","zip","parquet","feather","h5","hdf5","rds"}
RAW_EXT={"mzml","raw","mzxml","cdf","wiff","d","mgf"}


def scan_original_public_list(payload:bytes):
    if len(payload)>MAX:raise ValueError("Oversize public accession manifest")
    doc=json.loads(payload)
    if not isinstance(doc,dict) or not isinstance(doc.get("study"),list):
        raise ValueError("Official study file inventory schema missing")
    entries=doc["study"]
    extensions=Counter()
    categories=Counter()
    processed_candidates=[]
    total_bytes=0
    for source in entries:
        if not isinstance(source,dict):continue
        name=source.get("file")
        if not isinstance(name,str):continue
        path=PurePosixPath(name)
        if ".." in path.parts:
            raise ValueError("Unsafe source path, do not fetch")
        ext=(path.suffix.lstrip(".")).casefold() or "no-extension"
        extensions[ext]+=1
        classification=str(source.get("type","unknown"))
        categories[classification[:50]]+=1
        size=source.get("size")
        if isinstance(size,int) and size>0:total_bytes+=size
        if ext in CANDIDATE_EXT:
            # Never log the source filename: may contain subject/sample IDs.
            lower=name.lower()
            processed_candidates.append({
                "extension":ext,
                "is_author_isa_file":path.name.startswith(("a_","s_","m_","i_")),
                "filename_has_processed_keyword":any(x in lower for x in
                  ("processed","peak","intens","feature","aligned","quant")),
                "filename_has_pk_keyword":any(x in lower for x in
                  ("auc","cmax","pharmacokinetic","clinical","pk_","concentration")),
                "filename_has_rosuvastatin_keyword":"rosuvastatin" in lower,
                "source_original_filename_redacted":True,
            })
    # Important: 1833 file count from third-party OmicsDI is not evidence
    # that this endpoint exposes same list today. Report exactly returned.
    return {
        "study":STUDY,
        "n_original_public_file_entries":len(entries),
        "n_original_filename_entries":sum(extensions.values()),
        "extensions":dict(sorted(extensions.items())),
        "source_kind_counts":dict(sorted(categories.items())),
        "n_raw_MS_files":sum(extensions.get(e,0) for e in RAW_EXT),
        "n_possible_processed_structured_files":len(processed_candidates),
        "possible_processed_structured_files":processed_candidates[:70],
        "source_entries_truncated_for_display":max(0,len(processed_candidates)-70),
        "total_reported_bytes_when_available":total_bytes,
        "source_sample_IDs_or_PK_values_exported":False,
    }


def run():
    out={"doi":"10.3389/fphar.2021.752960", "study":STUDY,
         "public_original_omic_and_target_PK_40_human_join_qualified":False}
    for name,url in (
        ("official_complete_file_inventory",BASE+"/files?include_raw_data=true"),
        ("official_metadata_only_file_inventory",BASE+"/files?include_raw_data=false"),
    ):
        try:
            with urlopen(Request(url,headers={"User-Agent":"Omics2PK-public-MTBLS3128-original-file-inventory/1.1"}),timeout=45) as conn:
                blob=conn.read(MAX+1)
            out[name]={"access":"SOURCE_AVAILABLE",
                       "original_manifest_sha256":sha256(blob).hexdigest(),
                       "structure":scan_original_public_list(blob)}
        except (HTTPError,URLError,TimeoutError,OSError,ValueError,TypeError) as err:
            out[name]={"access":"SOURCE_UNAVAILABLE","error_type":type(err).__name__}
    return out


if __name__=="__main__":
    print(json.dumps(run(),indent=2,ensure_ascii=False,sort_keys=True))
