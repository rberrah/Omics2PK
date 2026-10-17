"""Read-only public MetaboLights MTBLS3128 FTP study data folder census.

Public MetaboLights ISA API enumerates only study ROOT; two directory entries
in include_raw_data=true do NOT mean only 8 files exist. Original deposited
mass-spectra files live in FTP /FILES/ subtree per official EMBL-EBI guide:
https://ebi-metabolights.github.io/guides/Files/

This script lists directory names anonymously, never downloads raw spectra,
never exports any participant/sample file names or human PK values.  No
model fit, no inferred crosswalk, no claim of 1833-file API completeness.
"""
from __future__ import annotations
from collections import Counter
from ftplib import FTP, error_perm, all_errors
import json
from pathlib import PurePosixPath
import re
import socket

HOST="ftp.ebi.ac.uk"
ROOT="/pub/databases/metabolights/studies/public/MTBLS3128/FILES"
CANDIDATE_EXT=frozenset({".csv",".tsv",".xlsx",".xls",".h5",".hdf5",".parquet",".rds",".zip"})
RAW_EXT=frozenset({".mzml",".mzxml",".raw",".wiff",".cdf",".mgf"})
MAX_DEPTH=4
MAX_ITEMS=20000
MAX_DIRS=500


def summarize_source_names(paths):
    exts=Counter()
    candidate=Counter()
    for relative in paths:
        p=PurePosixPath(relative)
        suffix=p.suffix.lower()
        exts[suffix if suffix else "no-extension"]+=1
        low=relative.lower()
        if suffix in CANDIDATE_EXT:
            candidate["structured_file"]+=1
            if any(w in low for w in ("processed","intens","aligned","peak","feature","quant")):
                candidate["structured_name_suggests_processed"]+=1
            if any(w in low for w in ("auc","cmax","concentration","pk_","clinical")):
                candidate["structured_name_suggests_PK_target"]+=1
    return {
        "original_listed_data_file_count":len(paths),
        "extension_counts":dict(sorted(exts.items())),
        "n_raw_MS_format_files":sum(exts[k] for k in RAW_EXT),
        "structured_candidate_counts":dict(sorted(candidate.items())),
        "individual_filename_or_subject_ID_exported":False,
        "raw_spectra_downloaded":False,
    }


def original_anonymous_ftp_census():
    ftp=FTP(timeout=18)
    n_dirs=0
    names=[]
    visited=set()
    try:
        ftp.connect(HOST,21,timeout=18)
        ftp.login()
        pending=[(ROOT,0)]
        while pending:
            directory,depth=pending.pop()
            if directory in visited:continue
            visited.add(directory)
            n_dirs+=1
            if n_dirs>MAX_DIRS:raise ValueError("Study exceeded original directory enumeration cap")
            ftp.cwd(directory)
            try:
                entries=list(ftp.mlsd())
            except (error_perm,NotImplementedError):
                # Support servers without MLSD; isdir inferred by attempting
                # CWD (only within the fixed original study subtree).
                entries=[]
                for n in ftp.nlst():
                    child=PurePosixPath(n).name
                    if child in (".","..") or not child:continue
                    try:
                        ftp.cwd(directory.rstrip("/")+"/"+child)
                        ftp.cwd(directory)
                        ty="dir"
                    except error_perm:
                        ty="file"
                    entries.append((child,{"type":ty}))
            for name,meta in entries:
                if name in (".","..") or "/" in name or "\\" in name:continue
                full=directory.rstrip("/")+"/"+name
                if not full.startswith(ROOT.rstrip("/")+"/"):continue
                ty=str(meta.get("type","")).lower()
                if ty=="dir":
                    if depth<MAX_DEPTH:pending.append((full,depth+1))
                elif ty in ("file",""):
                    names.append(full[len(ROOT)+1:])
                if len(names)>MAX_ITEMS:raise ValueError("Study exceeded bounded original file listing cap")
        return {"ftp_status":"ORIGINAL_DATA_FTP_CENSUS_AVAILABLE",
            "n_scanned_directories":n_dirs,
            "max_nested_directory_depth":MAX_DEPTH,
            "source":HOST+" [public /MTBLS3128/FILES; no subject paths exported]",
            "data_files":summarize_source_names(names)}
    finally:
        try:ftp.quit()
        except Exception:ftp.close()


def main():
    out={"study":"MTBLS3128","doi":"10.3389/fphar.2021.752960",
         "no_original_person_level_PK_omics_join_claim":True}
    try:
        out["official_public_FTP"]=original_anonymous_ftp_census()
    except (OSError,socket.timeout,EOFError,ValueError,error_perm) as err:
        out["official_public_FTP"]={"ftp_status":"PUBLIC_FTP_UNAVAILABLE_OR_INCOMPLETE",
                                    "error_type":type(err).__name__}
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=="__main__":main()
