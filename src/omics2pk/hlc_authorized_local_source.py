"""Private source inspection for authorized Synapse HLC files.

syn4629, syn88644, syn89614 may require registered or controlled access.
Do not distribute the original donor records, phenotype values or tokens.
Successful inspection is NOT proof of original GSM/donor/assay pairing.
"""
from __future__ import annotations

import csv
import gzip
from hashlib import sha256
from pathlib import Path
import itertools
import re
from zipfile import ZipFile, BadZipFile

SOURCES={"phenotype":"syn4629","expression":"syn88644","genotype":"syn89614"}
KEYWORDS={
    "cyp":re.compile(r"cyp\d|cytochrome|p450",re.I),
    "activity":re.compile(r"activity|substrate|hydroxyl|metabol|clearance|clint",re.I),
    "donor_key":re.compile(r"sample|donor|subject|patient|gsm|id",re.I),
    "rna":re.compile(r"rna|transcript|expression|probe|gene",re.I),
    "protein":re.compile(r"protein|proteom|abundance|pmol",re.I),
}
MAX_BYTES=650_000_000
MAX_ROWS=1_500_000
MAX_COLUMNS=75_000


class LocalHLCSourceError(ValueError):
    pass


def _source_digest(path:Path):
    digest=sha256()
    n=0
    with path.open("rb") as handle:
        for chunk in iter(lambda:handle.read(262144),b""):
            n+=len(chunk)
            if n>MAX_BYTES:raise LocalHLCSourceError("Source exceeds safe audit cap")
            digest.update(chunk)
    if n==0:raise LocalHLCSourceError("Empty source file")
    return n,digest.hexdigest()


def _scan_text(path:Path, gz:bool):
    fn=gzip.open if gz else open
    with fn(path,"rt",encoding="utf-8-sig",errors="replace",newline="") as f:
        first=f.readline()
        if not first or chr(0) in first or len(first)>4_000_000:
            return {"readable_tabular":False,"reason":"BINARY_OR_UNSUPPORTED_TEXT"}
        delimiter="\t" if first.count("\t")>=first.count(",") else ","
        rows=csv.reader(itertools.chain([first],f),delimiter=delimiter)
        heads=next(rows)
        if not 1<len(heads)<=MAX_COLUMNS:
            return {"readable_tabular":False,"reason":"INVALID_COLUMN_COUNT"}
        detected={role:any(pattern.search(h[:200]) for h in heads)
                  for role,pattern in KEYWORDS.items()}
        n=0
        malformed=0
        for row in rows:
            n+=1
            malformed+=len(row)!=len(heads)
            if n>MAX_ROWS:raise LocalHLCSourceError("Too many source rows")
        return {
            "readable_tabular":True,
            "format":"tsv" if delimiter=="\t" else "csv",
            "n_rows":n,
            "n_columns":len(heads),
            "inconsistent_width_rows":malformed,
            "source_header_semantics_present":detected,
            "no_actual_values_or_headers_printed":True,
            "source_GSM_to_donor_mapping_verified":False,
        }


def inspect_original_local_source(path:str|Path,*,role:str)->dict:
    if role not in SOURCES:raise LocalHLCSourceError("Unknown HLC source role")
    p=Path(path)
    if not p.is_file() or p.is_symlink():
        raise LocalHLCSourceError("Not an accessible regular original file")
    size,digest=_source_digest(p)
    with p.open("rb") as handle:magic=handle.read(12)
    if magic.startswith(bytes.fromhex("1f8b")):
        kind="gzip"
    elif magic.startswith(b"PK"):
        kind="zip"
    elif magic.startswith((b"RDX",b"RDA",b"RDS",b"\x89HDF")) or b"\x00" in magic:
        kind="binary"
    else:
        kind="text_candidate"
    info={
        "original_Synapse_id":SOURCES[role],
        "original_role":role,
        "bytes":size,
        "sha256":digest,
        "format":kind,
        "local_source_only":True,
        "source_person_IDs_and_assay_values_exported":False,
        "original_numeric_CYP_activity_qualified":False,
        "same_donor_GEO_crosswalk_qualified":False,
        "CYP_model_fitted":False,
    }
    if kind in ("gzip","text_candidate"):
        try:info["structure"]=_scan_text(p,kind=="gzip")
        except (UnicodeError,EOFError,OSError,csv.Error) as err:
            info["structure"]={"readable_tabular":False,
                               "reason":type(err).__name__}
    elif kind=="zip":
        try:
            with ZipFile(p) as archive:
                names=[x.filename.lower() for x in archive.infolist() if not x.is_dir()]
            info["structure"]={
                "archive_members":len(names),
                "n_spreadsheet_members":sum(s.endswith((".csv",".tsv",".xlsx")) for s in names),
                "n_R_serialization_members":sum(s.endswith((".rds",".rdata")) for s in names),
                "no_raw_archive_member_names_exported":True,
            }
        except BadZipFile as err:
            raise LocalHLCSourceError("Invalid original ZIP") from err
    else:
        info["structure"]={
            "requires_authorized_source_format_conversion":True,
            "not_an_empty_biological_dataset":True,
        }
    return info


def audit_local_hlc_sources(paths:dict[str,str|Path])->dict:
    if not paths or any(k not in SOURCES for k in paths):
        raise LocalHLCSourceError("Original HLC roles invalid")
    files={role:inspect_original_local_source(path,role=role)
           for role,path in paths.items()}
    return {
        "status":"AUTHORIZED_HLC_ORIGINAL_SOURCE_STRUCTURE_ONLY",
        "study":"syn4499 Human Liver Cohort",
        "sources":files,
        "n_source_files":len(files),
        "original_CYP_activity_per_donor_and_unit_verified":False,
        "GSE9588_GSM_to_original_Synapse_donor_ID_matched":False,
        "training_or_model_PK_validation_executed":False,
        "original_identifiers_values_credentials_never_logged":True,
    }
