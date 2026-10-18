"""Canonical ORIGINAL n33 patient PK versus predose metabolomic join test.

No row-order joins, no synthetic outcomes and NO source patient labels/values
in output. Official original Rotroff 2016 supplementary data only.

Source-specific original file layouts, as observed by independently executed
original-source schema audit:
- patient clinical 34-row worksheets: exactly 33 four-digit source person IDs
- metabolomic original sample worksheet: 33×A/B/C sample triplets, with
  sample label 4 letters + 4 digit same person ID + hyphen + A/B/C
Do not treat QC, blank, post-dose B/C, or group summaries as predictors.

This tool PROVES overlap of source person codes and classifies original
measured PK column headers. Fitting clinical PK requires confirming that
an outcome is ORIGINAL DRUG-exposure PK (not glucose AUC).
"""
from __future__ import annotations
from collections import Counter,defaultdict
from hashlib import sha256
from io import BytesIO
import json
import re
from urllib.request import Request,urlopen
from zipfile import ZipFile

from scripts.audit_rotroff33_metformin_original_pk_omics import (
    LIMIT,MEMBER_LIMIT,SUPPLEMENT,STUDY_DOI,_as_float
)
from scripts.audit_rotroff33_crosssheet_pk_omics import decode_workbook_sheets

SAMPLE=re.compile(r"^([A-Za-z]{4})([0-9]{4})-([ABC])$")
SOURCE_PERSON=re.compile(r"^[0-9]{4}$")
HEADER_WHITELIST=re.compile(r"^[A-Za-z0-9 \t_\-\/().%μµ]+$")
WORD_HINTS={
    "metformin":re.compile(r"metformin|(?:met)[ _-]*(?:concentration|auc|cmax)",re.I),
    "auc":re.compile(r"auc|area\s+under",re.I),
    "cmax":re.compile(r"c[_ ]?max|peak concentration|max(?:imum)?[ -]concentration",re.I),
    "glucose":re.compile(r"glucose|ogtt|glycemic|glycaemic",re.I),
    "pk_conc":re.compile(r"plasma|drug|serum|concentration|conc",re.I),
    "pk_cl":re.compile(r"clearance|^cl(?:r|sec|f)?$",re.I),
    "renal":re.compile(r"renal|urine|creatinine|gfr|crcl",re.I),
    "subject":re.compile(r"subject|id|patient|code",re.I),
    "time":re.compile(r"tmax|time|hr|hour",re.I),
    "age":re.compile(r"age",re.I),
    "sex":re.compile(r"sex|gender",re.I),
    "bmi":re.compile(r"bmi|weight",re.I),
}

def get_sheets(blob):
    if len(blob)>LIMIT or not blob.startswith(b"PK"):
        raise ValueError("Invalid official original publisher package")
    sheets=[]
    with ZipFile(BytesIO(blob)) as archive:
        for item in archive.infolist():
            if item.is_dir() or not item.filename.lower().endswith(".xlsx"):
                continue
            if item.file_size>MEMBER_LIMIT or ".." in item.filename.split("/"):
                raise ValueError("Unsafe original workbook")
            for rows in decode_workbook_sheets(archive.read(item)):
                sheets.append(rows)
    return sheets


def _clinical_ids(rows):
    if len(rows)!=34:return None
    result=[]
    for row in rows[1:]:
        text=row.get(0,"").strip()
        if not SOURCE_PERSON.fullmatch(text):return None
        result.append(text)
    if len(result)!=33 or len(set(result))!=33:
        raise ValueError("Ambiguous or duplicated clinical subject key")
    return set(result)


def _sample_index(rows):
    if not (98<=len(rows)<=145):return None
    groups=defaultdict(dict)
    for idx,row in enumerate(rows[1:],start=1):
        source=row.get(0,"").strip()
        m=SAMPLE.fullmatch(source)
        if not m:continue
        prefix,original_id,time=m.groups()
        if time in groups[original_id]:
            raise ValueError("Ambiguous repeated original omics person/timepoint")
        groups[original_id][time]=(idx,prefix)
    triples={k for k,v in groups.items() if set(v)=={"A","B","C"}}
    if len(triples)!=33 or len(groups)!=33:
        return {"study_subject_groups_original_match_count":len(groups),
                "n_complete_patient_triplets":len(triples),
                "source_molecular_subject_codes":None,
                "ambiguous":True}
    if len({prefix for v in groups.values() for _,prefix in v.values()})!=1:
        raise ValueError("Mixed original biological sample code prefixes")
    return {
        "study_subject_groups_original_match_count":len(groups),
        "n_complete_patient_triplets":len(triples),
        "source_molecular_subject_codes":triples,
        "ambiguous":False,
    }


def _safe_source_header_text(label):
    s=label.strip()
    if not 1<=len(s)<=58 or not HEADER_WHITELIST.fullmatch(s):
        return None
    # Names of source spreadsheet COLUMNS, never raw subject values.
    return s


def _clinical_column_structure(rows):
    roles={}
    original_headers=[]
    for k in sorted(rows[0]):
        header=rows[0][k]
        if _as_float(header) is not None:
            original_headers.append({"col_ordinal":k+1,"redacted":"untrusted_numeric_header"})
            continue
        raw=_safe_source_header_text(header)
        if raw is None:
            original_headers.append({"col_ordinal":k+1,"redacted":"non-structural_header"})
            continue
        n_measured=sum(_as_float(row.get(k)) is not None for row in rows[1:])
        kinds=[name for name,expr in WORD_HINTS.items() if expr.search(raw)]
        original_headers.append({"col_ordinal":k+1,"source_column_label":raw,
                                "semantic_categories":kinds,
                                "n_measured_original_donor_rows":n_measured})
    return original_headers


def qualify_person_key_overlap(sheets):
    with_omics=[]
    with_clinical=[]
    for idx,rows in enumerate(sheets,1):
        sample=_sample_index(rows)
        if sample:
            with_omics.append((idx,sample))
        ids=_clinical_ids(rows)
        if ids:
            with_clinical.append((idx,ids,rows))
    pairs=[]
    for i,sample in with_omics:
        if sample["source_molecular_subject_codes"] is None:continue
        for j,ids,rows in with_clinical:
            match=len(sample["source_molecular_subject_codes"]&ids)
            pairs.append({
                "original_omics_sheet_index":i,
                "original_clinical_sheet_index":j,
                "n_verified_numeric_patient_code_matches":match,
                "n_original_clinical_ids":len(ids),
                "n_complete_omics_ABC_triplet_ids":sample["n_complete_patient_triplets"],
                "one_to_one_original_33_person_code_coverage":match==len(ids)==33,
                "clinical_original_header_meaning":_clinical_column_structure(rows),
                "PK_endpoint_physical_units_still_require_source_codebook":True,
            })
    return {
        "n_candidate_molecular_sample_sheets":len(with_omics),
        "n_clinical_33_row_original_4digit_ID_tables":len(with_clinical),
        "n_unique_original_omics_ABC_triplets":max(
            (x["n_complete_patient_triplets"] for _,x in with_omics),default=0),
        "pairs":[*pairs],
        "exact_33_patient_source_omics_clinical_person_ID_pairing_qualified":any(
            p["one_to_one_original_33_person_code_coverage"] for p in pairs),
        "source_person_identifiers_or_molecular_numbers_exported":False,
        "not_a_measured_clinical_PK_model":True,
    }


def main():
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-validated-Rotroff33-original-patient-omics-clinical-join/1.0"}),timeout=65) as r:
        data=r.read(LIMIT+1)
    return {"DOI":STUDY_DOI,
            "official_original_supplement_archive_sha256":sha256(data).hexdigest(),
            **qualify_person_key_overlap(get_sheets(data))}

if __name__=="__main__":
    try:print(json.dumps(main(),sort_keys=True,indent=2))
    except (OSError,ValueError) as e:
        print(json.dumps({"status":"OFFICIAL_ORIGINAL_SOURCE_INCOMPLETE_ABSTAIN",
                          "error_type":type(e).__name__}))
