"""Source-only original Rotroff33 pseudonym layout diagnostic.

This does NOT infer sample pairing from row order. Shapes of original keys
are redacted (letters and digits replaced) and intersections are COUNTS only.
Distinct timepoint triplets and donor-level PK record counts are examined
inside original Europe PMC supplement, no sample ID or value logged.
"""
from __future__ import annotations
from collections import Counter
from io import BytesIO
import json
import re
from urllib.request import Request,urlopen
from zipfile import ZipFile

from scripts.audit_rotroff33_metformin_original_pk_omics import (
    SUPPLEMENT,LIMIT,MEMBER_LIMIT,_as_float
)
from scripts.audit_rotroff33_crosssheet_pk_omics import decode_workbook_sheets

def mask(s):
    return "".join("#" if c.isdigit() else "X" if c.isalpha() else c for c in s)[:45]

def _source_id_columns(rows):
    cols=[]
    for col in range(4):
        text=[r.get(col,"").strip() for r in rows[1:]]
        nonempty=[v for v in text if v]
        counter=Counter(nonempty)
        likely=[v for v in nonempty if 1<len(v)<=50 and
                re.fullmatch(r"[A-Za-z0-9_\- .]+",v)]
        cols.append({
            "column_position":col+1,
            "original_nonempty_values":len(nonempty),
            "n_unique_values":len(counter),
            "n_numeric_values":sum(_as_float(v) is not None for v in nonempty),
            "n_alphanumeric_short_values":len(likely),
            "three_most_common_anonymized_character_patterns":[
                {"token_shape":shape,"count":count}
                for shape,count in Counter(mask(v) for v in nonempty).most_common(3)
            ],
        })
    return cols

def _original_subject_suffix_keys(rows,col):
    # Exploratory candidate suffix A/B/C only; never proves original timepoint.
    by_parent={}
    for row in rows[1:]:
        t=row.get(col,"").strip()
        m=re.fullmatch(r"(.+?)[_\-. ]?([ABC])",t,re.I)
        if m and m.group(1) and not _as_float(m.group(1)):
            parent=m.group(1).casefold()
            by_parent.setdefault(parent,set()).add(m.group(2).upper())
    return {
        "n_source_distinct_provisional_ABC_parents":len(by_parent),
        "n_parents_with_all_three_A_B_C":sum(v=={"A","B","C"} for v in by_parent.values()),
        "n_parents_with_A":sum("A" in v for v in by_parent.values()),
        "source_subject_IDs_and_parent_keys_exported":False,
    }

def structure_from_rows(rows):
    cols=_source_id_columns(rows)
    return {
        "n_original_rows":len(rows),
        "first_four_columns_anonymized":cols,
        "suffix_candidates_first_four_columns":[
            {"column_position":i+1,**_original_subject_suffix_keys(rows,i)}
            for i in range(4)
        ],
        "never_pair_automatically":True,
    }

def run():
    with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-original-Rotroff33-privacy-safe-source-layout-v2/1.0"}),timeout=60) as stream:
        b=stream.read(LIMIT+1)
    if len(b)>LIMIT or not b.startswith(b"PK"):
        raise ValueError("Original official supplement unavailable")
    sheets=[]
    with ZipFile(BytesIO(b)) as z:
        for item in z.infolist():
            if item.is_dir() or not item.filename.lower().endswith(".xlsx"):
                continue
            if item.file_size>MEMBER_LIMIT:raise ValueError("Unsafe source workbook")
            for rows in decode_workbook_sheets(z.read(item)):
                # Source-wide layout only for original wide matrix, clinical n33.
                if len(rows) in range(30,43) or len(rows) in range(100,140):
                    sheets.append({"source_sheet_ordinal":len(sheets)+1,
                                   **structure_from_rows(rows)})
    return {"study":"Rotroff 2016 metformin 33 real humans original supplement",
            "n_target_shape_original_worksheets":len(sheets),
            "candidate_original_sample_label_layouts_only":sheets,
            "actual_same_human_pk_omics_ID_crosswalk_qualified":False,
            "individual_original_values_and_IDs_exported":False}


if __name__=="__main__":
    try:print(json.dumps(run(),indent=2,sort_keys=True))
    except (ValueError,OSError) as e:
        print(json.dumps({"status":"ORIGINAL_PUBLISHER_SOURCE_UNAVAILABLE_ABSTAIN","error_type":type(e).__name__,
            "actual_same_human_pk_omics_ID_crosswalk_qualified":False}))
