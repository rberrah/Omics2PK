"""Detect public ORIGINAL Rosuvastatin MTBLS3128 MAF processed intensity matrix.

Strictly source-derived: original file manifest, sample IDs, assay/MAF
column match; never log volunteer/sample IDs, molecular feature values or
original drug PK values. MAF may be just metabolite annotation, not
subject×feature intensities: require actual matched quantified numeric cells.
"""
from __future__ import annotations
from collections import Counter
from io import StringIO
import csv
import json
from urllib.error import HTTPError,URLError
from scripts.audit_rosuvastatin40_isa_patient_crosswalk import (
    STUDY, API, _get, find_metadata_filenames
)


def parse_csv(raw):
    txt=raw.decode("utf-8-sig")
    dialect=csv.excel_tab if "\t" in txt[:1500] else csv.excel
    reader=csv.reader(StringIO(txt),dialect=dialect)
    rows=[]
    for row in reader:
        rows.append(row)
        if len(rows)>1500:raise ValueError("Unexpectedly large original annotation table")
    return rows


def analyze_matrix(original,known_names):
    rows=parse_csv(original)
    if not rows:raise ValueError("Empty original metabolite annotation")
    header=rows[0]
    if len(header)>2500:raise ValueError("Unexpected original feature dimension")
    exact=[j for j,v in enumerate(header) if v in known_names]
    # Full matching of column labels is stronger than ambiguous substring
    # similarity; a partial match must not be treated as complete matrix.
    count=0
    nonempty=0
    for row in rows[1:]:
        for i in exact:
            if i>=len(row):continue
            v=row[i].strip()
            if v:
                nonempty+=1
                try:
                    f=float(v)
                    if f==f and f not in (float("inf"),float("-inf")):count+=1
                except ValueError:pass
    return {
        "n_rows_excluding_header":len(rows)-1,
        "n_columns":len(header),
        "n_assay_sample_columns_matched_exactly":len(exact),
        "n_numeric_values_in_matched_columns":count,
        "n_nonempty_values_in_matched_columns":nonempty,
        "has_original_processed_quantitative_sample_feature_matrix":(
            len(exact)>=15 and len(rows)>5 and count>=len(exact)*3),
        "original_metabolite_values_or_IDs_exported":False
    }


def audit():
    out={"study":STUDY,"status":"ORIGINAL_PROCESSED_MAF_MATRIX_GATE",
        "actual_human_auc_cmax_values_joined":False,"no_sample_identifiers_logged":True}
    try:listing=find_metadata_filenames(_get(API+"/files?include_raw_data=false"))
    except (ValueError,HTTPError,URLError,OSError) as err:
        out["source_status"]="SOURCE_UNAVAILABLE"
        out["source_error"]=type(err).__name__
        return out
    sample_names=set()
    for file,kind in listing:
        if not file.startswith("s_"):continue
        try:
            rows=parse_csv(_get(API+"/download?file="+file))
        except (ValueError,HTTPError,URLError,OSError) as err:
            continue
        if not rows:continue
        h=rows[0]
        for r in rows[1:]:
            if "Sample Name" in h:
                ix=h.index("Sample Name")
                if ix<len(r) and r[ix].strip():sample_names.add(r[ix])
    results=[]
    for file,kind in listing:
        if not file.startswith("m_"):continue
        try:results.append(analyze_matrix(_get(API+"/download?file="+file),sample_names))
        except (ValueError,HTTPError,URLError,OSError):
            results.append({"status":"ORIGINAL_SOURCE_UNAVAILABLE"})
    out.update({"n_exact_public_source_sample_names":len(sample_names),
                "n_original_metabolite_assignment_tables":len(results),
                "per_original_MAF_source":results,
                "has_quantified_MAF_source":any(x.get("has_original_processed_quantitative_sample_feature_matrix") for x in results)})
    return out

if __name__=="__main__":
    print(json.dumps(audit(),ensure_ascii=False,sort_keys=True,indent=2))
