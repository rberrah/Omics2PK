"""Patient-ID overlap diagnostics for original Rotroff33 public XLSX supplements.

This explicitly checks *candidate* 33-person labels across individual original
sheets, but NEVER exports sample name, ID, original PK value, metabolomic
intensity, or hashed join keys. The output is only aggregate dimensions,
occurrence of outcome/clinical metadata and exact overlap counts.

Any inferred ID match (especially against suffix-stripped sample names) is
diagnostic, NOT proof of paired predose omics↔PK without original codebook.
"""
from __future__ import annotations
from collections import Counter
from io import BytesIO
from hashlib import sha256
import json
import re
from urllib.request import Request,urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile,BadZipFile

from scripts.audit_rotroff33_metformin_original_pk_omics import (
    LIMIT, MEMBER_LIMIT, NS, PMCID, STUDY_DOI, SUPPLEMENT,
    _cell_str, _col_index, _as_float
)

PK_TERMS={
    "AUC_token":re.compile(r"(?<![a-z])(?:auc|area\s+under\s+the\s+curve)(?![a-z])",re.I),
    "Cmax_token":re.compile(r"(?<![a-z])(?:c[\s_\-(]*max|maximum\s+(?:drug|plasma)?\s*concentration)(?![a-z])",re.I),
    "metformin_token":re.compile(r"\bmetformin\b",re.I),
    "glucose_token":re.compile(r"\bglucose\b",re.I),
    "predose_token":re.compile(r"(?:predose|pre-dose|baseline|prior to|timepoint\s*A|time point\s*A)",re.I),
    "subject_token":re.compile(r"(?:subject|participant|patient|volunteer|study\s*ID|id|sample)",re.I),
}
SKIP_ID=re.compile(r"^(?:male|female|yes|no|na|n/a|control|treatment|baseline|total|mean|median|sd|se|sem|c|a|b|\d+(?:\.\d+)?)$",re.I)


def decode_workbook_sheets(data):
    if len(data)>MEMBER_LIMIT or not data.startswith(b"PK"):
        raise ValueError("Original workbook fails size/type guard")
    out=[]
    with ZipFile(BytesIO(data)) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            if z.getinfo("xl/sharedStrings.xml").file_size>MEMBER_LIMIT:
                raise ValueError("Oversize source strings")
            shared=["".join(x.itertext()) for x in ET.fromstring(
                z.read("xl/sharedStrings.xml")).iter(NS+"si")]
        for name in sorted(x for x in z.namelist() if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x)):
            if z.getinfo(name).file_size>MEMBER_LIMIT:
                raise ValueError("Oversize source sheet")
            root=ET.fromstring(z.read(name))
            rows=[]
            for node in root.iter(NS+"row"):
                row={}
                for cell in node.iter(NS+"c"):
                    col=_col_index(cell.attrib.get("r",""))
                    if col is not None and col<1500:
                        row[col]=_cell_str(cell,shared).strip()
                rows.append(row)
            out.append(rows)
    return out


def id_candidates(rows, *, last_col=3):
    """Candidate pseudonym tokens, never returned to user/report."""
    pool=set()
    for row in rows:
        for k in range(last_col):
            raw=row.get(k,"").strip()
            if (2<=len(raw)<=50 and not SKIP_ID.fullmatch(raw)
                    and _as_float(raw) is None and not any(
                        p.search(raw) for p in PK_TERMS.values()
                    )):
                pool.add(raw.casefold())
    return pool


def summarize_sheet(rows):
    max_numeric_per_row=max(
        (sum(_as_float(x) is not None for x in row.values()) for row in rows),
        default=0
    )
    head=[x for row in rows[:min(10,len(rows))] for x in row.values()]
    all_strings=(s for row in rows for s in row.values() if len(s)<90)
    all_role={k:0 for k in PK_TERMS}
    for token in all_strings:
        for role,pat in PK_TERMS.items():
            if pat.search(token):all_role[role]+=1
    return {
        "n_rows":len(rows),
        "max_numeric_per_row":max_numeric_per_row,
        "n_original_short_id_tokens_first_3_columns":len(id_candidates(rows)),
        "clinical_pk_role_tokens_found_anywhere":{k:n for k,n in all_role.items()},
        "header_pk_role_tokens_found_first_10_rows":{
            k:sum(bool(p.search(s)) for s in head)
            for k,p in PK_TERMS.items()},
        "possible_33person_rows_and_10plus_quantified_covariates":(
            30<=len(rows)<=42 and max_numeric_per_row>=10
        ),
        "possible_omics_33x3_or_transposed_matrix":(
            (90<=len(rows)<=145 and max_numeric_per_row>=100)
            or (350<=len(rows)<=420 and max_numeric_per_row>=90)
        ),
        "source_identifiers_and_values_exported":False,
    }


def analyze_package(raw):
    if len(raw)>LIMIT or not raw.startswith(b"PK"):
        raise ValueError("Unsafe original official article supplement")
    sheets=[]
    with ZipFile(BytesIO(raw)) as z:
        for item in z.infolist():
            if item.is_dir() or not item.filename.lower().endswith(".xlsx"):
                continue
            if item.file_size>MEMBER_LIMIT or ".." in item.filename.split("/"):
                raise ValueError("Unsafe source workbook")
            for rows in decode_workbook_sheets(z.read(item)):
                sheets.append((len(sheets)+1,rows,summarize_sheet(rows)))
    summaries=[s for _,_,s in sheets]
    candidates=[]
    for ia,rows_a,asum in sheets:
        for ib,rows_b,bsum in sheets:
            if ia>=ib:continue
            if not ((asum["possible_omics_33x3_or_transposed_matrix"] and
                    bsum["possible_33person_rows_and_10plus_quantified_covariates"]) or
                    (bsum["possible_omics_33x3_or_transposed_matrix"] and
                    asum["possible_33person_rows_and_10plus_quantified_covariates"])):
                continue
            idsA=id_candidates(rows_a)
            idsB=id_candidates(rows_b)
            candidates.append({
                "source_sheet_position_A":ia,
                "source_sheet_position_B":ib,
                "exact_casefolded_non_numeric_first_3_column_token_intersection":len(idsA&idsB),
                "n_original_token_candidates_A":len(idsA),
                "n_original_token_candidates_B":len(idsB),
                "no_source_key_or_omics_data_exported":True,
                "join_automatically_certified":False,
            })
    return {
        "study_DOI":STUDY_DOI,"PMC":PMCID,
        "original_source_archive_SHA256":sha256(raw).hexdigest(),
        "n_sheets_examined":len(sheets),
        "original_sheets_summary":summaries,
        "possible_sheet_pairs_with_exact_candidate_ID_overlap":candidates,
        "original_individual_patient_AUC_metabolome_join_verified":False,
        "n_original_patient_or_omics_measurements_exported":0,
    }


def main():
    try:
        with urlopen(Request(SUPPLEMENT,headers={"User-Agent":"Omics2PK-Rotroff33-crosssheet-pseudonym-test/1.0"}),timeout=60) as r:
            contents=r.read(LIMIT+1)
        return analyze_package(contents)
    except (ValueError,OSError,BadZipFile,ET.ParseError) as e:
        return {"source_status":"ORIGINAL_SOURCE_UNAVAILABLE_ABSTAIN","error_type":type(e).__name__,
                "original_individual_patient_AUC_metabolome_join_verified":False}


if __name__=="__main__":
    print(json.dumps(main(),indent=2,sort_keys=True))
