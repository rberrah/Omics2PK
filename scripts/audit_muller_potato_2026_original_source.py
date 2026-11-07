"""Original source 2026 time-resolved n20 controlled potato meal CYP2D6 plasma markers.

Müller et al Front Pharmacol 2026 DOI 10.3389/fphar.2026.1886427
https://github.com/JulianPeterMueller/Solanidine-standardized-potato-meal
source pinned 2e36a063cde740df0a803311831977fa344a17b6.
No drug PK was reported; never conflate genotypic-phenotypic agreement
with independently measured drug CL/F or systemic AUC.
"""
from __future__ import annotations
import json
from hashlib import sha256
from io import BytesIO
from zipfile import ZipFile
from urllib.parse import quote
from urllib.request import Request,urlopen
from scripts.audit_ariana_original_source import (
    _read_shared,_get_workbook_sheets,_read_sheet
)

ORIGINAL_REPO_SHA="2e36a063cde740df0a803311831977fa344a17b6"
ORIGINAL_URL=(
 "https://raw.githubusercontent.com/JulianPeterMueller/Solanidine-standardized-potato-meal/"
 + ORIGINAL_REPO_SHA + "/"+quote("Solanidine standardized potato meal.xlsx")
)
def original_potato_source_schema():
    with urlopen(Request(ORIGINAL_URL,headers={
        "User-Agent":"Omics2PK-scientific-original-public-source-audit/1.0",
    }),timeout=55) as response:
        b=response.read(300_000)
    if not 10_000<len(b)<300_000 or b[:2]!=b"PK":
        raise ValueError("Original public meal trial source not XLSX")
    with ZipFile(BytesIO(b)) as z:
        sh=_get_workbook_sheets(z)
        shared=_read_shared(z)
        rows={}
        for name,path in sh:
            r=_read_sheet(z,path,shared)
            header=r[0] if r else {}
            rows[name]={
                "original_source_n_nonempty_rows":len(r),
                "source_header_first_row_columns_only":{
                    str(k):str(v)[:130] for k,v in sorted(header.items())},
                "n_worksheet_columns":len(header),
                "source_patient_data_cells_never_logged_or_exported":True
            }
    return {
        "status":"ORIGINAL_MULLER_2026_CONTROLLED_POTATO_CYP2D6_TIME_SERIES_SOURCE_SCHEMA",
        "original_depositor_git_sha":ORIGINAL_REPO_SHA,
        "source_XLSX_SHA256":sha256(b).hexdigest(),
        "study_reported_n_independent_humans":20,
        "original_n_sampling_timepoints_reported":9,
        "plasma_specimen":"controlled potato diet and time-resolved EDTA plasma LC-MS/MS",
        "matched_CYP2D6_genotype":True,
        "identified_solanidine_OH_solanidine_SSDA":True,
        "no_actual_medicine_PK_AUC_or_CL_dataset":True,
        "schema_only_no_new_effect_estimate":True,
        "source_sheet_schemas":rows
    }

if __name__=="__main__":
    print(json.dumps(original_potato_source_schema(),indent=2,sort_keys=True))
