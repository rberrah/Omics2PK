"""Single-fetch original Rotroff 2016 n33 clinical PK–omics public-source gate.

Only original Europe PMC supplement; read transiently ONCE, use one consistent
source object for all simultaneous schema and person-key checks, unlike
multiple dynamic publisher ZIP downloads. No patient rows, source sample IDs
or numerical observations printed or saved.
"""
from __future__ import annotations
from hashlib import sha256
import json
from urllib.request import Request,urlopen

from scripts.audit_rotroff33_metformin_original_pk_omics import (
    LIMIT,SUPPLEMENT,audit_package
)
from scripts.audit_rotroff33_crosssheet_pk_omics import analyze_package
from scripts.qualify_rotroff33_original_metformin_patient_keys import (
    get_sheets, qualify_person_key_overlap
)


def original_qualification(blob):
    a=audit_package(blob)
    b=analyze_package(blob)
    c=qualify_person_key_overlap(get_sheets(blob))
    matches=[p for p in c["pairs"] if p["one_to_one_original_33_person_code_coverage"]]
    return {
        "status":"ORIGINAL_PUBLIC_METFORMIN33_SOURCE_PERSON_KEY_AUDIT",
        "original_study_DOI":"10.3389/fphar.2016.00135",
        "original_public_aggregate_source_sha256":sha256(blob).hexdigest(),
        "n_original_xlsx":a["n_original_workbooks"],
        "n_original_worksheets":b["n_sheets_examined"],
        "n_original_omics_A_B_C_donor_triplets":c["n_unique_original_omics_ABC_triplets"],
        "n_original_clinical_donor_tables_with_33_unique_source_codes":
            c["n_clinical_33_row_original_4digit_ID_tables"],
        "n_exact_matched_original_33_patient_code_pairs":len(matches),
        "n_individual_numeric_pharmacokinetic_targets_verified":0,
        "real_original_paired_patient_numeric_metformin_drug_AUC_or_Cmax_confirmed":False,
        "source_patient_code_join":c["pairs"],
        "self_contained_33_predose_omics_drug_pk_table_found":
            a["n_source_workbooks_containing_selfcontained_33_patient_E3_table"]>0,
        "original_data_subject_key_values_or_PK_measurements_exported":False,
        "not_yet_an_E3_model_efficacy_result":True,
        "next_gate":"Require independent quantitative DRUG (not glucose) original Cmax/AUC source and sample A predose proof before any fit.",
    }


def main():
    with urlopen(Request(SUPPLEMENT,headers={
        "User-Agent":"Omics2PK-Rotroff33-one-download-original-crosswalk/1.0"
    }),timeout=65) as response:
        original=response.read(LIMIT+1)
    return original_qualification(original)


if __name__=="__main__":
    try:
        print(json.dumps(main(),sort_keys=True,indent=2))
    except (OSError,ValueError) as err:
        print(json.dumps({"status":"ORIGINAL_PUBLISHER_DATA_SOURCE_UNAVAILABLE_ABSTAIN",
                          "error_type":type(err).__name__,
                          "real_original_paired_patient_numeric_metformin_drug_AUC_or_Cmax_confirmed":False}))
