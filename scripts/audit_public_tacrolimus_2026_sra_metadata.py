"""Read only ENA/INSDC index metadata for tacrolimus 2026 SRA clinical cohort.

Original article reports 36 fecal metagenomes under BioProject
PRJNA1416570, and n36 clinical PK C0/D rows in a separate PLOS XLSX.
Sequence read metadata are NOT an identity crosswalk and cannot be used to
infer original clinical PK results by aligning row number or clinical strata.
"""
from __future__ import annotations

import csv
import io
import json
from urllib.parse import urlencode
from urllib.request import Request,urlopen

SRA_PROJECT="PRJNA1416570"
API="https://www.ebi.ac.uk/ena/portal/api/filereport"
LIMIT=200_000

def audit():
    params=urlencode({
        "accession":SRA_PROJECT,"result":"read_run",
        "fields":"run_accession,sample_accession,sample_alias,study_accession",
        "format":"tsv"
    })
    req=Request(API+"?"+params,headers={"User-Agent":"Omics2PK-qualified-original-2026-clinical-metagenomics-source/1.0"})
    with urlopen(req,timeout=45) as r: blob=r.read(LIMIT+1)
    if len(blob)>LIMIT:raise ValueError("SRA metadata query unexpectedly large")
    response=blob.decode("utf-8-sig")
    rows=list(csv.DictReader(io.StringIO(response),delimiter="\t"))
    needed={"run_accession","sample_accession","sample_alias","study_accession"}
    if not rows or not needed.issubset(rows[0]):
        raise ValueError("Original SRA crosswalk source unavailable or has changed")
    unique_samples={r["sample_accession"] for r in rows if r["sample_accession"]}
    unique_alias={r["sample_alias"] for r in rows if r["sample_alias"]}
    if not unique_samples:
        raise ValueError("No original human metagenomic sequence sample IDs")
    return {
        "status":"ORIGINAL_2026_SRA_SEQUENCING_SAMPLE_INDEX_NO_CLINICAL_CROSSWALK",
        "original_bioproject":SRA_PROJECT,
        "ENA_index_endpoint":API,
        "n_original_sequence_runs":len(rows),
        "n_original_distinct_sequence_specimens":len(unique_samples),
        "n_original_distinct_sequence_alias_strings":len(unique_alias),
        "n_original_clinical_patient_rows_in_PLOS_supplement":36,
        "clinical_spreadsheet_contain_SRA_sample_accession":False,
        "clinical_spreadsheet_contain_registered_shared_stool_PK_id":False,
        "same_patient_PK_and_omics_row_match_verified":False,
        "needs_authors_crosswalk_with_original_patient_sample_ID":True,
        "all_early_PK_exposure_timing_information_verified":False,
        "sample_aliases_and_patient_source_rows_never_exported":True,
    }

if __name__=="__main__":
    print(json.dumps(audit(),sort_keys=True,indent=2))
