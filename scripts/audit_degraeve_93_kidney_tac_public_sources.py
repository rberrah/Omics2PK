"""Original 2024 93-human renal-transplant microbiome plus TAC AUC evidence audit.

DOI: 10.1002/cpt.3077; NCT04360031; BioProject PRJNA929904.
- 93 original CYP3A5 nonexpressers had original C0 and AUC_0_24
  ESTIMATED by a limited-sampling ML estimator, some missing C samples
  substituted by cohort median.
- AUC is not observed full NCA; should not be scored as such.
- 16S source at ENA is public, but clinical per-patient C0/AUC and
  PK↔16S source-ID crosswalk are not established as downloadable.
This original-source audit only verifies live ENA sample counts and source
archive availability metadata, never patient records or FASTQ.
"""
from __future__ import annotations
import csv
import io
import json
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from collections import Counter

PROJECT="PRJNA929904"
DOI="10.1002/cpt.3077"
EXPECTED_HUMAN_SUBJECTS=93
ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"
SUPPLEMENT_URLS=(
 "https://ascpt.onlinelibrary.wiley.com/doi/suppl/10.1002/cpt.3077/supinfo/cpt3077-sup-0001-Supinfo.zip",
 "https://onlinelibrary.wiley.com/doi/suppl/10.1002/cpt.3077/supinfo/cpt3077-sup-0001-Supinfo.zip",
)
CAP=400_000
ZIP_CAP=5_000_000


def ena_index():
    params=urlencode({
       "accession":PROJECT,"result":"read_run",
       "fields":"run_accession,sample_accession,sample_alias,study_accession",
       "format":"tsv"
    })
    req=Request(ENA+"?"+params,headers={"User-Agent":"Omics2PK-research-public-93-kidney-omics-PK-source-validation/1.0"})
    with urlopen(req,timeout=55) as response:
        payload=response.read(CAP+1)
    if len(payload)>CAP:raise ValueError("ENA source index too large")
    rows=list(csv.DictReader(io.StringIO(payload.decode("utf-8-sig")),delimiter="\t"))
    fields={"run_accession","sample_accession","sample_alias","study_accession"}
    if not rows or not fields.issubset(rows[0]):
        raise ValueError("Original ENA accession source lacks expected indexed fields")
    samples={r["sample_accession"] for r in rows if r["sample_accession"]}
    aliases={r["sample_alias"] for r in rows if r["sample_alias"]}
    studies={r["study_accession"] for r in rows if r["study_accession"]}
    if len(samples)!=EXPECTED_HUMAN_SUBJECTS:
        raise ValueError(f"Unexpected original 93-person sequence sample index; found {len(samples)}")
    return {
        "n_original_16s_runs":len(rows),
        "n_original_distinct_sequence_samples":len(samples),
        "n_original_distinct_16s_sample_aliases":len(aliases),
        "n_original_ENA_study_accessions":len(studies),
        "original_source_verified":True,
    }


def _archive_provenance():
    """Conservatively inspect official supplement without bypassing access controls.

    A 403 is not evidence of unavailable published data: only no automated
    access through this route. Do not forge an empty table, infer 93 PK values
    or redistribute any original clinical observations.
    """
    statuses=[]
    for url in SUPPLEMENT_URLS:
        try:
            req=Request(url,headers={"User-Agent":"Omics2PK-original-open-supplement-source-inventory/1.0"})
            with urlopen(req,timeout=18) as response:
                data=response.read(ZIP_CAP+1)
                mime=response.headers.get("Content-Type","unknown")
                status=response.status
            if len(data)>ZIP_CAP:
                statuses.append({"status":"oversized_or_source_changed","http_status":status})
            elif data[:2]==b"PK":
                import zipfile
                from hashlib import sha256
                with zipfile.ZipFile(io.BytesIO(data)) as z:
                    members=z.infolist()
                    safe=len(members)<=100 and sum(x.file_size for x in members)<20_000_000
                    names=[x.filename.rsplit("/",1)[-1] for x in members] if safe else []
                statuses.append({
                    "status":"ORIGINAL_ZIP_RETRIEVED_STRUCTURE_ONLY" if safe else "OVERSIZED_UNSAFE_ZIP",
                    "original_source_sha256":sha256(data).hexdigest(),
                    "n_members":len(members),"filenames_only":names,
                    "individual_clinical_AUC_x_16S_crosswalk_verified":False,
                })
            else:
                statuses.append({"status":"endpoint_did_not_return_original_zip",
                                 "http_status":status,"content_type":mime[:50]})
        except HTTPError as e:
            statuses.append({"status":"source_access_http_error","http_status":e.code})
        except (URLError,TimeoutError,OSError) as e:
            statuses.append({"status":"source_connection_or_access_blocked","error_type":type(e).__name__})
    return statuses


def run():
    try:
        source=ena_index()
        status="PUBLIC_93_HUMAN_TAC_16S_SEQUENCE_INDEX_NO_MATCHED_CLINICAL_PK_RECORDS_YET"
    except (HTTPError,URLError,TimeoutError,OSError) as exc:
        # An ENA 500 is transient unavailability, NOT missing experimental
        # data and NOT a successful current verification. Stay fail-closed.
        source={
            "original_source_verified":False,
            "status":"SOURCE_ENDPOINT_TEMPORARILY_UNAVAILABLE",
            "failure_type":type(exc).__name__,
            "n_original_distinct_sequence_samples":None,
        }
        status="ABSTAIN_TRANSIENT_ORIGINAL_ENA_SOURCE_UNAVAILABLE"
    return {
        "status":status,
        "original_article_doi":DOI,
        "original_bioproject":PROJECT,
        "original_clinical_trial":"NCT04360031",
        "original_pk_design":"cross_sectional_single_visit_C0_plus_limited_sampling_estimated_AUC0_24",
        "pk_AUC_original_status":"ALGORITHM_ESTIMATED_NOT_FULL_NONCOMPARTMENTAL_AUC",
        "reported_missing_PK_samples_source_handling":"COHORT_MEDIAN_IMPUTATION",
        "CYP3A5_genotype_inclusion":"ALL_NONEXPRESSERS",
        "n_original_human_recipients_reported":EXPECTED_HUMAN_SUBJECTS,
        "original_sample_source_audit":source,
        "publisher_supplement_access":_archive_provenance(),
        "original_source_subject_to_sequence_ID_join_verified":False,
        "unrestricted_individual_C0_x_estimated_AUC_and_molecular_matrix_verified":False,
        "zero_patient_identifiers_and_original_values_emitted":True,
        "not_a_clinical_longitudinal_PopPK_validation":True,
    }

if __name__=="__main__":
    print(json.dumps(run(),ensure_ascii=False,indent=2,sort_keys=True))
