"""Original publicly completed Jarmusch 2019 GNPS task feature matrix access audit.

Article DOI 10.1111/cts.12785 explicitly links job
https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=b3966445abad4658a3cdc63c8198a666
and shows input Feature Quantification Table:
alan_jarmusch/DrugMetabolism_Tsunoda/
MS1basedMolNets_gapfilled_filtered/GNPS_quant.csv

Only endpoints historically DOCUMENTED by the GNPS API are queried, without
credentials. No brute force, no bypass of permissions. Real feature intensities
and original sample filenames NEVER logged, uploaded or written to GitHub.
Even if a network "bucket" table is downloadable it is not guaranteed to be
the same assay/export as the original 2019 input. A verified exact biological
sample ID crosswalk is required before a molecular vs PGx clinical PK fit.
"""
from __future__ import annotations
from collections import Counter
from hashlib import sha256
from io import BytesIO
import json
from urllib.parse import urlencode
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from zipfile import ZipFile,BadZipFile

TASK="b3966445abad4658a3cdc63c8198a666"
BASE="https://gnps.ucsd.edu/ProteoSAFe/"
CAP=50_000_000
QUERIES={
    "documented_original_GNPS_job_parameters":{
        "action":"ManageParameters",
        "params":{"task":TASK},
        "cap":500_000,
    },
    "historical_GNPS_output_reformatted_peak_table":{
        "action":"DownloadResultFile",
        "params":{"task":TASK,"block":"main",
                  "file":"quantification_table_reformatted/"},
        "cap":CAP,
    },
    "historical_GNPS_completed_cytoscape_zip":{
        "action":"DownloadResult",
        "params":{"task":TASK,"view":"download_cytoscape_data"},
        "cap":CAP,
    },
}


def inspect_public_task_asset(descriptor):
    original_url=BASE+descriptor["action"]+"?"+urlencode(descriptor["params"])
    try:
        with urlopen(Request(original_url,headers={
            "User-Agent":"Omics2PK-public-GNPS-2019-Jarmusch-source-matrix-QC/1.0"
        }),timeout=40) as response:
            content_type=response.headers.get("Content-Type","").split(";")[0]
            original=response.read(descriptor["cap"]+1)
    except HTTPError as err:
        return {"source_status":"PUBLIC_GNPS_HTTP_DENIED_OR_UNAVAILABLE",
                "http_status":err.code}
    except (OSError,URLError,TimeoutError) as err:
        return {"source_status":"PUBLIC_GNPS_TRANSPORT_FAILED",
                "error_type":type(err).__name__}
    if len(original)>descriptor["cap"]:
        return {"source_status":"PUBLIC_GNPS_FILE_EXCEEDS_LIMIT_NO_ASSAY_CONTENT_PARSED",
                "downloaded_bytes_minimum":len(original)}
    result={
        "http_content_type":content_type,
        "source_response_bytes":len(original),
        "sha256_of_exact_public_endpoint_response":sha256(original).hexdigest(),
        "original_source_individual_feature_values_exported":False,
        "verified_per_person_predose_PK_crosswalk":False,
    }
    if original.startswith(b"PK\x03\x04"):
        try:
            with ZipFile(BytesIO(original)) as z:
                infos=[i for i in z.infolist() if not i.is_dir()]
                types=Counter(i.filename.lower().split(".")[-1] for i in infos)
                result["archive_n_files"]=len(infos)
                result["archive_file_extensions_no_identifiers"]=dict(types)
                result["archive_largest_member_bytes"]=max((i.file_size for i in infos),default=0)
                result["archive_candidate_feature_matrices"]=sum(
                    i.filename.lower().endswith((".csv",".tsv",".txt"))
                    and "quant" in i.filename.lower() for i in infos)
            result["source_status"]="PUBLIC_GNPS_ZIP_DOWNLOADABLE_NEEDS_ORIGINAL_ASSAY_COLUMN_QC"
        except (BadZipFile,ValueError):
            result["source_status"]="PUBLIC_GNPS_BAD_ARCHIVE_ABSTAIN"
        return result
    if (original[:200].lstrip().lower().startswith((b"<!doctype html",b"<html")) or
        content_type=="text/html"):
        result["source_status"]="GNPS_HTML_INTERSTITIAL_NOT_A_MOLECULAR_QUANTITATIVE_TABLE"
        return result
    if b"<?xml" in original[:100] or content_type in ("text/xml","application/xml"):
        result["source_status"]="GNPS_SOURCE_JOB_XML_METADATA_ONLY_NO_OMICS_MATRIX"
        return result
    if b"," in original[:300] or b"\t" in original[:300]:
        sample=original[:200000].decode("utf-8-sig",errors="replace")
        lines=sample.splitlines()
        delimiter="\t" if lines[0].count("\t")>lines[0].count(",") else ","
        headers=lines[0].split(delimiter)
        result["header_field_count_only"]=len(headers)
        result["some_columns_have_peak_area_label"]=any(
            "peak area" in name.lower() for name in headers)
        result["source_status"]="PUBLIC_GNPS_DELIMITED_TABLE_RETRIEVED_NEEDS_REAL_MATRIX_ID_PAIRING"
        return result
    result["source_status"]="GNPS_UNKNOWN_RESPONSE_NO_ORIGINAL_MATRIX_QUALIFIED"
    return result


def audit_public_original_GNPS():
    results={name:inspect_public_task_asset(conf) for name,conf in QUERIES.items()}
    return {
        "study":"Jarmusch et al. 2020 true PK and original MS sample features",
        "source_original_GNPS_completed_task":TASK,
        "original_published_quant_filename":"GNPS_quant.csv",
        "original_published_massIVE_dataset":"MSV000082493",
        "original_public_assay_table_access":results,
        "source_original_same_volunteer_sample_PK_mapping_verified":False,
        "genuine_prior_omics_vs_genotype_trial_executed":False,
        "original_source_raw_MS_peak_values_or_subjects_exported":False,
    }


if __name__=="__main__":
    print(json.dumps(audit_public_original_GNPS(),indent=2,sort_keys=True))
