"""Public Synapse HLC syn4499 access diagnostic, NO credentials and NO sample rows.

Historical Yang 2010 CYP phenotype donor crosswalk was once provided via
Sage Commons; HLC landing syn4499 lists several child resources. A Synapse
project landing record is NOT sufficient to certify unambiguous enzyme assay
activity × GSE9588 GSM original liver identity linkage. Probe only metadata
and count entity types; do NOT enumerate sample/data rows nor attempt private
authenticated files using user tokens from CI.

Original original CYP numerical activity table is required before training.
"""
from __future__ import annotations
from collections import Counter
import json
import re
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from urllib.parse import urlencode

SYNAPSE_ROOT="syn4499"
ROOT="https://repo-prod.prod.sagebase.org/repo/v1/entity/"
PUBLIC_SOURCE_IDS=("syn4499","syn3275753","syn88644","syn89614","syn4629")


def _json_get(uri):
    try:
        with urlopen(Request(uri,headers={
            "User-Agent":"Omics2PK-HLC-GSE9588-human-liver-CYP-source-metadata-audit/1.0",
            "Accept":"application/json"
        }),timeout=20) as r:
            data=r.read(350000)
            if len(data)>=350000:
                return {"status":"PUBLIC_METADATA_TOO_LARGE"}
        value=json.loads(data)
        if not isinstance(value,dict):return {"status":"NOT_JSON_OBJECT"}
        return {"status":"PUBLIC_JSON_RECEIVED","data":value}
    except HTTPError as e:
        return {"status":"HTTP_ACCESS_ERROR","http_status_code":e.code}
    except (OSError,URLError,ValueError,UnicodeError) as e:
        return {"status":"NOT_RETRIEVED","error_type":type(e).__name__}


def _probe_anonymous_file_metadata(synid,handle_id):
    """Never log presigned URLs or access tokens: metadata/anonymous gate only."""
    from urllib.parse import urlencode
    if not (re.fullmatch(r"syn\d+",str(synid)) and re.fullmatch(r"\d+",str(handle_id))):
        return "NO_SAFE_SOURCE_FILE_HANDLE"
    endpoint=("https://repo-prod.prod.sagebase.org/file/v1/file/"+str(handle_id)+"?"+
              urlencode({"redirect":"false","fileAssociateId":synid,
                         "fileAssociateType":"FileEntity"}))
    try:
        with urlopen(Request(endpoint,headers={
            "User-Agent":"Omics2PK-original-HLC-public-metadata-only/1.0"
        }),timeout=16) as response:
            # Do NOT print the response: may be a temporary presigned URL.
            response.read(128)
            return "ANONYMOUS_DOWNLOAD_URL_ENDPOINT_RETURNED_200_NEEDS_BINARY_QC"
    except HTTPError as e:
        return "ANONYMOUS_DOWNLOAD_HTTP_"+str(e.code)
    except (OSError,URLError):
        return "ANONYMOUS_DOWNLOAD_REQUEST_FAILED"


def audit_public_hlc_synapse():
    result={}
    for synid in PUBLIC_SOURCE_IDS:
        response=_json_get(ROOT+synid)
        row={"metadata_retrieval_status":response["status"]}
        if "http_status_code" in response:row["http_status_code"]=response["http_status_code"]
        if response["status"]=="PUBLIC_JSON_RECEIVED":
            source=response["data"]
            label=str(source.get("name",""))
            if (len(label)<120 and re.fullmatch(r"[A-Za-z0-9 _().-]+",label)
                    and not re.search(r"GSM\d+|PERSON_|DONOR_\d+",label,re.I)):
                row["public_original_asset_name"]=label
            row.update({
                "original_record_id_matches_requested_id":source.get("id")==synid,
                "entity_class":str(source.get("concreteType","")).split(".")[-1],
                "original_file_entity_has_data_file_handle":bool(source.get("dataFileHandleId")),
                "source_hlc_parent_folder_linked":source.get("parentId")==SYNAPSE_ROOT,
                "public_content_data_download_authorized_and_executed":False,
                "one_to_one_original_427_GSM_to_398_assay_ids_verified":False,
            })
            if source.get("dataFileHandleId"):
                row["anonymous_file_download_api_status"]=_probe_anonymous_file_metadata(
                    synid,source["dataFileHandleId"])
        result[synid]=row
    children=_json_get(ROOT+SYNAPSE_ROOT+"/children?"+urlencode({"limit":50,"offset":0}))
    root={"metadata_retrieval_status":children["status"]}
    if children["status"]=="PUBLIC_JSON_RECEIVED":
        items=children["data"].get("page",[])
        root["first_page_source_entity_types_no_ids"]=dict(Counter(
            str(x.get("type",x.get("concreteType","UNKNOWN"))).split(".")[-1] for x in items
        ))
        root["n_visible_first_page_entities"]=len(items)
        root["n_children_reported_by_API"]=children["data"].get("totalNumberOfResults")
        root["nextPageToken_present"]=bool(children["data"].get("nextPageToken"))
    if "http_status_code" in children:root["http_status_code"]=children["http_status_code"]
    return {
        "study":"Yang/Schadt HLC 9 CYP microsomal activities and liver GEO GSE9588",
        "original_Synapse_project_id":SYNAPSE_ROOT,
        "source_metadata":result,
        "source_hierarchy_summary":root,
        "verified_public_original_CYP_donor_activity_numeric_file_retrieved":False,
        "verified_original_GSE9588_GSM_to_CYP_activity_subject_ID_crosswalk":False,
        "original_same_donor_RNA_to_CYP_activity_heldout_benchmark_executed":False,
        "no_patient_and_private_Synapse_content_exported":True,
    }


if __name__=="__main__":
    print(json.dumps(audit_public_hlc_synapse(),indent=2,sort_keys=True))
