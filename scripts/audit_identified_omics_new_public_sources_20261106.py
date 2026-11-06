#!/usr/bin/env python3
"""Audit real public source repository metadata for new Omics2PK cohorts.

Do NOT label as source-acquired until donor-level genotype, omics AND
independent measured PK/function can all be crosswalked. This script
does not download human records and prints no raw biological rows.
All repositories are original public endpoints; no synthetic proof.
"""
from __future__ import annotations
import json
from collections import Counter
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError

SOURCES={
    "imatinib_genotype_target_protein":("PRIDE","PXD066482"),
    "tucatinib_CYP2C8_CYP3A4_CYP3A5":("PRIDE","PXD057282"),
    "liver_quant_proteomics_other_source":("PRIDE","PXD063438"),
    "pharmmet_liver_s9_drug_metabolite_library":("WORKBENCH","ST003991"),
}
BUDGET_BYTES=2_000_000

def _fetch_json(url):
    req=Request(url,headers={
        "Accept":"application/json",
        "User-Agent":"Omics2PK-public-source-original-access-audit/1.0",
    })
    with urlopen(req,timeout=35) as response:
        raw=response.read(BUDGET_BYTES+1)
        if len(raw)>BUDGET_BYTES:
            raise ValueError("Source public API output exceeded safe metadata budget")
        return json.loads(raw.decode("utf-8-sig"))

def _pride_file_list(accession):
    attempts=[]
    for version in ("v3","v2"):
        url=f"https://www.ebi.ac.uk/pride/ws/archive/{version}/projects/{accession}/files?pageSize=100"
        try:
            j=_fetch_json(url)
            # PRIDE v2 returns _embedded.files ; newer servers may return JSON list,
            # or a HAL-wrapped data list. Do not infer filenames from prose.
            if isinstance(j,list):items=j
            elif isinstance(j,dict):
                embedded=j.get("_embedded",{})
                items=(embedded.get("files")
                       or embedded.get("projectFiles")
                       or j.get("files")
                       or j.get("data")
                       or j.get("content")
                       or [])
            else:items=[]
            if not isinstance(items,list):
                raise ValueError("Unqualified PRIDE JSON file-list schema")
            if not items:
                attempts.append({"API_version":version,"status":"NO_STRUCTURED_FILE_RECORDS"})
                continue
            grouped=Counter()
            units=Counter()
            filetype=Counter()
            for f in items:
                if not isinstance(f,dict):
                    raise ValueError("PRIDE individual file listing unexpectedly nonobject")
                name=str(f.get("fileName") or f.get("filename") or f.get("name") or "")
                ext=name.rsplit(".",1)[-1].lower() if "." in name else "no_extension"
                grouped[ext]+=1
                category=str(f.get("fileCategory") or f.get("category") or f.get("fileType") or "unknown")
                filetype[category]+=1
                units[str(f.get("fileSizeBytes") is not None)]+=1
            return {
                "status":"REAL_PUBLIC_PRIDE_FILE_INDEX",
                "source_original_accession":accession,
                "API_version":version,
                "n_file_records_page1":len(items),
                "file_extension_counts":dict(sorted(grouped.items())),
                "source_file_category_counts":dict(sorted(filetype.items())),
                "vetted_source_numeric_donor_PGx_omics_function_join":False,
                "cannot_claim_patient_PK_from_raw_mass_spectra":True,
                "possible_pagination_unchecked":len(items)==100,
            }
        except (HTTPError,URLError,TimeoutError,ValueError,KeyError) as e:
            attempts.append({"API_version":version,"status":"SOURCE_METADATA_NOT_QUALIFIED",
                             "error_class":type(e).__name__})
    return {"status":"ABSTAIN_PRIDE_NO_QUALIFIED_PUBLIC_FILE_INDEX",
            "source_original_accession":accession,"attempts":attempts}

def _workbench_study(study):
    url=f"https://www.metabolomicsworkbench.org/rest/study/study_id/{study}/summary"
    try:
        j=_fetch_json(url)
        if not isinstance(j,(dict,list)):
            raise ValueError("Unqualified Metabolomics Workbench original summary")
        # Only whitelist abstract metadata keys: cannot emit clinical rows.
        out={}
        val=j if isinstance(j,dict) else {"records":len(j)}
        for key in ("study_id","study_title","title","study_summary","project_id",
                    "species","subject_type","analysis_type","records"):
            if key in val:
                out[key]=str(val[key])[:350]
        return {"status":"PUBLIC_ORIGINAL_WORKBENCH_STUDY_SUMMARY",
                "accession":study,
                "metadata":out,
                "source_donor_level_response_or_PK_not_verified":True,
                "not_clinical_PK":True}
    except (HTTPError,URLError,TimeoutError,ValueError,KeyError) as e:
        return {"status":"ABSTAIN_WORKBENCH_METADATA_FETCH",
                "accession":study,"error_class":type(e).__name__}

def _yee_supplement_probe():
    # Source supplement discovered in original PMC article. Probe only
    # content type + final redirect, no spreadsheet bytes or personal rows.
    base="https://pmc.ncbi.nlm.nih.gov/articles/PMC13601688/bin/"
    items={}
    for name in ("CPT-9999-0-s001.xlsx","CPT-9999-0-s003.xlsx"):
        u=base+name
        try:
            # This publisher may deny automated request; never say downloaded
            # unless HTTP response confirms correct MIME and file size.
            req=Request(u,method="HEAD",headers={
                "User-Agent":"Omics2PK-academic-source-audit/1.0",
            })
            with urlopen(req,timeout=25) as resp:
                items[name]={
                    "status":"ORIGINAL_PMC_SUPPLEMENT_ENDPOINT_RESPONDED",
                    "HTTP_status":resp.status,
                    "Content_Type":resp.headers.get("Content-Type"),
                    "Content_Length":resp.headers.get("Content-Length"),
                    "spreadsheet_contents_or_person_PK_crosswalk_examined":False,
                }
        except (HTTPError,URLError,TimeoutError,ValueError) as e:
            items[name]={"status":"SOURCE_SUPPLEMENT_HEAD_UNAVAILABLE",
                         "error_class":type(e).__name__}
    return items

def run():
    out={
        "status":"PUBLIC_IDENTIFIED_OMICS_SOURCE_ACCESS_TRIAGE_ONLY",
        "source_records_not_Pgx_PK_reproduction":True,
        "no_original_person_or_raw_sample_data_retrieved":True,
        "source_index":{}
    }
    for key,(kind,accession) in SOURCES.items():
        if kind=="PRIDE":
            out["source_index"][key]=_pride_file_list(accession)
        elif kind=="WORKBENCH":
            out["source_index"][key]=_workbench_study(accession)
    out["Yee_CYP3A_2026_original_PMC_supplements"]=_yee_supplement_probe()
    return out

if __name__=="__main__":
    print(json.dumps(run(),indent=2,sort_keys=True))
