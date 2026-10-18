"""Source-accurate, fail-closed audit of genuine human predose omics + measured PK.

Original manuscript records ONLY; public JATS XML and supplementary links.
Never turn an article's n=132, a plotted scatter, model R² or a list of 841
metabolites into a donor-level metabolome × measured BuCL matrix.

IMPORTANT McCune 2023: the 96 '2-week-preBU' samples are biological subjects
drawn FROM the 132 'preBU' training subjects, so original paper transfer
across timepoints is NOT a held-out independent-person validation.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import json
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError
from xml.etree import ElementTree as ET


PRIMARY={
    "busulfan_2023_predose_CL": {
        "pmcid":"PMC9888309","doi":"10.1002/cpt.2794",
        "reported_unique_humans":132,
        "reported_omics_features":841,
        "reported_immediate_predose_specimens":132,
        "reported_two_week_predose_specimens":96,
        "reported_biologically_shared_twoweek_people":96,
        "original_target":"Bayesian estimated IV busulfan clearance mL/min/kg normal fat mass",
    },
    "busulfan_2016_predose_CL": {
        "pmcid":"PMC5562150","doi":"10.1007/s11306-016-1106-6",
        "reported_unique_humans":108,
        "reported_predose_plasma_specimens":147,
        "reported_discovery_specimens":59,
        "reported_validation_specimens":88,
        "reported_omics_features":1885,
        "original_target":"Bayesian estimated IV busulfan clearance",
    },
    "busulfan_2021_metabolome_gvhd": {
        "pmcid":"PMC8214873","doi":"10.1021/acs.jproteome.0c00599",
        "reported_unique_humans":84,
        "reported_longitudinal_specimens":230,
        "reported_omics_features":763,
        "original_target":"GVHD and relapse with PK-guided busulfan therapy",
    },
    "jarmusch_2020_CYP_cocktail": {
        "pmcid":"PMC7485946","doi":"10.1111/cts.12785",
        "reported_unique_humans":14,
        "original_target":"measured CYP cocktail targeted PK and untargeted mass spectrometry",
    },
}


def official_jats_source_audit(record:dict)->dict:
    if not record["pmcid"].startswith("PMC") or not record["doi"]:
        raise ValueError("Original official source DOI/PMCID required")
    url=("https://www.ebi.ac.uk/europepmc/webservices/rest/"
         +record["pmcid"]+"/fullTextXML")
    try:
        with urlopen(Request(url,headers={
            "User-Agent":"Omics2PK-original-dose-timed-PK-source-availability/1.0"
        }),timeout=40) as resp:
            xml=resp.read(8_000_001)
        if len(xml)>8_000_000:
            raise ValueError("Original JATS source larger than validated bounded XML")
        root=ET.fromstring(xml)
    except (HTTPError,URLError,TimeoutError,OSError,ET.ParseError,ValueError) as err:
        return {
            "source_access_status":"ORIGINAL_PMC_JATS_NOT_RETRIEVED_ABSTAIN",
            "source_error_type":type(err).__name__,
            "real_original_person_matched_PK_omics_matrix_downloaded":False,
        }

    def tag(x):return x.tag.rsplit("}",1)[-1]
    # Require exact paper provenance, not a DOI cited as a background ref.
    front=next((x for x in root if tag(x)=="front"),None)
    article_ids=(
        [(n.attrib.get("pub-id-type"),"".join(n.itertext()).strip())
         for n in front.iter() if tag(n)=="article-id"]
        if front is not None else []
    )
    doi_valid=any(k=="doi" and v.lower()==record["doi"].lower()
                  for k,v in article_ids)
    pmc_valid=any(k in ("pmc","pmcid") and
                  v.lower().removeprefix("pmc")==record["pmcid"].lower().removeprefix("pmc")
                  for k,v in article_ids)
    verified=doi_valid and pmc_valid
    if not verified:
        return {"source_access_status":"OFFICIAL_JATS_DOI_PMCID_MISMATCH_ABSTAIN",
                "original_article_ids_crosschecked":False,
                "source_DOI_metadata_exact_match":doi_valid,
                "source_PMCID_metadata_match":pmc_valid,
                "article_id_schema_types_only":[k for k,_ in article_ids],
                "real_original_person_matched_PK_omics_matrix_downloaded":False}

    supplements=[]
    for node in root.iter():
        if tag(node)!="supplementary-material":continue
        links=[]
        for child in node.iter():
            if tag(child) not in {"media","graphic","ext-link"}:continue
            for attr,val in child.attrib.items():
                if attr.rsplit("}",1)[-1] in ("href","filename") and val:
                    links.append(val.split("/")[-1].lower())
        supplements.extend(links)
    categories=Counter()
    for filename in supplements:
        extension=filename.rsplit(".",1)[-1] if "." in filename else "unknown"
        categories[extension]+=1
    return {
        "source_access_status":"ORIGINAL_PMC_JATS_DOI_VERIFIED_AND_SUPPLEMENT_METADATA_INSPECTED",
        "original_article_ids_crosschecked":True,
        "original_fulltext_JATS_bytes":len(xml),
        "original_fulltext_JATS_SHA256":sha256(xml).hexdigest(),
        "n_official_supplement_file_links":len(supplements),
        "reported_original_supplement_extension_counts":dict(sorted(categories.items())),
        "individual_original_supplements_actually_downloaded_and_column_validated":False,
        "real_original_person_matched_PK_omics_matrix_downloaded":False,
        "original_CYP_or_GST_genotype_same_humans_qualified":False,
        "original_raw_human_records_or_filename_identifiers_exported":False,
    }


def audit_all_original_PK_sources():
    out={}
    for key,details in PRIMARY.items():
        observed=official_jats_source_audit(details)
        out[key]={
            "study":key,"doi":details["doi"],"pmcid":details["pmcid"],
            "published_sample_counts_are_NOT_proof_of_independent_subjects":True,
            "original_primary_paper_sample_metadata":{
                k:v for k,v in details.items() if k.startswith("reported_")
            },
            "published_original_endpoint":details["original_target"],
            "scientific_source_gate":observed,
            "source_patient_level_molecular_and_target_PK_model_eligible":False,
        }
    out["busulfan_2023_predose_CL"]["study_specific_validation_warning"]=(
        "96 two-week predose individuals are part of the original 132-person "
        "immediate predose training cohort. Published timepoint-transfer R2 is "
        "NOT independent patient-heldout or external validation."
    )
    return {"status":"REAL_ORIGINAL_PK_PUBLIC_SOURCE_GATE_ABSTAIN_UNTIL_ORIGINAL_INDIVIDUAL_ROWS",
            "original_studies_audited":len(out),
            "n_new_source_qualified_genotype_vs_omics_and_PK_datasets":0,
            "original_patient_rows_exposed":False,
            "sources":out}


if __name__=="__main__":
    print(json.dumps(audit_all_original_PK_sources(),sort_keys=True,indent=2))
