"""Original EPA 100-human-liver donor-matched RNA/enzymatic proof, no raw data output."""
import hashlib
import json

from scripts.audit_epa_liver_rna_activity import _get, SOURCE_DOI, PAPER_DOI
from scripts.inspect_epa_liver_crosswalk import parse_epa_workbook
from omics2pk.epa_liver_cyp_validation import (
    EXPECTED_GENES, EPALinkageError, read_paired_epa, validate_one_gene,
)

EXPECTED_SHA = "cfc9dc334064ea61e99ca6fd39c2c8ef7c621d000a37c4da2c95f26ca85e36dd"


def validate_public_epa_livers():
    data=_get()
    if hashlib.sha256(data).hexdigest()!=EXPECTED_SHA:
        raise EPALinkageError("US EPA original DOI workbook content changed; requalify before fitting")
    original=parse_epa_workbook(data)
    measurements,qc=read_paired_epa(original)
    endpoint_results={}
    for gene in EXPECTED_GENES:
        records=measurements[gene]
        try:
            endpoint_results[gene]={"status":"MEASURED_DONOR_HELDOUT",**validate_one_gene(records)}
        except EPALinkageError as exc:
            endpoint_results[gene]={
                "status":"ABSTAIN_INSUFFICIENT_LINKED_COMPLETE_ORIGINAL_DONORS",
                "n_eligible":len(records),
                "reason":str(exc)
            }
    return {
        "status":"ORIGINAL_HUMAN_LIVER_IN_VITRO_CYP_FUNCTION_NOT_PATIENT_PK",
        "source_repository_doi":SOURCE_DOI,
        "paper_doi":PAPER_DOI,
        "exact_document_sha256":EXPECTED_SHA,
        "unique_original_donors_linked_by_two_source_identifiers":100,
        "n_predeclared_cyp_gene_activity_pairs":len(EXPECTED_GENES),
        "n_fitted_pairs":sum(x["status"]=="MEASURED_DONOR_HELDOUT" for x in endpoint_results.values()),
        "n_rna_adds_to_age_sex":sum(x.get("RNA_reduces_adjusted_MAE_in_original_cohort",False) for x in endpoint_results.values()),
        "source_assay_qc":qc,
        "results":endpoint_results,
        "original_donor_and_raw_omics_values_exported":False,
        "limitations":[
            "Experimentally measured microsomal CYP probe substrate activity, not patient PK concentration, clearance or dose",
            "RNA is targeted qPCR, not global transcriptomics; gene panels prespecified, no SWEET or multiomics",
            "No PGx in this public matched workbook; comparison is beyond age/sex, NOT beyond clinical+PGx PopPK",
            "Same 100 donors contribute multiple CYP endpoints; eight outcomes are NOT eight independent cohorts",
            "Exploratory donor-level source reuse; bootstrap is paired-loss not nested refit, eight endpoint multiplicity not controlled",
            "RNA/protein/enzyme may be correlated and are not independent mechanistic intervention effects"
        ]
    }

if __name__=="__main__":
    print(json.dumps(validate_public_epa_livers(),sort_keys=True,indent=2,ensure_ascii=False))
