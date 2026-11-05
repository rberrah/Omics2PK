"""CYP3A4/CYP3A5 assay-role registry for functional-evidence qualification.

CYP3A5-selective is NOT mathematically exclusive. Metabolite formed during a
substrate-probe activity assay is the *functional response* and cannot be
recycled as baseline/pre-landmark omics. All source labels are fixed, documented
research categories, not a fitted probabilistic activity score.
"""
from __future__ import annotations


class CYP3AAssayError(ValueError):
    pass


ASSAY_ROLES = {
    "T5_N_OXIDATION": {
        "reference": "10.1124/dmd.113.054726",
        "class": "CYP3A5_SELECTIVE_PROBE",
        "readout": "functional_activity_assay",
        "supports_specific_CYP3A5_E2_target": True,
        "independent_baseline_omics": False,
        "note": "Strongly CYP3A5-selective in vitro; not fully CYP3A5-exclusive.",
    },
    "MIDAZOLAM_1_HYDROXYLATION": {
        "reference": "10.1124/dmd.119.088823",
        "class": "SHARED_CYP3A_PROBE",
        "readout": "functional_activity_assay",
        "supports_specific_CYP3A5_E2_target": False,
        "independent_baseline_omics": False,
        "note": "CYP3A4 and CYP3A5 can both catalyze this substrate reaction.",
    },
    "TESTOSTERONE_6B_HYDROXYLATION": {
        "reference": "10.1124/dmd.119.088823",
        "class": "SHARED_CYP3A_PROBE",
        "readout": "functional_activity_assay",
        "supports_specific_CYP3A5_E2_target": False,
        "independent_baseline_omics": False,
        "note": "CYP3A activity, not CYP3A5-specific functional phenotype.",
    },
    "FOUR_BETA_HYDROXYCHOLESTEROL_ABUNDANCE": {
        "reference": "10.1002/jcph.1391",
        "class": "ENDOGENOUS_SHARED_CYP3A_METABOLITE",
        "readout": "endogenous_molecular_biomarker",
        "supports_specific_CYP3A5_E2_target": False,
        "independent_baseline_omics": True,
        "note": "Slow-changing CYP3A marker; baseline independent assay required.",
    },
    "LAPATINIB_M1_FORMATION": {
        "reference": "10.1124/dmd.119.088823",
        "class": "DRUG_INCUBATION_METABOLISM_OUTCOME",
        "readout": "drug_metabolite_formation",
        "supports_specific_CYP3A5_E2_target": False,
        "independent_baseline_omics": False,
        "note": "Measured after drug addition, so not a pretreatment covariate.",
    },
    "VINORELBINE_PARENT_DEPLETION": {
        "reference": "10.1124/dmd.113.051094",
        "class": "DRUG_INCUBATION_METABOLISM_OUTCOME",
        "readout": "observed_drug_depletion",
        "supports_specific_CYP3A5_E2_target": False,
        "independent_baseline_omics": False,
        "note": "Total drug-specific metabolism with CYP3A4 and CYP3A5 roles.",
    },
}


def classify_cyp3a_assay(
    *,
    assay_code: str,
    assay_validated: bool,
    source_doi: str,
) -> dict:
    """Explicit assay identity; fail when missing, invented or unvalidated.

    Matching article DOI is mandatory for the *source of assay evidence*;
    sample origin and donor identities must be independently qualified.
    """
    if assay_code not in ASSAY_ROLES:
        raise CYP3AAssayError("Unknown CYP3A probe or biomarker identity")
    record = ASSAY_ROLES[assay_code]
    if assay_validated is not True:
        raise CYP3AAssayError("Assay validation not attested")
    if not isinstance(source_doi,str) or source_doi.strip()!=record["reference"]:
        raise CYP3AAssayError("Assay evidence DOI mismatch")
    return {
        "assay_code":assay_code,
        "source_reference":record["reference"],
        "assay_category":record["class"],
        "readout":record["readout"],
        "CYP3A5_selective_activity_endpoint_eligible":
            record["supports_specific_CYP3A5_E2_target"],
        "independent_baseline_molecular_measurement_possible":
            record["independent_baseline_omics"],
        "interpretation_note":record["note"],
        "no_fixed_genotype_activity_percentages":True,
        "not_a_patient_specific_clearance_or_dosing_model":True,
    }
