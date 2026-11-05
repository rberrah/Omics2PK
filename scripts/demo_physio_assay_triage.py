#!/usr/bin/env python3
"""Offline, metadata-only demonstration. NOT original patients or observed PK.

The CYP3A5 source allele *3 mechanism is literature qualified, but the
hypothetical assay options are NOT performed on a donor; none of the numerical
benefits, observed activity or PK outcomes are filled in.
"""
from __future__ import annotations
import json

from omics2pk.assay_aware_physio_triage import plan_physio_omic_acquisition

S="https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/"
T="https://pmc.ncbi.nlm.nih.gov/articles/PMC5709222/"
P="https://doi.org/10.1124/dmd.113.054825"

def assay(assay_id, method, tissue, gene=None, **kw):
    d={"assay_id":assay_id,"method":method,"tissue":tissue,
       "status":"proposed","collection_phase":"predose",**kw}
    if gene is not None:
        d["gene"]=gene
    return d

def demonstrate():
    tac=plan_physio_omic_acquisition(
        drug="tacrolimus",endpoint="AUC",
        variants=[{"gene":"CYP3A5","consequence":"splice_disrupting",
                   "star_allele":"*3","call_status":"complete",
                   "allele_state":"biallelic",
                   "source_variant_id":"NM_000777.5:c.219-237A>G / rs776746",
                   "transcript_id":"NM_000777.5","source_url":S}],
        assays=[
            assay("short_3prime_tag","RNA_3PRIME_TAG","liver","CYP3A5"),
            assay("total_RNA","RNA_GENE_BULK","liver","CYP3A5"),
            assay("isoform_RNA","RNA_LONG_READ_ISOFORM","liver","CYP3A5",
                  transcript_id="NM_000777.5"),
            assay("N_terminal_peptide","PROTEIN_TARGETED_PEPTIDE","liver","CYP3A5",
                  peptide_region="n_terminal"),
            assay("C_terminal_peptide","PROTEIN_TARGETED_PEPTIDE","liver","CYP3A5",
                  peptide_region="c_terminal"),
            assay("unresolved_bottomup","PROTEIN_BOTTOM_UP_SHARED_PEPTIDE","liver","CYP3A5"),
            assay("topdown_unverified","PROTEIN_TOP_DOWN_INTACT","liver","CYP3A5",
                  intact_identity_verified=False),
            assay("CYP3A5_selective_assay","ACTIVITY_CYP3A5_SELECTIVE_PROBE","liver","CYP3A5",
                  assay_validated=True,selectivity_evidence_url=P),
            assay("shared_midazolam_probe","ACTIVITY_CYP3A_SHARED_PROBE","liver","CYP3A5",
                  assay_validated=True),
            assay("plasma_4bHC_candidate","METABOLITE_MSMS_TARGETED","plasma",
                  candidate_genes=["CYP3A4","CYP3A5"],chemical_identity_grade="source_name_only"),
        ])
    met=plan_physio_omic_acquisition(
        drug="metformin",endpoint="CLr",
        assays=[
            assay("renal_OCT2_RNA","RNA_GENE_BULK","kidney","SLC22A2"),
            assay("renal_MATE1_protein","PROTEIN_TARGETED_PEPTIDE","kidney","SLC47A1",
                  peptide_region="unknown"),
            assay("OCT2_blood_RNA","RNA_GENE_BULK","whole_blood","SLC22A2"),
            assay("renal_OCT2_transport_function","ACTIVITY_TRANSPORT_SUBSTRATE",
                  "kidney","SLC22A2",assay_validated=True,substrate_drug="metformin"),
            assay("hepatic_OCT1_RNA_for_CLr","RNA_GENE_BULK","liver","SLC22A1"),
        ])
    return {
        "demo_type":"HYPOTHETICAL_ASSAY_METADATA_ONLY_NOT_PATIENT_PK",
        "is_original_patient_study":False,
        "data_driven_expected_performance_claim":False,
        "scenario_tacrolimus_CYP3A5_star3_splicing":tac,
        "scenario_metformin_renal_transport_not_Pgp":met,
    }

if __name__=="__main__":
    print(json.dumps(demonstrate(),ensure_ascii=False,sort_keys=True,indent=2))
