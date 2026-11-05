#!/usr/bin/env python3
"""Offline synthetic reference, NOT a real CYP3A5 transcript or donor.

Demonstrates why one nucleotide substitution can be synonymous, missense or
stop-gained and why the absence of a called variant is not proven activity.
No patient data, fitted model, observed omic intensity, PK or protein function.
"""
from __future__ import annotations
import json

from omics2pk.genetic_code_variant_effect import (
    annotate_genomic_snv, omic_implications, reconcile_curated_cyp3a5_allele,
)
from omics2pk.assay_aware_physio_triage import plan_physio_omic_acquisition

CDS="ATG"+"GAA"+"GCT"*28+"TAA"
REFERENCE={
    "gene":"CYP3A5","transcript_id":"NM_999999.1",  # SYNTHETIC
    "assembly":"GRCh38","chromosome":"7","strand":"+",
    "reference_source_url":"https://www.ncbi.nlm.nih.gov/refseq/",
    "spliced_cdna":CDS, "exons":[[100,184],[400,407]],
    "cds_start_offset":0,"cds_end_offset":len(CDS),
}
SOURCE="https://www.ncbi.nlm.nih.gov/refseq/"
DNA_ASSAY={"status":"partial_assay"}

def snv(pos,ref,alt):
    return {
        "gene":"CYP3A5","assembly":"GRCh38","chromosome":"7",
        "genomic_position_1based":pos,"genomic_ref_allele":ref,
        "genomic_alt_allele":alt,"variant_source_url":SOURCE,
    }

def assess(pos,ref,alt):
    source=snv(pos,ref,alt)
    genomic=annotate_genomic_snv(reference_transcript=REFERENCE,variant=source,
                                genome_assay=DNA_ASSAY)
    return {
        "genetic_code_effect":genomic,
        "hypotheses_for_omic_measurement":omic_implications(genomic),
    }

def run():
    examples={
        "missense_ARN_can_be_present_but_function_unproven":assess(104,"A","C"),
        "stop_gained_new_termination_and_possible_NMD":assess(103,"G","T"),
        "synonymous_same_amino_acid_but_splicing_not_excluded":assess(108,"T","C"),
        "intron_splice_donor_RNA_needed":assess(185,"A","G"),
    }
    # Separate real *published variant definition*, not a synthetic
    # patient's actual allele. This illustrates why translation-only
    # annotation for CYP3A5 *6 would miss the splice mechanism.
    cyp3a5_star6=reconcile_curated_cyp3a5_allele(
        sequence_effect={
            "gene":"CYP3A5","transcript_id":"NM_000777.5",
            "sequence_consequence":"synonymous_variant",
        },
        source_variant_id="rs10264272",hgvs_c="NM_000777.5:c.624G>A",
        star_allele="*6",
        source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/",
        variant_call_complete=True,
    )
    assays=[{"assay_id":"RNA_bulk","method":"RNA_GENE_BULK","gene":"CYP3A5",
             "tissue":"liver","status":"proposed","collection_phase":"predose"},
            {"assay_id":"RNA_assembled_isoform","method":"RNA_LONG_READ_ISOFORM",
             "gene":"CYP3A5","tissue":"liver","status":"proposed",
             "collection_phase":"predose","transcript_id":"NM_999999.1"},
            {"assay_id":"protein_C_terminal","method":"PROTEIN_TARGETED_PEPTIDE",
             "gene":"CYP3A5","tissue":"liver","status":"proposed",
             "collection_phase":"predose","peptide_region":"c_terminal"},
            {"assay_id":"enzyme_function_T5","method":"ACTIVITY_CYP3A5_SELECTIVE_PROBE",
             "gene":"CYP3A5","tissue":"liver","status":"proposed",
             "collection_phase":"predose","assay_validated":True,
             "probe_assay_code":"T5_N_OXIDATION",
             "selectivity_evidence_url":"https://doi.org/10.1124/dmd.113.054726"}]
    stop_design=plan_physio_omic_acquisition(
        drug="tacrolimus",endpoint="AUC",assays=assays,
        sequence_variants=[{"reference_transcript":REFERENCE,
                            "variant":snv(103,"G","T"),
                            "genome_assay":DNA_ASSAY}],
    )
    no_detected_variants=plan_physio_omic_acquisition(
        drug="tacrolimus",endpoint="AUC",assays=assays,
        genome_assay_background={
            "status":"complete_target_assay",
            "covered_transcript_ids":["NM_999999.1"],
            "assay_provenance_url":SOURCE})
    return {
        "provenance":"SYNTHETIC_REFERENCE_AND_GENOTYPES_NOT_REAL_CYP3A5_SEQUENCES_OR_PATIENTS",
        "not_original_genotyping_or_pharmacokinetics":True,
        "genetic_code_synthetic_examples":examples,
        "curated_literature_cyp3a5_star6_known_synonymous_splice_override":cyp3a5_star6,
        "stop_gained_physio_assay_design":stop_design,
        "no_variant_reported_does_not_prove_normal_protein":no_detected_variants,
    }

if __name__=="__main__":
    print(json.dumps(run(),sort_keys=True,indent=2,ensure_ascii=False))
