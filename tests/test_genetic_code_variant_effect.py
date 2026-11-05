"""Synthetic only. Translation and transcript-coordinate interpretation contracts.

Do not interpret synthetic exon coordinates or transcripts as real CYP3A5.
A real clinical VCF needs original assembly/REF, VEP/MANE/PharmVar and
whole-target validation before making any PGx assertion.
"""
import unittest
from omics2pk.genetic_code_variant_effect import (
    CODON_TO_AA, GeneticCodeEffectError, annotate_coding_edit,
    annotate_genomic_snv, omic_implications, reconcile_curated_cyp3a5_allele,
)

SOURCE="https://www.ncbi.nlm.nih.gov/refseq/"
CPIC="https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/"
SEQ="ATG"+"GAA"+"GCT"*28+"TAA"  # 93 bases; reference CDS ends in natural stop


def transcript(*, strand="+", exons=None, gene="CYP3A5"):
    return {
        "transcript_id":"NM_000777.5", "gene":gene,
        "assembly":"GRCh38","chromosome":"7","strand":strand,
        "reference_source_url":SOURCE,
        "spliced_cdna":SEQ,
        "exons":exons if exons is not None else (
            [[100,184],[400,407]] if strand=="+" else [[400,407],[100,184]]),
        "cds_start_offset":0,"cds_end_offset":len(SEQ),
    }


def variant(pos, ref, alt):
    return {
        "gene":"CYP3A5","assembly":"GRCh38","chromosome":"7",
        "genomic_position_1based":pos,"genomic_ref_allele":ref,
        "genomic_alt_allele":alt,"variant_source_url":SOURCE,
    }


ASSAY={"status":"partial_assay"}


class GeneticCodeConsequencesTest(unittest.TestCase):
    def test_standard_nuclear_code_reference(self):
        self.assertEqual(len(CODON_TO_AA),64)
        for codon,expected in [("ATG","M"),("GAA","E"),("TAA","*"),
                               ("TGA","*"),("TAG","*"),("GCT","A"),
                               ("GCC","A"),("TTC","F"),("TAC","Y")]:
            self.assertEqual(CODON_TO_AA[codon],expected)

    def test_exon_missense_can_have_normal_arn_and_impaired_protein(self):
        # Coding exon GAA -> GCA: glutamate -> alanine, function *UNKNOWN*.
        result=annotate_genomic_snv(
            reference_transcript=transcript(),variant=variant(104,"A","C"),
            genome_assay=ASSAY)
        self.assertEqual(result["sequence_consequence"],"missense_variant")
        self.assertEqual(result["reference_codon"],"GAA")
        self.assertEqual(result["observed_alternate_codon"],"GCA")
        self.assertEqual(result["predicted_protein_amino_acid_change"],"E2A")
        self.assertEqual(result["aa_position_1based"],2)
        self.assertFalse(result["clinical_protein_function_estimated"])
        implications=omic_implications(result)
        self.assertIn("validated_functional_activity",
                      implications["experimentally_informative_assay_classes"])
        self.assertTrue(implications["RNA_abundance_is_not_functional_protein"])
        self.assertFalse(implications["mutation_causes_loss_of_function_proven"])

    def test_exonic_synonymous_is_not_always_functionally_neutral(self):
        # GCT Ala -> GCC Ala, no amino-acid change; splicing remains unknown.
        result=annotate_genomic_snv(
            reference_transcript=transcript(),variant=variant(108,"T","C"),
            genome_assay=ASSAY)
        self.assertEqual(result["sequence_consequence"],"synonymous_variant")
        self.assertEqual(result["reference_aa"],"A")
        self.assertEqual(result["alternate_aa"],"A")
        self.assertFalse(result["splice_effect_verified"])
        self.assertIn("splice_junction_RNA",
                      omic_implications(result)["experimentally_informative_assay_classes"])

    def test_stop_gained_early_EJC_nmd_plausible_but_not_proven(self):
        # GAA -> TAA (stop) at cDNA 3, last exon junction at offset 85.
        result=annotate_genomic_snv(reference_transcript=transcript(),
                                   variant=variant(103,"G","T"),
                                   genome_assay=ASSAY)
        self.assertEqual(result["sequence_consequence"],"stop_gained")
        self.assertEqual(result["observed_alternate_codon"],"TAA")
        self.assertEqual(result["new_stop_cDNA_offset"],3)
        self.assertEqual(result["candidate_NMD_not_observed"],
                         "CLASSICAL_EJC_NMD_PLAUSIBLE_NOT_CONFIRMED")
        self.assertFalse(result["clinical_protein_function_estimated"])

    def test_native_stop_lost_not_naively_called_elongated_active_protein(self):
        result=annotate_genomic_snv(reference_transcript=transcript(),
                                   variant=variant(405,"T","C"),
                                   genome_assay=ASSAY)
        self.assertEqual(result["sequence_consequence"],"stop_lost")
        self.assertEqual(result["reference_aa"],"*")
        self.assertEqual(result["alternate_aa"],"Q")
        self.assertTrue(result["protein_misfolding_trafficking_and_catalytic_activity_unknown"])

    def test_in_frame_vs_frameshift_transcript_edit(self):
        source=transcript()
        fs=annotate_coding_edit(ref=source,transcript_offset=6,
                                ref_allele="G",alt_allele="GA")
        inframe=annotate_coding_edit(ref=source,transcript_offset=6,
                                     ref_allele="G",alt_allele="GAAA")
        self.assertEqual(fs["sequence_consequence"],"frameshift_variant")
        self.assertEqual(inframe["sequence_consequence"],
                         "inframe_insertion_or_deletion")
        self.assertFalse(fs["clinical_protein_function_estimated"])
        with self.assertRaises(GeneticCodeEffectError):
            annotate_coding_edit(ref=source,transcript_offset=6,
                                 ref_allele="A",alt_allele="GA")

    def test_intronic_donor_acceptor_and_deep_intron_separate(self):
        for pos,kind in [(185,"splice_donor_variant"),
                         (399,"splice_acceptor_variant"),
                         (190,"splice_region_variant"),
                         (250,"intron_variant")]:
            with self.subTest(pos=pos):
                result=annotate_genomic_snv(reference_transcript=transcript(),
                                            variant=variant(pos,"A","G"),
                                            genome_assay=ASSAY)
                self.assertEqual(result["sequence_consequence"],kind)
                self.assertIsNone(result["coding_aa_change"])
                self.assertFalse(result["splicing_observed_in_RNA"])
                self.assertFalse(result["genomic_reference_base_checked_against_original_genome"])

    def test_negative_strand_complements_forward_genomic_snv(self):
        # On negative strand: forward genomic T->G equals oriented A->C.
        result=annotate_genomic_snv(
            reference_transcript=transcript(strand="-"),
            variant=variant(403,"T","G"),genome_assay=ASSAY)
        self.assertEqual(result["sequence_consequence"],"missense_variant")
        self.assertEqual(result["transcript_reference_allele"],"A")
        self.assertEqual(result["transcript_alt_allele"],"C")
        self.assertEqual(result["predicted_protein_amino_acid_change"],"E2A")
        for pos,kind in [(399,"splice_donor_variant"),
                         (185,"splice_acceptor_variant")]:
            with self.subTest(pos=pos):
                result=annotate_genomic_snv(
                    reference_transcript=transcript(strand="-"),
                    variant=variant(pos,"A","G"),genome_assay=ASSAY)
                self.assertEqual(result["sequence_consequence"],kind)

    def test_no_variant_assayed_is_not_proof_of_normal_function(self):
        r=annotate_genomic_snv(reference_transcript=transcript(),
                               variant=variant(104,"A","C"),
                               genome_assay={"status":"complete_target_assay",
                                             "covered_transcript_ids":["NM_000777.5"],
                                             "assay_provenance_url":SOURCE})
        self.assertEqual(r["original_genetic_assay_status"],"complete_target_assay")
        self.assertTrue(r["no_variant_detected_is_not_proof_of_unmodified_function"])
        self.assertFalse(r["star_allele_diplotype_verified"])
        with self.assertRaises(GeneticCodeEffectError):
            annotate_genomic_snv(
                reference_transcript=transcript(),variant=variant(104,"A","C"),
                genome_assay={"status":"complete_target_assay"})
        with self.assertRaises(GeneticCodeEffectError):
            annotate_genomic_snv(
                reference_transcript=transcript(),variant=variant(104,"A","C"),
                genome_assay={"status":"wildtype"})
        with self.assertRaises(GeneticCodeEffectError):
            annotate_genomic_snv(
                reference_transcript=transcript(),variant=variant(104,"A","C"),
                genome_assay={"status":"not_assessed"})

    def test_reference_allele_and_assembly_drift_abstains(self):
        bad=[
            variant(104,"G","C"),
            {**variant(104,"A","C"),"assembly":"GRCh37"},
            {**variant(104,"A","C"),"gene":"ABCB1"},
            {**variant(104,"A","C"),"genomic_position_1based":False},
            {**variant(104,"A","C"),"chromosome":"6"},
        ]
        for v in bad:
            with self.subTest(v=v),self.assertRaises(GeneticCodeEffectError):
                annotate_genomic_snv(reference_transcript=transcript(),
                                    variant=v,genome_assay=ASSAY)
        bad_ref={**transcript(),"exons":[[100,185],[400,407]]}
        with self.assertRaises(GeneticCodeEffectError):
            annotate_genomic_snv(reference_transcript=bad_ref,
                                variant=variant(104,"A","C"),genome_assay=ASSAY)

    def test_synonymous_cyp3a5_star6_curated_splice_function_overrides_codon_only(self):
        effect={"gene":"CYP3A5","transcript_id":"NM_000777.5",
                "sequence_consequence":"synonymous_variant"}
        out=reconcile_curated_cyp3a5_allele(
            sequence_effect=effect,source_variant_id="rs10264272",
            hgvs_c="NM_000777.5:c.624G>A",
            star_allele="*6",source_url=CPIC,
            variant_call_complete=True)
        self.assertEqual(out["genetic_code_result"],"synonymous_variant")
        self.assertIn("ABERRANT_SPLICING",out["independent_literature_mechanism"])
        self.assertEqual(out["published_allele_function"],"no_function")
        self.assertFalse(out["in_this_subject_splicing_measured"])
        self.assertFalse(out["full_diplotype_or_variant_phasing_established"])
        with self.assertRaises(GeneticCodeEffectError):
            reconcile_curated_cyp3a5_allele(
                sequence_effect=effect,source_variant_id="rs10264272",
                hgvs_c="NM_000777.5:c.624G>A",
                star_allele="*6",source_url=CPIC,
                variant_call_complete=False)
        with self.assertRaises(GeneticCodeEffectError):
            reconcile_curated_cyp3a5_allele(
                sequence_effect={**effect,"sequence_consequence":"missense_variant"},
                source_variant_id="rs10264272",
                hgvs_c="NM_000777.5:c.624G>A",
                star_allele="*6",source_url=CPIC,
                variant_call_complete=True)


if __name__=="__main__":
    unittest.main()
