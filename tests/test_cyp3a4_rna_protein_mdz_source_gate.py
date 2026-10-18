"""Source qualification, held-human-CV and privacy synthetic tests (NOT original donor values)."""
from __future__ import annotations
import unittest

from omics2pk.cyp3a4_rna_protein_mdz_source_gate import (
    SOURCE_DOI, SourceNotQualified,source_contract,evaluate_original_source
)


def source_rows():
    result=[]
    for i in range(74):
        hist=("normal" if i<24 else "nafl" if i<50 else "nash")
        age=22+(i*7)%58
        mrna=0.8+((i*17)%87)/14
        protein=2.+((i*29)%71)/5
        activity=1.+2*protein+0.2*mrna+(i%5)*0.1
        result.append({
            "donor_id":f"SYNTHETIC_LIVER_{i:03d}",
            "source_doi":SOURCE_DOI,
            "row_kind":"individual",
            "age_years":str(age),"sex":"male" if i%2 else "female",
            "histology":hist,
            "CYP3A4_mRNA_relative_18S":str(mrna),
            "CYP3A4_protein_pmol_mg":str(protein),
            "midazolam_1hydroxylation_pmol_min_mg":str(activity),
        })
    return result


class TargetedSameLiverSourceTests(unittest.TestCase):
    def test_74_fictional_donor_model_comparisons_have_identical_folds(self):
        a=evaluate_original_source(source_rows())
        self.assertEqual(a["n_unique_independent_original_liver_donors"],74)
        self.assertEqual(a["n_original_donors_by_histology"],
                         {"normal":24,"nafl":26,"nash":24})
        self.assertIn("clinical_plus_RNA_and_protein",a["mae_by_model_original_activity_units"])
        self.assertLess(a["primary_delta_MAE_clinical_plus_RNA_protein_minus_clinical_plus_RNA"],0)
        self.assertTrue(a["not_measured_clinical_plasma_PK"])
        self.assertTrue(a["not_original_whole_omics_multiomics"])
        self.assertNotIn("SYNTHETIC_LIVER",str(a))
        self.assertTrue(a["original_individual_donor_rows_and_predictions_exported"] is False)

    def test_no_protein_per_person_no_model_fit(self):
        r=source_rows()
        r[5]["CYP3A4_protein_pmol_mg"]=""
        with self.assertRaises(SourceNotQualified):
            evaluate_original_source(r)

    def test_wrong_gene_assay_or_study_blocked(self):
        r=source_rows()
        r[4]["source_doi"]="10.9999/not_the_actual_study"
        with self.assertRaisesRegex(SourceNotQualified,"DOI"):
            source_contract(r)

    def test_no_source_replicates_pseudoreplication(self):
        r=source_rows()
        r[1]["donor_id"]=r[0]["donor_id"]
        with self.assertRaisesRegex(SourceNotQualified,"duplicated"):
            source_contract(r)

    def test_original_24_26_24_groups_required(self):
        r=source_rows()
        r[26]["histology"]="normal"
        with self.assertRaisesRegex(SourceNotQualified,"24 normal"):
            source_contract(r)

    def test_extra_simulated_PBPK_subject_invalid(self):
        r=source_rows()
        r[10]["row_kind"]="PBPK_simulated"
        with self.assertRaisesRegex(SourceNotQualified,"simulated"):
            source_contract(r)

    def test_same_people_cannot_be_joined_by_missing_code(self):
        r=source_rows()
        r[20]["donor_id"]=""
        with self.assertRaisesRegex(SourceNotQualified,"Missing"):
            source_contract(r)

    def test_original_gene_expression_not_whole_transcriptome(self):
        report=evaluate_original_source(source_rows())
        self.assertFalse(report["whole_transcriptome_assayed_in_this_model"])
        self.assertFalse(report["whole_proteome_assayed_in_this_model"])


if __name__=="__main__":
    unittest.main()
