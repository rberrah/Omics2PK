"""Actual NCI60 modality choice semantics, purely synthetic."""
import unittest
from omics2pk.modalities_information_budget import (
    ModalityBudgetError,choose_minimum_modalities
)
MODES={"tumor_only":(),"rna":("rna",),"protein":("prot",),
       "met":("met",),"rna_protein":("rna","prot"),
       "rna_met":("rna","met"),"prot_met":("prot","met"),
       "all":("rna","prot","met")}
class MultimodalCosts(unittest.TestCase):
    def test_cheapest_single_source_80_vs_full90(self):
        losses={"tumor_only":.35,"rna":.34,"protein":.34,"met":.328,
                "rna_protein":.330,"rna_met":.329,
                "prot_met":.328,"all":.325}
        eighty=choose_minimum_modalities(losses,MODES,fraction=.8)
        ninety=choose_minimum_modalities(losses,MODES,fraction=.9)
        self.assertEqual(eighty["arm"],"met")
        self.assertEqual(eighty["n_omic_modalities"],1)
        self.assertEqual(ninety["arm"],"all")
        self.assertEqual(ninety["n_omic_modalities"],3)
    def test_zero_gain_must_allow_no_biomarker(self):
        l={key:.35+len(modes)*.02 for key,modes in MODES.items()}
        for f in (.8,.9):
            c=choose_minimum_modalities(l,MODES,fraction=f)
            self.assertEqual(c["arm"],"tumor_only")
            self.assertEqual(c["status"],"NO_POSITIVE_INNER_CV_GAIN")
    def test_cost_weights_are_not_retrofitted_from_outer_test(self):
        losses={"tumor_only":.35,"rna":.33,"protein":.33,
                "met":.33,"rna_protein":.329,"rna_met":.329,
                "prot_met":.329,"all":.328}
        a=choose_minimum_modalities(losses,MODES,fraction=.8,
            source_assay_cost_weights={"rna":3.0,"prot":2.0,"met":1.0})
        self.assertEqual(a["arm"],"met")
        self.assertEqual(a["selected_weighted_assay_cost"],1)
    def test_mismatched_source_modality_schema_fails_closed(self):
        with self.assertRaises(ModalityBudgetError):
            choose_minimum_modalities({"tumor_only":.4},MODES,fraction=.8)
        with self.assertRaises(ModalityBudgetError):
            choose_minimum_modalities({k:.3 for k in MODES},
                   {**MODES,"all":("rna","rna")},fraction=.9)
if __name__=="__main__":
    unittest.main()
