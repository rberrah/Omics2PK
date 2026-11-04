"""STRING original-source interpretation: synthetic-only network evidence tests."""
import unittest
from scripts.audit_nci60_string_shap_mechanistic_network import (
    ALL,TOP_SHAP_SOURCE_RNA,TOP_SHAP_SOURCE_PROTEIN,
    LITERATURE_DOX_DRUG_MECHANISM,parse_string_output,
)

class STRINGNetworkEvidence(unittest.TestCase):
    def test_only_original_source_gene_edges_not_transporter_function(self):
        tsv="preferredName_A\tpreferredName_B\tscore\nIL18\tLGALS3\t0.55\nABCB1\tIL18\t0.91\n"
        d=parse_string_output(tsv)
        self.assertEqual(d["n_edges"],2)
        self.assertIn("IL18",d["original_high_shap_to_drug_prior_bridge_gene_symbols"])
        self.assertTrue(all(not x["is_drug_PK_mechanism_validated"] for x in d["edges"]))
        self.assertTrue(all(not x["is_direct_physical_interaction_confirmed"] for x in d["edges"]))
    def test_no_mapped_neighbors_or_non_original_symbols(self):
        tsv="preferredName_A\tpreferredName_B\tscore\nABCB1\tNOT_IN_MODEL\t0.98\n"
        self.assertEqual(parse_string_output(tsv)["n_edges"],0)
    def test_empty_graph_differs_from_technical_error(self):
        d=parse_string_output("")
        self.assertEqual(d["status"],"NO_INTERACTIONS_AT_FIXED_SCORE_CUTOFF")
        self.assertEqual(d["n_edges"],0)
        with self.assertRaises(ValueError):
            parse_string_output("source\ttarget\tscore\nIL18\tLGALS3\t.8\n")
    def test_original_gene_and_drug_mechanism_groups_kept_distinct(self):
        self.assertEqual(len(ALL),17)
        self.assertEqual(len(TOP_SHAP_SOURCE_RNA),5)
        self.assertEqual(len(TOP_SHAP_SOURCE_PROTEIN),5)
        self.assertIn("TOP2A",LITERATURE_DOX_DRUG_MECHANISM)
if __name__=="__main__":unittest.main()
