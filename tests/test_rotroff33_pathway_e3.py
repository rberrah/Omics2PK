"""Scientific schema guards, synthetic only. No source patient data printed."""
import unittest
import numpy as np
from omics2pk.rotroff33_pathway_e3 import (
   HYPOTHESES,PLANS,benchmark_one_target,qualify_original_feature_labels,E3SourceError
)
class SourcePathwayGates(unittest.TestCase):
    def test_names_are_partial_not_automated_KEGG_graph(self):
        self.assertEqual(HYPOTHESES['urea_cycle_partial'],('citrulline','ornithine','urea'))
        self.assertEqual(HYPOTHESES['purine_degradation_partial'],('inosine','hypoxanthine'))
        self.assertEqual(len(PLANS),7)
    def test_missing_original_binaries_fail_closed(self):
        with self.assertRaises((E3SourceError,ValueError)):
            qualify_original_feature_labels(b'fake')
    def test_do_not_invent_human_endpoint_or_original_384_feature_rows(self):
        rows=[{
          "clinical":[30.,1.,70.],"predose_omics":[.1]*384,
          "targets":{"metformin_AUClast_0_24_original_units":10.}
        } for _ in range(31)]
        annotation={"named":list(range(155)),"unknown":list(range(155,384)),
         "paths":{"urea_cycle_partial":[0,1,2],"purine_degradation_partial":[3,4]}}
        with self.assertRaises(E3SourceError):
            benchmark_one_target(rows,"metformin_AUClast_0_24_original_units",annotation)
if __name__=="__main__":unittest.main()
