"""NCI-60 source-specific biochemical label quality gates (synthetic-only)."""
import unittest
from omics2pk.nci60_biochemical_coverage import status_of_original_label,audit_labels,PATHWAY_PANELS
class TestNCI60BiochemicalEvidence(unittest.TestCase):
    def test_unknown_not_assigned_pathway_by_name_similarity(self):
        self.assertEqual(status_of_original_label("X-1595 - possible glutathione metabolite"),"unidentified_X_feature")
        self.assertEqual(status_of_original_label("Isobar 3 includes glucose 6-phosphate"),"unresolved_isobar_group")
        self.assertEqual(status_of_original_label("Possible Artifact-2336 (RI Ladder)"),"possible_artifact_or_retired")
        self.assertEqual(status_of_original_label("glutathione, reduced"),"named_candidate_identity_unverified")
    def test_exact_source_labels_only_not_similar_isoform_identity(self):
        p=PATHWAY_PANELS["glutathione_redox"]
        self.assertIn("glutathione, reduced",p["metabolites"])
        self.assertNotIn("glutathione, oxidized (GSSG)",p["metabolites"])
        self.assertIn("glutathione, oxidized (GSSG)",p["not_measured_key_metabolites"])
    def test_fail_closed_on_incomplete_source(self):
        with self.assertRaises(ValueError):
            audit_labels(("glutathione, reduced","X-1234"))
    def test_source_shaped_synthetic_only(self):
        core=set()
        for p in PATHWAY_PANELS.values():core.update(p["metabolites"])
        core=list(core)
        self.assertLessEqual(len(core),111)
        named=core+["synthetic_named_probe_"+str(i) for i in range(111-len(core))]
        mystery=["X-"+str(i+1000) for i in range(161)]
        isobar=["Isobar "+str(i)+" includes A, B" for i in range(7)]
        labels=named+mystery+isobar+["Possible Artifact-2336 (RI Ladder)"]
        report=audit_labels(labels)
        self.assertEqual(report["source_total_features"],280)
        self.assertEqual(report["candidate_named_features_not_structurally_confirmed"],111)
        self.assertEqual(report["source_unknown_isobar_artifact_count"],169)
        self.assertFalse(report["source_has_HMDB_ChEBI_InChIKey_or_MSMS_confidence_in_this_metadata"])
if __name__=="__main__":unittest.main()
