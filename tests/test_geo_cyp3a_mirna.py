"""Synthetic GEO-style SOFT fixtures; genuine sample values never in git."""
from __future__ import annotations

import gzip
from random import Random
import unittest

from omics2pk.geo_cyp3a_mirna import (
    GEO49012Error, FEATURES, FULL_PANEL, compare_geo_mirna_vs_cyp3a_activity,
    parse_original_geo_soft,
)


def synthetic_soft(*, wrong_series=False, wrong_pmid=False,
                   remove_206=False, reverse=False, paper_consistent=False):
    rng=Random(19)
    rows=list(range(34))
    if reverse:
        rows=list(reversed(rows))
    txt=[
        "^SERIES = GSE00000" if wrong_series else "^SERIES = GSE49012",
        "!Series_pubmed_id = 10000000" if wrong_pmid else "!Series_pubmed_id = 24058572",
    ]
    for i in rows:
        group="DISEASE" if i < 22 else "NORMAL"
        independent_state=rng.random() + i*.005
        y=(.16 if group=="DISEASE" else .42) + independent_state*.1
        if paper_consistent:
            y=.25 if group=="DISEASE" else .47
        mir27=100*(1+independent_state)
        txt.extend([
            f"^SAMPLE = GSM{1192100+i}",
            f"!Sample_title = RV-SYN-{i:03d}",
            "!Sample_platform_id = GPL17470",
            "!Sample_source_name_ch1 = Liver Tissue",
            "!Sample_organism_ch1 = Homo sapiens",
            f"!Sample_characteristics_ch1 = status: {group}",
            f"!Sample_characteristics_ch1 = disease state: {('NASH' if i < 22 else 'NORMAL')}",
            f"!Sample_characteristics_ch1 = hepatic cyp3a activity: {y:.8f}",
            "!sample_table_begin",
            "ID_REF\tVALUE",
        ])
        signals={
            "hsa-miR-27b":mir27,
            "hsa-miR-206":500+70*independent_state,
            "hsa-miR-21":310+i*.4,
            "hsa-miR-130a":240+i*.7,
            "hsa-miR-142":110+i*.6,
            "hsa-miR-103":125+2*rng.random(),
            "hsa-miR-107":105+2*rng.random(),
            "hsa-miR-142-5p":105+i*.6,
            "hsa-miR-142-3p":160+i*.5,
            "hsa-miR-27b*":60+i*.6,
            "hsa-miR-21*":70+i*.7,
            "hsa-miR-130a*":50+i*.8,
        }
        if remove_206:
            del signals["hsa-miR-206"]
        for k,v in signals.items():
            txt.append(f"{k}\t{v:.9f}")
        for j in range(405):
            txt.append(f"SYNTHETIC_FEATURE_{j:03d}\t{(10+rng.random()):.6f}")
        for j in range(65):
            txt.append(f"hsa-miR-SYNTEST_{j:03d}\t{(10+rng.random()):.6f}")
        txt.append("!sample_table_end")
    return "\n".join(txt).encode("utf8")


class ActualGeoContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original=parse_original_geo_soft(gzip.compress(synthetic_soft()))

    def test_parser_recovers_independent_source_activity_and_34_samples(self):
        rows=self.original
        self.assertEqual(len(rows),34)
        self.assertEqual(sum(x["status"]=="DISEASE" for x in rows),22)
        self.assertEqual(sum(x["status"]=="NORMAL" for x in rows),12)
        self.assertEqual(len({x["accession"] for x in rows}),34)
        self.assertIn("hsa-miR-27b",rows[0]["intensity"])

    def test_cyp3a_activity_not_cyp3a5_or_clinical_pk(self):
        r=compare_geo_mirna_vs_cyp3a_activity(self.original)
        self.assertTrue(r["not_CYP3A5_specific"])
        self.assertTrue(r["no_CYP3A5_genotype"])
        self.assertTrue(r["no_human_drug_PK"])
        self.assertFalse(r["independent_center_or_external_cohort_validation"])
        self.assertEqual(r["n_original_unique_liver_biological_samples"],34)
        self.assertEqual(r["n_outer_sample_heldout_folds"],5)
        self.assertNotIn("GSM1192",str(r))
        self.assertNotIn("RV-SYN",str(r))

    def test_frozen_literature_feature_arms_and_stratified_metrics(self):
        r=compare_geo_mirna_vs_cyp3a_activity(self.original)
        self.assertIn("mir27b_206_external_literature",r["miRNA_feature_sets_used"])
        self.assertEqual(r["miRNA_feature_sets_used"]["mir27b_206_external_literature"],
                         list(FEATURES["mir27b_206_external_literature"]))
        self.assertIn("disease_status_only",r["heldout_MAE_by_arm_source_activity_units"])
        self.assertIn("log_geometric_train_mean",r["heldout_MAE_by_arm_source_activity_units"])
        self.assertEqual(set(r["heldout_MAE_by_arm_and_disease_status"]["disease_status_only"]),
                         {"NORMAL","DISEASE"})

    def test_permuted_original_row_order_does_not_change_heldout_predictions(self):
        before=compare_geo_mirna_vs_cyp3a_activity(self.original)
        after=compare_geo_mirna_vs_cyp3a_activity(list(reversed(self.original)))
        self.assertEqual(before,after)

    def test_wrong_series_or_pmid_refused(self):
        with self.assertRaisesRegex(GEO49012Error,"Wrong GEO accession"):
            parse_original_geo_soft(synthetic_soft(wrong_series=True))
        with self.assertRaisesRegex(GEO49012Error,"PMID"):
            parse_original_geo_soft(synthetic_soft(wrong_pmid=True))

    def test_original_missing_outcome_rejected(self):
        raw=synthetic_soft().replace(b"hepatic cyp3a activity: 0.",b"hepatic CYP missing: 0.",1)
        with self.assertRaisesRegex(GEO49012Error,"measured CYP3A activity metadata absent"):
            parse_original_geo_soft(raw)

    def test_duplicate_original_sample_rejected(self):
        raw=synthetic_soft().replace(b"^SAMPLE = GSM1192101",b"^SAMPLE = GSM1192100")
        with self.assertRaisesRegex(GEO49012Error,"Duplicate"):
            parse_original_geo_soft(raw)

    def test_missing_prespecified_mirna_abstains_not_zero_imputes(self):
        rows=parse_original_geo_soft(synthetic_soft(remove_206=True))
        r=compare_geo_mirna_vs_cyp3a_activity(rows)
        self.assertEqual(r["miRNA_feature_arms_abstained"]["mir27b_206_external_literature"],
                         "PRESPECIFIED_MIRNA_FEATURE_ABSENT")
        self.assertNotIn("mir27b_206_external_literature",
                         r["heldout_MAE_by_arm_source_activity_units"])
        self.assertIn("mir27b_only",r["heldout_MAE_by_arm_source_activity_units"])

    def test_no_synthetic_source_can_claim_specific_or_in_vivo_pk(self):
        r=compare_geo_mirna_vs_cyp3a_activity(self.original)
        self.assertEqual(r["status"],"ABSTAIN_GEO_ACTIVITY_PUBLICATION_GROUP_PARITY_FAILURE")
        self.assertFalse(r["source_publication_parity_audit"]["source_endpoint_paper_parity_qualified"])
        self.assertFalse(r["source_publication_parity_audit"]["biological_superiority_inference_permitted"])
        self.assertEqual(r["n_original_individual_rows_or_predictions_exported"],0)
        self.assertTrue(r["source_upstream_normalization_not_train_frozen"])

    def test_possible_publication_aggregate_parity_is_not_validation(self):
        records=parse_original_geo_soft(synthetic_soft(paper_consistent=True))
        r=compare_geo_mirna_vs_cyp3a_activity(records)
        self.assertEqual(r["status"],"RESEARCH_E2_SUPPLIED_GEO_SOURCE_PARITY_POSSIBLE")
        audit=r["source_publication_parity_audit"]
        self.assertTrue(audit["source_endpoint_paper_parity_qualified"])
        self.assertFalse(audit["biological_superiority_inference_permitted"])
        self.assertTrue(audit["same_endpoint_and_original_donor_subset_assumed_for_parity_check"])
        self.assertEqual(r["n_original_individual_rows_or_predictions_exported"],0)

    def test_impossible_positive_disease_subset_forces_abstention(self):
        records=parse_original_geo_soft(synthetic_soft(paper_consistent=True))
        for sample in records:
            if sample["status"]=="DISEASE":
                sample["activity"]=0.6
        r=compare_geo_mirna_vs_cyp3a_activity(records)
        audit=r["source_publication_parity_audit"]
        self.assertEqual(r["status"],"ABSTAIN_GEO_ACTIVITY_PUBLICATION_GROUP_PARITY_FAILURE")
        self.assertFalse(audit["disease_subset_of_paper_nonnegative_activity_possible"])
        self.assertTrue(audit["paper_control_mean_consistent"])
        self.assertGreater(audit["GEO_cirrhosis_activity_sum"],
                           audit["maximum_compatible_source_total_with_mean_rounding"])

    def test_mature_and_star_strands_are_never_aliased(self):
        r=compare_geo_mirna_vs_cyp3a_activity(self.original)
        self.assertEqual(r["miRNA_feature_sets_used"]["mir142_5p_external_literature"],
                         ["hsa-miR-142-5p"])
        self.assertEqual(r["miRNA_feature_sets_used"]["mir27b_star_distinct_strand_control"],
                         ["hsa-miR-27b*"])
        self.assertNotEqual(r["miRNA_feature_sets_used"]["mir27b_only"],
                            r["miRNA_feature_sets_used"]["mir27b_star_distinct_strand_control"])
        self.assertTrue(r["mature_and_star_strands_not_interchangeable"])

    def test_biological_mirna_PCA_is_train_frozen_or_explicitly_abstains(self):
        r=compare_geo_mirna_vs_cyp3a_activity(self.original)
        self.assertGreaterEqual(r["source_n_biological_miRNA_feature_IDs"],50)
        if "disease_plus_all_miRNA_PCA2" in r["miRNA_feature_sets_used"]:
            self.assertIn("disease_plus_all_miRNA_PCA4",
                          r["heldout_MAE_by_arm_source_activity_units"])
            self.assertTrue(r["PCA_imputation_scaling_and_components_fit_within_outer_training_only"])
        else:
            self.assertEqual(r["miRNA_feature_arms_abstained"]["disease_plus_all_miRNA_PCA2"],
                             "OPTIONAL_NUMPY_SKLEARN_DEPENDENCIES_UNAVAILABLE")

    def test_wrong_tissue_or_source_platform_refused(self):
        with self.assertRaisesRegex(GEO49012Error,"tissue"):
            parse_original_geo_soft(synthetic_soft().replace(b"Liver Tissue",b"Muscle Tissue",1))
        with self.assertRaisesRegex(GEO49012Error,"platform"):
            parse_original_geo_soft(synthetic_soft().replace(b"GPL17470",b"GPL00000",1))


if __name__=="__main__":
    unittest.main()
