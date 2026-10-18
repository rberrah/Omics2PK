"""No clinical/person-level data. Synthetic prospective PK-PD contract tests."""
from __future__ import annotations
import unittest

from omics2pk.staged_pkpd_contract import (
    OmicsPDLinkError,qualify_temporal_pd_link,qualify_omics2pk_or_pd_staging
)


def fictional():
    omics=[{"sample_id":f"FAKE_SAMPLE_{i}","subject_id":f"FICTITIOUS_{i}",
            "omic":"metabolomics","collected_at":"2025-01-01T08:00:00+00:00",
            "tissue":"plasma"} for i in range(13)]
    treatment=[{"subject_id":f"FICTITIOUS_{i}","intervention":
                "L_carnitine" if i%2 else "placebo",
                "started_at":"2025-01-02T08:00:00+00:00",
                "allocation":"randomized"} for i in range(13)]
    outcome=[{"outcome_event_id":f"FAKE_OUTCOME_{i}","subject_id":f"FICTITIOUS_{i}",
              "endpoint":"90d_mortality","observed_value":str(i%2),
              "assessed_at":"2025-04-01T08:00:00+00:00"} for i in range(13)]
    return omics,treatment,outcome


class MatchedClinicalPDStagingTests(unittest.TestCase):
    def test_original_patient_prospectivity_and_biological_source_only(self):
        a,b,c=fictional()
        x=qualify_omics2pk_or_pd_staging(
            direction="PD",original_omics_samples=a,
            original_doses_or_treatments=b,original_observations=c,
            settings={"endpoint":"90d_mortality"})
        self.assertEqual(x["n_original_unique_subjects_with_PD_and_omics"],13)
        self.assertTrue(x["observational_association_not_causal_individual_PD_treatment_effect"])
        self.assertNotIn("FICTITIOUS",str(x))
        self.assertNotIn("FAKE_SAMPLE",str(x))

    def test_no_future_omics_post_treatment(self):
        a,b,c=fictional()
        a[2]["collected_at"]="2025-02-01T09:00:00+00:00"
        with self.assertRaisesRegex(OmicsPDLinkError,"predose"):
            qualify_temporal_pd_link(a,b,c,endpoint="90d_mortality")

    def test_no_fake_second_omic_layer(self):
        a,b,c=fictional()
        with self.assertRaisesRegex(OmicsPDLinkError,"predose"):
            qualify_temporal_pd_link(a,b,c,endpoint="90d_mortality",
                                       required_modalities=("metabolomics","proteomics"))

    def test_no_row_order_join_allowed(self):
        a,b,c=fictional()
        c[2]["subject_id"]="OTHER_PATIENT"
        with self.assertRaisesRegex(OmicsPDLinkError,"same-person"):
            qualify_temporal_pd_link(a,b,c,endpoint="90d_mortality")

    def test_original_outcome_cannot_precede_treatment(self):
        a,b,c=fictional()
        c[4]["assessed_at"]="2024-12-01T08:00:00Z"
        with self.assertRaisesRegex(OmicsPDLinkError,"before"):
            qualify_temporal_pd_link(a,b,c,endpoint="90d_mortality")

    def test_same_patient_repeated_assay_ambiguous(self):
        a,b,c=fictional()
        copy=dict(a[2]);copy["sample_id"]="OTHER_FAKE_ID";a.append(copy)
        with self.assertRaisesRegex(OmicsPDLinkError,"unambiguous"):
            qualify_temporal_pd_link(a,b,c,endpoint="90d_mortality")


if __name__=="__main__":unittest.main()
