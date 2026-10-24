"""Synthetic-only tests for BIOMARGIN -> ABIS aggregate feasibility."""
import unittest

from omics2pk.biomargin_abis_feasibility import (
    BiomarginAbisFeasibilityError,
    qualify_biomargin_abis_aggregate_feasibility,
)


def counts():
    return {
        "biomargin_molecular_patients":131,
        "exact_biomargin_abis_crosswalk_patients":70,
        "tacrolimus_tdm_patients":60,
        "tacrolimus_timed_dose_concentration_patients":48,
        "tacrolimus_preoutcome_omics_patients":35,
        "tacrolimus_cyp3a5_patients":28,
        "tacrolimus_repeated_omics_patients":12,
        "mpa_tdm_patients":45,
        "mpa_timed_dose_concentration_patients":30,
        "mpa_preoutcome_omics_patients":20,
        "independent_site_or_time_holdout_patients":18,
    }


def attestations():
    return {
        "deterministic_same_person_crosswalk_confirmed":True,
        "authorized_linkage_and_research_use_confirmed":True,
        "omics_collection_time_available":True,
        "no_direct_identifiers_in_aggregate_exchange":True,
    }


class BiomarginAbisFeasibilityTests(unittest.TestCase):
    def test_positive_aggregate_feasibility_without_patient_data(self):
        out=qualify_biomargin_abis_aggregate_feasibility(
            counts(),attestations())
        self.assertEqual(
            out["status"],
            "PAIRED_OMICS_TDM_FEASIBILITY_CONFIRMED_AT_AGGREGATE_LEVEL")
        self.assertTrue(out["tacrolimus_timed_pairing_exists"])
        self.assertTrue(out["tacrolimus_preoutcome_omics_pairing_exists"])
        self.assertTrue(out["tacrolimus_prospective_longitudinal_candidate_exists"])
        self.assertTrue(out["mpa_timed_pairing_exists"])
        self.assertTrue(out["cyp3a5_available_with_tacrolimus_timed_pk"])
        self.assertFalse(out["patient_level_data_requested_by_this_gate"])
        self.assertFalse(out["biomargin_itself_contains_pk"])
        self.assertFalse(out["abis_crosswalk_publicly_proven_by_this_gate"])

    def test_zero_crosswalk_returns_clean_stop(self):
        c=counts()
        for key in c:
            if key!="biomargin_molecular_patients":
                c[key]=0
        a=attestations()
        a["deterministic_same_person_crosswalk_confirmed"]=False
        a["authorized_linkage_and_research_use_confirmed"]=False
        a["omics_collection_time_available"]=False
        out=qualify_biomargin_abis_aggregate_feasibility(c,a)
        self.assertEqual(out["status"],"NO_BIOMARGIN_ABIS_PATIENT_CROSSWALK")
        self.assertFalse(out["tacrolimus_timed_pairing_exists"])
        self.assertFalse(out["mpa_timed_pairing_exists"])

    def test_nested_counts_fail_closed(self):
        c=counts()
        c["tacrolimus_timed_dose_concentration_patients"]=61
        with self.assertRaisesRegex(
            BiomarginAbisFeasibilityError,"Tacrolimus"):
            qualify_biomargin_abis_aggregate_feasibility(c,attestations())

        c=counts()
        c["mpa_preoutcome_omics_patients"]=31
        with self.assertRaisesRegex(BiomarginAbisFeasibilityError,"MPA"):
            qualify_biomargin_abis_aggregate_feasibility(c,attestations())

    def test_positive_match_requires_authorization_and_deterministic_link(self):
        a=attestations()
        a["deterministic_same_person_crosswalk_confirmed"]=False
        with self.assertRaisesRegex(
            BiomarginAbisFeasibilityError,"deterministic"):
            qualify_biomargin_abis_aggregate_feasibility(counts(),a)

        a=attestations()
        a["authorized_linkage_and_research_use_confirmed"]=False
        with self.assertRaisesRegex(
            BiomarginAbisFeasibilityError,"authorization"):
            qualify_biomargin_abis_aggregate_feasibility(counts(),a)

    def test_preoutcome_counts_require_collection_timing(self):
        a=attestations()
        a["omics_collection_time_available"]=False
        with self.assertRaisesRegex(
            BiomarginAbisFeasibilityError,"collection timing"):
            qualify_biomargin_abis_aggregate_feasibility(counts(),a)

    def test_direct_identifiers_never_allowed_in_aggregate_exchange(self):
        a=attestations()
        a["no_direct_identifiers_in_aggregate_exchange"]=False
        with self.assertRaisesRegex(
            BiomarginAbisFeasibilityError,"exclude direct identifiers"):
            qualify_biomargin_abis_aggregate_feasibility(counts(),a)

    def test_crosswalk_without_timed_pk_is_not_omics_pk_ready(self):
        c=counts()
        for key in (
            "tacrolimus_tdm_patients",
            "tacrolimus_timed_dose_concentration_patients",
            "tacrolimus_preoutcome_omics_patients",
            "tacrolimus_cyp3a5_patients",
            "tacrolimus_repeated_omics_patients",
            "mpa_tdm_patients",
            "mpa_timed_dose_concentration_patients",
            "mpa_preoutcome_omics_patients",
        ):
            c[key]=0
        out=qualify_biomargin_abis_aggregate_feasibility(c,attestations())
        self.assertEqual(
            out["status"],"CROSSWALK_EXISTS_BUT_NO_TIMED_DRUG_PK_PAIRING")
        self.assertFalse(out["tacrolimus_timed_pairing_exists"])

    def test_timed_tdm_without_preoutcome_omics_remains_intermediate(self):
        c=counts()
        c["tacrolimus_preoutcome_omics_patients"]=0
        c["tacrolimus_repeated_omics_patients"]=0
        c["mpa_preoutcome_omics_patients"]=0
        out=qualify_biomargin_abis_aggregate_feasibility(c,attestations())
        self.assertEqual(
            out["status"],
            "TIMED_TDM_PAIRING_EXISTS_BUT_PREOUTCOME_OMICS_NOT_CONFIRMED")

    def test_schema_is_frozen(self):
        c=counts(); c["unexpected"]=1
        with self.assertRaisesRegex(BiomarginAbisFeasibilityError,"frozen"):
            qualify_biomargin_abis_aggregate_feasibility(c,attestations())
        a=attestations(); a["unexpected"]=True
        with self.assertRaisesRegex(BiomarginAbisFeasibilityError,"frozen"):
            qualify_biomargin_abis_aggregate_feasibility(counts(),a)


if __name__=="__main__":
    unittest.main()
