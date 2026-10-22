"""Synthetic qualification for pre-dose longitudinal trough linkage."""
import copy
import unittest

from omics2pk.predose_trough_link_contract import (
    PredoseTroughLinkError,
    qualify_predose_trough_longitudinal_link,
)


def fixture():
    omics=[]
    doses=[]
    troughs=[]
    visits=("D7","M1","M3")
    for s in range(3):
        subject=f"SYNTH_{s:02d}"
        for j,visit in enumerate(visits):
            day=1+j*10
            omics.append({
                "sample_id":f"OM_{subject}_{visit}",
                "subject_id":subject,
                "visit":visit,
                "modality":"METABOLITE",
                "collected_at":f"2026-01-{day:02d}T07:45:00+01:00",
            })
            doses.extend([
                {
                    "dose_id":f"PREV_{subject}_{visit}",
                    "subject_id":subject,
                    "drug":"tacrolimus",
                    "administered_at":f"2026-01-{day-1 if day>1 else 1:02d}T20:00:00+01:00"
                        if day>1 else "2025-12-31T20:00:00+01:00",
                    "amount":2.0,
                    "dose_unit":"mg",
                    "route":"oral",
                },
                {
                    "dose_id":f"NEXT_{subject}_{visit}",
                    "subject_id":subject,
                    "drug":"tacrolimus",
                    "administered_at":f"2026-01-{day:02d}T08:00:00+01:00",
                    "amount":2.0,
                    "dose_unit":"mg",
                    "route":"oral",
                },
            ])
            troughs.append({
                "observation_id":f"C0_{subject}_{visit}",
                "subject_id":subject,
                "visit":visit,
                "drug":"tacrolimus",
                "previous_dose_id":f"PREV_{subject}_{visit}",
                "next_dose_id":f"NEXT_{subject}_{visit}",
                "sampled_at":f"2026-01-{day:02d}T07:50:00+01:00",
                "concentration":6.0+s+j/10,
                "concentration_unit":"ng/mL",
            })
    return omics,doses,troughs


class PredoseTroughLinkTests(unittest.TestCase):
    def test_three_visit_source_contract_and_prospective_pairs(self):
        out=qualify_predose_trough_longitudinal_link(
            *fixture(),drug="tacrolimus",visit_order=("D7","M1","M3"),
            max_omics_trough_gap_hours=1)
        self.assertEqual(out["n_unique_subjects_linked"],3)
        self.assertEqual(out["n_linked_same_visit_trough_omics_pairs"],9)
        self.assertEqual(out["linked_visit_counts"],{"D7":3,"M1":3,"M3":3})
        self.assertEqual(out["n_subjects_with_at_least_two_linked_visits"],3)
        self.assertEqual(out["n_earlier_omics_to_later_trough_pairs_available"],9)
        self.assertTrue(out["trough_bounded_by_exact_previous_and_next_dose_events"])
        self.assertTrue(out["same_visit_molecular_state_is_concurrent_exposure_sensing_not_prospective"])
        self.assertTrue(out["future_model_split_must_group_by_whole_subject"])
        self.assertFalse(out["individual_ids_times_doses_concentrations_exported"])
        self.assertNotIn("SYNTH_",str(out))

    def test_rejects_morning_dose_as_previous_dose(self):
        parts=list(fixture())
        parts[2][0]["previous_dose_id"]=parts[2][0]["next_dose_id"]
        with self.assertRaisesRegex(PredoseTroughLinkError,"must differ"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))

    def test_rejects_trough_after_next_morning_dose(self):
        parts=list(fixture())
        parts[2][0]["sampled_at"]="2026-01-01T08:05:00+01:00"
        with self.assertRaisesRegex(PredoseTroughLinkError,"strictly between"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))

    def test_rejects_unmatched_subject_or_unknown_dose(self):
        parts=list(fixture())
        parts[2][0]["subject_id"]="OTHER"
        with self.assertRaisesRegex(PredoseTroughLinkError,"disagree on subject"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))
        parts=list(fixture())
        parts[2][0]["previous_dose_id"]="UNKNOWN"
        with self.assertRaisesRegex(PredoseTroughLinkError,"unknown"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))

    def test_rejects_missing_or_ambiguous_same_visit_metabolomics(self):
        parts=list(fixture())
        parts[0]=parts[0][1:]
        with self.assertRaisesRegex(PredoseTroughLinkError,"exactly one"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))
        parts=list(fixture())
        duplicate=copy.deepcopy(parts[0][0])
        duplicate["sample_id"]="DUP"
        parts[0].append(duplicate)
        with self.assertRaisesRegex(PredoseTroughLinkError,"exactly one"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))

    def test_rejects_large_omics_trough_gap_and_mixed_units(self):
        parts=list(fixture())
        parts[0][0]["collected_at"]="2026-01-01T05:00:00+01:00"
        with self.assertRaisesRegex(PredoseTroughLinkError,"matching window"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"),
                max_omics_trough_gap_hours=1)
        parts=list(fixture())
        parts[2][0]["concentration_unit"]="ug/L"
        with self.assertRaisesRegex(PredoseTroughLinkError,"units"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))

    def test_partial_visits_are_allowed_without_inventing_327_complete_visits(self):
        omics,doses,troughs=fixture()
        troughs=[x for x in troughs if not (
            x["subject_id"]=="SYNTH_02" and x["visit"]=="M3")]
        out=qualify_predose_trough_longitudinal_link(
            omics,doses,troughs,drug="tacrolimus",
            visit_order=("D7","M1","M3"))
        self.assertEqual(out["n_linked_same_visit_trough_omics_pairs"],8)
        self.assertEqual(out["linked_visit_counts"]["M3"],2)
        self.assertEqual(out["n_earlier_omics_to_later_trough_pairs_available"],8)

    def test_timezone_and_positive_values_fail_closed(self):
        parts=list(fixture())
        parts[0][0]["collected_at"]="2026-01-01T07:45:00"
        with self.assertRaisesRegex(PredoseTroughLinkError,"time zone"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))
        parts=list(fixture())
        parts[2][0]["concentration"]=0
        with self.assertRaisesRegex(PredoseTroughLinkError,"positive"):
            qualify_predose_trough_longitudinal_link(
                *parts,drug="tacrolimus",visit_order=("D7","M1","M3"))


if __name__=="__main__":
    unittest.main()
