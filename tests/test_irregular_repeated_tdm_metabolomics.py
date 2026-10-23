"""Synthetic tests for irregular repeated pre-dose TDM/metabolomics linkage."""
import copy
import unittest

from omics2pk.irregular_repeated_tdm_metabolomics import (
    RepeatedTDMSourceError,
    qualify_irregular_repeated_tdm_metabolomics,
)


def fixture():
    omics=[]; doses=[]; troughs=[]
    # Irregular 3-5 visits per subject, intentionally not D7/M1/M3 labels.
    schedule={
        "SYN_A":[12,41,96,150],
        "SYN_B":[18,55,101],
        "SYN_C":[9,33,67,120,184],
    }
    for s_idx,(subject,days) in enumerate(schedule.items()):
        for j,day in enumerate(days):
            date_day=(j*3)+2
            biospecimen=f"BIO_{subject}_{j}"
            omics.append({
                "sample_id":f"OM_{subject}_{j}",
                "biospecimen_id":biospecimen,
                "subject_id":subject,
                "modality":"METABOLITE",
                "collected_at":f"2026-02-{date_day:02d}T07:45:00+01:00",
            })
            doses.extend([
                {
                    "dose_id":f"PREV_{subject}_{j}",
                    "subject_id":subject,
                    "drug":"tacrolimus",
                    "administered_at":f"2026-02-{date_day-1:02d}T20:00:00+01:00",
                    "amount":2.0,
                    "dose_unit":"mg",
                    "route":"oral",
                },
                {
                    "dose_id":f"NEXT_{subject}_{j}",
                    "subject_id":subject,
                    "drug":"tacrolimus",
                    "administered_at":f"2026-02-{date_day:02d}T08:00:00+01:00",
                    "amount":2.0,
                    "dose_unit":"mg",
                    "route":"oral",
                },
            ])
            troughs.append({
                "observation_id":f"C0_{subject}_{j}",
                "biospecimen_id":biospecimen,
                "subject_id":subject,
                "drug":"tacrolimus",
                "previous_dose_id":f"PREV_{subject}_{j}",
                "next_dose_id":f"NEXT_{subject}_{j}",
                "sampled_at":f"2026-02-{date_day:02d}T07:50:00+01:00",
                "concentration":6.0+s_idx+j/10,
                "concentration_unit":"ng/mL",
                "daily_dose":4.0,
                "daily_dose_unit":"mg/day",
                "weight_kg":60+s_idx*5,
                "postoperative_day":day,
            })
    return omics,doses,troughs


class IrregularRepeatedTDMTests(unittest.TestCase):
    def test_irregular_series_and_prospective_pair_count(self):
        out=qualify_irregular_repeated_tdm_metabolomics(*fixture())
        self.assertEqual(out["n_unique_subjects_linked"],3)
        self.assertEqual(out["n_linked_same_biospecimen_c0_omics_pairs"],12)
        self.assertEqual(out["n_subjects_with_at_least_two_linked_samples"],3)
        # 4C2 + 3C2 + 5C2 = 6 + 3 + 10
        self.assertEqual(out["n_earlier_omics_to_later_c0_pairs_available"],19)
        self.assertEqual(out["samples_per_subject_distribution"],{"3":1,"4":1,"5":1})
        self.assertFalse(out["fixed_nominal_visit_schedule_assumed"])
        self.assertTrue(out["daily_dose_and_weight_available_for_c0d_reconstruction"])
        self.assertTrue(out["future_model_split_must_group_by_whole_subject"])
        self.assertNotIn("SYN_A",str(out))

    def test_exact_biospecimen_link_required(self):
        parts=list(fixture())
        parts[2][0]["biospecimen_id"]="UNMATCHED"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"same-biospecimen"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)

    def test_biospecimen_subject_mismatch_rejected(self):
        parts=list(fixture())
        parts[0][0]["subject_id"]="OTHER"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"disagree on subject"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)

    def test_trough_must_be_between_exact_doses(self):
        parts=list(fixture())
        parts[2][0]["sampled_at"]="2026-02-02T08:05:00+01:00"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"strictly between"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)

    def test_missing_daily_dose_or_weight_rejected(self):
        parts=list(fixture())
        parts[2][0]["daily_dose"]=0
        with self.assertRaisesRegex(RepeatedTDMSourceError,"daily_dose"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)
        parts=list(fixture())
        parts[2][0]["weight_kg"]=None
        with self.assertRaisesRegex(RepeatedTDMSourceError,"weight_kg"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)

    def test_matching_window_and_units_fail_closed(self):
        parts=list(fixture())
        parts[0][0]["collected_at"]="2026-02-02T05:00:00+01:00"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"matching window"):
            qualify_irregular_repeated_tdm_metabolomics(*parts,max_omics_trough_gap_hours=1)
        parts=list(fixture())
        parts[2][0]["concentration_unit"]="ug/L"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"concentration units"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)
        parts=list(fixture())
        parts[2][0]["daily_dose_unit"]="mg"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"Daily dose units"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)

    def test_duplicate_chronology_rejected(self):
        parts=list(fixture())
        # Force a second event for same subject to the same trough timestamp while
        # keeping a distinct biospecimen and dose pair.
        subject="SYN_A"
        idx=[i for i,x in enumerate(parts[2]) if x["subject_id"]==subject][1]
        parts[2][idx]["sampled_at"]="2026-02-02T07:50:00+01:00"
        # Also move its molecular sample and bounding doses so only chronology
        # ambiguity, not the pre-dose gate, triggers.
        bio=parts[2][idx]["biospecimen_id"]
        for row in parts[0]:
            if row["biospecimen_id"]==bio:
                row["collected_at"]="2026-02-02T07:45:00+01:00"
        prev=parts[2][idx]["previous_dose_id"]; nxt=parts[2][idx]["next_dose_id"]
        for row in parts[1]:
            if row["dose_id"]==prev:
                row["administered_at"]="2026-02-01T20:00:00+01:00"
            if row["dose_id"]==nxt:
                row["administered_at"]="2026-02-02T08:00:00+01:00"
        with self.assertRaisesRegex(RepeatedTDMSourceError,"chronology ambiguous"):
            qualify_irregular_repeated_tdm_metabolomics(*parts)

    def test_non_target_drug_rows_do_not_create_links(self):
        parts=list(fixture())
        other=copy.deepcopy(parts[2][0])
        other["observation_id"]="CYCLO"
        other["drug"]="cyclosporine"
        parts[2].append(other)
        out=qualify_irregular_repeated_tdm_metabolomics(*parts)
        self.assertEqual(out["n_linked_same_biospecimen_c0_omics_pairs"],12)


if __name__=="__main__":
    unittest.main()
