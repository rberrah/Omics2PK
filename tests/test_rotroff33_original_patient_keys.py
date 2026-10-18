"""Original-clinical-person ID and A/B/C metabolomic sample gate, synthetic only."""
import unittest

from scripts.qualify_rotroff33_original_metformin_patient_keys import (
    _clinical_ids,_sample_index,qualify_person_key_overlap
)


def fixture(*, changed_one=False):
    omics=[{0:"Sample_ID"}]
    for i in range(33):
        s=f"SUBJ{i+1:04d}"
        for time in "ABC":
            omics.append({0:s+"-"+time,1:str(10+i),2:"metabolite"})
    for i in range(18): omics.append({0:f"QC_{i:02d}"})
    c1=[{0:"Original Subject ID",1:"Glucose AUC",2:"Metformin plasma Cmax",3:"Age"}]
    for i in range(33):
        v=f"{i+1:04d}" if not(changed_one and i==2) else "9999"
        c1.append({0:v,1:str(10+i),2:str(300+i),3:str(24+i)})
    return [omics,c1]


class SourceIndividualClinicalJoinTests(unittest.TestCase):
    def test_33_real_code_matching_rule_and_a_b_c(self):
        d=qualify_person_key_overlap(fixture())
        self.assertTrue(d["exact_33_patient_source_omics_clinical_person_ID_pairing_qualified"])
        self.assertEqual(d["n_unique_original_omics_ABC_triplets"],33)
        p=d["pairs"][0]
        self.assertEqual(p["n_verified_numeric_patient_code_matches"],33)
        self.assertTrue(p["one_to_one_original_33_person_code_coverage"])
        self.assertFalse(d["not_a_measured_clinical_PK_model"] is False)
        self.assertNotIn("SUBJ",str(d))

    def test_patient_code_not_joined_by_row_order(self):
        d=qualify_person_key_overlap(fixture(changed_one=True))
        self.assertFalse(d["exact_33_patient_source_omics_clinical_person_ID_pairing_qualified"])
        self.assertEqual(d["pairs"][0]["n_verified_numeric_patient_code_matches"],32)

    def test_qc_extra_sample_does_not_inflate_biological_patient_count(self):
        sample=_sample_index(fixture()[0])
        self.assertEqual(sample["n_complete_patient_triplets"],33)
        self.assertEqual(sample["study_subject_groups_original_match_count"],33)

    def test_duplicate_source_clinical_id_refused(self):
        c=fixture()[1]
        c[2][0]=c[1][0]
        with self.assertRaisesRegex(ValueError,"duplicated"):
            _clinical_ids(c)

    def test_malformed_source_ID_refused(self):
        c=fixture()[1]
        c[1][0]="19"
        self.assertIsNone(_clinical_ids(c))


if __name__=="__main__":unittest.main()
