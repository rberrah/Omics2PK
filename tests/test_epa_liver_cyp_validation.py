"""Synthetic-only multi-endpoint 100-liver alignment and anti-leakage tests."""
import unittest
from omics2pk.epa_liver_cyp_validation import (
 EXPECTED_GENES,EPALinkageError,read_paired_epa,validate_one_gene
)


def fake_workbook(n=100):
    acols=["KUMC Sample ID#","XenoTech Sample ID#","Gender","Age (years)",*EXPECTED_GENES]
    rcols=["KUMC Sample ID#","XenoTech Sample ID#","Gender","Age (years)",*EXPECTED_GENES,"CYP7A1"]
    activity=[]
    rna=[]
    for i in range(n):
        base=[f"KUMC_SYNTH_{i}",f"XT_SYNTH_{i}","Female" if i%2 else "Male",float(i%70+4)]
        x=[float(20+i*(g+1)/3) for g in range(len(EXPECTED_GENES))]
        y=[max(0.1,(1.2*v+(i%4)*2)) for v in x]
        activity.append(base+y)
        rna.append(base+x+[4])
    return {"Figure 1":(["ID"],[[str(i)] for i in range(n)]),
            "Figure 2":(acols,activity),
            "Table 3":(rcols,rna)}


class OriginalSourceEvidenceTests(unittest.TestCase):
    def test_exact_two_identifier_100_source_subjects(self):
        table,qc=read_paired_epa(fake_workbook())
        self.assertEqual(len(table),8)
        self.assertEqual(len(table["CYP3A4"]),100)
        self.assertEqual(qc["CYP2B6"]["n_complete_independent_donors"],100)

    def test_wrong_one_subject_crosswalk_is_rejected(self):
        wb=fake_workbook()
        wb["Table 3"][1][1][1]="XT_UNKNOWN"
        with self.assertRaisesRegex(EPALinkageError,"do not match"):
            read_paired_epa(wb)

    def test_repeat_subject_in_activity_data_is_rejected(self):
        wb=fake_workbook()
        wb["Figure 2"][1][1][:2]=wb["Figure 2"][1][0][:2]
        with self.assertRaisesRegex(EPALinkageError,"Duplicated"):
            read_paired_epa(wb)

    def test_fails_on_cross_modal_demographic_identity_conflict(self):
        wb=fake_workbook()
        wb["Figure 2"][1][1][3]=999
        with self.assertRaisesRegex(EPALinkageError,"contradiction"):
            read_paired_epa(wb)

    def test_heldout_model_shares_no_original_donor_rows(self):
        rows,_=read_paired_epa(fake_workbook())
        result=validate_one_gene(rows["CYP3A4"])
        self.assertEqual(result["n_independent_original_donors"],100)
        self.assertEqual(result["n_donor_heldout_folds"],100)
        self.assertNotIn("XT_SYNTH",repr(result))
        self.assertIn("age_sex_plus_rna",result["mae_original_units"])

    def test_at_least_sixty_independent_donors(self):
        rows,_=read_paired_epa(fake_workbook())
        with self.assertRaisesRegex(EPALinkageError,">=60"):
            validate_one_gene(rows["CYP3A4"][:59])

if __name__=="__main__":
    unittest.main()
