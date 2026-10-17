"""Original human predose LC-MS metabolite assignment files: honest matrix gate."""
import unittest
from scripts.audit_rosuvastatin40_original_maf_intensity import analyze_matrix,parse_csv


class SourceMAFMatrixTests(unittest.TestCase):
    def test_original_exact_sample_columns_have_numeric_intensities(self):
        ids={f"PRIVATE_PATIENT_{k}" for k in range(15)}
        header="metabolite_identification\t"+" \t".join(sorted(ids)).replace(" \t","\t")
        lines=[header]+["feature_"+str(i)+"\t"+"\t".join(str(10*i+k+1) for k in range(15))
                        for i in range(6)]
        d=analyze_matrix(("\n".join(lines)+"\n").encode(),ids)
        self.assertEqual(d["n_assay_sample_columns_matched_exactly"],15)
        self.assertEqual(d["n_numeric_values_in_matched_columns"],90)
        self.assertTrue(d["has_original_processed_quantitative_sample_feature_matrix"])
        self.assertNotIn("PRIVATE_PATIENT_",repr(d))

    def test_metabolite_annotation_without_individual_abundance_is_not_matrix(self):
        orig=b"metabolite_identification\tmz\trt\tconfidence\nMETAB_X\t123.4\t1.1\tmedium\n"
        d=analyze_matrix(orig,{"PRIVATE_PATIENT_0"})
        self.assertFalse(d["has_original_processed_quantitative_sample_feature_matrix"])
        self.assertEqual(d["n_assay_sample_columns_matched_exactly"],0)
        self.assertNotIn("METAB_X",repr(d))

    def test_partial_label_not_falsely_joined(self):
        orig=b"name\tPRIVATE_PATIENT_1-technicalrep\nitem\t42\n"
        d=analyze_matrix(orig,{"PRIVATE_PATIENT_1"})
        self.assertFalse(d["has_original_processed_quantitative_sample_feature_matrix"])

    def test_nonnumeric_unmeasured_cells_rejected_as_quantitation(self):
        ids={f"PERSON_PRIVATE{i}" for i in range(15)}
        header="item\t"+"\t".join(sorted(ids))
        orig=("\n".join([header]+["f"+str(i)+"\t"+"\t".join("NA" for _ in range(15)) for i in range(7)])+"\n").encode()
        d=analyze_matrix(orig,ids)
        self.assertEqual(d["n_numeric_values_in_matched_columns"],0)
        self.assertFalse(d["has_original_processed_quantitative_sample_feature_matrix"])


if __name__=="__main__":
    unittest.main()
