"""Synthetic-only original rosuvastatin 40-person local source contract tests.

No publisher original biological records or target values are included.
"""
from __future__ import annotations
import csv
import json
import tempfile
import unittest
from pathlib import Path

from omics2pk.rosuvastatin_original_local import (
    SourceGateError, original_csv, qualify_exact_subject_link, qualify_combined_source, summarize_local_sources
)


def _csv(path,headers,rows):
    with open(path,"w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


class RealSourceGateTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.omics=self.root/"omic.csv"
        self.pk=self.root/"pk.csv"
        self.people=[f"PRIVATE_ORIGINAL_HUMAN_{i:03d}" for i in range(40)]
        self.omicsrows=[{
            "original_subject":p,
            "metabolite_001":str(10.0+i/10),
            "metabolite_002":str(50.0+i/5),
            "metabolite_003":str(1000.0+i),
        } for i,p in enumerate(self.people)]
        self.pkrows=[{
            "human_subject":p,
            "rosu_auc_original":str(75.0+i),
            "rosu_cmax_original":str(9.0+i/4),
        } for i,p in enumerate(self.people)]
        self.write()

    def write(self):
        _csv(self.omics,["original_subject","metabolite_001","metabolite_002","metabolite_003"],
             self.omicsrows)
        _csv(self.pk,["human_subject","rosu_auc_original","rosu_cmax_original"],
             self.pkrows)

    def qualify(self, **kw):
        return qualify_exact_subject_link(self.omics,self.pk,
            omics_id_col="original_subject",pk_id_col="human_subject",
            auc_col="rosu_auc_original",cmax_col="rosu_cmax_original",**kw)

    def test_actual_40_matching_only_means_structural_proof(self):
        d=self.qualify()
        self.assertEqual(d["n_biological_subjects"],40)
        self.assertEqual(d["n_usable_quantitative_omic_features"],3)
        self.assertEqual(d["n_original_PK_auc_and_cmax_both_numeric"],40)
        self.assertFalse(d["model_fit_performed"])
        self.assertFalse(d["original_source_provenance_independently_authenticated"])
        self.assertFalse(d["predose_timing_independently_verified"])
        self.assertFalse(d["quantitative_PK_units_and_assay_independently_verified"])
        self.assertNotIn("PRIVATE_ORIGINAL_HUMAN",json.dumps(d))
        self.assertNotIn("75.0",json.dumps(d))
        self.assertNotIn("metabolite_001",json.dumps(d))

    def test_inventory_never_leaks_columns_or_rows(self):
        d=summarize_local_sources([self.omics,self.pk])
        self.assertEqual(d["n_original_source_files"],2)
        self.assertTrue(d["files"][1]["headers_suggest_auc"])
        self.assertTrue(d["files"][1]["headers_suggest_cmax"])
        self.assertFalse(d["qualified_original_PK_omic_person_join"])
        self.assertNotIn("human_subject",repr(d))
        self.assertNotIn("PRIVATE_ORIGINAL_HUMAN",repr(d))

    def test_mismatched_ids_must_not_join_by_row_order(self):
        self.omicsrows[0]["original_subject"]="DISTINCT_PRIVATE_PERSON"
        self.write()
        with self.assertRaisesRegex(SourceGateError,"exact strict 40-person"):
            self.qualify()

    def test_duplicate_measurement_not_an_independent_person(self):
        self.omicsrows[-1]["original_subject"]=self.omicsrows[0]["original_subject"]
        self.write()
        with self.assertRaisesRegex(SourceGateError,"Duplicate original omics"):
            self.qualify()

    def test_missing_PK_endpoints_hard_fail(self):
        self.pkrows[0]["rosu_auc_original"]="NA"
        self.write()
        with self.assertRaisesRegex(SourceGateError,"Missing original source quantitative AUC"):
            self.qualify()

    def test_cannot_train_with_target_in_features(self):
        for row in self.omicsrows:
            row["AUC"]=row["metabolite_001"]
        _csv(self.omics,["original_subject","metabolite_001","metabolite_002","metabolite_003","AUC"],
             self.omicsrows)
        with self.assertRaisesRegex(SourceGateError,"PK outcome or postdose"):
            self.qualify()

    def test_wrong_original_human_count_rejected(self):
        self.pkrows.pop()
        self.write()
        with self.assertRaisesRegex(SourceGateError,"biological-person count"):
            self.qualify()

    def test_source_data_annotated_but_not_numeric_rejected(self):
        for row in self.omicsrows:
            for key in ("metabolite_001","metabolite_002","metabolite_003"):
                row[key]=""
        self.write()
        with self.assertRaisesRegex(SourceGateError,"not contain enough"):
            self.qualify()

    def test_author_origins_can_join_only_via_explicit_one_to_one_mapping(self):
        for i,row in enumerate(self.omicsrows):
            row["original_subject"]=f"LAB_SAMPLE_{i:03d}"
        self.write()
        cross=self.root/"author_crosswalk.csv"
        _csv(cross,["omics_sample_id","pk_subject_id"],[
            {"omics_sample_id":f"LAB_SAMPLE_{i:03d}","pk_subject_id":p}
            for i,p in enumerate(self.people)])
        d=self.qualify(author_crosswalk_path=cross)
        self.assertEqual(d["n_biological_subjects"],40)
        self.assertTrue(d["source_SHA256_author_crosswalk_if_supplied"])
        with self.assertRaisesRegex(SourceGateError,"exact strict 40-person"):
            self.qualify()

    def test_original_authors_combined_csv_is_only_structural_evidence(self):
        original=self.root/"author_combined.csv"
        rows=[]
        for i,p in enumerate(self.people):
            rows.append({"subject":p,"AUC":self.pkrows[i]["rosu_auc_original"],
                         "Cmax":self.pkrows[i]["rosu_cmax_original"],
                         "molecular_feature_a":self.omicsrows[i]["metabolite_001"],
                         "molecular_feature_b":self.omicsrows[i]["metabolite_002"]})
        _csv(original,["subject","AUC","Cmax","molecular_feature_a","molecular_feature_b"],rows)
        got=qualify_combined_source(original,person_id_col="subject",auc_col="AUC",cmax_col="Cmax")
        self.assertEqual(got["n_unique_original_biological_humans"],40)
        self.assertEqual(got["n_candidate_numeric_features_excluding_explicit_PK_targets"],2)
        self.assertTrue(got["target_columns_explicitly_excluded_from_feature_candidates"])
        self.assertFalse(got["model_fit_or_predictive_performance_claim"])
        self.assertNotIn("PRIVATE_ORIGINAL_HUMAN",str(got))
        self.assertNotIn("molecular_feature_a",str(got))

    def test_combined_extra_PK_target_column_rejected(self):
        original=self.root/"author_combined.csv"
        rows=[{"subject":p,"AUC":str(100+i),"Cmax":str(10+i),
               "rosuvastatin_concentration_postdose":str(i+1),
               "molecular_a":str(i+2),"molecular_b":str(i+3)} for i,p in enumerate(self.people)]
        _csv(original,list(rows[0]),rows)
        with self.assertRaisesRegex(SourceGateError,"extra potential postdose PK"):
            qualify_combined_source(original,person_id_col="subject",auc_col="AUC",cmax_col="Cmax")

    def test_html_gdrive_is_not_data(self):
        self.omics.write_text("<!doctype html><body>Sign in to Google Drive</body>")
        with self.assertRaisesRegex(SourceGateError,"Google Drive returned HTML"):
            self.qualify()

    def test_bad_csv_has_no_implicit_inferred_column(self):
        with self.assertRaisesRegex(SourceGateError,"Exact original ID column missing"):
            qualify_exact_subject_link(self.omics,self.pk,
                omics_id_col="not_actually_present",pk_id_col="human_subject",
                auc_col="rosu_auc_original",cmax_col="rosu_cmax_original")

    def test_semi_colon_original_csv_supported(self):
        f=self.root/"semicolon.csv"
        f.write_text("original_subject;feature\nPRIVATE_A;1.5\nPRIVATE_B;2.5\n")
        h,r,_=original_csv(f)
        self.assertEqual(len(h),2)
        self.assertEqual(len(r),2)


if __name__=="__main__":
    unittest.main()
