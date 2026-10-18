"""No original patients/values in these entirely simulated assay examples."""
from __future__ import annotations
import csv
from io import StringIO
import importlib.util
import unittest

from omics2pk.jarmusch14_public_clinical_pk import (
    ENDPOINTS,JarmuschOriginalPKError,
    read_author_metadata,benchmark_original_author_Jarmusch14
)


def fictional_source():
    primary=["ATTRIBUTE_Subject","ATTRIBUTE_Study_Day","ATTRIBUTE_Sample_Type",
             "MS_SampleType","Enrollment_Subject_Age","ATTRIBUTE_sex",
             "Day0_Vitals_weight_kg","PharmGeno_CYP2C19",
             "PharmGeno_CYP3A5","PharmGeno_CYP3A4",
             "PharmKinetics_CYP3A4","PharmKinetics_CYP2C19",
             "PharmKinetics_CYP1A2"]
    fields=[*primary,*[f"UNUSED_PUBLISHED_SOURCE_COLUMN_{k}" for k in range(202)]]
    h=StringIO();wr=csv.writer(h,delimiter="\t")
    wr.writerow(fields)
    for i in range(14):
        for day in (1,9):
            for replicate in range(5):
                obs={
                    "ATTRIBUTE_Subject":f"FAKE_HUMAN_{i:02d}",
                    "ATTRIBUTE_Study_Day":str(day),
                    "ATTRIBUTE_Sample_Type":"blood_plasma",
                    "MS_SampleType":"sample",
                    "Enrollment_Subject_Age":str(23+i*2),
                    "ATTRIBUTE_sex":"male" if i%2 else "female",
                    "Day0_Vitals_weight_kg":str(64+i),
                    "PharmGeno_CYP2C19":("1_1","1_17","1_2")[i%3],
                    "PharmGeno_CYP3A5":("3_3","1_3")[i%2],
                    "PharmGeno_CYP3A4":"1_1",
                    "PharmKinetics_CYP3A4":str(5.+i*.31+(2 if day==9 else 0)),
                    "PharmKinetics_CYP2C19":str(0.8+i*.1+(.2 if day==9 else 0)),
                    "PharmKinetics_CYP1A2":str(17.+i*.4+(1 if day==9 else 0))
                }
                if i==13 and day==9:
                    obs["PharmKinetics_CYP3A4"]=""
                    obs["PharmKinetics_CYP1A2"]=""
                wr.writerow([obs.get(k,"") for k in fields])
    # Source analytical blank/QC rows have no volunteer ID or PK target.
    for _ in range(1423-14*2*5):
        wr.writerow(["" for _ in fields])
    return h.getvalue()


class JarmuschHumanMatchedOriginalContractTests(unittest.TestCase):
    def test_original_14_people_repeated_same_person_MS_records_not_independent(self):
        source=read_author_metadata(fictional_source())
        self.assertEqual(source["n_independent_humans"],14)
        self.assertEqual(source["n_original_MS_metadata_rows"],1423)

    @unittest.skipUnless(importlib.util.find_spec("numpy") and
                         importlib.util.find_spec("sklearn"),"Optional scientific packages")
    def test_exact_LOPO_three_real_PK_endpoint_comparators(self):
        result=benchmark_original_author_Jarmusch14(fictional_source())
        self.assertEqual(result["n_unique_original_volunteers"],14)
        self.assertEqual(result["endpoint_analysis"]["CYP3A_midazolam_oral_CLF"]["n_independent_donors_with_day1_and_day9_actual_NCA_PK"],13)
        self.assertEqual(result["endpoint_analysis"]["CYP2C19_omeprazole_metabolite_AUC_ratio"]["n_independent_donors_with_day1_and_day9_actual_NCA_PK"],14)
        self.assertIn("clinical_plus_CYP_genotype_plus_PRIOR_drug_probe_PK",
                      result["endpoint_analysis"]["CYP3A_midazolam_oral_CLF"]["MAE_comparison"])
        self.assertEqual(len(result["endpoint_analysis"]["CYP1A2_caffeine_oral_CLF"]["MAE_comparison"]),2)
        self.assertNotIn("FAKE_HUMAN",str(result))

    def test_replicated_conflicting_PK_source_must_fail(self):
        source=f is not None if False else fictional_source()
        lines=source.splitlines()
        # Change one of the replicated same human/day values, no biological
        # subject information included in public fixture.
        lines[2]=lines[2].replace("5.0\t0.8\t17.0","900\t0.8\t17.0")
        with self.assertRaises(JarmuschOriginalPKError):
            read_author_metadata("\n".join(lines))


if __name__=="__main__":
    unittest.main()
