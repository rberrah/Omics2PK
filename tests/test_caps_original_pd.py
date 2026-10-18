"""Fictional-only original CaPS RACE source contract: no actual human observations."""
from __future__ import annotations

from io import StringIO
import csv
import importlib.util
import math
import unittest


def fake_caps_files(n=215):
    clinical=[("studyid","dose","carntx","survivaldays","Age","Sex","Calculated.T0.SOFA.Value")]
    nmr=[("studyid","t0nmronly",*[f"MET_NMR_{i:02d}" for i in range(27)])]
    ac=[("timepoint","studyid",*[f"AC_LCMS_{i:02d}" for i in range(24)])]
    for i in range(n):
        s=f"FICTIONAL{i:04d}"
        dose=(0,6,12,18)[i%4]
        died=(i*23)%43<21
        age=25+(i*13)%58
        sofa=5+(i*7)%15
        clinical.append((s,dose,"Placebo" if dose==0 else "CarnTx",
                         35 if died else 365,age,"Male" if i%2 else "Female",sofa))
        nmr.append((s,1,*[round(2.+(i*j%17)/3.,4) for j in range(1,28)]))
        ac.append(("T0",s,*[round(1.+(i*(j+3)%19)/7.,4) for j in range(24)]))
    def serial(rows):
        h=StringIO()
        writer=csv.writer(h);writer.writerows(rows)
        return h.getvalue()
    return {"RACE_Demographics.csv":serial(clinical),
            "raw_nmr.csv":serial(nmr),
            "raw_ac.csv":serial(ac)}


class PublicCapsStudyPureContractTests(unittest.TestCase):
    @unittest.skipUnless(importlib.util.find_spec("numpy") and
                         importlib.util.find_spec("sklearn"),
                         "Research-only numpy/sklearn dependencies not installed")
    def test_215_same_original_fake_people_and_same_cv_folds(self):
        from omics2pk.caps_original_pd import load_original_caps,evaluate_original_caps
        dat=load_original_caps(files=fake_caps_files())
        self.assertEqual(dat["source_counts"]["n_true_same_original_people_both_omics_and_outcome"],215)
        result=evaluate_original_caps(dat)
        self.assertEqual(result["n_original_independent_people"],215)
        self.assertTrue(result["same_original_people_and_outer_split_all_models"])
        self.assertTrue(result["predictive_PD_outcome_not_individual_causal_treatment_benefit"])
        self.assertTrue(result["not_transcriptome_proteome_multiomics"])
        self.assertEqual(set(result["metrics"]),{
            "clinical_only","clinical_plus_NMR","clinical_plus_LCMS",
            "clinical_plus_both","clinical_plus_both_treatment_interactions"})
        self.assertNotIn("FICTIONAL",str(result))

    @unittest.skipUnless(importlib.util.find_spec("numpy"),"numpy optional")
    def test_post_treatment_omics_rejected(self):
        from omics2pk.caps_original_pd import load_original_caps,OriginalPDSourceError
        parts=fake_caps_files()
        parts["raw_ac.csv"]=parts["raw_ac.csv"].replace("T0,FICTIONAL0002","T12,FICTIONAL0002")
        with self.assertRaisesRegex(OriginalPDSourceError,"predose"):
            load_original_caps(files=parts)

    @unittest.skipUnless(importlib.util.find_spec("numpy"),"numpy optional")
    def test_mismatched_randomized_treatment_and_dose_rejected(self):
        from omics2pk.caps_original_pd import load_original_caps,OriginalPDSourceError
        parts=fake_caps_files()
        parts["RACE_Demographics.csv"]=parts["RACE_Demographics.csv"].replace(
            "FICTIONAL0000,0,Placebo,","FICTIONAL0000,0,CarnTx,")
        with self.assertRaisesRegex(OriginalPDSourceError,"dose"):
            load_original_caps(files=parts)

    @unittest.skipUnless(importlib.util.find_spec("numpy"),"numpy optional")
    def test_original_source_single_person_duplicate_never_pseudoreplicated(self):
        from omics2pk.caps_original_pd import load_original_caps,OriginalPDSourceError
        parts=fake_caps_files()
        lines=parts["raw_nmr.csv"].splitlines()
        lines.append(lines[1])
        parts["raw_nmr.csv"]="\n".join(lines)+"\n"
        with self.assertRaisesRegex(OriginalPDSourceError,"Duplicate"):
            load_original_caps(files=parts)


if __name__=="__main__":unittest.main()
