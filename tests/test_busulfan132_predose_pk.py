"""ALL fake donors and plasma values; no original McCune patient source here."""
from __future__ import annotations
import importlib.util
import unittest

from omics2pk.busulfan132_predose_pk import (
    SOURCE_DOI,FEATURES,BusulfanSourceNotQualified,
    qualify_authorized_original_busulfan,audit_original_busulfan_study
)


def fake_original_shape(with_pgx=False):
    # Strictly fabricated 132 patients and 96 repeated premolecular samples.
    rows=[]
    for i in range(132):
        identity=f"FICTIONAL_ORIGINAL_BUSULFAN_{i:03d}"
        # target determined by target-independent chemical synthetic signature
        cl=2.0+(i%19)*.07
        clinical={"age_years":str(22+i%65),
                  "normal_fat_mass_kg":str(30+i%38),
                  "sex":"male" if i%2 else "female"}
        if with_pgx:
            clinical.update({"PGx_original_assay_verified":True,
                             "GSTA1_original_function_class":i%3,
                             "GSTM1_original_copy_number":i%2})
        for period in (("immediate_predose","2026-01-17T10:00:00+00:00"),
                       ("two_week_predose","2026-01-04T09:00:00+00:00")):
            if period[0]=="two_week_predose" and i>=96:continue
            vals=[round(2+j*.01+(i%13)*.02+(j%3)*(.2 if period[0]=="immediate_predose" else .12),4)
                  for j in range(FEATURES)]
            rows.append(dict(source_doi=SOURCE_DOI,subject_id=identity,
                             sample_id=f"FAKE_SAMPLE_{i:03d}_{period[0]}",
                             time_window=period[0],collected_at=period[1],
                             first_busulfan_dose_at="2026-01-18T10:00:00+00:00",
                             busulfan_CL_per_NFM=str(cl),
                             CL_unit="mL/min/kg_NFM",
                             metabolite_abundances=vals,**clinical))
    return rows


class OriginalBusulfanMethodSafetyTests(unittest.TestCase):
    def test_132_original_source_structural_patient_identity_not_proof_data_is_available(self):
        fake=fake_original_shape()
        x=audit_original_busulfan_study(fake)
        self.assertEqual(x["n_original_unique_biological_humans"],132)
        self.assertEqual(x["n_original_two_week_samples"],96)
        self.assertEqual(x["n_original_pre_dose_biological_samples"],228)
        self.assertFalse(x["genuine_source_doc_authorized_and_independently_reviewed"])
        self.assertFalse(x["original_genotype_available"])
        self.assertNotIn("FICTIONAL_ORIGINAL",str(x))

    def test_same_donor_two_periods_cannot_have_two_clearances(self):
        fake=fake_original_shape()
        fake[1]["busulfan_CL_per_NFM"]="1000"
        with self.assertRaisesRegex(BusulfanSourceNotQualified,"Inconsistent"):
            qualify_authorized_original_busulfan(fake)

    def test_no_postdose_original_metabolomic_assay_masquerade(self):
        fake=fake_original_shape()
        fake[0]["collected_at"]="2026-01-20T10:00:00+00:00"
        with self.assertRaisesRegex(BusulfanSourceNotQualified,"PRE-dose"):
            qualify_authorized_original_busulfan(fake)

    def test_no_original_shortcut_from_missing_841_source_values(self):
        fake=fake_original_shape()
        fake[0]["metabolite_abundances"]=fake[0]["metabolite_abundances"][:840]
        with self.assertRaisesRegex(BusulfanSourceNotQualified,"841"):
            qualify_authorized_original_busulfan(fake)

    def test_optional_true_study_PGx_needs_every_observation(self):
        fake=fake_original_shape(with_pgx=True)
        del fake[-1]["PGx_original_assay_verified"]
        with self.assertRaisesRegex(BusulfanSourceNotQualified,"provenance"):
            qualify_authorized_original_busulfan(fake)

    def test_source_cannot_use_unmatched_liver_ID_approximations(self):
        fake=fake_original_shape()
        fake[0]["subject_id"]="INCORRECT_OTHER_CLINICAL_HUMAN"
        with self.assertRaises(BusulfanSourceNotQualified):
            qualify_authorized_original_busulfan(fake)

    @unittest.skipUnless(importlib.util.find_spec("numpy") and
                         importlib.util.find_spec("sklearn"),
                         "Optional scientific numpy/scikit-learn")
    def test_real_method_fits_only_on_other_people_not_the_other_timepoint(self):
        from omics2pk.busulfan132_predose_pk import benchmark_authorized_original_busulfan
        x=benchmark_authorized_original_busulfan(fake_original_shape())
        self.assertEqual(x["n_true_person_heldout_two_week_assays"],96)
        self.assertTrue(x["each_two_week_predose_heldout_person_excluded_from_all_training_timepoints"])
        self.assertFalse(x["source_original_specific_two_loci_PGx_present"])
        self.assertFalse(x["genuine_source_documents_and_crosswalk_independently_reviewed"])
        self.assertNotIn("FICTIONAL_ORIGINAL",str(x))
        self.assertEqual(len(x["results"]),2)


if __name__=="__main__":unittest.main()
