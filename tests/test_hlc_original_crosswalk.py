"""Fictional-only 427 HLC samples and observed-enzyme label crosswalk gates."""
from __future__ import annotations
import unittest

from omics2pk.hlc_original_crosswalk import (
    audit_explicit_hlc_crosswalk,HumanLiverCrosswalkError
)

SHA="a"*64


def fake_source():
    gsm=[f"GSM{i+100000}" for i in range(427)]
    code=[{"gsm_id":g,"donor_id":f"FICTIONAL_LIVER_{i:03d}",
           "codebook_sha256":SHA} for i,g in enumerate(gsm[:420])]
    pheno=[{"donor_id":f"FICTIONAL_LIVER_{i:03d}",
            "assay_id":"CYP3A4_midazolam_hydroxylation",
            "measured_value":str(i*3+1),
            "assay_units":"pmol/min/mg",
            "original_entity_id":"syn4629"} for i in range(9,398)]
    return gsm,code,pheno


class HLCMatchedDonorSyntheticTests(unittest.TestCase):
    def test_exact_matched_pheno_and_GSM_source_only(self):
        ids,book,pheno=fake_source()
        r=audit_explicit_hlc_crosswalk(
            ids,book,pheno,assay_id="CYP3A4_midazolam_hydroxylation"
        )
        self.assertEqual(r["source_original_GSM_donor_matches_with_assay"],389)
        self.assertEqual(r["source_crosswalk_unique_original_GSM_donor_pairs"],420)
        self.assertFalse(r["qualified_independent_biological_donor_numeric_model"])
        self.assertFalse(r["original_source_codebook_proven_authentic_by_this_function"])
        self.assertNotIn("FICTIONAL_LIVER",str(r))
        self.assertNotIn("GSM10000",str(r))

    def test_nonunique_biological_donor_crosswalk_not_allowed(self):
        ids,book,pheno=fake_source()
        book[3]["donor_id"]=book[2]["donor_id"]
        with self.assertRaisesRegex(HumanLiverCrosswalkError,"Repeated"):
            audit_explicit_hlc_crosswalk(ids,book,pheno,assay_id="CYP3A4_midazolam_hydroxylation")

    def test_phenotype_relabeling_must_be_rejected(self):
        ids,book,pheno=fake_source()
        pheno[4]["original_entity_id"]="unrelated_study"
        with self.assertRaisesRegex(HumanLiverCrosswalkError,"syn4629"):
            audit_explicit_hlc_crosswalk(ids,book,pheno,assay_id="CYP3A4_midazolam_hydroxylation")

    def test_physical_unit_mismatch_rejected(self):
        ids,book,pheno=fake_source()
        pheno[3]["assay_units"]="mg/L"
        with self.assertRaisesRegex(HumanLiverCrosswalkError,"units inconsistent"):
            audit_explicit_hlc_crosswalk(ids,book,pheno,assay_id="CYP3A4_midazolam_hydroxylation")

    def test_unsourced_crosswalk_digest_rejected(self):
        ids,book,pheno=fake_source()
        book[0]["codebook_sha256"]="unsourced"
        with self.assertRaisesRegex(HumanLiverCrosswalkError,"hash"):
            audit_explicit_hlc_crosswalk(ids,book,pheno,assay_id="CYP3A4_midazolam_hydroxylation")

    def test_original_398_activity_n_is_not_equal_427_expression_n_by_magic(self):
        ids,book,pheno=fake_source()
        r=audit_explicit_hlc_crosswalk(ids,book,pheno,assay_id="CYP3A4_midazolam_hydroxylation")
        self.assertNotEqual(r["source_original_GSM_donor_matches_with_assay"],427)


if __name__=="__main__":
    unittest.main()
