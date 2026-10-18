"""Synthetic controls only: original EPA records NEVER included in public tests."""
from __future__ import annotations
from random import Random
import unittest

from omics2pk.epa_liver_cyp_validation import EXPECTED_GENES, EPALinkageError
from omics2pk.epa_rna_negative_control import (
    _permute_rna_within_sex, _check_records,
    evaluate_epa_permutation_negative_control,
)


def fixture():
    base=[]
    for i in range(100):
        age=25.0 + (i*3)%45
        sex=float(i%2)
        rna=0.75+float((i*7)%71)/5
        base.append(((f"FAKE_SOURCE_K_{i:03}",f"FAKE_SOURCE_X_{i:03}"),
                     rna,20.+rna*2.75+(i%7)*0.02,age,sex))
    return {gene:list(base) if j!=2 else list(base[:-6])
            for j,gene in enumerate(EXPECTED_GENES)}


class FakeDonorMatchingTests(unittest.TestCase):
    def test_rna_permutations_preserve_original_outcome_and_age_sex(self):
        rows=fixture()["CYP1A2"]
        perm=_permute_rna_within_sex(rows,Random(42))
        self.assertEqual([r[:1]+r[2:] for r in perm],[r[:1]+r[2:] for r in rows])
        for sex in (0.,1.):
            self.assertEqual(
                sorted(r[1] for r in rows if r[4]==sex),
                sorted(r[1] for r in perm if r[4]==sex),
            )
        self.assertNotEqual([r[1] for r in rows],[r[1] for r in perm])

    def test_real_source_heldout_permutation_contract(self):
        report=evaluate_epa_permutation_negative_control(
            fixture(),n_permutations=19,seed=20261018
        )
        self.assertEqual(report["n_unique_original_biological_donors"],100)
        self.assertEqual(report["n_prespecified_enzyme_endpoints"],8)
        self.assertTrue(report["not_patient_PK_or_PopPK"])
        self.assertEqual(report["endpoints"]["CYP2C9"]["n_unique_original_human_liver_donors"],94)
        self.assertLess(
            report["endpoints"]["CYP1A2"]["real_matched_RNA_delta_mae"],0
        )
        self.assertNotIn("FAKE_SOURCE",repr(report))
        for gene in EXPECTED_GENES:
            e=report["endpoints"][gene]
            self.assertGreaterEqual(e["one_sided_empirical_permutation_p"],0.05)
            self.assertLessEqual(e["one_sided_empirical_permutation_p"],1.0)
            self.assertFalse(e["source_rows_or_individual_predictions_exported"])

    def test_wrong_source_gene_count_rejected(self):
        d=fixture();d.pop("CYP1A2")
        with self.assertRaises(EPALinkageError):
            evaluate_epa_permutation_negative_control(d,n_permutations=19)

    def test_duplicate_donor_rejected(self):
        d=fixture()
        original=d["CYP2C19"][0]
        r=d["CYP2C19"][1]
        d["CYP2C19"][1]=(original[0],*r[1:])
        with self.assertRaisesRegex(EPALinkageError,"Duplicate"):
            evaluate_epa_permutation_negative_control(d,n_permutations=19)

    def test_invalid_permutations_rejected(self):
        for val in (True,18,400):
            with self.assertRaisesRegex(EPALinkageError,"19..399"):
                evaluate_epa_permutation_negative_control(fixture(),n_permutations=val)

    def test_sex_stratum_too_small_rejected(self):
        rows=fixture()["CYP1A2"]
        rows[:]=[(*row[:4],float(i<2)) for i,row in enumerate(rows)]
        with self.assertRaisesRegex(EPALinkageError,"per sex"):
            _check_records(rows)


if __name__=="__main__":
    unittest.main()
