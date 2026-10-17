"""Synthetic-only multiphenotype clustered bootstrap, provenance and multiplicity."""
import unittest

from omics2pk.epa_liver_cyp_validation import EXPECTED_GENES, EPALinkageError
from omics2pk.epa_eight_cyp_sensitivity import (
    calculate_eight_gene_sensitivity, _quantile, _point_loss_contrast,
)


def source():
    data={}
    for k,gene in enumerate(EXPECTED_GENES):
        records=[]
        for i in range(100):
            if k == 2 and i >= 94:
                continue
            key=(f"ORIG_FAKE{i:03d}",f"XT_FAKE{i:03d}")
            age=20.+(i*7)%60
            male=float(i%2)
            mrna=0.1+(i*13+k*17)%73/3.
            effect=0.5*mrna if k != 4 else 0.
            activity=max(.01,20.+age*.05+effect+(i%11)*1.6)
            records.append((key,mrna,activity,age,male))
        data[gene]=records
    return data


class ClusteredOriginalSourceTests(unittest.TestCase):
    def test_eight_endpoints_one_shared_biological_sample(self):
        s=source()
        result=calculate_eight_gene_sensitivity(s,n_resamples=110)
        self.assertEqual(result["original_biological_donors"],100)
        self.assertEqual(result["n_related_CYP_endpoints"],8)
        self.assertEqual(result["n_prespecified_shared_donor_bootstrap_draws"],110)
        self.assertEqual(result["results"]["CYP2C9"]["n_original_donors"],94)
        self.assertTrue(result["original_individual_donor_records_exported"] is False)
        self.assertNotIn("ORIG_FAKE",str(result))
        self.assertEqual(len(result["results"]),8)

    def test_wrong_gene_universe_rejected(self):
        s=source()
        del s["CYP2C9"]
        with self.assertRaisesRegex(EPALinkageError,"Eight prespecified"):
            calculate_eight_gene_sensitivity(s,n_resamples=100)

    def test_duplicate_original_donor_rejected(self):
        s=source()
        rows=s["CYP1A2"]
        rows[1]=(rows[0][0],*rows[1][1:])
        with self.assertRaisesRegex(EPALinkageError,"Duplicate"):
            calculate_eight_gene_sensitivity(s,n_resamples=100)

    def test_fewer_than_original_100_donors_rejected(self):
        s=source()
        for gene in s:
            s[gene]=s[gene][:-1]
        with self.assertRaisesRegex(EPALinkageError,"100 total"):
            calculate_eight_gene_sensitivity(s,n_resamples=100)

    def test_fixed_random_seed_and_adjusted_interval_wider(self):
        s=source()
        first=calculate_eight_gene_sensitivity(s,n_resamples=101)
        second=calculate_eight_gene_sensitivity(s,n_resamples=101)
        self.assertEqual(first,second)
        for g in EXPECTED_GENES:
            row=first["results"][g]
            ci95=row["paired_percentile_CI_95_unadjusted"]
            cibon=row["paired_percentile_CI_99_375_Bonferroni_8"]
            self.assertLessEqual(cibon[0],ci95[0])
            self.assertGreaterEqual(cibon[1],ci95[1])

    def test_percentile_is_monotone(self):
        self.assertAlmostEqual(_quantile([10,20,30,40],0.5),25.)
        self.assertLess(_quantile([3,1,8],.01),_quantile([3,1,8],.99))


if __name__=="__main__":
    unittest.main()
