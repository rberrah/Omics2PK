"""Synthetic-only tests; never pass fixtures as original human cohort evidence."""
import unittest
import numpy as np
from omics2pk.epa_liver_cyp_validation import EPALinkageError
from scripts.benchmark_original_epa_cyp_8rna_saturation import _metrics_cv, GENES


def synthetic_rows(seed=20261105):
    rng=np.random.default_rng(seed)
    n=73
    ages=rng.uniform(20,75,size=n)
    sexes=rng.integers(0,2,size=n)
    logrna=rng.normal(1.3,.5,size=(n,len(GENES)))
    # CYP3A4 is deliberately predictive of *synthetic* outcome.
    y=np.maximum(.01,np.expm1(2.0+.6*logrna[:,GENES.index("CYP3A4")]
                             +rng.normal(0,.1,size=n)))
    return [(f"FAKE_DONOR_{i:03d}",float(ages[i]),float(sexes[i]),
             float(y[i]),*[float(v) for v in logrna[i]])
            for i in range(n)]


class OriginalEPASaturationSyntheticContracts(unittest.TestCase):
    def test_all_arms_complete_and_not_original_PK(self):
        r=_metrics_cv(synthetic_rows(),"CYP3A4",repeats=1)
        self.assertEqual(r["n_original_independent_donors"],73)
        self.assertEqual(len(r["MAE_in_original_activity_units"]),8)
        self.assertEqual(len(r["exploratory_paired_person_contrasts"]),6)
        self.assertEqual(r["predictor_tissue"],"same human liver bulk qPCR")
        self.assertFalse(r["per_donor_values_and_IDs_emitted"])
        self.assertNotIn("FAKE_DONOR",str(r))
        self.assertEqual(r,_metrics_cv(synthetic_rows(),"CYP3A4",repeats=1))
        self.assertLess(r["MAE_in_original_activity_units"]["clinical_plus_target_RNA|alpha_1"],
                        r["MAE_in_original_activity_units"]["clinical|alpha_1"])

    def test_reordering_is_not_permitted_as_accidental_split_change(self):
        rows=synthetic_rows()
        result=_metrics_cv(rows,"CYP3A4",repeats=1)
        self.assertEqual(result,_metrics_cv(rows[::-1],"CYP3A4",repeats=1))

    def test_rejects_wrong_gene_and_invalid_source(self):
        rows=synthetic_rows()
        with self.assertRaises(EPALinkageError):
            _metrics_cv(rows,"CYP3A5",repeats=1)
        with self.assertRaises(EPALinkageError):
            _metrics_cv(rows,"CYP3A4",repeats=9)
        bad=list(rows)
        bad[0]=(bad[0][0],bad[0][1],bad[0][2],float("nan"),*bad[0][4:])
        with self.assertRaises(EPALinkageError):
            _metrics_cv(bad,"CYP3A4",repeats=1)
        bad=list(rows)
        bad[1]=(bad[0][0],*bad[1][1:])
        with self.assertRaises(EPALinkageError):
            _metrics_cv(bad,"CYP3A4",repeats=1)


if __name__=="__main__":
    unittest.main()
