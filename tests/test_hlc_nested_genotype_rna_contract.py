"""Artificial only, authorized HLC genotype/RNA nested source contract."""
import unittest
try:
    import numpy as np
    from scripts.benchmark_authorized_hlc_nested_pgx_rna import (
        _fold_fit,_source_repeated_nested_cv,HLCNestedOriginalError)
except ModuleNotFoundError as exc:
    if exc.name not in ("numpy","scipy","sklearn","pandas"):
        raise
    np=None
    _source_repeated_nested_cv=None


@unittest.skipUnless(_source_repeated_nested_cv is not None,
    "Optional source science packages missing")
class AuthorizedHLCNestedSourceContract(unittest.TestCase):
    def make_synthetic(self):
        rng=np.random.default_rng(20261106)
        n=61
        age=rng.uniform(20,80,n)
        sex=rng.integers(0,2,n)
        rs=rng.integers(0,3,n).astype(float)
        x1=rng.normal(size=(n,2))
        x2=rng.normal(size=(n,2))
        source_activity=np.maximum(.1,np.expm1(4+.7*x1[:,0]+rng.normal(0,.1,n)))
        clin=np.column_stack((age,sex,rs))
        return source_activity,{
            "clinical_PGx":clin,
            "clinical_PGx_plus_target_RNA":np.column_stack((clin,x1)),
            "clinical_PGx_plus_target_RNA_POR_ABCB1":
                np.column_stack((clin,x1,x2)),
        }

    def test_same_donor_fixed_split_reproducible_and_leakage_free(self):
        y,designs=self.make_synthetic()
        r=_source_repeated_nested_cv(y,designs,seeds=(25,))
        self.assertEqual(r,_source_repeated_nested_cv(y,designs,seeds=(25,)))
        self.assertLess(
            r["clinical_PGx_plus_target_RNA"]["observed_human_liver_CYP_activity_MAE"],
            r["clinical_PGx"]["observed_human_liver_CYP_activity_MAE"])
        self.assertTrue(r["clinical_PGx_plus_target_RNA"]["no_original_donor_values_returned"])
        self.assertNotIn("_losses_private",str(r))
        self.assertEqual(r["clinical_PGx"]["n_same_original_genotyped_livers"],61)
        self.assertEqual(sum(r["clinical_PGx"]["nested_4fold_alpha_selected_counts_across_3x5_outer_splits"].values()),5)

    def test_reject_holdout_leakage_and_bad_source(self):
        y,designs=self.make_synthetic()
        with self.assertRaises(HLCNestedOriginalError):
            _fold_fit(designs["clinical_PGx"],y,[0,1,2,3],[2,4],1.)
        with self.assertRaises(HLCNestedOriginalError):
            _source_repeated_nested_cv(y[:2],{k:v[:2] for k,v in designs.items()},seeds=(25,))
        with self.assertRaises(HLCNestedOriginalError):
            _source_repeated_nested_cv(y,{"clinical_PGx":designs["clinical_PGx"][:-1]},seeds=(25,))


if __name__=="__main__":
    unittest.main()
