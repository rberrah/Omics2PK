"""Synthetic tests for the authorized HLC -> inductive SWEET bridge."""
import unittest

try:
    import numpy as np
    from omics2pk.hlc_authorized_sweet_bridge import (
        prepare_authorized_hlc_sweet_input,
        source_center_from_hlc_donor_id,
    )
    from omics2pk.hlc_real_cyp_validation import OriginalHLCQualificationError
except ModuleNotFoundError as exc:
    if exc.name != "numpy":
        raise
    np = None


@unittest.skipIf(np is None, "NumPy optional extra not installed")
class HLCAuthorizedSweetBridgeTests(unittest.TestCase):
    def test_public_center_codebook_is_strict(self):
        self.assertEqual(source_center_from_hlc_donor_id("1110005"), "Vanderbilt")
        self.assertEqual(source_center_from_hlc_donor_id("2220048"), "Pittsburgh")
        self.assertEqual(source_center_from_hlc_donor_id("3330053"), "Merck")
        with self.assertRaisesRegex(OriginalHLCQualificationError, "prefix"):
            source_center_from_hlc_donor_id("4440001")
        with self.assertRaisesRegex(OriginalHLCQualificationError, "seven-digit"):
            source_center_from_hlc_donor_id("111005")

    def test_prepare_reuses_exact_source_ids_and_all_target_probes(self):
        rng=np.random.default_rng(20261022)
        n_each=105
        ids=tuple(
            [f"111{i:04d}" for i in range(n_each)]
            + [f"222{i:04d}" for i in range(n_each)]
            + [f"333{i:04d}" for i in range(n_each)]
        )
        n=len(ids)
        features=tuple(f"probe_{i:03d}" for i in range(24))
        clinical=np.column_stack((
            rng.uniform(18,85,n),
            np.asarray([i%2 for i in range(n)],float),
        ))
        rna=rng.normal(size=(n,len(features)))
        y=np.exp(1+0.3*rna[:,2]-0.15*rna[:,3]+0.01*clinical[:,0])
        parsed={
            "sample_ids":list(ids),
            "clinical":clinical,
            "rna":rna,
            "feature_ids":list(features),
            "gene_rows":{"CYP3A4":np.asarray([2,3,4],int)},
            "outcomes":{"CYP3A4M_activity":y},
        }
        out=prepare_authorized_hlc_sweet_input(
            parsed,activity_key="CYP3A4M_activity")
        self.assertEqual(out["n_complete"],315)
        self.assertEqual(out["expression"].shape,(24,315))
        self.assertEqual(out["direct_cyp_features"],(
            "probe_002","probe_003","probe_004"))
        self.assertEqual(out["center_counts"],{
            "Merck":105,"Pittsburgh":105,"Vanderbilt":105})
        self.assertEqual(out["activity_unit"],"nmol/min/mg protein")

    def test_prepare_rejects_unknown_center_and_nonfinite_rna(self):
        rng=np.random.default_rng(88)
        ids=tuple(
            [f"111{i:04d}" for i in range(101)]
            + [f"222{i:04d}" for i in range(101)]
            + [f"333{i:04d}" for i in range(101)]
        )
        features=tuple(f"p{i}" for i in range(10))
        clinical=np.column_stack((np.full(303,50.0),np.arange(303)%2))
        rna=rng.normal(size=(303,10))
        base={
            "sample_ids":list(ids),
            "clinical":clinical,
            "rna":rna,
            "feature_ids":list(features),
            "gene_rows":{"CYP3A4":np.asarray([0,1],int)},
            "outcomes":{"CYP3A4M_activity":np.ones(303)},
        }
        bad=dict(base)
        bad["sample_ids"]=list(ids)
        bad["sample_ids"][0]="4440001"
        with self.assertRaisesRegex(OriginalHLCQualificationError,"prefix"):
            prepare_authorized_hlc_sweet_input(
                bad,activity_key="CYP3A4M_activity")
        bad2=dict(base)
        bad2["rna"]=rna.copy()
        bad2["rna"][0,0]=np.nan
        with self.assertRaisesRegex(OriginalHLCQualificationError,"missing/nonfinite"):
            prepare_authorized_hlc_sweet_input(
                bad2,activity_key="CYP3A4M_activity")


if __name__=="__main__":
    unittest.main()
