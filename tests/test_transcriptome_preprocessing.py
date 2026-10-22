"""Synthetic tests for train-only transcriptome preprocessing."""
import unittest

try:
    import numpy as np
    from omics2pk.transcriptome_preprocessing import (
        MODE_COUNTS_CPM_LOG1P_ZSCORE,
        MODE_PASSTHROUGH,
        MODE_TRAIN_ZSCORE,
        TranscriptomePreprocessingError,
        fit_transcriptome_preprocessor,
        preprocess_outer_split,
        transform_transcriptome,
    )
except ModuleNotFoundError as exc:
    if exc.name!="numpy":
        raise
    np=None


@unittest.skipIf(np is None,"NumPy optional extra not installed")
class TranscriptomePreprocessingTests(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(20261022)
        self.genes=tuple(f"G{i}" for i in range(10))
        self.ids=tuple(f"S{i:02}" for i in range(14))
        self.train=self.ids[:10]
        self.test=self.ids[10:]
        self.x=rng.normal(size=(10,14))

    def split(self,x=None,mode=MODE_TRAIN_ZSCORE):
        return preprocess_outer_split(
            self.x if x is None else x,
            genes=self.genes,subject_ids=self.ids,
            train_ids=self.train,test_ids=self.test,mode=mode)

    def test_test_values_cannot_change_training_preprocessor_or_training_transform(self):
        a=self.split()
        changed=self.x.copy()
        changed[:,10:]+=np.linspace(-100,100,10)[:,None]
        b=self.split(changed)
        self.assertEqual(a["preprocessing_fingerprint"],b["preprocessing_fingerprint"])
        np.testing.assert_array_equal(a["training"],b["training"])
        self.assertFalse(a["uses_test_subjects_to_fit_preprocessing"])
        self.assertFalse(a["pk_pd_target_used"])
        self.assertFalse(a["batch_correction_executed"])

    def test_train_zscore_uses_training_statistics_and_train_median_imputation(self):
        x=self.x.copy()
        x[2,1]=np.nan
        x[2,11]=np.nan
        out=self.split(x)
        train=out["training"]
        self.assertAlmostEqual(float(train[2].mean()),0.0,places=12)
        self.assertAlmostEqual(float(train[2].std()),1.0,places=12)
        frozen=out["frozen"]
        expected=(frozen.medians[2]-frozen.means[2])/frozen.scales[2]
        self.assertAlmostEqual(float(out["test"][2,1]),float(expected),places=12)

    def test_passthrough_is_exact_and_rejects_missing(self):
        out=self.split(mode=MODE_PASSTHROUGH)
        np.testing.assert_array_equal(out["training"],self.x[:,:10])
        np.testing.assert_array_equal(out["test"],self.x[:,10:])
        bad=self.x.copy(); bad[0,12]=np.nan
        with self.assertRaises(TranscriptomePreprocessingError):
            self.split(bad,mode=MODE_PASSTHROUGH)

    def test_raw_counts_use_per_sample_cpm_before_train_zscore(self):
        rng=np.random.default_rng(77)
        counts=rng.poisson(lam=50,size=(10,14)).astype(float)+1
        a=preprocess_outer_split(
            counts,genes=self.genes,subject_ids=self.ids,
            train_ids=self.train,test_ids=self.test,
            mode=MODE_COUNTS_CPM_LOG1P_ZSCORE)
        scaled=counts.copy()
        scaled[:,12]*=7
        b=preprocess_outer_split(
            scaled,genes=self.genes,subject_ids=self.ids,
            train_ids=self.train,test_ids=self.test,
            mode=MODE_COUNTS_CPM_LOG1P_ZSCORE)
        self.assertEqual(a["preprocessing_fingerprint"],b["preprocessing_fingerprint"])
        np.testing.assert_allclose(a["test"][:,2],b["test"][:,2],rtol=0,atol=1e-12)
        np.testing.assert_array_equal(a["training"],b["training"])

    def test_constant_training_gene_fails_closed(self):
        x=self.x.copy(); x[0,:10]=2.0
        with self.assertRaisesRegex(TranscriptomePreprocessingError,"Constant"):
            self.split(x)

    def test_gene_order_and_subject_partition_are_strict(self):
        frozen=fit_transcriptome_preprocessor(
            self.x[:,:10],genes=self.genes,subject_ids=self.train,
            mode=MODE_TRAIN_ZSCORE)
        with self.assertRaisesRegex(TranscriptomePreprocessingError,"Gene identity"):
            transform_transcriptome(
                frozen,self.x[:,10:],genes=tuple(reversed(self.genes)))
        with self.assertRaisesRegex(TranscriptomePreprocessingError,"partition"):
            preprocess_outer_split(
                self.x,genes=self.genes,subject_ids=self.ids,
                train_ids=self.train,test_ids=self.test[:-1],
                mode=MODE_TRAIN_ZSCORE)

    def test_raw_counts_reject_negative_and_empty_library(self):
        counts=np.ones((10,14))
        counts[0,0]=-1
        with self.assertRaisesRegex(TranscriptomePreprocessingError,"nonnegative"):
            self.split(counts,mode=MODE_COUNTS_CPM_LOG1P_ZSCORE)
        counts=np.ones((10,14)); counts[:,0]=0
        with self.assertRaisesRegex(TranscriptomePreprocessingError,"library"):
            self.split(counts,mode=MODE_COUNTS_CPM_LOG1P_ZSCORE)


if __name__=="__main__":
    unittest.main()
