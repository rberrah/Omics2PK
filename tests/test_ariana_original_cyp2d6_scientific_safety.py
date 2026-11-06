"""No genotype leakage into pure clinical/null arms; no drug target-derived ratio."""
from __future__ import annotations
import unittest
import numpy as np

from scripts.benchmark_ariana_original_cyp2d6_solanidine_vs_pgx_20261106 import (
    _fit_predict,ARMS
)

class ARIANAScientificSourceSafety(unittest.TestCase):
    def source_shape_fixture(self):
        rng=np.random.default_rng(20261106)
        rows=[]
        for i in range(30):
            rows.append({
                "plasma_metoprolol_ng_ml":float(np.exp(
                    4+.2*(i%3)+rng.normal(0,.12))),
                "age":55+i*.3,"weight":74+(i%5),
                "log_daily_dose":float(np.log(95 if i%2 else 190)),
                "activity_score":float(i%3),
                "ln_oh_solanidine_ratio":float(rng.normal()),
                "ln_ssda_solanidine_ratio":float(rng.normal()),
                "sex":"male" if i%2 else "female",
                "CYP2D6_inhibitor":"1" if i%5==0 else "0",
                "copy_number_variation":"2" if i%3 else "3",
            })
        return rows

    def test_null_model_must_not_see_copy_number(self):
        rows=self.source_shape_fixture()
        modified=[dict(r) for r in rows]
        for j,r in enumerate(modified):
            r["copy_number_variation"]=str((j%4)+4)
        train=np.arange(28)
        test=np.array([28,29])
        for arm in ("clinical_only","clinical_plus_OH_solanidine_without_PGx"):
            before=_fit_predict(rows,ARMS[arm],train,test,alpha=1.)
            after=_fit_predict(modified,ARMS[arm],train,test,alpha=1.)
            np.testing.assert_array_equal(before,after)

    def test_features_prespecified_without_target_drug_metabolite(self):
        for name,feature_set in ARMS.items():
            self.assertNotIn("ln_oh_metoprolol_ratio",feature_set)
            self.assertNotIn("plasma_metoprolol_ng_ml",feature_set)
            self.assertNotIn("weighted_concentration",feature_set)
            if name.startswith("clinical_PGx"):
                self.assertIn("activity_score",feature_set)

    def test_never_use_overlap_or_invalid_shrinkage(self):
        rows=self.source_shape_fixture()
        with self.assertRaises(Exception):
            _fit_predict(rows,ARMS["clinical_PGx_CYP2D6_AS"],
                         np.arange(20),np.array([19,20]),1.)


if __name__=="__main__":
    unittest.main()
