"""Synthetic source-contract tests only; no invented data used as original proof."""
import unittest
from math import log

try:
    import numpy as np
    from scripts.benchmark_original_wegler15_target_vs_five_cyps import (
        _source_panel, _fit_heldout, benchmark_records,
        WeglerPanelAbstention, DRUG_TERMS, SYSTEMS)
except ModuleNotFoundError as exc:
    if exc.name not in {"numpy","sklearn","scipy","xgboost","pypdf"}:
        raise
    np=None
    _source_panel=_fit_heldout=benchmark_records=None


def make_original_shape_fixture(shuffle=False, n=15):
    rng=np.random.default_rng(20261105)
    result=[]
    for j in range(n):
        values=rng.uniform(.8,5.,size=len(DRUG_TERMS))
        for k,drug in enumerate(DRUG_TERMS):
            c=float(np.exp(.5 + .7*np.log(values[k])
                           +rng.normal(0,.08)))
            result.append({
                "donor_id":f"ARTIFICIAL_{j:03d}",
                "drug":drug,
                "cyp_hlm_pmol_per_mg":float(values[k]),
                "cyp_hh_pmol_per_mg":float(values[k]),
                "source_eq1":c*2,
                "source_eq4":c,
            })
    if shuffle:
        rng.shuffle(result)
    return result


@unittest.skipUnless(_fit_heldout is not None,"optional scientific numerical extras unavailable")
class WeglerTargetFiveOriginalShapeContracts(unittest.TestCase):
    def test_source_panel_index_donor_invariance(self):
        r=make_original_shape_fixture()
        ids,x,y,j=_source_panel(r,"HH","midazolam")
        ids2,x2,y2,j2=_source_panel(make_original_shape_fixture(shuffle=True),"HH","midazolam")
        self.assertEqual(ids,ids2)
        self.assertEqual(j,j2)
        np.testing.assert_equal(x,x2)
        np.testing.assert_equal(y,y2)
        d=_fit_heldout(x,y,j)
        self.assertEqual(d["n_independent_original_donors"],15)
        self.assertEqual(len(d["MAE_log_original_intrinsic_clearance"]),6)
        self.assertTrue(d["same_complete_case_donor_set_all_models"])
        self.assertNotIn("ARTIFICIAL_",str(d))
        self.assertEqual(d,_fit_heldout(x2,y2,j2))

    def test_missing_endpoint_is_complete_case_across_arms(self):
        rows=make_original_shape_fixture()
        rows[0]["source_eq4"]=None
        ids,x,y,j=_source_panel(rows,"HH",rows[0]["drug"])
        self.assertEqual(len(ids),14)
        self.assertEqual(_fit_heldout(x,y,j)["n_independent_original_donors"],14)
        rows=make_original_shape_fixture()
        rows[0]["cyp_hh_pmol_per_mg"]=None
        ids,x,y,j=_source_panel(rows,"HH","midazolam")
        self.assertEqual(len(ids),14)

    def test_invalid_source_identity_and_nonpositive_protein_rejected_or_excluded(self):
        rows=make_original_shape_fixture()
        with self.assertRaises(WeglerPanelAbstention):
            _source_panel(rows,"PLASMA","midazolam")
        with self.assertRaises(WeglerPanelAbstention):
            _source_panel(rows,"HH","not_drug")
        with self.assertRaises(WeglerPanelAbstention):
            _source_panel(rows[:-1],"HH","midazolam")
        with self.assertRaises(WeglerPanelAbstention):
            _source_panel(rows+rows[:1],"HH","midazolam")
        with self.assertRaises(WeglerPanelAbstention):
            _fit_heldout(np.zeros((15,4)),np.ones(15),0)
        with self.assertRaises(WeglerPanelAbstention):
            _fit_heldout(np.zeros((15,5)),np.ones(15),6)


if __name__=="__main__":
    unittest.main()
