"""Nested omic gain-retention design: synthetic-only and no data leakage tests."""
import unittest
from omics2pk.information_budget import (
    InformationBudgetError, select_gain_retention_budget,
    cumulative_shap_coverage, pareto_frontier, validate_budgets
)

class InformationBudgetSemantics(unittest.TestCase):
    def test_90pct_gain_preservation_not_90pct_shap(self):
        # Clinical loss=10. Best internal molecular MAE=6 (gain 4).
        # k=3 retains gain 3.2 -> 80%; k=8 retains >=3.6 ->90%.
        cv={0:10.0,1:9.0,2:7.5,3:6.8,5:6.5,8:6.3,12:6.0}
        r=select_gain_retention_budget(cv,.8)
        self.assertEqual(r["k"],3)
        r=select_gain_retention_budget(cv,.9)
        self.assertEqual(r["k"],8)
        self.assertEqual(r["best_positive_inner_gain"],4.0)
    def test_no_positive_gain_selects_no_omics(self):
        for cv in ({0:10,1:11,3:11.5,10:10.1},{0:10.0,1:10,5:10}):
            for fraction in (.8,.9):
                selected=select_gain_retention_budget(cv,fraction)
                self.assertEqual(selected["k"],0)
                self.assertEqual(selected["status"],"NO_POSITIVE_INNER_VALIDATED_GAIN")
    def test_shap_90pct_not_verified_pk_gain(self):
        info=cumulative_shap_coverage([9.,1.,0.],[0])
        self.assertAlmostEqual(info["fraction_of_absolute_SHAP"],.9)
        self.assertIsNone(info["fraction_of_predictive_gain"])
        self.assertIn("NOT_PREDICTIVE_GAIN",info["status"])
    def test_pareto_baseline_can_dominate_every_biomarker(self):
        curve={"0":{"k":0,"loss":10},"1":{"k":1,"loss":11},
               "3":{"k":3,"loss":10.1}}
        frontier=pareto_frontier(curve)["non_dominated_models"]
        self.assertEqual([item["model"] for item in frontier],["0"])
    def test_pareto_decreasing_loss_has_multiple_cost_benefit_points(self):
        curve={"0":{"k":0,"loss":10},"1":{"k":1,"loss":8},
               "2":{"k":2,"loss":9},"5":{"k":5,"loss":7}}
        self.assertEqual([item["model"] for item in pareto_frontier(curve)["non_dominated_models"]],
                         ["0","1","5"])
    def test_wrong_budget_schedule_rejected(self):
        with self.assertRaises(InformationBudgetError):
            validate_budgets([0,1,10],9)
        with self.assertRaises(InformationBudgetError):
            select_gain_retention_budget({1:3},.8)
        with self.assertRaises(InformationBudgetError):
            cumulative_shap_coverage([1,2],[0,0])

if __name__=="__main__":
    unittest.main()
