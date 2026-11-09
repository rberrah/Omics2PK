"""Contract tests for naming the destroyed-signal replacement.

The diagnostic's number supports three incompatible readings until the
replacement is named: a non-human duplicate of the destroyed transcript, an
admissible paralogue that plausibly measures the same transcript, or a different
gene whose association may be real. This module guards the naming and the
refusal to call any of them a mechanism.
"""
from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

try:
    import pytest  # noqa: F401
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

try:
    import numpy  # noqa: F401  (import guard only)
    import sklearn  # noqa: F401  (import guard only)
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"numpy and scikit-learn are required: {exc}") from exc

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "name_destroyed_signal_replacement.py"
POOLS = REPO_ROOT / "registry" / "QUALIFIED_PROBE_POOLS.json"


def load():
    spec = importlib.util.spec_from_file_location("namer", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NAMER = load()


class TestTheStratumLookup(unittest.TestCase):
    def pools(self):
        return {
            "admissible_features": {"A_qualified_human": ["f1", "f2"]},
            "undetermined_features": ["f3"],
            "conflicting_features": ["f4"],
            "symbol_names_another_gene_features": ["f5"],
            "platform_self_contradictory_features": ["f6"],
        }

    def test_every_stratum_is_distinguished(self):
        index = NAMER.verdict_index(self.pools())
        self.assertEqual(index["f1"], "admissible_pool_a")
        self.assertEqual(index["f3"], "undetermined")
        self.assertEqual(index["f4"], "conflicting")
        self.assertEqual(index["f5"], "symbol_names_another_gene")
        self.assertEqual(index["f6"], "platform_self_contradictory")

    def test_an_unknown_feature_is_not_silently_called_admissible(self):
        index = NAMER.verdict_index(self.pools())
        self.assertNotIn("f99", index)

    def test_admissibility_wins_over_a_refused_listing(self):
        # A feature cannot be in both, but if the file ever said so, the
        # stronger claim must not be invented from the weaker listing.
        pools = self.pools()
        pools["undetermined_features"].append("f1")
        self.assertEqual(NAMER.verdict_index(pools)["f1"], "admissible_pool_a")

    def test_the_committed_pool_file_feeds_the_lookup(self):
        index = NAMER.verdict_index(json.loads(POOLS.read_text(encoding="utf-8")))
        # The two features behind the finding that produced this script.
        self.assertEqual(index["10033668425"], "conflicting")   # CYP2D17
        self.assertEqual(index["10023822237"], "admissible_pool_a")  # CYP2A7


class TestTheScriptRefusesToOverclaim(unittest.TestCase):
    def source(self):
        # Flattened: a phrase that a line break splits is still the same
        # statement, and a test that fails on reflowed text teaches nothing.
        return " ".join(SCRIPT.read_text(encoding="utf-8").split())

    def test_it_states_that_a_name_is_not_a_mechanism(self):
        source = self.source()
        self.assertIn("what_this_does_not_establish", source)
        self.assertIn("not a demonstrated mechanism", source)
        self.assertIn("never a confirmed one", source)

    def test_it_says_cross_hybridisation_is_not_checked(self):
        self.assertIn("no source distributes the", self.source())

    def test_it_names_the_three_outcomes_the_cohort_contains(self):
        source = self.source()
        for phrase in ("non-human duplicate", "admissible paralogue",
                       "different gene"):
            self.assertIn(phrase, source)

    def test_it_changes_no_arm_loss_or_margin(self):
        source = self.source()
        self.assertIn("changes no arm, loss, margin or contrast", source)
        for forbidden in ("NON_INFERIORITY", "paired_contrast", "loss_summary"):
            self.assertNotIn(forbidden, source)

    def test_the_target_gene_comes_from_the_frozen_rule(self):
        source = self.source()
        self.assertIn("seeds_from_endpoint_label", source)
        self.assertNotIn('split("_")[0]', source)

    def test_it_reuses_the_runner_rather_than_reimplementing_it(self):
        source = self.source()
        for call in ("H2A.repeated_folds", "H2A.fold_matrix",
                     "H2A._spearman_strength", "H2A.build_universe",
                     "H2A.restrict_universe_to_qualified_pool"):
            self.assertIn(call, source)

    def test_a_missing_target_is_reported_as_not_assessable(self):
        self.assertIn("not_assessable", self.source())


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
