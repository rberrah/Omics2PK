"""Tests for the F-series functional-omics layer.

The suite is organised around the four things that could make the experiment
dishonest rather than merely wrong:

* a cytochrome re-entering the regulon it is predicted from (circularity);
* the frozen table moving between the document, the module and the runner;
* the "sign inversion" control being reported as a control when it is an
  algebraic identity;
* a number being returned where `NOT_EVALUABLE` is the truthful answer.

Everything here runs on synthetic matrices and a synthetic network. No test
reads the Human Liver Cohort, so the suite stays runnable without the archives
and no test can be made to pass by looking at a result.
"""
from __future__ import annotations

import importlib.util
import re
import sys
import unittest
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk import functional_omics as FO  # noqa: E402

PRESPEC = REPO_ROOT / "docs" / (
    "OMICS2PK_F_SERIES_FUNCTIONAL_OMICS_PRESPECIFICATION_20261109.md")
RUNNER = REPO_ROOT / "scripts" / "run_f_series_hlc_functional_omics.py"


def load_runner():
    spec = importlib.util.spec_from_file_location("_f_series_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FakeFrame:
    """The three columns the sign resolver reads, without requiring pandas."""

    def __init__(self, rows):
        self._rows = list(rows)
        self.columns = ["source_genesymbol", "target_genesymbol",
                        "consensus_stimulation", "consensus_inhibition"]

    def __len__(self):
        return len(self._rows)

    def __getitem__(self, key):
        index = self.columns.index(key)
        return _FakeColumn([row[index] for row in self._rows])


class _FakeColumn:
    def __init__(self, values):
        self._values = list(values)

    def astype(self, _kind):
        return self

    @property
    def str(self):
        return self

    def strip(self):
        return _FakeColumn([str(v).strip() for v in self._values])

    def lower(self):
        return _FakeColumn([str(v).lower() for v in self._values])

    def eq(self, other):
        return _FakeColumn([str(v) == other for v in self._values])

    def to_numpy(self, dtype=None):
        return np.asarray(self._values, dtype=dtype if dtype else None)

    def contains(self, *_args, **_kwargs):
        return _FakeColumn([False] * len(self._values))


def network(rows):
    """rows: (source, target, stimulation, inhibition) with python booleans."""
    return FakeFrame([(s, t, "True" if up else "False", "True" if down else "False")
                      for s, t, up, down in rows])


# --------------------------------------------------------------------------
class TestCircularityGuards(unittest.TestCase):
    """The defect that made the first draft of this experiment worthless."""

    def test_the_predicted_enzyme_leaves_every_regulon(self):
        signed = FO.resolve_signed_network(network([
            ("NR1I2", "CYP3A4", True, False),
            ("NR1I2", "ABCB1", True, False),
            ("NFKB1", "CYP3A4", False, True),
            ("NFKB1", "IL6", True, False),
        ]))
        for regulator, regulon in signed.regulons().items():
            self.assertNotIn("CYP3A4", regulon,
                             f"{regulator} still carries the predicted enzyme")

    def test_every_cytochrome_leaves_not_only_the_study_enzymes(self):
        """A paralogue measures nearly the same transcript; it goes too."""
        signed = FO.resolve_signed_network(network([
            ("AHR", "CYP1A1", True, False),
            ("AHR", "CYP2A7", True, False),
            ("AHR", "CYP3A4P2", True, False),
            ("AHR", "ALDH1A1", True, False),
        ]))
        self.assertEqual(set(signed.regulons()["AHR"]), {"ALDH1A1"})

    def test_the_guard_counts_what_it_removed(self):
        signed = FO.resolve_signed_network(network([
            ("AHR", "CYP1A2", True, False),
            ("AHR", "NQO1", True, False),
        ]))
        self.assertEqual(signed.audit["dropped_cytochrome_target"], 1)
        self.assertTrue(signed.audit["cytochrome_targets_removed"])

    def test_candidate_regulators_are_read_before_the_removal(self):
        """Selecting regulators by the enzyme is a prior; keeping it is circular."""
        frame = network([
            ("NR1I2", "CYP3A4", True, False),
            ("NFKB1", "CYP3A4", False, True),
            ("AHR", "NQO1", True, False),
        ])
        self.assertEqual(FO.candidate_regulators(frame, "CYP3A4"),
                         ["NFKB1", "NR1I2"])
        self.assertNotIn("CYP3A4",
                         FO.resolve_signed_network(frame).regulons().get("NR1I2", {}))

    def test_an_alias_probe_is_refused_even_with_a_clean_array_symbol(self):
        audit = FO.cytochrome_guard_audit(
            array_symbols=["ALB", "FLJ20359", "CYP3A4"],
            feature_ids=["1", "2", "3"],
            resolved_symbols={"1": "ALB", "2": "CYP2W1", "3": "CYP3A4"},
            admissible_features={"1", "2", "3"})
        self.assertEqual(audit["refused_by_array_symbol"], 1)
        self.assertEqual(audit["refused_by_resolved_alias"], 1)
        self.assertEqual(audit["resolved_alias_examples"][0]["array_symbol"],
                         "FLJ20359")
        self.assertEqual(set(audit["refused_feature_ids"]), {"2", "3"})

    def test_the_audit_reports_zero_rather_than_asserting_it(self):
        audit = FO.cytochrome_guard_audit(
            array_symbols=["ALB"], feature_ids=["1"],
            resolved_symbols={"1": "ALB"}, admissible_features={"1"})
        self.assertEqual(audit["refused_by_resolved_alias"], 0)
        self.assertIn("how_to_read", audit)

    def test_a_cytochrome_in_the_scoring_universe_stops_the_run(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.assert_no_cytochrome_reaches_scoring(["ALB", "CYP2D6", "NQO1"])
        self.assertTrue(FO.assert_no_cytochrome_reaches_scoring(["ALB", "NQO1"]))

    def test_the_alias_path_also_stops_the_run(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.assert_no_cytochrome_reaches_scoring(
                ["ALB", "FLJ20359"], resolved_symbols={"FLJ20359": "CYP2W1"})

    def test_the_pattern_recognises_the_shapes_a_list_would_miss(self):
        for symbol in ("CYP3A4", "CYP2D6", "CYP2A7", "CYP1A1", "CYP4F22",
                       "CYP27C1", "CYP3A4P2", "cyp3a5"):
            self.assertTrue(FO.is_cytochrome(symbol), symbol)
        for symbol in ("ALB", "NR1I2", "CYPRIN", "ACYP1", "", None):
            self.assertFalse(FO.is_cytochrome(symbol), symbol)

    def test_the_three_guards_are_all_still_in_the_module(self):
        source = (REPO_ROOT / "src" / "omics2pk" / "functional_omics.py").read_text(
            encoding="utf-8")
        self.assertIn("drop_cytochrome_targets", source)
        self.assertIn("refused_by_resolved_alias", source)
        self.assertIn("assert_no_cytochrome_reaches_scoring", source)


# --------------------------------------------------------------------------
class TestSignResolution(unittest.TestCase):
    def test_one_sided_flags_become_plus_and_minus_one(self):
        signed = FO.resolve_signed_network(network([
            ("A", "X", True, False), ("B", "X", False, True)]))
        self.assertEqual(signed.edges[("A", "X")], FO.STIMULATION)
        self.assertEqual(signed.edges[("B", "X")], FO.INHIBITION)

    def test_an_indeterminate_edge_leaves_rather_than_weighing_zero(self):
        """A zero-weight member would still count towards the evaluability floor."""
        signed = FO.resolve_signed_network(network([
            ("A", "X", True, True), ("A", "Y", False, False),
            ("A", "Z", True, False)]))
        self.assertEqual(set(signed.regulons()["A"]), {"Z"})
        self.assertEqual(signed.audit["dropped_sign_indeterminate"], 2)

    def test_identical_duplicates_collapse_to_one_edge(self):
        signed = FO.resolve_signed_network(network([
            ("RELA", "X", False, True), ("RELA", "X", False, True)]))
        self.assertEqual(signed.audit["duplicate_rows_collapsed"], 1)
        self.assertEqual(signed.edges[("RELA", "X")], FO.INHIBITION)

    def test_a_contradictory_duplicate_is_dropped_not_arbitrated(self):
        signed = FO.resolve_signed_network(network([
            ("A", "X", True, False), ("A", "X", False, True),
            ("A", "Y", True, False)]))
        self.assertNotIn(("A", "X"), signed.edges)
        self.assertEqual(signed.audit["dropped_conflicting_duplicate"], 1)

    def test_an_empty_resolution_raises_instead_of_returning_nothing(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.resolve_signed_network(network([("A", "X", True, True)]))

    def test_a_missing_column_is_named(self):
        class Bare:
            columns = ["source_genesymbol"]
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.resolve_signed_network(Bare())


class TestComplexRegulators(unittest.TestCase):
    def test_a_complex_expands_to_its_members(self):
        self.assertEqual(FO.complex_members("NFKB1_RELA"), ["NFKB1", "RELA"])

    def test_a_plain_regulator_is_its_own_member(self):
        self.assertEqual(FO.complex_members("NR1I2"), ["NR1I2"])

    def test_the_separator_is_the_one_the_resource_uses(self):
        self.assertEqual(FO.COMPLEX_SEPARATOR, "_")


# --------------------------------------------------------------------------
class TestUnivariateLinearModel(unittest.TestCase):
    """The score is the published statistic, checked against its closed form."""

    def setUp(self):
        rng = np.random.default_rng(11)
        self.genes, self.donors = 40, 7
        self.values = rng.normal(size=(self.genes, self.donors))
        self.weights = np.zeros((self.genes, 2))
        self.weights[:10, 0] = 1.0
        self.weights[10:16, 0] = -1.0
        self.weights[5:20, 1] = 1.0

    def test_it_matches_an_explicit_least_squares_fit_per_donor(self):
        scores = FO.ulm_activities(self.values, self.weights)
        for column in range(self.weights.shape[1]):
            design = np.column_stack(
                [np.ones(self.genes), self.weights[:, column]])
            for donor in range(self.donors):
                y = self.values[:, donor]
                beta, *_ = np.linalg.lstsq(design, y, rcond=None)
                residual = y - design @ beta
                sigma2 = (residual @ residual) / (self.genes - 2)
                covariance = sigma2 * np.linalg.inv(design.T @ design)
                expected = beta[1] / np.sqrt(covariance[1, 1])
                self.assertAlmostEqual(scores[donor, column], expected, places=8)

    def test_a_weight_vector_with_no_spread_yields_no_score_not_zero(self):
        weights = np.ones((self.genes, 1))
        scores = FO.ulm_activities(self.values, weights)
        self.assertTrue(np.all(np.isnan(scores)))

    def test_mismatched_shapes_are_refused(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.ulm_activities(self.values, np.ones((self.genes + 1, 1)))

    def test_two_genes_cannot_carry_a_t_statistic(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.ulm_activities(np.ones((2, 3)), np.array([[1.0], [-1.0]]))

    def test_the_weight_matrix_places_signs_at_the_right_genes(self):
        matrix = FO.weight_matrix(
            {"A": {"g1": 1, "g3": -1}}, ["A"], ["g1", "g2", "g3"])
        self.assertEqual(list(matrix[:, 0]), [1.0, 0.0, -1.0])

    def test_a_regulon_gene_absent_from_the_universe_is_simply_absent(self):
        matrix = FO.weight_matrix(
            {"A": {"g1": 1, "offarray": -1}}, ["A"], ["g1", "g2"])
        self.assertEqual(list(matrix[:, 0]), [1.0, 0.0])


class TestFoldSafeStandardisation(unittest.TestCase):
    def test_the_mean_and_spread_come_from_the_training_donors_only(self):
        values = np.array([[1.0, 2.0, 30.0], [0.0, 0.0, 0.0]])
        out = FO.standardise_on_training(values, train=[0, 1])
        self.assertAlmostEqual(float(out[0, :2].mean()), 0.0, places=12)
        self.assertGreater(out[0, 2], 10.0)

    def test_a_gene_with_no_training_spread_becomes_a_constant_zero(self):
        values = np.array([[5.0, 5.0, 9.0], [1.0, 2.0, 3.0]])
        out = FO.standardise_on_training(values, train=[0, 1])
        self.assertTrue(np.all(out[0] == 0.0))

    def test_one_training_donor_is_refused(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.standardise_on_training(np.ones((3, 3)), train=[0])

    def test_a_one_dimensional_input_is_refused(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.standardise_on_training(np.ones(5), train=[0, 1])


# --------------------------------------------------------------------------
class TestFullInversionIsAnIdentity(unittest.TestCase):
    """The control a reader expects, shown to be worthless rather than omitted."""

    def setUp(self):
        rng = np.random.default_rng(3)
        self.values = rng.normal(size=(60, 12))
        self.regulons = {"A": {f"g{i}": (1 if i % 3 else -1) for i in range(20)}}
        self.genes = [f"g{i}" for i in range(60)]

    def test_inverting_every_sign_negates_every_activity_exactly(self):
        flipped = FO.perturb_regulon_signs(
            self.regulons, 1.0, np.random.default_rng(0))
        base = FO.ulm_activities(
            self.values, FO.weight_matrix(self.regulons, ["A"], self.genes))
        other = FO.ulm_activities(
            self.values, FO.weight_matrix(flipped, ["A"], self.genes))
        np.testing.assert_allclose(other, -base, atol=1e-10)

    def test_a_ridge_on_the_negated_column_predicts_identically(self):
        runner = load_runner()
        rng = np.random.default_rng(5)
        x = rng.normal(size=(40, 3))
        y = x @ np.array([1.0, -2.0, 0.5]) + rng.normal(scale=0.1, size=40)
        negated = x.copy()
        negated[:, 2] *= -1.0
        first, _ = runner.H2A.ridge_prediction(
            x[:30], y[:30], x[30:], inner_folds=3, seed=1)
        second, _ = runner.H2A.ridge_prediction(
            negated[:30], y[:30], negated[30:], inner_folds=3, seed=1)
        np.testing.assert_allclose(first, second, atol=1e-10)

    def test_a_partial_flip_genuinely_changes_the_score(self):
        partial = FO.perturb_regulon_signs(
            self.regulons, 0.5, np.random.default_rng(7))
        base = FO.ulm_activities(
            self.values, FO.weight_matrix(self.regulons, ["A"], self.genes))
        other = FO.ulm_activities(
            self.values, FO.weight_matrix(partial, ["A"], self.genes))
        self.assertFalse(np.allclose(other, base, atol=1e-6))
        self.assertFalse(np.allclose(other, -base, atol=1e-6))

    def test_the_document_says_the_full_inversion_is_not_a_control(self):
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertIn("algebraic_identity_check", text)
        self.assertIn("never `control`", text)

    def test_the_runner_labels_it_as_the_identity_check(self):
        runner = load_runner()
        names = [s["name"] for s in runner.control_settings(2)]
        self.assertIn("algebraic_identity_check", names)
        self.assertEqual(sum(1 for n in names if n == "algebraic_identity_check"), 1)


class TestPartialPerturbation(unittest.TestCase):
    def setUp(self):
        self.regulons = {"A": {f"g{i}": 1 for i in range(8)},
                         "B": {f"h{i}": -1 for i in range(4)}}

    def test_the_flipped_count_is_the_requested_fraction(self):
        out = FO.perturb_regulon_signs(
            self.regulons, 0.5, np.random.default_rng(0))
        flipped = sum(1 for gene, weight in out["A"].items()
                      if weight != self.regulons["A"][gene])
        self.assertEqual(flipped, 4)

    def test_a_zero_fraction_changes_nothing(self):
        out = FO.perturb_regulon_signs(
            self.regulons, 0.0, np.random.default_rng(0))
        self.assertEqual(out, self.regulons)

    def test_the_membership_never_changes_only_the_signs(self):
        out = FO.perturb_regulon_signs(
            self.regulons, 0.75, np.random.default_rng(1))
        for regulator in self.regulons:
            self.assertEqual(set(out[regulator]), set(self.regulons[regulator]))

    def test_a_fraction_outside_the_unit_interval_is_refused(self):
        for bad in (-0.1, 1.1):
            with self.assertRaises(FO.FunctionalOmicsError):
                FO.perturb_regulon_signs(
                    self.regulons, bad, np.random.default_rng(0))

    def test_the_same_seed_reproduces_the_same_perturbation(self):
        first = FO.perturb_regulon_signs(
            self.regulons, 0.5, np.random.default_rng(42))
        second = FO.perturb_regulon_signs(
            self.regulons, 0.5, np.random.default_rng(42))
        self.assertEqual(first, second)


class TestStructurallyComparableNull(unittest.TestCase):
    def setUp(self):
        self.regulons = {"A": {"g1": 1, "g2": 1, "g3": -1},
                         "B": {"g2": -1, "g4": -1}}
        self.universe = [f"g{i}" for i in range(1, 30)]

    def test_the_null_keeps_the_size_and_the_signed_composition(self):
        out = FO.rewire_regulons(
            self.regulons, self.universe, np.random.default_rng(0))
        for regulator, regulon in self.regulons.items():
            self.assertEqual(len(out[regulator]), len(regulon))
            self.assertEqual(
                sum(1 for w in out[regulator].values() if w > 0),
                sum(1 for w in regulon.values() if w > 0))
            self.assertEqual(
                sum(1 for w in out[regulator].values() if w < 0),
                sum(1 for w in regulon.values() if w < 0))

    def test_it_changes_which_genes_are_named(self):
        out = FO.rewire_regulons(
            self.regulons, self.universe, np.random.default_rng(2))
        self.assertNotEqual(set(out["A"]), set(self.regulons["A"]))

    def test_no_cytochrome_can_be_drawn_into_a_null_regulon(self):
        out = FO.rewire_regulons(
            self.regulons, self.universe + ["CYP3A4", "CYP2D6"],
            np.random.default_rng(0))
        for regulon in out.values():
            for gene in regulon:
                self.assertFalse(FO.is_cytochrome(gene))

    def test_members_are_drawn_without_replacement(self):
        out = FO.rewire_regulons(
            {"A": {f"g{i}": 1 for i in range(1, 25)}}, self.universe,
            np.random.default_rng(0))
        self.assertEqual(len(set(out["A"])), 24)

    def test_a_regulon_larger_than_the_universe_raises(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.rewire_regulons({"A": {f"g{i}": 1 for i in range(50)}},
                               ["g1", "g2"], np.random.default_rng(0))

    def test_an_empty_universe_raises(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            FO.rewire_regulons(self.regulons, ["CYP3A4"],
                               np.random.default_rng(0))


class TestDestroyedLink(unittest.TestCase):
    def test_the_donor_permutation_keeps_the_values_and_moves_the_donors(self):
        values = np.arange(12, dtype=float).reshape(3, 4)
        permuted, order = FO.permute_donors(values, np.random.default_rng(0))
        self.assertEqual(sorted(order.tolist()), [0, 1, 2, 3])
        np.testing.assert_allclose(np.sort(permuted, axis=1),
                                   np.sort(values, axis=1))


# --------------------------------------------------------------------------
class TestEvaluability(unittest.TestCase):
    def test_a_regulator_under_the_floor_is_refused_with_a_reason(self):
        regulons = {"A": {f"g{i}": 1 for i in range(10)}, "B": {"g1": 1, "g2": -1}}
        kept, refused = FO.evaluable_regulators(
            regulons, ["A", "B"], [f"g{i}" for i in range(10)])
        self.assertEqual(kept, ["A"])
        self.assertIn("B", refused)
        self.assertIn("fewer than 5", refused["B"]["reason"])
        self.assertEqual(refused["B"]["admissible_members"], 2)

    def test_inadmissible_members_do_not_count_towards_the_floor(self):
        regulons = {"A": {f"g{i}": 1 for i in range(10)}}
        kept, refused = FO.evaluable_regulators(regulons, ["A"], ["g0", "g1"])
        self.assertEqual(kept, [])
        self.assertEqual(refused["A"]["members"], 10)
        self.assertEqual(refused["A"]["admissible_members"], 2)

    def test_a_regulator_with_no_surviving_regulon_is_named(self):
        kept, refused = FO.evaluable_regulators({}, ["A"], ["g1"])
        self.assertEqual(kept, [])
        self.assertIn("no regulon survived", refused["A"]["reason"])

    def test_a_cytochrome_member_is_not_counted_even_if_admissible(self):
        regulons = {"A": {"CYP3A4": 1, "CYP2D6": 1, "g1": 1, "g2": 1, "g3": 1}}
        kept, refused = FO.evaluable_regulators(
            regulons, ["A"], ["CYP3A4", "CYP2D6", "g1", "g2", "g3"])
        self.assertEqual(kept, [])
        self.assertEqual(refused["A"]["admissible_members"], 3)

    def test_no_admissible_target_probe_makes_only_the_transcript_arms_refuse(self):
        verdict = FO.endpoint_evaluability(
            endpoint="CYP1A2_activity", target_admissible_probes=0,
            evaluable_regulator_count=30, donors_with_both=398)
        for arm in ("F1", "F2", "F5"):
            self.assertEqual(verdict["verdict"][arm], "NOT_EVALUABLE")
        for arm in ("F0", "F4p"):
            self.assertEqual(verdict["verdict"][arm], "EVALUABLE")
        self.assertIn("not a negative result", verdict["reasons"]["F1/F2/F5"])

    def test_too_few_regulators_refuses_only_the_activity_arms(self):
        verdict = FO.endpoint_evaluability(
            endpoint="X_activity", target_admissible_probes=2,
            evaluable_regulator_count=2, donors_with_both=300)
        self.assertEqual(verdict["verdict"]["F4p"], "NOT_EVALUABLE")
        self.assertEqual(verdict["verdict"]["F5"], "NOT_EVALUABLE")
        self.assertEqual(verdict["verdict"]["F1"], "EVALUABLE")

    def test_too_few_donors_refuses_every_arm(self):
        verdict = FO.endpoint_evaluability(
            endpoint="X_activity", target_admissible_probes=2,
            evaluable_regulator_count=40, donors_with_both=10)
        self.assertTrue(all(v == "NOT_EVALUABLE"
                            for v in verdict["verdict"].values()))
        self.assertIn("below the frozen floor", verdict["reasons"]["all_arms"])

    def test_a_fully_evaluable_endpoint_says_so(self):
        verdict = FO.endpoint_evaluability(
            endpoint="CYP3A4M_activity", target_admissible_probes=3,
            evaluable_regulator_count=50, donors_with_both=397)
        self.assertTrue(all(v == "EVALUABLE"
                            for v in verdict["verdict"].values()))
        self.assertEqual(verdict["reasons"], {})


# --------------------------------------------------------------------------
class TestFreezeAgreesWithTheDocument(unittest.TestCase):
    """The freeze lives in two places and must say the same thing in both."""

    @classmethod
    def setUpClass(cls):
        cls.text = PRESPEC.read_text(encoding="utf-8")

    def row(self, label):
        match = re.search(
            r"^\|\s*" + re.escape(label) + r"\s*\|\s*(.+?)\s*\|\s*$",
            self.text, re.MULTILINE)
        self.assertIsNotNone(match, f"the document no longer states {label!r}")
        return match.group(1)

    def test_the_fold_seed_matches(self):
        self.assertEqual(int(self.row("Fold seed").strip("`")), FO.FROZEN["seed"])

    def test_the_bootstrap_count_matches(self):
        self.assertEqual(int(self.row("Bootstrap draws").replace(",", "")),
                         FO.FROZEN["n_bootstrap"])

    def test_the_relative_threshold_matches(self):
        self.assertEqual(float(self.row("Relative MAE threshold")),
                         FO.FROZEN["relative_mae_threshold"])

    def test_the_trim_fraction_matches(self):
        self.assertEqual(float(self.row("Trim fraction")),
                         FO.FROZEN["trim_fraction"])

    def test_the_three_floors_match(self):
        self.assertEqual(int(self.row("Regulon floor (admissible members)")),
                         FO.FROZEN["regulon_floor"])
        self.assertEqual(int(self.row("Regulator floor per endpoint")),
                         FO.FROZEN["regulator_floor"])
        self.assertEqual(int(self.row("Donor floor per endpoint")),
                         FO.FROZEN["donor_floor"])

    def test_the_control_replicate_count_matches(self):
        self.assertEqual(int(self.row("Control replicates").split()[0]),
                         FO.FROZEN["control_replicates"])

    def test_the_status_ceiling_matches_and_is_exploratory(self):
        self.assertEqual(self.row("Status ceiling").strip("`"),
                         FO.FROZEN["status_ceiling"])
        self.assertEqual(FO.FROZEN["status_ceiling"], "exploratory")

    def test_the_runner_checks_the_document_itself(self):
        runner = load_runner()
        checked = runner.check_freeze_against_document(PRESPEC)
        self.assertEqual(len(checked["checked"]), 9)
        self.assertEqual(len(checked["sha256"]), 64)

    def test_the_runner_refuses_a_document_that_disagrees(self):
        import tempfile

        runner = load_runner()
        edited = self.text.replace("| Fold seed | `20261109` |",
                                   "| Fold seed | `1` |")
        self.assertNotEqual(edited, self.text)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "edited.md"
            path.write_text(edited, encoding="utf-8")
            with self.assertRaises(FO.FunctionalOmicsError):
                runner.check_freeze_against_document(path)

    def test_the_runner_refuses_a_document_missing_a_row(self):
        import tempfile

        runner = load_runner()
        edited = re.sub(r"^\| Status ceiling .*$", "", self.text,
                        flags=re.MULTILINE)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "edited.md"
            path.write_text(edited, encoding="utf-8")
            with self.assertRaises(FO.FunctionalOmicsError):
                runner.check_freeze_against_document(path)


class TestPrimaryContrastIsTheCommissionedOne(unittest.TestCase):
    """The commissioning instruction put F5 against F1, not F4p against F1."""

    def test_the_primary_contrast_is_f5_against_f1(self):
        self.assertEqual(FO.PRIMARY_ARM, "F5")
        self.assertEqual(FO.PRIMARY_COMPARATOR, "F1")

    def test_the_primary_endpoint_is_the_named_demonstrator(self):
        self.assertEqual(FO.PRIMARY_ENDPOINT, "CYP3A4M_activity")

    def test_both_key_secondaries_are_present(self):
        self.assertIn(("F4p", "F1"), FO.KEY_SECONDARY_CONTRASTS)
        self.assertIn(("F5", "F2"), FO.KEY_SECONDARY_CONTRASTS)

    def test_the_document_states_the_primary_contrast_in_words(self):
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertIn("`F5` versus `F1` on `CYP3A4M_activity`", text)

    def test_the_same_enzyme_replicate_is_not_an_independent_validation(self):
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertEqual(FO.SAME_ENZYME_REPLICATE, "CYP3A4T_activity")
        self.assertIn("never counted as an independent validation", text)

    def test_f3_is_reserved_with_a_stated_reason_not_silently_absent(self):
        self.assertIn("F3", FO.RESERVED_UNIMPLEMENTED_ARMS)
        self.assertIn("xenobiotic", FO.RESERVED_UNIMPLEMENTED_ARMS["F3"])
        self.assertNotIn("F3", FO.ARMS)

    def test_the_arms_are_exactly_the_five_specified(self):
        self.assertEqual(FO.ARMS, ("F0", "F1", "F2", "F4p", "F5"))

    def test_the_f_labels_do_not_collide_with_the_protocol_arms(self):
        for label in FO.ARMS:
            self.assertFalse(label.startswith("M"))
            self.assertFalse(label.startswith("P"))
            self.assertNotEqual(label, "HY")


class TestRunnerGovernance(unittest.TestCase):
    def setUp(self):
        self.runner = load_runner()

    def test_the_estimator_is_h2a_s_own(self):
        self.assertTrue(self.runner.check_estimator_machinery_unchanged())

    def test_the_fold_seed_differs_from_h2a_s(self):
        self.assertNotEqual(self.runner.H2A.SEED, FO.FROZEN["seed"])

    def test_the_ridge_grid_is_the_frozen_one(self):
        self.assertEqual(tuple(self.runner.H2A.RIDGE_ALPHAS),
                         tuple(FO.FROZEN["ridge_alphas"]))

    def test_a_drifting_estimator_stops_the_run(self):
        original = self.runner.H2A.INNER_FOLDS
        try:
            self.runner.H2A.INNER_FOLDS = 7
            with self.assertRaises(FO.FunctionalOmicsError):
                self.runner.check_estimator_machinery_unchanged()
        finally:
            self.runner.H2A.INNER_FOLDS = original

    def test_the_control_settings_are_enumerated_before_any_of_them_runs(self):
        settings = self.runner.control_settings(FO.FROZEN["control_replicates"])
        kinds = {}
        for setting in settings:
            kinds[setting["kind"]] = kinds.get(setting["kind"], 0) + 1
        self.assertEqual(kinds["real"], 1)
        self.assertEqual(kinds["rewired_null"], FO.FROZEN["control_replicates"])
        self.assertEqual(kinds["donor_permuted"], FO.FROZEN["control_replicates"])
        self.assertEqual(
            kinds["sign_perturbation"],
            FO.FROZEN["control_replicates"] * len(FO.FROZEN["perturbation_fractions"])
            + 1)

    def test_every_control_setting_carries_its_own_seed(self):
        settings = self.runner.control_settings(3)
        seeds = [s["seed"] for s in settings if s["kind"] != "real"]
        self.assertTrue(all(seed is not None for seed in seeds))
        by_kind: dict[str, set] = {}
        for setting in settings:
            if setting["kind"] in ("rewired_null", "donor_permuted"):
                by_kind.setdefault(setting["kind"], set()).add(setting["seed"])
        self.assertEqual(len(by_kind["rewired_null"] & by_kind["donor_permuted"]), 0)

    def test_the_setting_names_are_unique(self):
        names = [s["name"] for s in self.runner.control_settings(5)]
        self.assertEqual(len(names), len(set(names)))

    def test_an_unknown_control_kind_raises(self):
        with self.assertRaises(FO.FunctionalOmicsError):
            self.runner.activities_for_setting(
                {"kind": "invented", "seed": 1, "fraction": None},
                regulons={}, regulators=[], genes=[],
                standardised=np.ones((5, 3)), universe_genes=[])

    def test_a_repeat_is_never_counted_as_a_second_donor(self):
        block = np.array([[1.0, 3.0, np.nan], [3.0, 5.0, np.nan]])
        values, mask = self.runner.per_donor(block)
        self.assertEqual(values.size, 2)
        np.testing.assert_allclose(values, [2.0, 4.0])
        self.assertEqual(mask.tolist(), [True, True, False])

    def test_the_bonferroni_interval_widens_with_the_family(self):
        rng = np.random.default_rng(0)
        difference = rng.normal(loc=-0.3, scale=1.0, size=200)
        reference = np.abs(rng.normal(loc=5.0, scale=1.0, size=200))
        single = self.runner.bonferroni_interval(
            difference, reference, seed=1, family_size=1)
        family = self.runner.bonferroni_interval(
            difference, reference, seed=1, family_size=9)
        self.assertLess(family["ci_low"], single["ci_low"])
        self.assertGreater(family["ci_high"], single["ci_high"])
        self.assertAlmostEqual(family["alpha"], 0.05 / 9)

    def test_the_two_sided_binomial_is_right_on_a_known_case(self):
        self.assertAlmostEqual(self.runner._binomial_p(10, 10), 2 / 1024)
        self.assertAlmostEqual(self.runner._binomial_p(5, 10), 1.0)

    def test_a_dry_run_flag_exists_so_the_machinery_can_be_checked_unfrozen(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn("--dry-run", source)
        self.assertIn("stopping before any model is fitted", source)


class TestNoClinicalClaim(unittest.TestCase):
    """Nothing here may be presented as a transfer to clinical pharmacokinetics."""

    def test_the_runner_writes_the_transfer_as_not_evaluated(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn('"clinical_pk_transfer": "NOT_EVALUATED"', source)

    def test_no_module_function_converts_an_activity_into_an_exposure(self):
        for forbidden in ("clearance", "auc", "concentration", "trough", "cmin",
                          "dose"):
            self.assertNotIn(
                forbidden, [name.lower() for name in FO.__all__],
                f"the module exposes a {forbidden} function")

    def test_the_endpoint_class_is_stated_as_in_vitro(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn('"endpoint_class": "functional_in_vitro"', source)

    def test_the_document_names_the_missing_data_rather_than_scheduling_the_step(self):
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertIn("ABSTAIN", text)
        self.assertIn("concentration-time pharmacokinetics", text)

    def test_the_document_forbids_presenting_a_simulation_as_a_validation(self):
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertIn("no F-series artefact contains a simulated concentration",
                      text)


class TestReportingSeparation(unittest.TestCase):
    def test_the_runner_emits_three_blocks_and_names_them(self):
        source = RUNNER.read_text(encoding="utf-8")
        for block in ('"predictive_performance"', '"biological_relevance"',
                      '"analytical_validity"'):
            self.assertIn(block, source)

    def test_the_biology_block_states_what_it_is_not(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn("what_it_is_not", source)
        self.assertIn("never upgrades the predictive block", source)

    def test_the_document_forbids_borrowing_between_blocks(self):
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertIn("never upgrades block 1", text)

    def test_not_evaluable_is_reported_rather_than_dropped(self):
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn('"not_evaluable"', source)
        text = PRESPEC.read_text(encoding="utf-8")
        self.assertIn("never counted as a negative result", text)

    def test_the_status_ceiling_travels_with_every_contrast(self):
        runner = load_runner()
        source = RUNNER.read_text(encoding="utf-8")
        self.assertIn('"status_ceiling": FO.FROZEN["status_ceiling"]', source)
        self.assertEqual(FO.FROZEN["status_ceiling"], "exploratory")
        self.assertTrue(callable(runner.contrast))


class TestPrespecificationContent(unittest.TestCase):
    """The document has to carry the commitments that make the run honest."""

    @classmethod
    def setUpClass(cls):
        cls.text = PRESPEC.read_text(encoding="utf-8")

    def test_it_records_the_unfavourable_prior_in_advance(self):
        self.assertIn("The honest prior is unfavourable", self.text)
        self.assertIn("0 of\n   80 comparisons", self.text)

    def test_it_says_a_null_result_is_a_result(self):
        self.assertIn("A null result is therefore the expected outcome", self.text)

    def test_it_states_the_three_part_success_criterion(self):
        self.assertIn("relative MAE reduction **≥ 10 %**", self.text)
        self.assertIn("excludes zero", self.text)
        self.assertIn("no reversal on **both**", self.text)

    def test_it_names_cyp1a2_as_not_evaluable_before_the_run(self):
        self.assertIn("`CYP1A2_activity` is `NOT_EVALUABLE`", self.text)

    def test_it_states_the_three_place_cytochrome_removal(self):
        self.assertIn("Removing every cytochrome, in three places", self.text)
        self.assertIn("822 of 64,516", self.text)

    def test_it_states_the_donor_count_that_matters(self):
        self.assertIn("**397**", self.text)
        self.assertIn("not 487", self.text)

    def test_it_lists_what_was_deliberately_not_built(self):
        self.assertIn("Deliberately not built", self.text)
        for item in ("SIGNOR", "PROGENy", "COSMOS"):
            self.assertIn(item, self.text)

    def test_it_forbids_moving_the_rules_after_the_results(self):
        self.assertIn("What may not be done after seeing the results", self.text)

    def test_it_declares_the_f_labels_do_not_replace_the_protocol_arms(self):
        self.assertIn("do\n**not** replace `M0`–`M5`", self.text)

    def test_it_carries_a_change_log(self):
        self.assertIn("## 13. Change log", self.text)
        self.assertIn("| 1.0.0 | 2026-11-09 |", self.text)

    def test_the_sign_rule_names_the_consensus_columns_not_the_flags(self):
        """The two disagree on 4,806 edges, so which one is read is not cosmetic."""
        self.assertIn("the consensus columns decide", self.text)
        self.assertIn("2,590 (4.01 %)", self.text)
        self.assertIn("are **not** a restatement of the flags", self.text)

    def test_the_candidate_rule_requires_a_resolvable_sign(self):
        self.assertIn("**sign-resolvable** CollecTRI edge", self.text)
        self.assertIn("sign-resolvable candidates", self.text)

    def test_it_records_the_two_pre_freeze_corrections_rather_than_hiding_them(self):
        self.assertIn("before the freeze and before any model was fitted",
                      self.text)
        self.assertIn("No threshold, floor, arm, endpoint, loss or primary "
                      "contrast changed.", self.text)

    def test_the_resolved_sign_composition_is_stated(self):
        self.assertIn("26 are activating and 24 repressing", self.text)

    def test_the_f2_arm_is_declared_narrower_than_f5(self):
        self.assertIn("`F2` is therefore the narrower arm", self.text)


if __name__ == "__main__":
    unittest.main()
