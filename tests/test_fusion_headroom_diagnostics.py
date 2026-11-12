"""Tests for the corrected fusion-headroom diagnostics.

The module exists because an earlier ceiling was mathematically wrong: a
constant weight vector cannot express a gate, so the best constant mixture does
not bound contextual fusion. The first test class is that counterexample, and
it is the one that must never be allowed to pass by accident.

The second thing under test is the verdict rule. Treating an interval that
contains zero as a refutation confuses non-detection with absence, which is the
error this project corrects everywhere else and had committed here.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PROTOCOL = ROOT / "docs" / "OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md"
RESULTS = ROOT / "docs" / "OMICS2PK_FUSION_HEADROOM_RESULTS_20261111.md"
ARTEFACT = ROOT / "reports" / "data" / "fusion_headroom_diagnostics.json"


@pytest.fixture(scope="module")
def mod():
    from omics2pk import fusion_headroom_diagnostics
    return fusion_headroom_diagnostics


@pytest.fixture(scope="module")
def payload():
    if not ARTEFACT.is_file():
        pytest.skip("the diagnostic has not been run in this checkout")
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


def _cells(payload):
    for name, block in payload["cohorts"].items():
        if "per_repeat" in block:
            for label, inner in block["per_repeat"].items():
                yield f"{name}::{label}", inner
        else:
            yield name, block


# --------------------------------------------------------------------------
# 1. the counterexample — the reason the module exists
# --------------------------------------------------------------------------
class TestConstantWeightsCannotSeeAGate:
    def test_the_counterexample_has_the_shape_it_claims(self, mod):
        predictions, truth, context = \
            mod.counterexample_constant_weights_cannot_see_a_gate(200)
        assert predictions.shape == (200, 2)
        assert np.all(truth == 0.0)
        # The two arms bracket the truth nowhere strictly, but each is exactly
        # right on one half, which is what makes a gate perfect.
        assert np.isclose(np.abs(predictions[:, 0] - truth).mean(), 0.5)
        assert np.isclose(np.abs(predictions[:, 1] - truth).mean(), 0.5)
        assert context.shape == (200,)

    def test_every_constant_weight_keeps_the_same_loss(self, mod):
        predictions, truth, _ = \
            mod.counterexample_constant_weights_cannot_see_a_gate(200)
        for weight in np.linspace(0.0, 1.0, 21):
            mixed = predictions @ np.array([weight, 1.0 - weight])
            assert np.isclose(np.abs(mixed - truth).mean(), 0.5), (
                f"weight {weight} should keep the loss at 0.5; if it does not "
                "the counterexample no longer demonstrates the point")

    def test_the_constant_ceiling_shows_no_gain(self, mod):
        predictions, truth, _ = \
            mod.counterexample_constant_weights_cannot_see_a_gate(200)
        assert mod.constant_weight_oracle(
            predictions, truth)["mean_loss"] == pytest.approx(0.5, abs=1e-6)

    def test_the_per_individual_ceilings_reach_zero(self, mod):
        predictions, truth, _ = \
            mod.counterexample_constant_weights_cannot_see_a_gate(200)
        assert mod.individual_selection_oracle(
            predictions, truth)["mean_loss"] == pytest.approx(0.0, abs=1e-12)
        assert mod.individual_convex_oracle(
            predictions, truth)["mean_loss"] == pytest.approx(0.0, abs=1e-12)

    def test_the_context_predicts_the_winner_perfectly_here(self, mod):
        predictions, truth, context = \
            mod.counterexample_constant_weights_cannot_see_a_gate(200)
        probe = mod.context_informativeness(predictions, truth, context)
        assert probe["sign_agreement"] == pytest.approx(1.0)

    def test_the_self_test_asserts_all_of_it(self, mod):
        out = mod.self_test()
        assert out["constant_weight_oracle"] == pytest.approx(0.5, abs=1e-6)
        assert out["individual_convex_oracle"] == pytest.approx(0.0, abs=1e-12)
        assert "would have discarded the hypothesis it was meant to test" in \
            out["what_this_proves"]

    def test_the_self_test_fails_if_the_convex_ceiling_is_broken(self, mod,
                                                                monkeypatch):
        # If someone replaces the per-individual ceiling with the constant one,
        # the self-test must notice. It is the only guard against the original
        # error returning.
        monkeypatch.setattr(mod, "individual_convex_oracle",
                            mod.constant_weight_oracle)
        with pytest.raises(mod.FusionDiagnosticsError):
            mod.self_test()


# --------------------------------------------------------------------------
# 2. the ceilings are ordered as the mathematics requires
# --------------------------------------------------------------------------
class TestCeilingOrdering:
    @pytest.fixture
    def cohort(self):
        rng = np.random.default_rng(4)
        truth = rng.normal(size=120)
        arms = np.column_stack([
            truth + rng.normal(scale=0.5, size=120),
            truth + rng.normal(scale=0.7, size=120),
            truth * 0.5 + rng.normal(scale=0.6, size=120),
        ])
        return arms, truth

    def test_the_convex_ceiling_is_the_tightest(self, mod, cohort):
        predictions, truth = cohort
        convex = mod.individual_convex_oracle(predictions, truth)["mean_loss"]
        for other in (mod.best_single_oracle(predictions, truth)["mean_loss"],
                      mod.constant_weight_oracle(predictions,
                                                 truth)["mean_loss"],
                      mod.individual_selection_oracle(predictions,
                                                      truth)["mean_loss"]):
            assert convex <= other + 1e-9

    def test_the_selection_ceiling_is_below_the_best_single(self, mod, cohort):
        predictions, truth = cohort
        assert (mod.individual_selection_oracle(predictions,
                                                truth)["mean_loss"]
                <= mod.best_single_oracle(predictions, truth)["mean_loss"]
                + 1e-9)

    def test_the_constant_ceiling_is_never_worse_than_the_best_single(
            self, mod, cohort):
        predictions, truth = cohort
        assert (mod.constant_weight_oracle(predictions, truth)["mean_loss"]
                <= mod.best_single_oracle(predictions, truth)["mean_loss"]
                + 1e-9)

    def test_the_honest_choice_can_be_worse_than_the_oracle_choice(
            self, mod, cohort):
        predictions, truth = cohort
        assert (mod.best_single_honest(predictions, truth)["mean_loss"]
                >= mod.best_single_oracle(predictions, truth)["mean_loss"]
                - 1e-9)

    def test_the_bracketing_rate_is_a_proportion(self, mod, cohort):
        predictions, truth = cohort
        rate = mod.individual_convex_oracle(
            predictions, truth)["bracketing_rate"]
        assert 0.0 <= rate <= 1.0

    def test_identical_arms_give_no_headroom_anywhere(self, mod):
        rng = np.random.default_rng(5)
        truth = rng.normal(size=80)
        one = truth + rng.normal(scale=0.4, size=80)
        predictions = np.column_stack([one, one, one])
        reference = np.abs(one - truth).mean()
        for name in ("constant_weight_oracle", "individual_selection_oracle",
                     "individual_convex_oracle"):
            assert getattr(mod, name)(
                predictions, truth)["mean_loss"] == pytest.approx(
                    reference, abs=1e-9), name

    def test_each_ceiling_declares_what_it_bounds(self, mod, cohort):
        predictions, truth = cohort
        assert "constant-weight mixtures only" in mod.constant_weight_oracle(
            predictions, truth)["bounds"]
        assert "EVERY method in the family" in mod.individual_convex_oracle(
            predictions, truth)["bounds"]

    def test_at_least_two_arms_are_required(self, mod):
        with pytest.raises(mod.FusionDiagnosticsError,
                           match="at least two arms"):
            mod.individual_convex_oracle(np.zeros((5, 1)), np.zeros(5))


# --------------------------------------------------------------------------
# 3. the verdict rule does not confuse non-detection with absence
# --------------------------------------------------------------------------
class TestVerdictRule:
    def test_an_interval_excluding_zero_is_headroom(self, mod):
        gap = {"gap": 0.5, "interval": [0.2, 0.8]}
        assert mod.verdict(gap, 1.0)["verdict"] == "HEADROOM_PRESENT"

    def test_an_upper_bound_below_the_margin_is_the_refutation(self, mod):
        gap = {"gap": 0.001, "interval": [-0.01, 0.02]}
        out = mod.verdict(gap, 1.0, margin=0.05)
        assert out["verdict"] == "REFUTED_CEILING_BELOW_MARGIN"
        assert "IS a refutation" in out["how_to_read"]

    def test_an_interval_containing_zero_and_reaching_past_the_margin_is_not(
            self, mod):
        gap = {"gap": 0.02, "interval": [-0.05, 0.40]}
        out = mod.verdict(gap, 1.0, margin=0.05)
        assert out["verdict"] == "NO_DETECTABLE_HEADROOM_UNDERPOWERED"
        assert "nothing is refuted" in out["how_to_read"]

    def test_the_rule_states_the_principle_explicitly(self, mod):
        out = mod.verdict({"gap": 0.0, "interval": [-1.0, 1.0]}, 1.0)
        assert out["an_interval_containing_zero_is_not_a_refutation"] is True

    def test_a_non_positive_reference_loss_is_refused(self, mod):
        with pytest.raises(mod.FusionDiagnosticsError,
                           match="reference loss must be positive"):
            mod.verdict({"gap": 0.0, "interval": [0.0, 0.0]}, 0.0)

    def test_the_margin_scales_with_the_reference(self, mod):
        gap = {"gap": 0.0, "interval": [-0.01, 0.03]}
        assert mod.verdict(gap, 1.0, margin=0.05)["verdict"] == \
            "REFUTED_CEILING_BELOW_MARGIN"
        assert mod.verdict(gap, 0.1, margin=0.05)["verdict"] == \
            "NO_DETECTABLE_HEADROOM_UNDERPOWERED"


# --------------------------------------------------------------------------
# 4. the probe, and that it too is three-valued
# --------------------------------------------------------------------------
class TestContextProbe:
    def test_an_informative_context_is_detected(self, mod):
        predictions, truth, context = \
            mod.counterexample_constant_weights_cannot_see_a_gate(300)
        probe = mod.context_informativeness(predictions, truth, context)
        assert probe["verdict"] == "SIGNAL_REACHES_MARGIN", (
            "a perfect gate must clear the declared 0.30 margin, not merely "
            "exclude zero")

    def test_a_pure_noise_context_at_a_large_n_is_refuted(self, mod):
        rng = np.random.default_rng(6)
        n = 600
        truth = rng.normal(size=n)
        predictions = np.column_stack([
            truth + rng.normal(scale=0.5, size=n),
            truth + rng.normal(scale=0.5, size=n)])
        context = rng.normal(size=(n, 3))
        probe = mod.context_informativeness(predictions, truth, context)
        assert probe["verdict"] == "REFUTED_CONTEXT_UNINFORMATIVE"

    def test_a_pure_noise_context_at_a_small_n_is_not_refuted(self, mod):
        rng = np.random.default_rng(7)
        n = 20
        truth = rng.normal(size=n)
        predictions = np.column_stack([
            truth + rng.normal(scale=0.5, size=n),
            truth + rng.normal(scale=0.5, size=n)])
        context = rng.normal(size=(n, 3))
        probe = mod.context_informativeness(predictions, truth, context)
        assert probe["verdict"] == "NO_DETECTABLE_SIGNAL_UNDERPOWERED", (
            "at n = 20 a null probe must not be called a refutation")

    def test_the_bonferroni_interval_is_wider(self, mod):
        rng = np.random.default_rng(8)
        n = 300
        truth = rng.normal(size=n)
        predictions = np.column_stack([
            truth + rng.normal(scale=0.5, size=n),
            truth + rng.normal(scale=0.8, size=n)])
        context = rng.normal(size=(n, 2))
        narrow = mod.context_informativeness(predictions, truth, context,
                                             family_size=1)
        wide = mod.context_informativeness(predictions, truth, context,
                                           family_size=20)
        assert (wide["spearman_interval_bonferroni"][1]
                - wide["spearman_interval_bonferroni"][0]
                >= narrow["spearman_interval_bonferroni"][1]
                - narrow["spearman_interval_bonferroni"][0])

    def test_it_declares_the_minimum_interesting_association(self, mod):
        predictions, truth, context = \
            mod.counterexample_constant_weights_cannot_see_a_gate(100)
        probe = mod.context_informativeness(predictions, truth, context)
        assert probe["minimum_interesting_rank_association"] == 0.30

    def test_a_misaligned_context_is_refused(self, mod):
        predictions, truth, _ = \
            mod.counterexample_constant_weights_cannot_see_a_gate(100)
        with pytest.raises(mod.FusionDiagnosticsError, match="disagree on n"):
            mod.context_informativeness(predictions, truth, np.zeros(50))


# --------------------------------------------------------------------------
# 5. weight learning is evaluated out of sample
# --------------------------------------------------------------------------
class TestNestedStacking:
    def test_it_refits_the_weights_per_partition(self, mod):
        rng = np.random.default_rng(9)
        truth = rng.normal(size=150)
        predictions = np.column_stack([
            truth + rng.normal(scale=0.5, size=150),
            truth + rng.normal(scale=0.9, size=150)])
        out = mod.nested_constant_stacking(predictions, truth)
        assert out["weights_refitted_per_partition"] is True
        assert "as opposed to the ceiling" in out["what_this_is"]

    def test_it_can_be_worse_than_the_hindsight_ceiling(self, mod):
        # The point of the nested figure: having to learn the weights costs
        # something, so it must not be allowed to equal the ceiling by
        # construction.
        rng = np.random.default_rng(10)
        truth = rng.normal(size=60)
        predictions = np.column_stack([
            truth + rng.normal(scale=1.0, size=60),
            rng.normal(scale=1.0, size=60)])
        ceiling = mod.constant_weight_oracle(predictions, truth)["mean_loss"]
        nested = mod.nested_constant_stacking(predictions,
                                              truth)["mean_loss"]
        assert nested >= ceiling - 1e-9


# --------------------------------------------------------------------------
# 6. the module refuses to be read as a performance
# --------------------------------------------------------------------------
class TestFraming:
    def test_it_states_it_implements_no_fusion(self, mod):
        out = mod.what_this_module_does_not_do()
        assert out["implements_no_fusion_method"] is True
        assert out["fits_no_gate"] is True

    def test_it_says_ceilings_are_not_performances(self, mod):
        out = mod.what_this_module_does_not_do()
        assert "None may be quoted as a result" in out[
            "ceilings_are_not_performances"]

    def test_it_says_a_permissive_ceiling_licenses_nothing(self, mod):
        out = mod.what_this_module_does_not_do()
        assert "cannot promise" in out[
            "a_permissive_ceiling_licenses_nothing"]

    def test_the_module_docstring_records_the_correction(self, mod):
        text = mod.__doc__
        assert "That rule is wrong" in text
        assert "does not bound contextual fusion" in text


# --------------------------------------------------------------------------
# 7. the executed artefact, and the documents that report it
# --------------------------------------------------------------------------
class TestArtefact:
    def test_it_carries_the_correction_and_its_proof(self, payload):
        block = payload["corrects"]
        assert "does not bound contextual fusion" in block["what_was_wrong"]
        assert block["proof_that_the_correction_matters"][
            "constant_weight_oracle"] == pytest.approx(0.5, abs=1e-6)
        assert block["proof_that_the_correction_matters"][
            "individual_convex_oracle"] == pytest.approx(0.0, abs=1e-12)

    def test_the_convex_ceiling_is_never_above_the_constant_one(self, payload):
        for name, block in _cells(payload):
            ceilings = block["ceilings"]
            assert (ceilings["individual_convex_oracle"]["mean_loss"]
                    <= ceilings["constant_weight_oracle"]["mean_loss"]
                    + 1e-6), name

    def test_every_cell_has_both_verdicts_and_two_probes(self, payload):
        for name, block in _cells(payload):
            assert set(block["verdicts"]) == {"stacking_family",
                                              "whole_fusion_family"}
            assert set(block["context_probes"]) == {
                "clinical_covariates",
                "arm_predictions_and_their_disagreement"}, name

    def test_no_cell_was_refuted_on_the_governing_ceiling(self, payload):
        # The executed finding. If this ever changes, the results document and
        # the roadmap must both be revised.
        verdicts = {block["verdicts"]["whole_fusion_family"]["verdict"]
                    for _, block in _cells(payload)}
        assert verdicts == {"HEADROOM_PRESENT"}

    def test_the_constant_ceiling_was_refuted_somewhere(self, payload):
        # The proof that the old rule would have stopped the programme.
        refuted = [name for name, block in _cells(payload)
                   if block["verdicts"]["stacking_family"]["verdict"]
                   == "REFUTED_CEILING_BELOW_MARGIN"]
        assert len(refuted) >= 1, (
            "if no cell refutes the constant-weight ceiling, the claim that "
            "the old rule would have stopped the programme is unsupported")

    def test_it_exports_no_individual_level_data(self, payload):
        assert payload["no_individual_level_data_in_this_artefact"] is True

        def arrays(obj):
            found = []
            if isinstance(obj, dict):
                for value in obj.values():
                    found += arrays(value)
            elif isinstance(obj, list):
                numeric = [x for x in obj if isinstance(x, (int, float))]
                if len(numeric) >= 10:
                    found.append(len(numeric))
                for value in obj:
                    if not isinstance(value, (int, float)):
                        found += arrays(value)
            return found

        assert not arrays(payload)


class TestDocuments:
    def test_the_protocol_is_at_version_three(self):
        text = PROTOCOL.read_text(encoding="utf-8")
        assert "**Version 3.0.0 — 2026-11-11.**" in text

    def test_the_protocol_keeps_every_earlier_change_log_entry(self):
        """A version bump must not quietly drop the record of what was wrong."""
        text = PROTOCOL.read_text(encoding="utf-8")
        for version in ("| 1.0.0 |", "| **2.0.0** |", "| **3.0.0** |"):
            assert version in text

    def test_the_protocol_requires_a_nested_evaluation(self):
        flat = re.sub(r"[\s*_`>|]+", " ", PROTOCOL.read_text(encoding="utf-8"))
        assert "genuinely nested evaluation" in flat
        assert "A null contrast is only readable if the harness can detect a "             "gain that exists" in flat

    def test_the_protocol_applies_each_margin_to_the_lower_bound(self):
        flat = re.sub(r"[\s*_`>|]+", " ", PROTOCOL.read_text(encoding="utf-8"))
        assert "applies its margin to both sides" in flat
        # the flattening regex strips underscores along with markup
        assert "GAIN REACHES MARGIN" in flat
        assert "SIGNAL REACHES MARGIN" in flat

    def test_the_protocol_records_that_q2_is_not_authorised(self):
        flat = re.sub(r"[\s*_`>|]+", " ", PROTOCOL.read_text(encoding="utf-8"))
        assert "Q1 was answered negatively" in flat
        assert "NOT authorised" in flat

    def test_the_protocol_records_the_error_rather_than_hiding_it(self):
        flat = re.sub(r"[\s*_`>|]+", " ", PROTOCOL.read_text(encoding="utf-8"))
        assert "version 1.0.0 had it wrong" in flat.lower()
        assert "wrong in the direction that would have discarded the " \
            "hypothesis it was meant to test" in flat

    def test_the_protocol_splits_the_two_questions(self):
        flat = re.sub(r"[\s*_`>|]+", " ", PROTOCOL.read_text(encoding="utf-8"))
        assert "Two co-primary questions" in flat
        assert "the first can succeed while the second fails" in flat

    def test_the_protocol_no_longer_contains_the_old_verdict_name(self):
        assert "REFUTED_NO_HEADROOM" not in PROTOCOL.read_text(
            encoding="utf-8")

    def test_the_results_document_matches_the_artefact(self, payload):
        text = RESULTS.read_text(encoding="utf-8")
        for name, block in _cells(payload):
            value = block["ceilings"]["individual_convex_oracle"]["mean_loss"]
            forms = {f"{value:.3f}", f"{value:.2f}", f"{value:.1f}"}
            assert any(form in text for form in forms), (
                f"{name}'s convex ceiling is in the artefact and not in the "
                "results document")

    def test_the_results_document_states_what_is_untested(self):
        flat = re.sub(r"[\s*_`>|]+", " ", RESULTS.read_text(encoding="utf-8"))
        assert "Q2 — does mechanistic knowledge help exploit that context?" \
            in flat or "Untested" in flat

    def test_the_results_document_flags_the_cyp3a4_contradiction(self):
        flat = re.sub(r"[\s*_`>|]+", " ", RESULTS.read_text(encoding="utf-8"))
        assert "The two CYP3A4 substrates disagree" in flat
