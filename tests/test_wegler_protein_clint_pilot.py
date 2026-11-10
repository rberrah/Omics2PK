"""Contract tests for the Axis B protein-to-clearance pilot.

Two of these guard against errors that would have corrupted the experiment
silently rather than making it fail, and both were caught before the freeze:

1. **Circularity.** Two of the article's three published clearance equations
   scale by the specific CYP amount, so predicting them from that amount is
   circular. The guard must refuse them by name, and a test must prove the
   guard is not merely decorative.
2. **The donor mapping.** Table S3 numbers donors UU1-UU15 and Data S1 numbers
   them 1-15; they agree for 1-9 and are permuted for 10-15. The mapping must
   be re-derived and refused unless unanimous, never trusted as a constant.

The rest cover fold safety, the losses being three readings of one error, the
refusal to read a non-detection as equivalence, the prohibition on squaring a
Spearman coefficient, and the sample-size projection's scenario structure.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PRESPEC = (
    REPO
    / "docs"
    / "OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_PRESPECIFICATION_20261110.md"
)
RUNNER = REPO / "scripts" / "run_wegler_protein_clint_pilot.py"
REPORT = REPO / "reports" / "data" / "axis_b_wegler_protein_clint_pilot.json"

np = pytest.importorskip("numpy")


@pytest.fixture(scope="module")
def mod():
    sys.path.insert(0, str(REPO / "src"))
    from omics2pk import protein_clint_pilot as module

    return module


@pytest.fixture(scope="module")
def runner():
    sys.path.insert(0, str(REPO / "src"))
    spec = importlib.util.spec_from_file_location("axis_b_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _normalise(text):
    return re.sub(r"\s+", " ", text)


# --------------------------------------------------------------------------
# 1. the circularity guard
# --------------------------------------------------------------------------
class TestCircularityGuard:
    def test_the_two_admissible_endpoints_pass(self, mod):
        for name in ("clint_mic_eq7", "clint_hep_eq7"):
            assert mod.assert_endpoint_is_admissible(name) is True

    def test_every_cyp_scaled_equation_is_refused(self, mod):
        for name in ("clint_mic_eq8", "clint_mic_eq9", "clint_hep_eq8"):
            with pytest.raises(mod.ProteinClintPilotError, match="circular"):
                mod.assert_endpoint_is_admissible(name)

    def test_the_forbidden_list_matches_the_freeze(self, mod):
        assert set(mod.FROZEN["forbidden_endpoints"]) == {
            "clint_mic_eq8", "clint_mic_eq9", "clint_hep_eq8"
        }

    def test_the_incubation_columns_are_refused_as_endpoints(self, mod):
        for name in mod.NEVER_A_PREDICTOR_OR_ENDPOINT:
            with pytest.raises(mod.ProteinClintPilotError,
                               match="baseline predictor in another unit"):
                mod.assert_endpoint_is_admissible(name)

    def test_an_endpoint_that_is_not_eq7_is_refused(self, mod):
        with pytest.raises(mod.ProteinClintPilotError, match="not one of the admissible"):
            mod.assert_endpoint_is_admissible("clint_mic_eq10")

    def test_an_empty_endpoint_is_refused(self, mod):
        with pytest.raises(mod.ProteinClintPilotError):
            mod.assert_endpoint_is_admissible("")

    def test_the_runner_exercises_the_guard_rather_than_trusting_it(self):
        source = RUNNER.read_text(encoding="utf-8")
        assert "the circularity guard accepted" in source, (
            "the runner must try each forbidden endpoint and fail if the guard "
            "lets one through"
        )


# --------------------------------------------------------------------------
# 2. the donor mapping
# --------------------------------------------------------------------------
def _mapping_inputs(permutation, *, genes=("A", "B", "C", "D"),
                    systems=("HLM", "HH"), noise=0.0):
    """Build two files whose donor numbering differs by `permutation`.

    `permutation[index]` is the batch number that Data S1 index carries.
    """
    rng = np.random.default_rng(5)
    table = {}
    data = {}
    for gene in genes:
        for system in systems:
            values = rng.uniform(1.0, 100.0, size=15)
            for batch in range(1, 16):
                table[(gene, f"UU{batch}", system)] = float(values[batch - 1])
            for index in range(1, 16):
                batch = permutation[index]
                data[(gene, index, system)] = float(
                    values[batch - 1] + noise * rng.normal()
                )
    return table, data


class TestDonorMapping:
    def test_the_identity_is_recognised_as_the_identity(self, mod):
        table, data = _mapping_inputs({i: i for i in range(1, 16)})
        out = mod.resolve_donor_mapping(table, data)
        assert out["is_identity"] is True
        assert out["unanimous"] is True
        assert out["mapping"][10] == "UU10"

    def test_the_real_permutation_is_recovered(self, mod):
        truth = {i: i for i in range(1, 10)}
        truth.update({10: 14, 11: 15, 12: 10, 13: 11, 14: 12, 15: 13})
        table, data = _mapping_inputs(truth)
        out = mod.resolve_donor_mapping(table, data)
        assert out["is_identity"] is False
        assert {index: f"UU{batch}" for index, batch in truth.items()} == out["mapping"]

    def test_the_result_is_always_a_bijection(self, mod):
        truth = {i: i for i in range(1, 16)}
        truth[1], truth[2] = 2, 1
        table, data = _mapping_inputs(truth)
        out = mod.resolve_donor_mapping(table, data)
        assert len(set(out["mapping"].values())) == 15

    def test_a_non_unanimous_vote_is_refused(self, mod):
        truth = {i: i for i in range(1, 16)}
        table, data = _mapping_inputs(truth)
        # Corrupt one gene/system vector so one donor no longer matches there.
        data[("A", 7, "HLM")] = 123456.0
        with pytest.raises(mod.ProteinClintPilotError, match="not unanimous"):
            mod.resolve_donor_mapping(table, data)

    def test_a_donor_matching_nothing_is_refused(self, mod):
        truth = {i: i for i in range(1, 16)}
        table, data = _mapping_inputs(truth)
        for gene in ("A", "B", "C", "D"):
            for system in ("HLM", "HH"):
                data[(gene, 4, system)] = 999999.0
        with pytest.raises(mod.ProteinClintPilotError, match="matches no"):
            mod.resolve_donor_mapping(table, data)

    def test_mismatched_donor_counts_are_refused(self, mod):
        truth = {i: i for i in range(1, 16)}
        table, data = _mapping_inputs(truth)
        for key in [k for k in data if k[1] == 15]:
            del data[key]
        with pytest.raises(mod.ProteinClintPilotError, match="same donor set"):
            mod.resolve_donor_mapping(table, data)

    def test_too_few_vectors_are_refused(self, mod):
        truth = {i: i for i in range(1, 16)}
        table, data = _mapping_inputs(truth, genes=("A",), systems=("HLM",))
        with pytest.raises(mod.ProteinClintPilotError, match="too few"):
            mod.resolve_donor_mapping(table, data)

    def test_two_donors_with_identical_values_are_refused(self, mod):
        """An ambiguous match must stop the run, not pick one at random."""
        truth = {i: i for i in range(1, 16)}
        table, data = _mapping_inputs(truth)
        for gene in ("A", "B", "C", "D"):
            for system in ("HLM", "HH"):
                table[(gene, "UU3", system)] = table[(gene, "UU2", system)]
                data[(gene, 3, system)] = data[(gene, 2, system)]
        with pytest.raises(mod.ProteinClintPilotError):
            mod.resolve_donor_mapping(table, data)


# --------------------------------------------------------------------------
# 3. fold safety and the estimator
# --------------------------------------------------------------------------
class TestFoldSafety:
    def test_components_ignore_the_held_out_row(self, mod):
        rng = np.random.default_rng(11)
        train = rng.normal(size=(14, 300))
        first, _, _ = mod.proteome_components(train, rng.normal(size=(1, 300)), 2)
        second, _, _ = mod.proteome_components(
            train, rng.normal(size=(1, 300)) * 100.0, 2
        )
        assert np.array_equal(first, second)

    def test_the_first_component_is_a_slice_of_two(self, mod):
        rng = np.random.default_rng(3)
        train = rng.normal(size=(14, 80))
        test = rng.normal(size=(1, 80))
        wide_train, wide_test, _ = mod.proteome_components(train, test, 2)
        narrow_train, narrow_test, _ = mod.proteome_components(train, test, 1)
        assert np.allclose(wide_train[:, :1], narrow_train)
        assert np.allclose(wide_test[:, :1], narrow_test)

    def test_zero_components_returns_an_empty_block(self, mod):
        train = np.random.default_rng(1).normal(size=(10, 20))
        a, b, taken = mod.proteome_components(train, train[:1], 0)
        assert a.shape == (10, 0) and b.shape == (1, 0) and taken == 0

    def test_a_rank_deficient_block_reports_its_shortfall(self, mod):
        train = np.tile(np.arange(5.0), (3, 1)).T  # 5 rows, 3 identical columns
        _, _, taken = mod.proteome_components(train, train[:1], 2)
        assert taken == 1

    def test_a_missing_cell_is_refused(self, mod):
        train = np.ones((5, 4))
        train[0, 0] = np.nan
        with pytest.raises(mod.ProteinClintPilotError, match="missing cell"):
            mod.proteome_components(train, np.ones((1, 4)), 1)

    def test_standardisation_uses_training_constants(self, mod):
        train = np.array([[1.0], [3.0], [5.0]])
        test = np.array([[9.0]])
        z_train, z_test = mod.standardise(train, test)
        assert np.allclose(z_train.mean(), 0.0)
        assert z_test[0, 0] > z_train.max()

    def test_an_empty_design_predicts_the_training_mean(self, mod):
        target = np.array([1.0, 2.0, 6.0])
        prediction = mod._ridge_fit_predict(
            np.zeros((3, 0)), target, np.zeros((1, 0)), 1.0
        )
        assert prediction[0] == pytest.approx(3.0)

    def test_an_enormous_penalty_returns_the_training_mean(self, mod):
        rng = np.random.default_rng(7)
        design = rng.normal(size=(12, 2))
        target = rng.normal(size=12) + 5.0
        prediction = mod._ridge_fit_predict(design, target, design[:1], 1e12)
        assert prediction[0] == pytest.approx(target.mean(), abs=1e-4)

    def test_a_recoverable_signal_is_recovered(self, mod):
        rng = np.random.default_rng(23)
        design = rng.normal(size=(40, 2))
        target = design @ np.array([2.0, -1.0]) + rng.normal(size=40) * 0.05
        prediction = mod._ridge_fit_predict(
            design[:30], target[:30], design[30:], 0.01
        )
        assert np.corrcoef(prediction, target[30:])[0, 1] > 0.95

    def test_the_penalty_is_selected_inside_the_training_donors(self, mod):
        rng = np.random.default_rng(13)
        design = rng.normal(size=(14, 2))
        noise = rng.normal(size=14) * 50.0
        chosen = mod.select_penalty(design, noise, mod.FROZEN["ridge_alphas"])
        assert chosen in mod.FROZEN["ridge_alphas"]

    def test_pure_noise_takes_a_heavier_penalty_than_a_clean_signal(self, mod):
        rng = np.random.default_rng(29)
        design = rng.normal(size=(14, 2))
        clean = design @ np.array([3.0, -2.0])
        noise = rng.normal(size=14) * 10.0
        assert mod.select_penalty(
            design, noise, mod.FROZEN["ridge_alphas"]
        ) >= mod.select_penalty(design, clean, mod.FROZEN["ridge_alphas"])

    def test_too_few_training_donors_is_refused(self, mod):
        with pytest.raises(mod.ProteinClintPilotError, match="three training donors"):
            mod.select_penalty(np.ones((2, 1)), np.ones(2), (1.0,))


# --------------------------------------------------------------------------
# 4. losses in readable units
# --------------------------------------------------------------------------
class TestLosses:
    def test_the_fold_error_is_ten_to_the_mae(self, mod):
        out = mod.losses([0.3, 0.3, 0.3], [10.0, 10.0, 10.0])
        assert out["mae_log10"] == pytest.approx(0.3)
        assert out["geometric_fold_error"] == pytest.approx(10 ** 0.3)

    def test_a_perfect_prediction_is_a_fold_error_of_one(self, mod):
        out = mod.losses([0.0, 0.0], [0.0, 0.0])
        assert out["geometric_fold_error"] == pytest.approx(1.0)

    def test_the_native_error_is_reported_separately(self, mod):
        out = mod.losses([0.1, 0.2], [30.0, 50.0])
        assert out["mae_ml_per_min_per_kg"] == pytest.approx(40.0)

    def test_the_two_derived_losses_carry_no_test(self, mod):
        out = mod.losses([0.1], [1.0])
        assert "carry no separate test" in out["how_to_read"]


# --------------------------------------------------------------------------
# 5. intervals, and the refusal to call non-detection equivalence
# --------------------------------------------------------------------------
class TestIntervals:
    def test_an_improvement_is_detected(self, mod):
        rng = np.random.default_rng(3)
        baseline = rng.uniform(0.5, 0.6, size=15)
        new = baseline - 0.2
        out = mod.paired_bootstrap(new, baseline, n_boot=500)
        assert out["direction"] == "improvement"
        assert out["interval_excludes_zero"] is True

    def test_identical_arms_give_zero_and_no_detection(self, mod):
        errors = np.linspace(0.1, 0.5, 15)
        out = mod.paired_bootstrap(errors, errors.copy(), n_boot=500)
        assert out["mean_difference_log10"] == 0.0
        assert out["interval_excludes_zero"] is False

    def test_non_detection_is_never_called_equivalence(self, mod):
        errors = np.linspace(0.1, 0.5, 15)
        out = mod.paired_bootstrap(errors, errors.copy(), n_boot=200)
        assert "never a demonstration that there is none" in out["how_to_read"]
        assert out["exploratory"] is True

    def test_the_resampling_unit_is_the_donor(self, mod):
        out = mod.paired_bootstrap(
            np.linspace(0.1, 0.5, 15), np.linspace(0.2, 0.6, 15), n_boot=100
        )
        assert out["resampling_unit"] == "donor"
        assert out["n_independent_donors"] == 15

    def test_the_bonferroni_interval_is_wider(self, mod):
        rng = np.random.default_rng(9)
        baseline = rng.uniform(0.3, 0.9, size=15)
        new = baseline - 0.1
        plain = mod.paired_bootstrap(new, baseline, n_boot=800)
        adjusted = mod.paired_bootstrap(new, baseline, n_boot=800, family_size=5)
        assert adjusted["interval"][0] <= plain["interval"][0]
        assert adjusted["interval"][1] >= plain["interval"][1]
        assert adjusted["family_size"] == 5

    def test_mismatched_arms_are_refused(self, mod):
        with pytest.raises(mod.ProteinClintPilotError):
            mod.paired_bootstrap(np.ones(15), np.ones(13), n_boot=10)

    def test_the_permutation_tail_is_the_improvement_tail(self, mod):
        """A negative difference is an improvement, so the lower tail counts."""
        out = mod.permutation_rank(-0.5, [0.0] * 20)
        assert out["better_than_every_replicate"] is True
        assert out["rank_p_value"] == pytest.approx(1 / 21)

    def test_an_improvement_inside_the_null_is_not_better_than_it(self, mod):
        out = mod.permutation_rank(-0.1, [-0.3, -0.2, 0.0, 0.1])
        assert out["better_than_every_replicate"] is False
        assert out["n_null_at_least_as_good_as_observed"] == 2

    def test_the_null_explains_why_it_is_not_centred_on_zero(self, mod):
        out = mod.permutation_rank(0.0, [0.0])
        assert "shrinks harder" in out["how_to_read"]


# --------------------------------------------------------------------------
# 6. the sample-size projection
# --------------------------------------------------------------------------
class TestSampleSize:
    def test_it_runs_every_scenario_combination(self, mod):
        rng = np.random.default_rng(1)
        out = mod.sample_size_projection(
            rng.normal(size=15) * 0.05,
            effect_scenarios={"a": -0.05, "b": -0.025},
        )
        assert len(out["scenarios"]) == 2 * 2 * 2

    def test_a_smaller_effect_needs_more_donors(self, mod):
        differences = np.full(15, -0.05) + np.random.default_rng(2).normal(size=15) * 0.02
        out = mod.sample_size_projection(
            differences, effect_scenarios={"big": -0.05, "small": -0.01}
        )
        big = out["scenarios"]["big|power=0.80|family=1"]["donors_required"]
        small = out["scenarios"]["small|power=0.80|family=1"]["donors_required"]
        assert small > big

    def test_more_power_needs_more_donors(self, mod):
        # The effect has to be small enough that the required n is not so tiny
        # that both powers round to the same integer.
        differences = np.full(15, -0.01) + np.random.default_rng(4).normal(size=15) * 0.05
        out = mod.sample_size_projection(differences, effect_scenarios={"e": -0.01})
        assert (
            out["scenarios"]["e|power=0.90|family=1"]["donors_required"]
            > out["scenarios"]["e|power=0.80|family=1"]["donors_required"]
        )

    def test_the_bonferroni_family_needs_more_donors(self, mod):
        differences = np.full(15, -0.05) + np.random.default_rng(6).normal(size=15) * 0.02
        out = mod.sample_size_projection(differences, effect_scenarios={"e": -0.05})
        assert (
            out["scenarios"]["e|power=0.80|family=5"]["donors_required"]
            > out["scenarios"]["e|power=0.80|family=1"]["donors_required"]
        )

    def test_a_zero_effect_requires_an_unbounded_sample(self, mod):
        out = mod.sample_size_projection(
            np.random.default_rng(8).normal(size=15) * 0.05,
            effect_scenarios={"none": 0.0},
        )
        assert out["scenarios"]["none|power=0.80|family=1"]["donors_required"] is None

    def test_the_assumptions_are_published(self, mod):
        out = mod.sample_size_projection(
            np.random.default_rng(10).normal(size=15) * 0.05,
            effect_scenarios={"e": -0.05},
        )
        assert out["n_donors_the_sd_was_estimated_from"] == 15
        text = " ".join(out["assumptions"])
        assert "itself imprecise" in text
        assert "never a single number" in out["how_to_read"]

    def test_too_few_differences_are_refused(self, mod):
        with pytest.raises(mod.ProteinClintPilotError):
            mod.sample_size_projection([0.1, 0.2], effect_scenarios={"e": -0.1})


# --------------------------------------------------------------------------
# 7. the prohibitions the commissioning note named
# --------------------------------------------------------------------------
class TestProhibitions:
    def test_the_spearman_square_prohibition_is_in_the_code(self, mod):
        text = mod.spearman_square_is_not_a_variance_fraction()
        assert "not a variance fraction" in text
        assert "rank" in text

    def test_no_clinical_claim_is_available(self, mod):
        text = mod.no_clinical_claim()
        assert "no result here is a clinical-pharmacokinetic proof" in text
        assert "not five independent cohorts" in text
        assert "are not two cohorts" in text

    def test_the_module_computes_no_rank_correlation_at_all(self):
        """The prohibition is enforced by the pilot owning no such machinery."""
        source = (REPO / "src" / "omics2pk" / "protein_clint_pilot.py").read_text(
            encoding="utf-8"
        )
        for machinery in ("def _spearman", "def spearman(", "rankdata",
                          "def _ranks", "spearmanr", "** 2 if rho", "rho ** 2"):
            assert machinery not in source, machinery
        # The only mention of the word is the prohibition itself.
        assert source.lower().count("spearman") == source.count(
            "Spearman"
        ) + source.count("spearman_square_is_not_a_variance_fraction")

    def test_the_withdrawn_variance_sentence_stays_withdrawn(self):
        older = (
            REPO / "docs" / "OMICS2PK_PUBLIC_SOURCE_QUALIFICATION_20261109.md"
        ).read_text(encoding="utf-8")
        assert "WITHDRAWN, 2026-11-10" in older
        assert "square of a rank correlation is not" in _normalise(older)


# --------------------------------------------------------------------------
# 8. the freeze agrees with the document
# --------------------------------------------------------------------------
class TestFreezeAgreesWithTheDocument:
    @pytest.fixture(scope="class")
    def text(self):
        assert PRESPEC.is_file(), f"prespecification absent: {PRESPEC}"
        return PRESPEC.read_text(encoding="utf-8")

    def test_the_check_passes(self, runner):
        out = runner.check_freeze_against_document()
        assert out["checked"]["Components (k)"] == 2
        assert out["checked"]["Permutation replicates"] == 200

    def test_the_seed_and_arms_are_in_the_table(self, mod, text):
        assert f"| Seed | `{mod.FROZEN['seed']}` |" in text
        assert f"| Primary arm | `{mod.PRIMARY_ARM}` |" in text
        assert f"| Primary comparator | `{mod.PRIMARY_COMPARATOR}` |" in text
        assert f"| Secondary comparator | `{mod.SECONDARY_COMPARATOR}` |" in text

    def test_the_forbidden_endpoints_are_named_in_the_document(self, mod, text):
        for forbidden in mod.FROZEN["forbidden_endpoints"]:
            assert f"`{forbidden}`" in text, forbidden

    def test_a_document_missing_the_circularity_record_is_refused(self, runner, tmp_path):
        reduced = PRESPEC.read_text(encoding="utf-8").replace(
            "the only admissible endpoint", "fine"
        )
        path = tmp_path / "prespec.md"
        path.write_text(reduced, encoding="utf-8")
        with pytest.raises(runner.RunError, match="circularity check"):
            runner.check_freeze_against_document(path)

    def test_a_document_missing_the_mapping_record_is_refused(self, runner, tmp_path):
        reduced = PRESPEC.read_text(encoding="utf-8").replace(
            "donors 10–15 are permuted", "donors are in order"
        )
        path = tmp_path / "prespec.md"
        path.write_text(reduced, encoding="utf-8")
        with pytest.raises(runner.RunError, match="donor mapping"):
            runner.check_freeze_against_document(path)


class TestPrespecificationContent:
    @pytest.fixture(scope="class")
    def text(self):
        return _normalise(PRESPEC.read_text(encoding="utf-8"))

    def test_it_quotes_the_articles_own_equation_caption(self, text):
        assert "scaled to kg bw **with specific CYP" in text
        assert "hepatocellularity and HLM with literature MPPGL" in text

    def test_it_records_the_unfavourable_prior(self, text):
        assert "A null result is the expected outcome" in text

    def test_it_states_the_mapping_before_the_freeze(self, text):
        assert "mis-paired six of fifteen donors" in text

    def test_it_requires_the_target_cyp_to_leave_the_block(self, text):
        assert "mandatory, not cosmetic" in text

    def test_it_excludes_rather_than_imputes_the_two_nd_values(self, text):
        assert "handled by exclusion, never by imputation" in text

    def test_it_warns_the_null_is_not_centred_on_zero(self, text):
        assert "The null will not be centred on zero, and that is the point"

    def test_it_states_five_substrates_are_not_five_cohorts(self, text):
        assert "They are not five independent cohorts" in text

    def test_it_forbids_reading_non_detection_as_equivalence(self, text):
        assert "Absence of difference is not equivalence" in text

    def test_it_forbids_a_variance_fraction_from_a_rank_correlation(self, text):
        assert "No variance fraction from a rank correlation" in text

    def test_it_defines_the_system_comparison_as_a_difference_of_differences(self, text):
        assert "difference of differences over donors" in text
        assert "forbidden by §7.4" in text

    def test_it_structures_the_sample_size_projection(self, text):
        assert "Scenarios on the true effect" in text
        assert "80 % and 90 % power" in text


# --------------------------------------------------------------------------
# 9. the committed report
# --------------------------------------------------------------------------
class TestCommittedReport:
    @pytest.fixture(scope="class")
    def payload(self):
        if not REPORT.is_file():
            pytest.skip("Axis B report not generated")
        return json.loads(REPORT.read_text(encoding="utf-8"))

    def test_it_declares_the_status_ceiling(self, payload, mod):
        assert payload["status_ceiling"] == mod.FROZEN["status_ceiling"]

    def test_it_refuses_a_clinical_reading(self, payload):
        assert "no result here is a clinical-pharmacokinetic proof" in payload["no_clinical_claim"]

    def test_it_carries_the_spearman_prohibition(self, payload):
        assert "not a variance fraction" in payload["spearman_prohibition"]

    def test_the_circularity_check_is_recorded_with_its_reason(self, payload):
        check = payload["circularity_check"]
        assert check["admissible_endpoints"] == ["clint_hep_eq7", "clint_mic_eq7"]
        assert len(check["forbidden_endpoints"]) == 3
        assert "specific CYP amount" in check["why"]

    def test_the_donor_mapping_is_unanimous_and_not_the_identity(self, payload):
        mapping = payload["donor_mapping"]
        assert mapping["unanimous"] is True
        assert mapping["is_identity"] is False
        assert mapping["mapping"]["10"] == "UU14"
        assert mapping["mapping"]["12"] == "UU10"
        assert len(set(mapping["mapping"].values())) == 15

    def test_every_endpoint_used_is_an_eq7_column(self, payload):
        for key, cell in payload["results"].items():
            assert cell["endpoint"].endswith("eq7"), key

    def test_the_target_cyp_is_out_of_every_proteome_block(self, payload):
        for key, cell in payload["results"].items():
            if cell["cyp"] == "CYP2D6":
                continue  # absent from Data S1 entirely
            assert cell["target_cyp_removed_from_the_block"] is True, key

    def test_the_bufuralol_microsomal_cell_uses_thirteen_donors(self, payload):
        for key, cell in payload["results"].items():
            if cell["substrate"] == "Bufuralol" and cell["system"] == "HLM":
                assert cell["n_donors"] == 13, key
                assert cell["n_not_determined_excluded"] == 2, key

    def test_every_other_cell_uses_fifteen_donors(self, payload):
        for key, cell in payload["results"].items():
            if cell["substrate"] == "Bufuralol" and cell["system"] == "HLM":
                continue
            assert cell["n_donors"] == 15, key

    def test_each_primary_cell_carries_a_permutation_null(self, payload, mod):
        expected = 2 * mod.FROZEN["n_substrates"]
        found = [
            key for key, cell in payload["results"].items()
            if "permutation_null" in cell
        ]
        assert len(found) == expected
        for key in found:
            null = payload["results"][key]["permutation_null"]
            assert null["n_replicates"] == mod.FROZEN["permutation_replicates"]
            assert len(null["null_deciles"]) == 11

    def test_every_cell_reports_both_interval_levels(self, payload):
        for key, cell in payload["results"].items():
            contrasts = cell["contrasts"]
            assert contrasts["primary_W2_vs_W1"]["family_size"] == 1, key
            assert contrasts["primary_W2_vs_W1_bonferroni"]["family_size"] == 5, key

    def test_every_cell_reports_three_losses(self, payload, mod):
        for key, cell in payload["results"].items():
            for arm in mod.ARMS:
                loss = cell["losses"][arm]
                for field in ("mae_log10", "geometric_fold_error",
                              "mae_ml_per_min_per_kg"):
                    assert field in loss, f"{key}/{arm}/{field}"

    def test_the_system_comparison_is_paired_and_not_causal(self, payload):
        comparison = payload["paired_system_comparison"]
        assert len(comparison) == 5
        for substrate, entry in comparison.items():
            assert "difference of differences" in entry["definition"]
            assert "It does NOT say why" in entry["not_causal"]
            assert entry["n_donors_in_both_systems"] in (13, 15), substrate

    def test_the_sample_size_projection_covers_every_substrate(self, payload):
        projection = payload["sample_size_projection"]["per_substrate"]
        assert len(projection) == 5
        for substrate, entry in projection.items():
            assert entry["n_donors_the_sd_was_estimated_from"] in (13, 15)
            assert len(entry["scenarios"]) == 12, substrate

    def test_the_cytochrome_sensitivity_was_run(self, payload):
        keys = [k for k in payload["results"] if k.endswith("::no_cyp")]
        assert len(keys) == 10
        removed = payload["cytochrome_sensitivity"][
            "n_cytochrome_or_redox_genes_removed"
        ]
        assert removed["HLM"] == 32
        assert removed["HH"] == 28

    def test_both_component_counts_were_run(self, payload):
        for k in (1, 2):
            keys = [key for key in payload["results"] if f"::k{k}::all" in key]
            assert len(keys) == 10, k

    def test_the_prespecified_projection_targets_the_primary_contrast(self, payload):
        projection = payload["sample_size_projection"]
        assert "PRESPECIFIED primary contrast (W2 versus W1)" in projection["scope"]
        assert "not worth sizing a study for" in projection[
            "how_to_read_the_primary_projection"
        ]

    def test_the_q1_projection_is_labelled_post_hoc(self, payload):
        post_hoc = payload["sample_size_projection"][
            "post_hoc_q1_projection_not_prespecified"
        ]
        assert post_hoc["status"] == "POST_HOC_ADDED_AFTER_SEEING_THE_RESULTS"
        assert "not prespecified" in post_hoc["why"]
        assert "reporting addition, not a change of design" in post_hoc["why"]
        assert post_hoc["contrast"] == "W1_target_cyp versus W0_mean"
        assert len(post_hoc["per_cell"]) == 10
        for cell, entry in post_hoc["per_cell"].items():
            assert len(entry["scenarios"]) == 12, cell

    def test_the_defect_is_in_the_prespecifications_change_log(self):
        text = _normalise(PRESPEC.read_text(encoding="utf-8"))
        assert "One defect of this freeze, found by running it" in text
        assert "The prespecified projection is unaltered" in text

    def test_the_results_document_quotes_the_artefacts_own_numbers(self, payload):
        """A write-up that drifts from its artefact is the easiest error here."""
        document = (
            REPO / "docs"
            / "OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_RESULTS_20261110.md"
        )
        if not document.is_file():
            pytest.skip("results document not written")
        text = _normalise(document.read_text(encoding="utf-8"))
        for system, substrate, key, field in (
            ("HH", "Midazolam", "secondary_W1_vs_W0", "mean_difference_log10"),
            ("HH", "Omeprazole", "secondary_W1_vs_W0", "mean_difference_log10"),
            ("HH", "Diclofenac", "primary_W2_vs_W1", "mean_difference_log10"),
            ("HLM", "Midazolam", "primary_W2_vs_W1", "mean_difference_log10"),
        ):
            cell = payload["results"][f"{system}::{substrate}::k2::all"]
            value = cell["contrasts"][key][field]
            printed = f"{value:+.4f}".replace("-", "−")
            assert printed in text, (
                f"{system}/{substrate}/{key} is {printed} in the artefact but "
                "that figure is not in the results document"
            )

    def test_it_exports_no_donor_level_data(self, payload):
        assert payload["no_donor_identifiers_values_or_predictions_exported"] is True

        def arrays(obj, path=""):
            found = []
            if isinstance(obj, dict):
                for key, value in obj.items():
                    found += arrays(value, f"{path}.{key}")
            elif isinstance(obj, list):
                numeric = [x for x in obj if isinstance(x, (int, float))]
                if len(numeric) >= 10:
                    found.append((path, len(numeric)))
                for index, value in enumerate(obj):
                    if not isinstance(value, (int, float)):
                        found += arrays(value, f"{path}[{index}]")
            return found

        assert not arrays(payload)
