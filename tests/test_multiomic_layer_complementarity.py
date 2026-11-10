"""Contract tests for the Axis A masked-modality reconstruction.

What these tests exist to prevent, in order of how much damage each would do:

1. a prediction that saw the held-out individual's own target or its own layer
   statistics, which would make every number meaningless;
2. a complementarity claim that is really a second assay of the same molecules,
   which is the specific error the commissioning note named;
3. a parsimony selection written into the code after looking at results rather
   than returned by the frozen rule;
4. a drift between the prespecification document and the module's freeze;
5. a pharmacological sentence attached to a cohort that has no drug in it.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PRESPEC = (
    REPO
    / "docs"
    / "OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_PRESPECIFICATION_20261110.md"
)
REPORT = REPO / "reports" / "data" / "axis_a_eatris_masked_modality_reconstruction.json"

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")


@pytest.fixture(scope="module")
def mod():
    import sys

    sys.path.insert(0, str(REPO / "src"))
    from omics2pk import multiomic_layer_complementarity as module

    return module


def _normalise(text):
    return re.sub(r"\s+", " ", text)


# --------------------------------------------------------------------------
# 1. the family partition
# --------------------------------------------------------------------------
class TestFamilyPartition:
    def test_the_twelve_layers_are_partitioned_exactly(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        assert len(layers) == 12
        assert len(set(layers)) == 12
        assert mod.assert_exact_partition(layers) is True

    def test_a_missing_layer_is_refused(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        with pytest.raises(mod.LayerComplementarityError):
            mod.assert_exact_partition(layers[:-1])

    def test_an_undeclared_layer_is_refused(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        with pytest.raises(mod.LayerComplementarityError):
            mod.assert_exact_partition(layers + ["spatial-transcriptomics"])

    def test_every_family_is_reachable_from_every_layer(self, mod):
        for family, layers in mod.FAMILIES.items():
            for layer in layers:
                assert mod.family_of(layer) == family

    def test_an_unknown_layer_has_no_family(self, mod):
        with pytest.raises(mod.LayerComplementarityError):
            mod.family_of("not-a-layer")

    def test_adjacency_is_symmetric(self, mod):
        for family, others in mod.BIOCHEMICALLY_ADJACENT.items():
            for other in others:
                assert family in mod.BIOCHEMICALLY_ADJACENT[other]


# --------------------------------------------------------------------------
# 2. the redundancy trap -- the error the commissioning note named
# --------------------------------------------------------------------------
class TestSameFamilySiblingIsKeptSeparable:
    def test_the_primary_target_has_a_sibling_and_it_is_named(self, mod):
        assert mod.siblings_of("proteins") == ("proteomics-ipp",)

    def test_the_primary_arm_excludes_the_sibling(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        arms = mod.arms_for("proteins", layers)
        assert "proteomics-ipp" not in arms[mod.PRIMARY_ARM]
        assert "proteins" not in arms[mod.PRIMARY_ARM]

    def test_the_sibling_arm_is_exactly_the_sibling(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        arms = mod.arms_for("proteins", layers)
        assert arms[mod.COMPLEMENTARITY_COMPARATOR] == ("proteomics-ipp",)

    def test_the_upper_reference_includes_the_sibling_but_never_the_target(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        arms = mod.arms_for("proteins", layers)
        assert "proteomics-ipp" in arms["A3_all_but_target"]
        assert "proteins" not in arms["A3_all_but_target"]

    def test_no_arm_ever_contains_the_target(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        for target in ("proteins", "acylcarnitines", "mRNA", "EM-seq"):
            for name, arm in mod.arms_for(target, layers).items():
                assert target not in arm, f"{name} contains its own target"

    def test_the_secondary_target_declares_its_adjacent_family(self, mod):
        family = mod.family_of("acylcarnitines")
        assert "lipidome" in mod.BIOCHEMICALLY_ADJACENT[family]

    def test_a_target_absent_from_the_configuration_is_refused(self, mod):
        with pytest.raises(mod.LayerComplementarityError):
            mod.arms_for("mRNA", ("proteins", "EM-seq"))

    def test_leave_one_family_out_removes_exactly_that_family(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        arms = mod.arms_for("proteins", layers)
        full = set(arms[mod.PRIMARY_ARM])
        for family in ("transcriptome", "lipidome", "metabolome_targeted"):
            removed = set(arms[f"L_{family}"])
            single = set(arms[f"B_{family}"])
            assert removed == full - single
            assert single and single < full


# --------------------------------------------------------------------------
# 3. fold safety -- nothing crosses from held-out to training
# --------------------------------------------------------------------------
class TestFoldSafety:
    def test_components_ignore_the_held_out_rows_entirely(self, mod):
        rng = np.random.default_rng(7)
        train = rng.normal(size=(40, 25))
        test_a = rng.normal(size=(6, 25))
        test_b = rng.normal(size=(6, 25)) * 50.0 + 100.0
        first, _, _ = mod.layer_components(train, test_a, 4)
        second, _, _ = mod.layer_components(train, test_b, 4)
        # Replacing the held-out block entirely must leave the TRAINING scores
        # bit-for-bit identical: if it does not, the basis saw the test rows.
        assert np.array_equal(first, second)

    def test_held_out_scores_match_a_projection_computed_the_long_way(self, mod):
        rng = np.random.default_rng(11)
        train = rng.normal(size=(30, 12))
        test = rng.normal(size=(5, 12))
        scores_train, scores_test, taken = mod.layer_components(
            train, test, 3, clip=1e9
        )
        mean = train.mean(axis=0)
        sd = train.std(axis=0, ddof=0)
        z_train = (train - mean) / sd
        z_test = (test - mean) / sd
        _, _, right = np.linalg.svd(z_train - z_train.mean(axis=0), full_matrices=False)
        expected_test = (z_test - z_train.mean(axis=0)) @ right[:taken].T
        # Sign of a principal direction is arbitrary, so compare up to sign.
        for column in range(taken):
            assert np.allclose(
                np.abs(scores_test[:, column]), np.abs(expected_test[:, column]),
                atol=1e-8,
            )

    def test_the_first_j_components_are_a_slice_of_a_k_component_extraction(self, mod):
        rng = np.random.default_rng(3)
        train = rng.normal(size=(30, 40))
        test = rng.normal(size=(5, 40))
        wide_train, wide_test, _ = mod.layer_components(train, test, 10)
        narrow_train, narrow_test, _ = mod.layer_components(train, test, 3)
        assert np.allclose(wide_train[:, :3], narrow_train)
        assert np.allclose(wide_test[:, :3], narrow_test)

    def test_a_constant_training_column_is_dropped(self, mod):
        train = np.ones((20, 4))
        train[:, 1] = np.arange(20)
        test = np.ones((4, 4))
        test[:, 1] = np.arange(4)
        test[:, 2] = 99.0  # constant in train, different in test: must not enter
        scores_train, scores_test, taken = mod.layer_components(train, test, 3)
        assert taken == 1
        assert scores_train.shape[1] == 1 and scores_test.shape[1] == 1

    def test_a_wholly_constant_layer_contributes_nothing(self, mod):
        train = np.full((15, 6), 2.5)
        test = np.full((3, 6), 9.0)
        scores_train, scores_test, taken = mod.layer_components(train, test, 4)
        assert taken == 0
        assert scores_train.shape == (15, 0)
        assert scores_test.shape == (3, 0)

    def test_a_missing_cell_in_a_layer_block_is_refused(self, mod):
        train = np.ones((10, 3))
        train[0, 0] = np.nan
        with pytest.raises(mod.LayerComplementarityError):
            mod.layer_components(train, np.ones((2, 3)), 2)

    def test_the_clip_uses_training_units_only(self, mod):
        train = np.concatenate([np.zeros((20, 1)), np.ones((20, 1))], axis=0)
        train = np.column_stack([train.ravel(), train.ravel()])
        test = np.array([[1000.0, 1000.0]])
        _, scores_test, _ = mod.layer_components(train, test, 1, clip=4.0)
        # A held-out value a thousand SDs out cannot produce a score that
        # scales with it: the clip bounds the standardised block first.
        assert np.all(np.abs(scores_test) < 100.0)

    def test_design_standardisation_uses_training_constants(self, mod):
        train = np.array([[1.0, 10.0], [3.0, 30.0], [5.0, 50.0]])
        test = np.array([[7.0, 70.0]])
        z_train, z_test = mod.standardise_design(train, test)
        assert np.allclose(z_train.mean(axis=0), 0.0)
        assert np.allclose(z_train.std(axis=0, ddof=0), 1.0)
        # The held-out row is outside the training range and must stay so.
        assert np.all(z_test > z_train.max(axis=0))

    def test_a_constant_design_column_does_not_divide_by_zero(self, mod):
        train = np.column_stack([np.ones(5), np.arange(5.0)])
        test = np.column_stack([np.ones(2), np.arange(2.0)])
        z_train, z_test = mod.standardise_design(train, test)
        assert np.isfinite(z_train).all() and np.isfinite(z_test).all()


# --------------------------------------------------------------------------
# 4. the estimator
# --------------------------------------------------------------------------
class TestRidgePerTarget:
    def test_each_target_gets_its_own_penalty(self, mod):
        rng = np.random.default_rng(5)
        design = rng.normal(size=(60, 6))
        signal = design @ rng.normal(size=(6, 1))
        noise = rng.normal(size=(60, 1)) * 30.0
        targets = np.column_stack([signal.ravel(), noise.ravel()])
        _, chosen, _ = mod.ridge_per_target(
            design, targets, design[:5], alphas=(0.01, 10000.0), inner_folds=3
        )
        assert chosen.shape == (2,)
        # The pure-noise target should take the heavier penalty.
        assert chosen[1] >= chosen[0]

    def test_predictions_are_on_the_original_target_scale(self, mod):
        rng = np.random.default_rng(13)
        design = rng.normal(size=(50, 4))
        targets = 1000.0 + design @ rng.normal(size=(4, 3)) * 5.0
        predictions, _, mean = mod.ridge_per_target(design, targets, design[:8])
        assert predictions.shape == (8, 3)
        assert np.all(predictions > 500.0)
        assert np.allclose(mean, targets.mean(axis=0))

    def test_the_returned_mean_is_the_training_mean(self, mod):
        rng = np.random.default_rng(17)
        design = rng.normal(size=(30, 3))
        targets = rng.normal(size=(30, 2)) + 7.0
        _, _, mean = mod.ridge_per_target(design, targets, design[:4])
        assert np.allclose(mean, targets.mean(axis=0))

    def test_an_enormous_penalty_returns_the_training_mean(self, mod):
        rng = np.random.default_rng(19)
        design = rng.normal(size=(40, 5))
        targets = rng.normal(size=(40, 3)) + 2.0
        predictions, _, mean = mod.ridge_per_target(
            design, targets, design[:6], alphas=(1e12,), inner_folds=3
        )
        assert np.allclose(predictions, mean[None, :], atol=1e-3)

    def test_mismatched_rows_are_refused(self, mod):
        with pytest.raises(mod.LayerComplementarityError):
            mod.ridge_per_target(np.ones((10, 2)), np.ones((9, 1)), np.ones((2, 2)))

    def test_a_missing_training_target_is_refused(self, mod):
        targets = np.ones((10, 1))
        targets[3, 0] = np.nan
        with pytest.raises(mod.LayerComplementarityError):
            mod.ridge_per_target(np.ones((10, 2)), targets, np.ones((2, 2)))

    def test_a_recoverable_signal_is_actually_recovered(self, mod):
        """A sanity floor: if this fails the estimator cannot learn at all."""
        rng = np.random.default_rng(23)
        design = rng.normal(size=(80, 4))
        weights = np.array([[2.0], [-1.0], [0.5], [0.0]])
        targets = design @ weights + rng.normal(size=(80, 1)) * 0.1
        predictions, _, _ = mod.ridge_per_target(
            design[:60], targets[:60], design[60:]
        )
        correlation = np.corrcoef(predictions.ravel(), targets[60:].ravel())[0, 1]
        assert correlation > 0.9


# --------------------------------------------------------------------------
# 5. scoring and intervals
# --------------------------------------------------------------------------
class TestScoring:
    def test_perfect_prediction_gives_r2_of_one(self, mod):
        error = np.zeros((3, 10))
        total = np.ones((3, 10))
        assert np.allclose(mod.r2_per_feature(error, total), 1.0)

    def test_predicting_the_training_mean_gives_r2_of_zero(self, mod):
        total = np.full((4, 12), 2.0)
        assert np.allclose(mod.r2_per_feature(total.copy(), total), 0.0)

    def test_worse_than_the_mean_gives_a_negative_r2(self, mod):
        total = np.ones((2, 8))
        error = np.full((2, 8), 3.0)
        assert np.all(mod.r2_per_feature(error, total) < 0.0)

    def test_a_zero_variance_feature_returns_not_a_number(self, mod):
        total = np.zeros((1, 5))
        error = np.zeros((1, 5))
        assert np.isnan(mod.r2_per_feature(error, total)).all()

    def test_shape_disagreement_is_refused(self, mod):
        with pytest.raises(mod.LayerComplementarityError):
            mod.r2_per_feature(np.ones((2, 5)), np.ones((2, 4)))

    def test_the_bootstrap_resamples_individuals_not_features(self, mod):
        rng = np.random.default_rng(29)
        total = np.full((6, 40), 1.0)
        better = rng.uniform(0.1, 0.4, size=(6, 40))
        worse = rng.uniform(0.6, 0.9, size=(6, 40))
        out = mod.paired_bootstrap_median_delta(better, worse, total, n_boot=200)
        assert out["resampling_unit"] == "individual"
        assert out["n_independent_individuals"] == 40
        assert out["median_delta_r2"] > 0.0
        assert out["interval_excludes_zero"] is True

    def test_identical_arms_give_a_zero_difference(self, mod):
        total = np.full((5, 30), 2.0)
        error = np.full((5, 30), 1.0)
        out = mod.paired_bootstrap_median_delta(error, error.copy(), total, n_boot=100)
        assert out["median_delta_r2"] == 0.0
        assert out["interval_excludes_zero"] is False

    def test_the_bootstrap_is_deterministic_under_the_frozen_seed(self, mod):
        total = np.full((4, 25), 1.5)
        a = np.random.default_rng(1).uniform(size=(4, 25))
        b = np.random.default_rng(2).uniform(size=(4, 25))
        first = mod.paired_bootstrap_median_delta(a, b, total, n_boot=150)
        second = mod.paired_bootstrap_median_delta(a, b, total, n_boot=150)
        assert first["descriptive_95_interval"] == second["descriptive_95_interval"]

    def test_a_non_significant_interval_is_not_called_equivalence(self, mod):
        total = np.full((3, 20), 1.0)
        error = np.full((3, 20), 0.5)
        out = mod.paired_bootstrap_median_delta(error, error.copy(), total, n_boot=100)
        assert "never that there is none" in out["how_to_read"]

    def test_mismatched_blocks_are_refused(self, mod):
        with pytest.raises(mod.LayerComplementarityError):
            mod.paired_bootstrap_median_delta(
                np.ones((2, 5)), np.ones((2, 5)), np.ones((2, 6)), n_boot=10
            )


class TestPermutationRank:
    def test_exceeding_every_replicate_gives_the_smallest_possible_p(self, mod):
        out = mod.permutation_rank(1.0, [0.1] * 20)
        assert out["exceeds_every_replicate"] is True
        assert out["rank_p_value"] == pytest.approx(1 / 21)

    def test_landing_inside_the_null_is_reported_as_such(self, mod):
        out = mod.permutation_rank(0.5, [0.0, 0.4, 0.6, 1.0])
        assert out["exceeds_every_replicate"] is False
        assert out["n_null_at_least_observed"] == 2
        assert out["rank_p_value"] == pytest.approx(3 / 5)

    def test_the_null_range_is_reported(self, mod):
        out = mod.permutation_rank(0.0, [-1.0, 0.0, 2.0])
        assert out["null_min"] == -1.0
        assert out["null_max"] == 2.0


# --------------------------------------------------------------------------
# 6. parsimony must be produced by the rule, not written down
# --------------------------------------------------------------------------
class TestParsimony:
    def test_the_enumeration_covers_every_subset(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        subsets = mod.family_subsets("proteins", layers)
        # Five non-target families in the full twelve-layer set.
        assert len(subsets) == 2 ** 5
        assert "S_none" in subsets
        assert subsets["S_none"] == ((), ())

    def test_every_subset_excludes_the_target_and_its_family(self, mod):
        layers = sorted(x for group in mod.FAMILIES.values() for x in group)
        for families, members in mod.family_subsets("proteins", layers).values():
            assert "proteins" not in members
            assert "proteomics-ipp" not in members
            assert "proteome" not in families

    def test_the_fewest_families_wins_within_the_tolerance(self, mod):
        table = {
            "S_a": {"families": ("lipidome",), "n_features": 600,
                    "median_delta_r2": 0.090},
            "S_ab": {"families": ("lipidome", "transcriptome"),
                     "n_features": 17000, "median_delta_r2": 0.095},
        }
        out = mod.parsimonious_subset(table, tolerance=0.01)
        assert out["selected"] == "S_a"
        assert out["n_families"] == 1

    def test_a_gain_beyond_the_tolerance_is_not_given_up(self, mod):
        table = {
            "S_a": {"families": ("lipidome",), "n_features": 600,
                    "median_delta_r2": 0.05},
            "S_ab": {"families": ("lipidome", "transcriptome"),
                     "n_features": 17000, "median_delta_r2": 0.20},
        }
        out = mod.parsimonious_subset(table, tolerance=0.01)
        assert out["selected"] == "S_ab"

    def test_ties_on_family_count_are_broken_by_features(self, mod):
        table = {
            "S_big": {"families": ("transcriptome",), "n_features": 16386,
                      "median_delta_r2": 0.10},
            "S_small": {"families": ("metabolome_targeted",), "n_features": 79,
                        "median_delta_r2": 0.10},
        }
        out = mod.parsimonious_subset(table, tolerance=0.01)
        assert out["selected"] == "S_small"

    def test_the_empty_subset_can_win_when_nothing_helps(self, mod):
        table = {
            "S_none": {"families": (), "n_features": 0, "median_delta_r2": 0.0},
            "S_a": {"families": ("lipidome",), "n_features": 600,
                    "median_delta_r2": 0.002},
        }
        out = mod.parsimonious_subset(table, tolerance=0.01)
        assert out["selected"] == "S_none"
        assert out["n_families"] == 0

    def test_the_rule_and_its_optimism_are_reported(self, mod):
        table = {"S_none": {"families": (), "n_features": 0,
                            "median_delta_r2": 0.0}}
        out = mod.parsimonious_subset(table)
        assert out["performance_of_the_selected_subset_is_optimistically_biased"] is True
        assert "fewest families" in out["rule"]
        assert "unbiased" in out["why"]

    def test_an_empty_enumeration_is_refused(self, mod):
        with pytest.raises(mod.LayerComplementarityError):
            mod.parsimonious_subset({})

    def test_no_subset_name_is_hard_coded_as_the_winner(self):
        """The selection must come out of the rule, not out of the source."""
        source = (REPO / "src" / "omics2pk" / "multiomic_layer_complementarity.py").read_text(
            encoding="utf-8"
        )
        runner = (
            REPO / "scripts" / "run_eatris_plus_masked_modality_reconstruction.py"
        ).read_text(encoding="utf-8")
        for text in (source, runner):
            assert '"selected": "S_' not in text
            assert "'selected': 'S_" not in text


# --------------------------------------------------------------------------
# 7. scale and missingness rules
# --------------------------------------------------------------------------
class TestScaleAndMissingness:
    def test_the_log_rule_fires_on_a_wide_positive_layer(self, mod):
        assert mod.log_rule_applies(np.array([[1.0, 2.0], [5000.0, 9000.0]])) is True

    def test_the_log_rule_does_not_fire_on_a_narrow_layer(self, mod):
        assert mod.log_rule_applies(np.array([[11.9, 20.0], [30.0, 33.6]])) is False

    def test_the_log_rule_never_fires_on_a_layer_with_negatives(self, mod):
        wide = np.array([[-4.2, 1.0], [2.0, 24000.0]])
        assert mod.log_rule_applies(wide) is False

    def test_an_all_missing_layer_does_not_fire_the_log_rule(self, mod):
        assert mod.log_rule_applies(np.full((3, 2), np.nan)) is False

    def test_columns_with_any_missing_cell_are_dropped(self, mod):
        block = np.ones((5, 4))
        block[2, 1] = np.nan
        block[0, 3] = np.nan
        mask = mod.complete_columns(block)
        assert mask.tolist() == [True, False, True, False]

    def test_the_screen_reads_missingness_and_not_values(self, mod):
        first = np.ones((4, 3))
        second = np.arange(12.0).reshape(4, 3)
        first[1, 2] = np.nan
        second[1, 2] = np.nan
        assert mod.complete_columns(first).tolist() == mod.complete_columns(second).tolist()


# --------------------------------------------------------------------------
# 8. the freeze must agree with the document
# --------------------------------------------------------------------------
class TestFreezeAgreesWithTheDocument:
    @pytest.fixture(scope="class")
    def text(self):
        assert PRESPEC.is_file(), f"prespecification absent: {PRESPEC}"
        return PRESPEC.read_text(encoding="utf-8")

    def test_the_document_exists_and_declares_its_freeze_date(self, text):
        assert "frozen 2026-11-10, before any model was fitted" in _normalise(text)

    def test_the_seed_matches(self, mod, text):
        assert f"| Fold seed | `{mod.FROZEN['seed']}` |" in text

    def test_the_fold_structure_matches(self, mod, text):
        expected = (
            f"| Outer folds × repeats × inner | {mod.FROZEN['outer_folds']} × "
            f"{mod.FROZEN['outer_repeats']} × {mod.FROZEN['inner_folds']} |"
        )
        assert expected in text

    def test_the_component_count_matches(self, mod, text):
        assert f"| Components per layer (k) | {mod.FROZEN['components_per_layer']} |" in text

    def test_the_permutation_count_matches(self, mod, text):
        assert (
            f"| Permutation replicates | {mod.FROZEN['permutation_replicates']} |"
            in text
        )

    def test_the_parsimony_tolerance_matches(self, mod, text):
        assert f"| Parsimony tolerance (median ΔR²) | {mod.FROZEN['parsimony_tolerance']} |" in text

    def test_the_targets_and_arms_match(self, mod, text):
        assert f"| Primary target | `{mod.FROZEN['primary_target']}` |" in text
        assert f"| Secondary target | `{mod.FROZEN['secondary_target']}` |" in text
        assert f"| Primary arm | `{mod.PRIMARY_ARM}` |" in text
        assert f"| Primary comparator | `{mod.PRIMARY_COMPARATOR}` |" in text
        assert (
            f"| Complementarity comparator | `{mod.COMPLEMENTARITY_COMPARATOR}` |"
            in text
        )

    def test_the_status_ceiling_matches(self, mod, text):
        assert f"| Status ceiling | `{mod.FROZEN['status_ceiling']}` |" in text

    def test_every_configuration_checksum_is_recorded(self, mod, text):
        for name, spec in mod.CONFIGURATIONS.items():
            assert spec["individual_set_sha256_prefix"] in text, name

    def test_the_commissioned_counts_are_recorded(self, mod, text):
        normalised = _normalise(text)
        for n in (79, 55, 114):
            assert f"**{n}**" in normalised

    def test_the_primary_contrast_is_stated_as_prose(self, mod, text):
        assert f"`{mod.PRIMARY_ARM}` versus `{mod.PRIMARY_COMPARATOR}`" in text

    def test_the_runner_reads_the_document_back(self):
        runner = (
            REPO / "scripts" / "run_eatris_plus_masked_modality_reconstruction.py"
        ).read_text(encoding="utf-8")
        assert "check_freeze_against_document" in runner
        assert "the document says" in runner


# --------------------------------------------------------------------------
# 9. the prespecification's own content
# --------------------------------------------------------------------------
class TestPrespecificationContent:
    @pytest.fixture(scope="class")
    def text(self):
        return _normalise(PRESPEC.read_text(encoding="utf-8"))

    def test_it_records_the_unfavourable_prior(self, text):
        assert "The prior is unfavourable and is recorded as such" in text

    def test_it_forbids_a_pharmacological_statement(self, text):
        assert "**Not pharmacological.**" in text
        assert "No pharmacological, enzymatic or pharmacokinetic statement" in text

    def test_it_fixes_the_failure_sentence_before_the_result(self, text):
        assert (
            "attributable to a second assay of the same molecules, not to "
            "inter-family information" in text
        )

    def test_it_forbids_the_cross_configuration_comparison(self, text):
        assert "No paired contrast crosses a configuration" in text
        assert (
            "performance differences between distinct sample sizes are not "
            "modality effects" in text.lower()
        )

    def test_it_declares_the_parsimony_optimism(self, text):
        assert "optimistically biased" in text

    def test_it_declares_the_no_imputation_rule_for_the_primary_contrast(self, text):
        assert "there is no imputation inside a layer" in text

    def test_it_records_why_haematology_is_out_of_the_baseline(self, text):
        assert "complete, so the primary contrast needs no imputation" in text

    def test_it_declares_the_measured_log_rule_targets(self, text):
        assert "`lipidomics-neg` and `lipidomics-pos`" in text

    def test_it_states_the_three_separate_reports(self, text):
        assert "Predictive performance" in text
        assert "Biological relevance" in text
        assert "Analytical validity" in text

    def test_it_lists_what_may_not_be_done_afterwards(self, text):
        assert "What may not be done after seeing the results" in text

    def test_its_change_log_records_the_pre_freeze_correction(self, text):
        assert "before the freeze and before any model was fitted" in text


# --------------------------------------------------------------------------
# 10. no clinical or pharmacological claim may attach to this cohort
# --------------------------------------------------------------------------
class TestNoPharmacologicalClaim:
    def test_the_module_carries_the_disclaimer(self, mod):
        sentence = mod.no_pharmacological_claim()
        for phrase in ("no drug", "pharmacokinetics", "no hepatic interpretation"):
            assert phrase in sentence

    def test_the_module_mentions_no_pharmacokinetic_target(self, mod):
        source = (REPO / "src" / "omics2pk" / "multiomic_layer_complementarity.py").read_text(
            encoding="utf-8"
        )
        for forbidden in ("clint", "auc", "clearance_target", "trough", "cmin"):
            assert forbidden not in source.lower().replace(
                "intrinsic clearance", ""
            ).replace("drug metabolism", "")


# --------------------------------------------------------------------------
# 11. the committed report, when it exists
# --------------------------------------------------------------------------
class TestCommittedReport:
    @pytest.fixture(scope="class")
    def payload(self):
        if not REPORT.is_file():
            pytest.skip("Axis A report not generated")
        return json.loads(REPORT.read_text(encoding="utf-8"))

    def test_it_declares_the_status_ceiling(self, payload, mod):
        assert payload["status_ceiling"] == mod.FROZEN["status_ceiling"]

    def test_it_carries_the_no_pharmacology_sentence(self, payload):
        assert "pharmacokinetics" in payload["no_pharmacological_claim"]

    def test_it_exports_no_identifier_value_or_prediction(self, payload):
        assert payload["no_identifier_value_or_prediction_exported"] is True

    def test_every_configuration_reached_its_commissioned_size(self, payload, mod):
        for name, spec in mod.CONFIGURATIONS.items():
            assert payload["configurations"][name]["n_individuals"] == spec["n"]
            assert payload["configurations"][name][
                "individual_set_sha256"
            ].startswith(spec["individual_set_sha256_prefix"])

    def test_the_log_rule_fired_only_on_the_two_declared_layers(self, payload):
        assert payload["log10_applied_to"] == ["lipidomics-neg", "lipidomics-pos"]

    def test_the_primary_cell_carries_a_permutation_control(self, payload, mod):
        key = f"{mod.FROZEN['primary_configuration']}::{mod.FROZEN['primary_target']}"
        control = payload["results"][key]["permutation_control"]
        assert control["n_replicates"] == mod.FROZEN["permutation_replicates"]
        assert control["rank_p_value"] is not None

    def test_every_cell_carries_a_permutation_control(self, payload):
        for key, cell in payload["results"].items():
            assert "permutation_control" in cell, key

    def test_the_verdict_is_one_of_the_three_prespecified_strings(self, payload):
        allowed = {
            "SUPPORTED_WITH_COMPLEMENTARITY",
            "SUPPORTED_BUT_ATTRIBUTABLE_TO_THE_SAME_FAMILY_SIBLING",
            "NOT_SUPPORTED",
        }
        for key, cell in payload["results"].items():
            assert cell["prespecified_decision"]["verdict"] in allowed, key

    def test_no_arm_in_the_report_contains_its_own_target(self, payload):
        for key, cell in payload["results"].items():
            target = cell["target_layer"]
            for name, arm in cell["arms"].items():
                assert target not in arm["layers"], f"{key}:{name}"

    def test_the_primary_arm_in_the_report_excludes_the_sibling(self, payload, mod):
        for key, cell in payload["results"].items():
            siblings = set(cell["target_same_family_sibling"])
            arm = set(cell["arms"][mod.PRIMARY_ARM]["layers"])
            assert not (arm & siblings), key

    def test_the_forbidden_comparison_is_written_down(self, payload):
        text = payload["results"]["C9_no_rna::proteins"]["forbidden_comparison"]
        assert "not modality effects" in text

    def test_a_decision_is_never_claimed_without_its_permutation_criterion(self, payload):
        for key, cell in payload["results"].items():
            decision = cell["prespecified_decision"]
            if decision["primary_claim_supported"]:
                assert decision["criteria"]["exceeds_every_permutation_replicate"] is True, key
