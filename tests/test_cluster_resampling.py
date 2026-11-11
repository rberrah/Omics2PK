"""Contracts for the resampling unit, and for the direct ablation contrast.

A reviewer pointed out that this project's folds are grouped by protein
position while its bootstrap resamples individual **substitutions**, of which
up to nineteen share one residue. Grouped folds do not make a grouped
interval, and only the interval decides a verdict — so the published intervals
were too narrow in the one direction that manufactures a verdict. The same
reading found that no contrast measured what the newest knowledge block
contributed; `M4` beating `M1` cannot, because `M4` already contained every
earlier block.

What this file holds, in order of how much damage its absence would do:

1. **A cluster bootstrap widens by the right factor.** On a statistic that is
   constant inside a cluster of `k` rows, the width must grow by about
   `sqrt(k)`, because `k` correlated rows carry the information of one.
   -> `TestClusterBootstrapIsCorrect`.
2. **It is still able to detect a real effect.** Widening every interval until
   nothing is ever significant is not rigour, it is the opposite error.
   -> `TestClusterBootstrapStillDetects`.
3. **The row bootstrap is untouched.** The frozen HLC analyses were published
   with it, and a row there is a donor, which *is* the unit of independence.
   -> `TestRowBootstrapUnchanged`.
4. **The ablation refuses to pair runs that scored different rows.**
   -> `TestAblationPairsOrRefuses`.
5. **One verdict ladder, shared.** The last time a second contrast grew its
   own copy, the copy called +0.62 % a gain against a 5 % margin.
   -> `TestLadderIsShared`.
"""
from __future__ import annotations

import pathlib

import pytest

np = pytest.importorskip("numpy")

ROOT = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def fhd():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import fusion_headroom_diagnostics
    return fusion_headroom_diagnostics


@pytest.fixture(scope="module")
def vms():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import variant_mechanism_omic_selection
    return variant_mechanism_omic_selection


@pytest.fixture(scope="module")
def sa():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import structural_annotation
    return structural_annotation


def _clustered_gap(k, n_clusters=150, seed=0):
    """A gap that is CONSTANT inside each cluster of `k` rows."""
    rng = np.random.default_rng(seed)
    labels = np.repeat(np.arange(n_clusters), k)
    per_cluster = rng.normal(0.02, 0.10, size=n_clusters)
    difference = per_cluster[labels]
    return labels, np.zeros(difference.size), -difference


class TestClusterBootstrapIsCorrect:
    @pytest.mark.parametrize("k", [2, 4, 9])
    def test_width_grows_like_the_square_root_of_the_cluster_size(self, fhd, k):
        labels, reference, ceiling = _clustered_gap(k)
        row = fhd.paired_gap(reference, ceiling)
        clustered = fhd.paired_gap(reference, ceiling, groups=labels)
        width_row = row["interval"][1] - row["interval"][0]
        width_cluster = clustered["interval"][1] - clustered["interval"][0]
        assert width_cluster / width_row == pytest.approx(np.sqrt(k), rel=0.15)

    def test_the_point_estimate_is_untouched(self, fhd):
        labels, reference, ceiling = _clustered_gap(4)
        row = fhd.paired_gap(reference, ceiling)
        clustered = fhd.paired_gap(reference, ceiling, groups=labels)
        assert clustered["gap"] == pytest.approx(row["gap"])

    def test_it_reports_the_unit_and_how_many_there_were(self, fhd):
        labels, reference, ceiling = _clustered_gap(4, n_clusters=150)
        clustered = fhd.paired_gap(reference, ceiling, groups=labels)
        assert "cluster" in clustered["resampling_unit"]
        assert clustered["n_resampling_units"] == 150
        row = fhd.paired_gap(reference, ceiling)
        assert row["resampling_unit"] == "individual"
        assert row["n_resampling_units"] == 600

    def test_three_clusters_is_the_floor(self, fhd):
        difference = np.linspace(0.0, 1.0, 10)
        with pytest.raises(Exception):
            fhd.paired_gap(np.zeros(10), -difference,
                           groups=np.zeros(10, dtype=int))

    def test_a_label_per_row_is_required(self, fhd):
        with pytest.raises(Exception):
            fhd.paired_gap(np.zeros(10), np.zeros(10),
                           groups=np.arange(4))

    def test_singleton_clusters_reproduce_the_row_bootstrap_distribution(
            self, fhd):
        # One row per cluster IS the row bootstrap, so the widths must agree.
        rng = np.random.default_rng(1)
        difference = rng.normal(0.03, 0.2, size=400)
        row = fhd.paired_gap(np.zeros(400), -difference)
        singleton = fhd.paired_gap(np.zeros(400), -difference,
                                   groups=np.arange(400))
        width_row = row["interval"][1] - row["interval"][0]
        width_one = singleton["interval"][1] - singleton["interval"][0]
        assert width_one == pytest.approx(width_row, rel=0.10)


class TestClusterBootstrapStillDetects:
    def test_a_large_real_effect_survives_clustering(self, fhd):
        rng = np.random.default_rng(2)
        labels = np.repeat(np.arange(200), 4)
        per_cluster = rng.normal(0.30, 0.05, size=200)
        difference = per_cluster[labels]
        clustered = fhd.paired_gap(np.zeros(800), -difference, groups=labels)
        assert clustered["excludes_zero"]

    def test_the_positive_control_holds_under_the_clustered_unit(self, vms):
        report = vms.self_test()
        block = report["position_clustered_bootstrap"]
        assert block["M5_vs_M4"] == "GAIN_REACHES_MARGIN"
        assert block["n_clusters"] >= 3

    def test_the_direction_of_the_width_change_is_not_asserted(self, vms):
        # It was, once, and the assertion fired on its own positive control:
        # in that cohort only the mechanism indicator is shared within a
        # position, so the gap has little intra-cluster correlation and the
        # clustered interval comes out slightly narrower. Direction is an
        # empirical matter here, exactly as it was for the arm-selection leak.
        source = (ROOT / "src" / "omics2pk"
                  / "variant_mechanism_omic_selection.py").read_text(
                      encoding="utf-8")
        assert "No assertion on the DIRECTION of the width change" in source


class TestRowBootstrapUnchanged:
    def test_omitting_groups_reproduces_the_published_path_exactly(self, fhd):
        rng = np.random.default_rng(3)
        reference = rng.normal(1.0, 0.3, size=250)
        ceiling = reference - rng.normal(0.05, 0.1, size=250)
        first = fhd.paired_gap(reference, ceiling, seed=11, n_boot=500)
        second = fhd.paired_gap(reference, ceiling, seed=11, n_boot=500,
                                groups=None)
        assert first["interval"] == second["interval"]
        assert first["resampling_unit"] == "individual"


class TestAblationPairsOrRefuses:
    @pytest.fixture(scope="class")
    def two_runs(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            vms.conditionally_gateable_cohort()
        folds = vms.position_groups(positions, 5, 0)
        plain = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                             positions)
        extra = np.column_stack(
            [mechanism, np.random.default_rng(0).normal(size=(truth.size, 2))])
        richer = vms.evaluate(truth, abundance, substitution, extra, folds,
                              positions)
        return plain, richer, positions

    def test_it_returns_a_verdict_from_the_shared_ladder(self, vms, two_runs):
        plain, richer, positions = two_runs
        out = vms.ablation(richer, plain, "M4_fixed_union", label="noise",
                           groups=positions)
        assert out["verdict"] in {
            "GAIN_REACHES_MARGIN", "GAIN_PRESENT_BUT_BELOW_MARGIN",
            "GAIN_PRESENT_MARGIN_UNDECIDED", "DETECTABLY_WORSE",
            "REFUTED_GAIN_BELOW_MARGIN", "NO_DETECTABLE_GAIN_UNDERPOWERED"}
        assert out["added_layer"] == "noise"
        assert out["arm"] == "M4_fixed_union"

    def test_two_columns_of_noise_do_not_earn_a_margin(self, vms, two_runs):
        plain, richer, positions = two_runs
        out = vms.ablation(richer, plain, "M4_fixed_union", label="noise",
                           groups=positions)
        assert out["verdict"] != "GAIN_REACHES_MARGIN"

    def test_it_refuses_runs_that_scored_different_rows(self, vms, two_runs):
        plain, richer, positions = two_runs
        broken = dict(richer)
        mask = np.array(richer["scored_rows"], dtype=bool).copy()
        mask[np.flatnonzero(mask)[0]] = False
        broken["scored_rows"] = mask
        with pytest.raises(Exception) as caught:
            vms.ablation(broken, plain, "M4_fixed_union", label="x")
        assert "different rows" in str(caught.value)

    def test_it_refuses_a_result_without_the_row_mask(self, vms, two_runs):
        plain, richer, _ = two_runs
        stripped = {k: v for k, v in richer.items() if k != "scored_rows"}
        with pytest.raises(Exception):
            vms.ablation(stripped, plain, "M4_fixed_union", label="x")

    def test_it_refuses_an_arm_that_does_not_exist(self, vms, two_runs):
        plain, richer, _ = two_runs
        with pytest.raises(Exception):
            vms.ablation(richer, plain, "M9_invented", label="x")

    def test_ablating_an_identical_block_gives_a_gap_of_zero(self, vms,
                                                             two_runs):
        plain, _, positions = two_runs
        out = vms.ablation(plain, plain, "M4_fixed_union", label="nothing",
                           groups=positions)
        assert out["gap"] == pytest.approx(0.0)
        assert out["verdict"] in {"REFUTED_GAIN_BELOW_MARGIN",
                                  "NO_DETECTABLE_GAIN_UNDERPOWERED"}


class TestLadderIsShared:
    def test_contrast_and_ablation_call_the_same_ladder(self, vms):
        source = (ROOT / "src" / "omics2pk"
                  / "variant_mechanism_omic_selection.py").read_text(
                      encoding="utf-8")
        assert source.count("def _ladder(") == 1
        # Neither contrast nor ablation may spell a rung out for itself.
        assert source.count('"GAIN_REACHES_MARGIN"') == 1 + source.count(
            'verdict"] == "GAIN_REACHES_MARGIN"')

    def test_a_gain_below_the_margin_is_never_called_present(self, vms):
        assert vms._ladder([0.001, 0.002], 0.05) == \
            "GAIN_PRESENT_BUT_BELOW_MARGIN"
        assert vms._ladder([0.06, 0.09], 0.05) == "GAIN_REACHES_MARGIN"
        assert vms._ladder([-0.01, 0.09], 0.05) == \
            "NO_DETECTABLE_GAIN_UNDERPOWERED"
        assert vms._ladder([-0.05, -0.01], 0.05) == "DETECTABLY_WORSE"

    def test_the_margin_is_applied_to_the_lower_bound(self, vms):
        # An upper bound above the margin earns nothing on its own.
        assert vms._ladder([0.001, 0.40], 0.05) != "GAIN_REACHES_MARGIN"


class TestContiguousBlocks:
    def test_neighbouring_residues_share_a_block(self, vms):
        labels = vms.contiguous_blocks([1, 2, 9, 10, 11, 25], 10)
        assert labels[0] == labels[1] == labels[2]
        assert labels[3] == labels[4]
        assert labels[5] not in labels[:5]

    def test_width_one_makes_every_position_its_own_block(self, vms):
        labels = vms.contiguous_blocks([4, 5, 6], 1)
        assert len(set(labels.tolist())) == 3

    def test_a_block_must_be_at_least_one_residue_wide(self, vms):
        with pytest.raises(Exception):
            vms.contiguous_blocks([1, 2, 3], 0)


class TestRatioIntervalSharesOneImplementation:
    def test_the_sequence_and_structural_rulers_use_the_same_routine(self, sa):
        script = (ROOT / "scripts"
                  / "run_variant_mechanism_benchmark.py").read_text(
                      encoding="utf-8")
        assert "SA.discordance_ratio_interval" in script
        source = (ROOT / "src" / "omics2pk"
                  / "structural_annotation.py").read_text(encoding="utf-8")
        assert source.count("def discordance_ratio_interval(") == 1

    def test_clustering_a_position_level_predictor_widens_the_interval(self, sa):
        # `reach` is a function of the POSITION, so every substitution at one
        # residue carries an identical value of the explanatory variable. This
        # is the case where a row bootstrap is most badly wrong.
        rng = np.random.default_rng(4)
        n_pos, per_pos = 120, 8
        labels = np.repeat(np.arange(n_pos), per_pos)
        reach_by_pos = rng.random(n_pos) < 0.25
        reach = reach_by_pos[labels]
        # The outcome has to be correlated WITHIN a position too, or there is
        # nothing for clustering to correct: each residue gets its own
        # propensity, which is what a real residue does.
        propensity = np.where(reach_by_pos, rng.beta(2.0, 2.4, size=n_pos),
                              rng.beta(0.4, 5.0, size=n_pos))
        discordant = rng.random(reach.size) < propensity[labels]
        row = sa.discordance_ratio_interval(reach, discordant, seed=5)
        clustered = sa.discordance_ratio_interval(reach, discordant,
                                                  groups=labels, seed=5)
        assert (clustered[1] - clustered[0]) > (row[1] - row[0])
        assert clustered[0] > 1.0  # a real enrichment is still detected

    def test_too_few_clusters_returns_no_interval_rather_than_a_wrong_one(
            self, sa):
        reach = np.array([True] * 10 + [False] * 10)
        discordant = np.array([True] * 5 + [False] * 15)
        assert all(np.isnan(v) for v in sa.discordance_ratio_interval(
            reach, discordant, groups=np.zeros(20, dtype=int)))

    def test_mismatched_lengths_are_refused(self, sa):
        with pytest.raises(ValueError):
            sa.discordance_ratio_interval(np.array([True, False]),
                                          np.array([True]))


class TestTheBenchmarkDeclaresItsUnit:
    @pytest.fixture(scope="class")
    def script(self):
        return (ROOT / "scripts"
                / "run_variant_mechanism_benchmark.py").read_text(
                    encoding="utf-8")

    def test_every_contrast_in_the_run_loop_passes_a_cluster(self, script):
        body = script[script.index("        contrasts[label] = {"):]
        body = body[:body.index("        }")]
        # The design-effect block deliberately calls the row bootstrap once,
        # to measure what it was doing, so it is excluded from the count.
        reported = body[:body.index("design_effect_on_M4_vs_M1")]
        assert reported.count("VMS.contrast(") >= 6
        assert reported.count("VMS.contrast(") == reported.count(
            "groups=cluster")

    def test_the_design_effect_block_is_the_only_row_wise_call(self, script):
        body = script[script.index("        contrasts[label] = {"):]
        body = body[:body.index("        }")]
        effect = body[body.index("design_effect_on_M4_vs_M1"):]
        assert effect.count("VMS.contrast(") == 2
        assert effect.count("groups=cluster") == 1

    def test_the_design_effect_is_reported_not_asserted(self, script):
        assert "design_effect_on_M4_vs_M1" in script
        assert "verdict_changed" in script

    def test_a_threshold_free_block_is_run(self, script):
        assert "ablation_of_a_threshold_free_block" in script
        assert "CONTINUOUS_COLUMNS = [0, 1, 5]" in script

    def test_the_per_variant_mask_stays_out_of_the_artefact(self, script):
        assert 'if k not in ("losses", "scored_rows")' in script


class TestBrafPreRegistration:
    """The pre-registration is a promise, so its text is a contract too.

    A registration that can be quietly edited after the result is in is worth
    nothing. These assertions pin the parts that make it one: a singular
    primary contrast, refutation conditions, and the record that the primary
    came back refuted rather than being swapped for a secondary.
    """

    @pytest.fixture(scope="class")
    def prereg(self):
        return (ROOT / "docs"
                / "OMICS2PK_BRAF_PREREGISTRATION_20261111.md").read_text(
                    encoding="utf-8")

    def test_it_records_that_the_primary_was_refuted(self, prereg):
        assert "THE PRIMARY CONTRAST IS REFUTED" in prereg
        assert "REFUTED_GAIN_BELOW_MARGIN" in prereg

    def test_the_primary_contrast_is_singular_and_unchanged(self, prereg):
        assert "primary contrast, declared singular" in prereg.lower()
        assert "do not substitute a\n  secondary" in prereg
        assert "the primary contrast was not substituted" in prereg

    def test_it_declared_the_failure_mode_before_the_run(self, prereg):
        # Section 2.1 exists so that "BRAF variants act allosterically" cannot
        # become a post-hoc excuse for a null.
        assert "A reason this may fail that is not about the method" in prereg
        assert "cannot become a post-hoc excuse" in prereg

    def test_every_quantitative_prediction_is_scored(self, prereg):
        assert "Every quantitative prediction in \u00a76.1 failed" in prereg
        for measured in ("0.655", "21.2 %", "\u22120.290", "+1.420"):
            assert measured in prereg

    def test_the_enrichment_is_recorded_as_excluding_one_from_below(
            self, prereg):
        assert "excludes 1 from below" in prereg

    def test_the_post_hoc_control_is_labelled_post_hoc(self, prereg):
        assert "One deviation: a control added after the scores were read" \
            in prereg
        assert "labelled post-hoc throughout" in prereg
        assert "attack an unexpected positive" in prereg

    def test_the_surviving_claim_is_scoped_not_general(self, prereg):
        assert "it does remove the right to generalise" in prereg.lower()
        assert "detects that mismatch automatically" \
            in prereg.lower()

    def test_the_numbering_guard_is_not_the_salt_bridge(self, prereg):
        assert "12.69" in prereg
        assert "2.875" in prereg

    def test_the_results_document_points_at_the_refutation(self):
        results = (ROOT / "docs"
                   / "OMICS2PK_MAVEDB_VARIANT_MECHANISM_RESULTS_20261111.md"
                   ).read_text(encoding="utf-8")
        assert "it refuted the\n  rule" in results
        assert "0.655 [0.42, 0.96]" in results
