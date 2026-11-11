"""Contracts for the variant -> mechanism -> measurement-selection benchmark.

The guarantees this module exists to hold, in order of how much damage their
absence would do:

1. **A protein position is never split across folds.** Four substitutions at
   one residue are not four independent observations, and a random split would
   train on `p.Arg34Ala` and test on `p.Arg34Gly`. -> `TestGroupedSplits`.
2. **The harness can see a conditional rule that exists**, so a null means
   something. -> `TestPositiveControl`.
3. **Mechanism hypotheses are ranked candidates, never labels.**
   -> `TestMechanismRulesPropose`.
4. **The conditional arm collapses exactly to the fixed union** when the
   conditional part carries nothing, so `M5` vs `M4` isolates that part.
   -> `TestConditionalArmDegenerates`.
5. The verdict ladder applies its margin to the lower bound, as everywhere
   else in this project. -> `TestVerdictLadder`.
"""
from __future__ import annotations

import pathlib

import pytest

np = pytest.importorskip("numpy")

ROOT = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def vms():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import variant_mechanism_omic_selection
    return variant_mechanism_omic_selection


@pytest.fixture(scope="module")
def annotation():
    return {"length": 164, "functional_sites": (47, 49, 63, 67),
            "motifs": ((48, 69),)}


# --------------------------------------------------------------------------
# 1. parsing
# --------------------------------------------------------------------------
class TestParsing:
    def test_a_missense_variant(self, vms):
        assert vms.parse_protein_variant("p.Arg34Ala") == ("R", 34, "A")

    def test_nonsense_keeps_the_position_and_marks_the_stop(self, vms):
        assert vms.parse_protein_variant("p.Gly17Ter") == ("G", 17, None)

    def test_synonymous_is_marked_rather_than_silently_parsed(self, vms):
        assert vms.parse_protein_variant("p.Leu10=") == ("L", 10, "=")

    @pytest.mark.parametrize("bad", ["c.101A>G", "p.Xyz12Ala", "Arg34Ala",
                                     "", None, "p.Arg34delins", "NA"])
    def test_anything_else_is_refused_not_guessed(self, vms, bad):
        assert vms.parse_protein_variant(bad) is None

    def test_a_three_digit_position_parses(self, vms):
        assert vms.parse_protein_variant("p.Asn117Arg") == ("N", 117, "R")


# --------------------------------------------------------------------------
# 2. the mechanism rules PROPOSE
# --------------------------------------------------------------------------
class TestMechanismRulesPropose:
    def test_a_ranked_list_is_returned_never_a_single_label(self, vms,
                                                            annotation):
        out = vms.mechanism_hypotheses("R", 47, "A", annotation)
        assert isinstance(out, list) and out
        assert {"mechanism", "prior_weight", "argument"} <= set(out[0])

    def test_hypotheses_are_sorted_by_prior_weight(self, vms, annotation):
        out = vms.mechanism_hypotheses("L", 15, "P", annotation)
        weights = [row["prior_weight"] for row in out]
        assert weights == sorted(weights, reverse=True)

    def test_a_variant_at_a_binding_site_argues_for_activity(self, vms,
                                                             annotation):
        out = vms.mechanism_hypotheses("R", 47, "A", annotation)
        assert out[0]["mechanism"] == "functional_site_proximal"
        assert "ACTIVITY" in out[0]["argument"]

    def test_a_nonsense_variant_argues_that_neither_assay_is_better(
            self, vms, annotation):
        out = vms.mechanism_hypotheses("G", 17, None, annotation)
        assert out[0]["mechanism"] == "nonsense_truncation"
        assert "neither" in out[0]["argument"]

    def test_a_destabilising_change_argues_for_abundance(self, vms,
                                                         annotation):
        out = vms.mechanism_hypotheses("I", 120, "D", annotation)
        names = [row["mechanism"] for row in out]
        assert "buried_core_candidate" in names
        argument = next(r["argument"] for r in out
                        if r["mechanism"] == "buried_core_candidate")
        assert "ABUNDANCE" in argument

    def test_an_unannotatable_variant_says_so_rather_than_inventing(
            self, vms):
        sparse = {"length": 490, "functional_sites": (435,), "motifs": ()}
        out = vms.mechanism_hypotheses("A", 250, "S", sparse)
        assert out
        assert any("does not support" in row["argument"]
                   or row["prior_weight"] <= 0.35 for row in out)

    def test_every_emitted_mechanism_is_declared(self, vms, annotation):
        seen = set()
        for position in range(1, 165, 7):
            for replacement in ("A", "D", "P", "W", None):
                for row in vms.mechanism_hypotheses("L", position, replacement,
                                                    annotation):
                    seen.add(row["mechanism"])
        assert seen <= set(vms.MECHANISM_CLASSES), (
            f"undeclared mechanism emitted: {seen - set(vms.MECHANISM_CLASSES)}")


# --------------------------------------------------------------------------
# 3. the grouped splits — the guarantee that matters most
# --------------------------------------------------------------------------
class TestGroupedSplits:
    def test_a_position_never_appears_in_two_folds(self, vms):
        positions = np.repeat(np.arange(1, 61), 4)
        folds = vms.position_groups(positions, 5, 7)
        for position in np.unique(positions):
            assigned = np.unique(folds[positions == position])
            assert assigned.size == 1, (
                f"position {position} was split across folds {assigned}; four "
                "substitutions at one residue are not independent")

    def test_every_fold_is_used(self, vms):
        positions = np.repeat(np.arange(1, 61), 3)
        folds = vms.position_groups(positions, 5, 7)
        assert np.unique(folds).size == 5

    def test_region_blocks_are_contiguous(self, vms):
        positions = np.arange(1, 165)
        folds = vms.region_groups(positions, 8, 164)
        for block in np.unique(folds):
            inside = positions[folds == block]
            assert inside.max() - inside.min() + 1 == inside.size, (
                "a region block must be a contiguous stretch of sequence")

    def test_region_blocks_cover_the_sequence_in_order(self, vms):
        positions = np.arange(1, 165)
        folds = vms.region_groups(positions, 8, 164)
        assert list(folds) == sorted(folds), (
            "region blocks must increase along the sequence")

    def test_the_region_split_is_the_harder_one(self, vms):
        """A region split puts whole neighbourhoods out of sample."""
        positions = np.arange(1, 165)
        region = vms.region_groups(positions, 8, 164)
        position = vms.position_groups(positions, 8, 7)
        held_region = positions[region == 0]
        held_position = positions[position == 0]
        assert held_region.std() < held_position.std(), (
            "region-held positions should be clustered, position-held ones "
            "scattered; otherwise the two splits are not different tests")


# --------------------------------------------------------------------------
# 4. features
# --------------------------------------------------------------------------
class TestFeatures:
    def test_substitution_features_have_the_declared_width(self, vms):
        row = vms.substitution_features("R", 34, "A")
        assert len(row) == len(vms.SUBSTITUTION_NAMES)

    def test_a_nonsense_variant_is_flagged_in_the_substitution_block(self, vms):
        row = vms.substitution_features("R", 34, None)
        assert row[vms.SUBSTITUTION_NAMES.index("is_nonsense")] == 1.0

    def test_mechanism_features_have_the_declared_width(self, vms,
                                                        annotation):
        row = vms.mechanism_features("R", 47, "A", annotation)
        assert len(row) == len(vms.MECHANISM_NAMES)

    def test_a_site_adjacent_variant_is_marked(self, vms, annotation):
        row = vms.mechanism_features("R", 48, "A", annotation)
        assert row[vms.MECHANISM_NAMES.index("within_3_of_site")] == 1.0
        assert row[vms.MECHANISM_NAMES.index("inside_motif")] == 1.0

    def test_a_distant_variant_is_not_marked(self, vms, annotation):
        row = vms.mechanism_features("L", 130, "A", annotation)
        assert row[vms.MECHANISM_NAMES.index("within_3_of_site")] == 0.0
        assert row[vms.MECHANISM_NAMES.index("inside_motif")] == 0.0

    def test_build_design_returns_two_aligned_blocks(self, vms, annotation):
        variants = [("R", 34, "A"), ("G", 17, None), ("L", 130, "P")]
        substitution, mechanism = vms.build_design(variants, annotation)
        assert substitution.shape == (3, len(vms.SUBSTITUTION_NAMES))
        assert mechanism.shape == (3, len(vms.MECHANISM_NAMES))

    def test_the_blocks_are_finite_everywhere(self, vms, annotation):
        variants = [("L", p, r) for p in range(1, 165, 3)
                    for r in ("A", "D", "W", None)]
        substitution, mechanism = vms.build_design(variants, annotation)
        assert np.isfinite(substitution).all() and np.isfinite(mechanism).all()


# --------------------------------------------------------------------------
# 5. the positive control
# --------------------------------------------------------------------------
class TestPositiveControl:
    def test_the_self_test_passes(self, vms):
        report = vms.self_test()
        assert report["with_real_mechanism"]["M5_vs_M4"] == "GAIN_REACHES_MARGIN"
        assert not report["with_permuted_mechanism"]["M5_vs_M4"].startswith(
            "GAIN_REACHES")

    def test_the_interpretable_gate_is_recorded_as_the_weaker_form(self, vms):
        """It blends predictors already pooled across regimes, and loses."""
        report = vms.self_test()
        assert (report["with_real_mechanism"]["relative_gain"]
                > report["with_real_mechanism"]["M6_relative_gain"])

    def test_the_trust_gate_actually_varies_its_weight(self, vms):
        report = vms.self_test()
        assert report["with_real_mechanism"]["blend_spread"] > 0.01


# --------------------------------------------------------------------------
# 6. the conditional arm collapses to the fixed union when it should
# --------------------------------------------------------------------------
class TestConditionalArmDegenerates:
    @staticmethod
    def _cohort(seed, conditional):
        rng = np.random.default_rng(seed)
        positions = np.repeat(np.arange(1, 121), 4)
        n = positions.size
        truth = rng.normal(size=n)
        flag = (positions <= 60).astype(float)
        if conditional:
            abundance = np.where(flag > 0, rng.normal(size=n),
                                 truth + rng.normal(scale=0.2, size=n))
        else:
            abundance = truth + rng.normal(scale=0.2, size=n)
        substitution = rng.normal(size=(n, 3))
        mechanism = np.column_stack([flag, rng.normal(size=(n, 2))])
        return truth, abundance, substitution, mechanism, positions

    def test_with_no_conditional_structure_M5_does_not_beat_M4(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(11, conditional=False)
        folds = vms.position_groups(positions, 5, vms.FROZEN["seed"])
        result = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                              positions)
        cell = vms.contrast(result, "M5_conditional_interaction",
                            "M4_fixed_union")
        assert not cell["verdict"].startswith("GAIN_REACHES"), (
            "M5 must not win where abundance is uniformly informative; it "
            "would mean the extra columns are buying flexibility, not "
            "conditional knowledge")

    def test_with_conditional_structure_M5_does_beat_M4(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(11, conditional=True)
        folds = vms.position_groups(positions, 5, vms.FROZEN["seed"])
        result = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                              positions)
        cell = vms.contrast(result, "M5_conditional_interaction",
                            "M4_fixed_union")
        assert cell["relative_gain"] > 0.0

    def test_permuting_the_mechanism_removes_the_advantage(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(11, conditional=True)
        folds = vms.position_groups(positions, 5, vms.FROZEN["seed"])
        result = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                              positions, permute_mechanism=True)
        cell = vms.contrast(result, "M5_conditional_interaction",
                            "M4_fixed_union")
        assert not cell["verdict"].startswith("GAIN_REACHES")

    def test_permuting_abundance_degrades_the_abundance_arm(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(12, conditional=False)
        folds = vms.position_groups(positions, 5, vms.FROZEN["seed"])
        real = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                            positions)
        dead = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                            positions, permute_abundance=True)
        assert (dead["mean_loss"]["M1_abundance_only"]
                > real["mean_loss"]["M1_abundance_only"])

    def test_all_declared_comparators_are_scored(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(13, conditional=True)
        folds = vms.position_groups(positions, 5, vms.FROZEN["seed"])
        result = vms.evaluate(truth, abundance, substitution, mechanism, folds,
                              positions)
        assert set(result["mean_loss"]) == set(vms.COMPARATORS)

    def test_a_misaligned_block_is_refused(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(14, conditional=True)
        folds = vms.position_groups(positions, 5, vms.FROZEN["seed"])
        with pytest.raises(vms.VariantSelectionError, match="disagrees on n"):
            vms.evaluate(truth, abundance[:-1], substitution, mechanism, folds,
                         positions)

    def test_too_few_groups_is_refused(self, vms):
        truth, abundance, substitution, mechanism, positions = \
            self._cohort(15, conditional=True)
        folds = np.zeros(truth.size, dtype=int)
        with pytest.raises(vms.VariantSelectionError, match="grouped split"):
            vms.evaluate(truth, abundance, substitution, mechanism, folds,
                         positions)


# --------------------------------------------------------------------------
# 7. the ladder, identical to the rest of the project
# --------------------------------------------------------------------------
class TestVerdictLadder:
    @staticmethod
    def _result(better, worse):
        return {"losses": {"a": np.asarray(better, float),
                           "b": np.asarray(worse, float)}}

    def test_a_tiny_real_gain_is_not_called_a_gain(self, vms):
        rng = np.random.default_rng(31)
        worse = 1.0 + rng.normal(scale=0.0005, size=600)
        better = worse - 0.005
        assert vms.contrast(self._result(better, worse), "a", "b")["verdict"] \
            == "GAIN_PRESENT_BUT_BELOW_MARGIN"

    def test_a_large_gain_reaches_the_margin(self, vms):
        rng = np.random.default_rng(32)
        worse = 1.0 + rng.normal(scale=0.01, size=600)
        better = worse - 0.3
        assert vms.contrast(self._result(better, worse), "a", "b")["verdict"] \
            == "GAIN_REACHES_MARGIN"

    def test_being_worse_is_named(self, vms):
        rng = np.random.default_rng(33)
        worse = 1.0 + rng.normal(scale=0.01, size=600)
        better = worse + 0.3
        assert vms.contrast(self._result(better, worse), "a", "b")["verdict"] \
            == "DETECTABLY_WORSE"

    def test_a_wide_null_is_not_a_refutation(self, vms):
        rng = np.random.default_rng(34)
        worse = 1.0 + rng.normal(scale=2.0, size=20)
        better = worse + rng.normal(scale=2.0, size=20)
        cell = vms.contrast(self._result(better, worse), "a", "b")
        assert cell["verdict"] == "NO_DETECTABLE_GAIN_UNDERPOWERED"
        assert cell["an_interval_containing_zero_is_not_a_refutation"] is True

    def test_the_margin_matches_the_rest_of_the_project(self, vms):
        from omics2pk import fusion_headroom_diagnostics as fhd
        assert (vms.FROZEN["minimum_interesting_relative_gain"]
                == fhd.FROZEN["minimum_interesting_relative_gain"])

    def test_an_unknown_comparator_is_refused(self, vms):
        with pytest.raises(vms.VariantSelectionError, match="unknown"):
            vms.contrast(self._result([1.0], [1.0]), "a", "zzz")


# --------------------------------------------------------------------------
# 8. what the module refuses to claim
# --------------------------------------------------------------------------
class TestScopeDisclaimers:
    def test_it_denies_being_a_pharmacokinetic_result(self, vms):
        out = vms.what_this_module_does_not_do()
        assert "not a clearance" in out["is_not_a_pharmacokinetic_result"].lower()

    def test_it_denies_being_a_multi_omic_measurement_in_individuals(self, vms):
        out = vms.what_this_module_does_not_do()
        assert "engineered single variants" in \
            out["is_not_a_multi_omic_measurement_in_individuals"]

    def test_it_separates_measurement_type_from_biomarker_choice(self, vms):
        out = vms.what_this_module_does_not_do()
        assert "separate question" in out["measurement_type_versus_biomarker"]

    def test_it_refuses_to_read_an_annotation_as_causal(self, vms):
        out = vms.what_this_module_does_not_do()
        assert "not a causal claim" in out["does_not_treat_an_annotation_as_causal"]

    def test_the_docstring_names_the_assay_system_confound(self, vms):
        source = (ROOT / "src" / "omics2pk"
                  / "variant_mechanism_omic_selection.py").read_text(
                      encoding="utf-8")
        assert "humanised yeast" in source
        assert "NUDT15 is the primary gene" in source


# --------------------------------------------------------------------------
# 9. the runner's declarations
# --------------------------------------------------------------------------
class TestRunnerDeclarations:
    @pytest.fixture(scope="class")
    def source(self):
        return (ROOT / "scripts"
                / "run_variant_mechanism_benchmark.py").read_text(
                    encoding="utf-8")

    def test_the_four_score_sets_are_named_by_urn(self, source):
        for urn in ("00000055-a-1", "00000055-b-1",
                    "00000095-a-1", "00000095-b-1"):
            assert urn in source

    def test_the_embedded_annotation_can_be_verified_against_uniprot(self,
                                                                     source):
        assert "verify_uniprot" in source
        assert "embedded annotation is stale" in source

    def test_nudt15_is_primary_and_cyp2c9_is_a_sensitivity(self, source):
        assert '"role": "primary"' in source
        assert '"role": "sensitivity"' in source
        assert "assay_confound" in source

    def test_the_annotation_sparsity_is_called_a_finding(self, source):
        # the sentence wraps in the comment, so match a non-wrapping fragment
        assert "asymmetric between proteins" in source
        assert "exactly one annotated binding residue" in source

    def test_both_groupings_are_run(self, source):
        assert "grouped_by_position" in source
        assert "grouped_by_contiguous_region" in source

    def test_the_specificity_control_is_run(self, source):
        assert "predicting_abundance_instead_specificity_control" in source

    def test_the_annotation_coverage_is_measured_not_assumed(self, source):
        assert "share_of_discordant_variants_the_annotation_reaches" in source

    def test_the_downloaded_files_must_stay_outside_the_repository(self, source):
        assert "must stay outside the" in source

    def test_the_artefact_declares_it_holds_no_individual_data(self, source):
        assert "no_individual_level_data_in_this_artefact" in source


# --------------------------------------------------------------------------
# 10. the two documents must keep saying what was measured
# --------------------------------------------------------------------------
ARCHITECTURE = (ROOT / "docs"
                / "OMICS2PK_VARIANT_MECHANISM_SELECTION_ARCHITECTURE_20261111.md")
RESULTS = (ROOT / "docs"
           / "OMICS2PK_MAVEDB_VARIANT_MECHANISM_RESULTS_20261111.md")


class TestArchitectureDocument:
    @pytest.fixture(scope="class")
    def text(self):
        return ARCHITECTURE.read_text(encoding="utf-8")

    def test_it_is_versioned(self, text):
        assert "**Version 1.2.0 \u2014 2026-11-11.**" in text

    def test_it_answers_all_six_commissioned_items(self, text):
        for item in range(1, 7):
            assert f"Deliverable {item}" in text

    def test_it_states_the_headline_both_ways(self, text):
        assert "confirmed — once the knowledge layer is adequate" in text
        assert "6.74" in text
        assert "+15.83 % / +9.47 %" in text

    def test_it_names_annotation_coverage_as_the_constraint(self, text):
        assert "binding constraint is annotation coverage, not the" in text

    def test_it_preserves_h1_h5_and_the_fusing_record(self, text):
        assert "H1\u2013H5 and every frozen result are untouched" in text
        assert "suspended, not deleted" in text

    def test_it_forbids_reading_an_interaction_as_causal(self, text):
        assert "An interaction is not a cause" in text
        assert "not that the partner is a useful" in text

    def test_it_keeps_the_degree_stratified_null_mandatory(self, text):
        assert "degree-stratified" in text
        assert "network proximity measures how well studied a gene is" in text

    def test_it_records_the_pharmvar_access_ladder(self, text):
        assert "401" in text and "download file for the same data returns" in text

    def test_it_explains_why_the_pgx_layer_is_not_a_feature(self, text):
        assert "cannot be a feature" in text
        assert "97 % missing" in text

    def test_it_requires_no_gnn_first(self, text):
        assert "No GNN" in text

    def test_it_requires_a_positive_control(self, text):
        assert "positive control is not optional" in text

    def test_it_lists_what_is_absent_and_needed(self, text):
        assert "choosing the **selector itself** inside training folds" in text

    def test_it_declares_independence_from_thesis(self, text):
        assert "Nothing is\nimported" in text or "Nothing is imported" in text

    def test_it_keeps_clinical_pk_not_evaluable(self, text):
        assert "`NOT_EVALUABLE`" in text


class TestResultsDocument:
    @pytest.fixture(scope="class")
    def text(self):
        return RESULTS.read_text(encoding="utf-8")

    def test_it_is_versioned(self, text):
        assert "**Version 3.0.0 \u2014 2026-11-11.**" in text

    def test_it_reports_the_enrichment_with_its_interval(self, text):
        assert "6.74\u00d7 [4.82, 9.42]" in text

    def test_it_reports_the_coverage_that_limits_the_rule(self, text):
        assert "48.8 %" in text and "12.5 %" in text

    def test_it_reports_that_no_conditional_cell_reaches_the_margin(self, text):
        assert "reaches the declared 5 % margin in **none**" in text

    def test_it_reports_the_grouping_that_changed_a_conclusion(self, text):
        # Version 1.0.0 read the two splits as a sign reversal and 2.0.0
        # repeated it. Under a position-clustered interval the sequence layer
        # adds nothing measurable in EITHER direction, so that claim is
        # withdrawn and the withdrawal has to be visible in its place.
        assert "in either direction, under either split" in text
        assert "2.8" in text
        assert "load-bearing" in text

    def test_it_reports_the_positive_control(self, text):
        assert "+26.6 %" in text
        assert "a null below is about the data, not the harness" in text

    def test_it_admits_the_first_conditional_form_failed_its_own_control(
            self, text):
        assert "failed its own positive control at\n+3.2 %" in text \
            or "failed its own positive control at +3.2 %" in text

    def test_it_flags_the_cyp2c9_assay_system_confound(self, text):
        assert "humanised\nyeast" in text or "humanised yeast" in text
        assert "difference of *system*" in text

    def test_it_states_cyp2c9_is_untestable_rather_than_refuted(self, text):
        assert "untestable here, not refuted" in text

    def test_it_refuses_to_claim_the_two_genes_replicate(self, text):
        assert "not independent evidence of one another" in text

    def test_it_denies_a_pharmacokinetic_reading(self, text):
        assert "not a clearance, not an AUC, not a concentration" in text

    def test_it_keeps_the_premise_honest_about_abundance(self, text):
        assert "the single most useful thing measured here" in text

    def test_it_records_that_the_pgx_layer_was_not_tested(self, text):
        assert "it was not" in text.lower()

    def test_the_status_ceiling_is_exploratory(self, text):
        assert "Status ceiling: `exploratory`" in text


# --------------------------------------------------------------------------
# 11. the structural layer: the documents must keep reporting what it changed
# --------------------------------------------------------------------------
class TestStructuralLayerIsReported:
    @pytest.fixture(scope="class")
    def results(self):
        return RESULTS.read_text(encoding="utf-8")

    @pytest.fixture(scope="class")
    def architecture(self):
        return ARCHITECTURE.read_text(encoding="utf-8")

    def test_the_title_no_longer_calls_the_rule_useless(self, results):
        first = results.splitlines()[0]
        assert "becomes useful" in first, (
            "the conclusion changed when the structural layer was added; a "
            "title left saying otherwise is a stale claim")

    def test_the_diagnosis_is_recorded_as_tested_not_asserted(self, results):
        assert "diagnosis was refutable, and it was tested" in results

    def test_the_coverage_gain_is_reported(self, results):
        for figure in ("48.8 %", "95.9 %", "12.5 %", "87.5 %"):
            assert figure in results

    def test_the_enrichment_gain_is_reported_with_intervals(self, results):
        assert "55.35\u00d7 [21.09, 277.53]" in results
        assert "29.70\u00d7 [10.78, 105.08]" in results

    def test_cyp2c9_is_recorded_as_moving_from_untestable_to_supported(
            self, results):
        # Not "confirmed". The structures were chosen after these
        # proteins' discordant residues had been seen, which no internal
        # validation undoes, so the strongest admissible word is "supported"
        # until the pre-registered third protein has run.
        assert "from *untestable* to *supported*" in results
        assert "PENDING INDEPENDENT" in results

    def test_the_split_that_no_longer_reverses_is_reported(self, results):
        assert "+15.83 %" in results and "+9.47 %" in results

    def test_the_conditional_arm_is_still_reported_as_failing(self, results):
        assert "reaches the margin nowhere" in results

    def test_the_numbering_guard_measurement_is_reported(self, results):
        assert "2.378 \u00c5" in results
        assert "2.11 \u00c5" in results

    def test_why_a_sequence_ruler_had_to_fail_is_shown(self, results):
        assert "134 residues" in results
        assert "I-helix" in results

    def test_the_structural_layer_limits_are_declared(self, results):
        assert "crystallised ligand complex" in results \
            or "crystallised\nligand complex" in results

    def test_the_architecture_records_the_module(self, architecture):
        assert "`structural_annotation.py`" in architecture

    def test_the_architecture_marks_item_one_done_with_its_residue(
            self, architecture):
        assert "DONE 2026-11-11 and it" in architecture
        assert "Still outstanding from this item" in architecture

    def test_the_architecture_records_why_alphafold_was_not_used(
            self, architecture):
        assert "no ligands" in architecture

    def test_the_architecture_separates_selecting_from_weighting(
            self, architecture):
        assert "helps **select** information" in architecture
        assert "does not help **weight**" in architecture

    def test_both_documents_keep_their_earlier_change_log_entries(
            self, results, architecture):
        assert "| 1.0.0 | 2026-11-11 | First issue." in results
        assert "| 1.0.0 | 2026-11-11 | First issue." in architecture
