"""Tests for the Axis B2 validation: the design, the guards, and the documents.

The point of these tests is not that the code runs. It is that the design
cannot drift into testing something other than what the protocol froze, and
that the two things this run exists to establish -- a replication at risk and
a refitting interval -- cannot be quietly weakened.
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

PRESPEC = (ROOT / "docs"
           / "OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_PRESPECIFICATION_20261110.md")
RESULTS = (ROOT / "docs"
           / "OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_RESULTS_20261110.md")
ARTEFACT = ROOT / "reports" / "data" / "plasma_biomarker_validation.json"
QUALIFICATION = (ROOT / "docs"
                 / "OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_20261110.md")


@pytest.fixture(scope="module")
def mod():
    from omics2pk import plasma_biomarker_validation
    return plasma_biomarker_validation


@pytest.fixture(scope="module")
def payload():
    if not ARTEFACT.is_file():
        pytest.skip("the run has not been executed in this checkout")
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


@pytest.fixture
def cohort():
    """A small synthetic cohort with a planted abundance-activity relation."""
    rng = np.random.default_rng(7)
    rows, order = {}, []
    for index in range(15):
        latent = rng.normal()
        donor = f"D{index:02d}"
        rows[donor] = {
            "cyp3a4_abundance": float(10.0 ** (1.0 + 0.5 * latent)),
            "cyp3a5_abundance": float(10.0 ** (0.3 * rng.normal())),
            "plasma_4b_hc_cholesterol_ratio": float(
                10.0 ** (-0.7 + 0.2 * rng.normal())),
            "midazolam_1_hydroxylation": float(
                10.0 ** (1.0 + 0.6 * latent + 0.2 * rng.normal())),
            "ibrutinib_m37": float(10.0 ** (1.0 + 0.4 * latent
                                            + 0.3 * rng.normal())),
        }
        rows[donor]["cyp3a_total_abundance"] = (
            rows[donor]["cyp3a4_abundance"] + rows[donor]["cyp3a5_abundance"])
        order.append(donor)
    return rows, order


# --------------------------------------------------------------------------
# 1. the design cannot drift
# --------------------------------------------------------------------------
class TestDesignIntegrity:
    def test_every_contrast_permutes_exactly_what_its_arm_adds(self, mod):
        for name, spec in mod.CONTRASTS.items():
            added = (set(mod.ARMS[spec["arm"]])
                     - set(mod.ARMS[spec["comparator"]]))
            assert set(spec["permuted_block"]) == added, (
                f"{name} would permute something other than the information "
                "the contrast is about, so its null would not be its null")

    def test_the_baseline_arm_has_no_predictors(self, mod):
        assert mod.ARMS["V0_mean"] == ()

    def test_the_primary_contrast_is_the_incremental_one(self, mod):
        assert mod.FROZEN["primary_arm"] == "V3_protein_plus_plasma"
        assert mod.FROZEN["primary_comparator"] == "V1_hepatic_protein"

    def test_the_family_size_matches_the_design(self, mod):
        endpoints = {mod.FROZEN["primary_endpoint"],
                     mod.FROZEN["secondary_endpoint"]}
        assert (mod.FROZEN["declared_family_size"]
                == len(mod.CONTRASTS) * len(endpoints))

    def test_at_most_two_predictors_anywhere(self, mod):
        for arm, columns in {**mod.ARMS, **mod.SENSITIVITY_ARMS}.items():
            assert len(columns) <= 2, f"{arm} is wider than the freeze allows"

    def test_the_seed_and_grid_come_from_the_frozen_pilot(self, mod):
        from omics2pk import protein_clint_pilot as pilot
        assert mod.FROZEN["seed"] == pilot.FROZEN["seed"]
        assert mod.FROZEN["ridge_alphas"] == pilot.FROZEN["ridge_alphas"]


# --------------------------------------------------------------------------
# 2. the circularity guard
# --------------------------------------------------------------------------
class TestCircularityGuard:
    @pytest.mark.parametrize("name", [
        "cyp3a4_abundance", "cyp3a5_abundance", "cyp3a_total_abundance",
        "plasma_4b_hc_cholesterol_ratio",
    ])
    def test_a_predictor_may_never_be_the_endpoint(self, mod, name):
        with pytest.raises(mod.PlasmaBiomarkerValidationError,
                           match="predicting its own input"):
            mod.assert_endpoint_is_admissible(name)

    def test_the_real_endpoints_are_admitted(self, mod):
        for name in (mod.FROZEN["primary_endpoint"],
                     mod.FROZEN["secondary_endpoint"]):
            assert mod.assert_endpoint_is_admissible(name) == name

    def test_every_predictor_column_is_in_the_forbidden_list(self, mod):
        used = set()
        for columns in {**mod.ARMS, **mod.SENSITIVITY_ARMS}.values():
            used |= set(columns)
        assert used <= set(mod.FORBIDDEN_ENDPOINTS), (
            "a predictor exists that the endpoint guard would not refuse")


# --------------------------------------------------------------------------
# 3. the fitting is honestly out of sample
# --------------------------------------------------------------------------
class TestLeaveOneDonorOut:
    def test_it_predicts_every_donor_once(self, mod, cohort):
        rows, order = cohort
        target = mod.log_target(rows, order, "midazolam_1_hydroxylation")
        block = mod.design(rows, order, ("cyp3a4_abundance",))
        predicted, errors, penalties = mod.leave_one_donor_out(block, target)
        assert predicted.shape == errors.shape == penalties.shape == (15,)

    def test_a_planted_relation_beats_the_mean(self, mod, cohort):
        rows, order = cohort
        arm, _ = mod.score_arm(rows, order, "midazolam_1_hydroxylation",
                               ("cyp3a4_abundance",))
        base, _ = mod.score_arm(rows, order, "midazolam_1_hydroxylation", ())
        assert arm["mae_log10"] < base["mae_log10"]

    def test_an_empty_arm_is_the_training_mean(self, mod, cohort):
        rows, order = cohort
        score, errors = mod.score_arm(rows, order, "ibrutinib_m37", ())
        assert score["n_predictors"] == 0
        assert np.all(errors >= 0.0)

    def test_too_few_donors_is_refused(self, mod):
        with pytest.raises(mod.PlasmaBiomarkerValidationError,
                           match="at least four donors"):
            mod.leave_one_donor_out(np.ones((3, 1)), np.ones(3))

    def test_a_non_positive_endpoint_raises_rather_than_dropping_a_donor(
            self, mod, cohort):
        rows, order = cohort
        rows[order[0]]["midazolam_1_hydroxylation"] = 0.0
        with pytest.raises(mod.PlasmaBiomarkerValidationError,
                           match="never by silently dropping"):
            mod.log_target(rows, order, "midazolam_1_hydroxylation")

    def test_a_missing_predictor_cell_raises(self, mod, cohort):
        rows, order = cohort
        rows[order[2]]["cyp3a4_abundance"] = float("nan")
        with pytest.raises(mod.PlasmaBiomarkerValidationError,
                           match="missing or non-finite"):
            mod.design(rows, order, ("cyp3a4_abundance",))


# --------------------------------------------------------------------------
# 4. the refitting interval, which is why this module exists
# --------------------------------------------------------------------------
class TestRefittingInterval:
    def test_it_declares_that_it_refits(self, mod, cohort):
        rows, order = cohort
        out = mod.out_of_bag_bootstrap(
            rows, order, "midazolam_1_hydroxylation",
            ("cyp3a4_abundance",), (), n_boot=200)
        assert out["refits_the_model_inside_each_draw"] is True

    def test_it_is_wider_than_the_stored_error_interval(self, mod, cohort):
        from omics2pk import protein_clint_pilot as pilot
        rows, order = cohort
        _, arm_errors = mod.score_arm(rows, order,
                                      "midazolam_1_hydroxylation",
                                      ("cyp3a4_abundance",))
        _, base_errors = mod.score_arm(rows, order,
                                       "midazolam_1_hydroxylation", ())
        stored = pilot.paired_bootstrap(arm_errors, base_errors, n_boot=800)
        refit = mod.out_of_bag_bootstrap(
            rows, order, "midazolam_1_hydroxylation",
            ("cyp3a4_abundance",), (), n_boot=800)
        width = lambda block: block["interval"][1] - block["interval"][0]
        assert width(refit) > width(stored), (
            "an interval that refits cannot be narrower than one that does "
            "not; if it is, the resampling is leaking")

    def test_it_scores_only_donors_left_out_of_the_draw(self, mod, cohort):
        rows, order = cohort
        out = mod.out_of_bag_bootstrap(
            rows, order, "midazolam_1_hydroxylation",
            ("cyp3a4_abundance",), (), n_boot=300)
        assert "out-of-bag" in out["scoring_donors_per_draw"]
        assert out["n_draws_used"] + out["n_draws_skipped_too_small"] == 300

    def test_it_refuses_to_report_an_interval_from_too_few_usable_draws(
            self, mod, cohort):
        rows, order = cohort
        with pytest.raises(mod.PlasmaBiomarkerValidationError,
                           match="would not mean anything"):
            mod.out_of_bag_bootstrap(
                rows, order, "midazolam_1_hydroxylation",
                ("cyp3a4_abundance",), (), n_boot=20)

    def test_it_names_what_it_still_cannot_price(self, mod, cohort):
        rows, order = cohort
        out = mod.out_of_bag_bootstrap(
            rows, order, "ibrutinib_m37", ("cyp3a4_abundance",), (),
            n_boot=200)
        assert "fifteen purchased hepatocyte lots" in \
            out["what_it_still_does_not_price"]

    def test_the_bonferroni_level_widens_the_interval(self, mod, cohort):
        rows, order = cohort
        narrow = mod.out_of_bag_bootstrap(
            rows, order, "ibrutinib_m37", ("cyp3a4_abundance",), (),
            n_boot=600, family_size=1)
        wide = mod.out_of_bag_bootstrap(
            rows, order, "ibrutinib_m37", ("cyp3a4_abundance",), (),
            n_boot=600, family_size=6)
        assert (wide["interval"][1] - wide["interval"][0]
                >= narrow["interval"][1] - narrow["interval"][0])


# --------------------------------------------------------------------------
# 5. the null
# --------------------------------------------------------------------------
class TestPermutationNull:
    def test_it_destroys_the_added_information_only(self, mod, cohort):
        rows, order = cohort
        values = mod.permutation_null(rows, order,
                                      "midazolam_1_hydroxylation",
                                      "replication_V1_vs_V0", replicates=8)
        assert len(values) == 8
        assert all(np.isfinite(v) for v in values)

    def test_an_unknown_contrast_is_refused(self, mod, cohort):
        rows, order = cohort
        with pytest.raises(mod.PlasmaBiomarkerValidationError,
                           match="not a frozen contrast"):
            mod.permutation_null(rows, order, "ibrutinib_m37", "made_up")

    def test_deciles_replace_the_replicate_list(self, mod):
        out = mod.deciles(list(range(200)))
        assert set(out) == {f"p{q}" for q in range(0, 101, 10)}
        assert out["p0"] <= out["p50"] <= out["p100"]


# --------------------------------------------------------------------------
# 6. the refusals
# --------------------------------------------------------------------------
class TestRefusals:
    def test_it_names_the_five_things_it_cannot_claim(self, mod):
        out = mod.what_this_cannot_claim()
        for key in ("not_an_intrinsic_clearance", "not_a_pharmacokinetic_result",
                    "not_a_clinical_biomarker_validation",
                    "two_substrates_are_not_two_cohorts",
                    "absence_of_a_gain_is_not_equivalence"):
            assert key in out and len(out[key]) > 40

    def test_absence_of_a_gain_is_not_called_equivalence(self, mod):
        text = mod.what_this_cannot_claim()["absence_of_a_gain_is_not_equivalence"]
        assert "not detected" in text or "no difference was detected" in text


# --------------------------------------------------------------------------
# 7. the artefact exports nothing about a donor
# --------------------------------------------------------------------------
class TestArtefact:
    def test_it_declares_no_donor_level_export(self, payload):
        assert payload[
            "no_donor_identifiers_values_or_predictions_exported"] is True

    def test_it_holds_no_long_numeric_list(self, payload):
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

    def test_the_per_donor_table_is_outside_the_repository(self, payload):
        where = pathlib.Path(payload[
            "per_donor_table_is_outside_this_repository"])
        assert ROOT not in where.parents and ROOT != where.parent

    def test_every_contrast_carries_both_intervals_and_a_null(self, payload):
        for block in payload["results"].values():
            for name, cell in block.items():
                assert cell["stored_error_bootstrap"][
                    "refits_the_model_inside_each_draw"] is False
                assert cell["refitting_out_of_bag_bootstrap"][
                    "refits_the_model_inside_each_draw"] is True
                assert "rank_p_value" in cell["permutation_control"]
                assert cell["verdict"] in {"SUPPORTED", "NOT_SUPPORTED"}

    def test_the_refitting_interval_governs_the_verdict(self, payload):
        for block in payload["results"].values():
            for name, cell in block.items():
                if cell["verdict"] == "SUPPORTED":
                    assert cell["refitting_out_of_bag_bootstrap"][
                        "interval_excludes_zero"] is True, (
                        f"{name} is SUPPORTED on an interval that includes "
                        "zero, which the frozen verdict rule forbids")

    def test_the_permutation_null_is_computed_for_every_contrast(self, payload):
        expected = len(payload["results"]) * 3
        got = sum(len(block) for block in payload["results"].values())
        assert got == expected
        for block in payload["results"].values():
            for cell in block.values():
                assert cell["permutation_control"]["n_replicates"] == 200

    def test_the_wegler_comparator_refuses_to_pool(self, payload):
        text = payload["wegler_comparator"]["why_this_is_not_a_pooled_estimate"]
        assert "not one estimate measured twice" in text


# --------------------------------------------------------------------------
# 8. the documents say what the artefact says
# --------------------------------------------------------------------------
def _normalise(text):
    return re.sub(r"[\s*_`>|]+", " ", text)


class TestDocuments:
    def test_both_documents_exist(self):
        assert PRESPEC.is_file() and RESULTS.is_file()
        assert QUALIFICATION.is_file()

    def test_the_prespecification_admits_it_is_not_blind(self):
        text = _normalise(PRESPEC.read_text(encoding="utf-8"))
        assert "not a blind prespecification" in text

    def test_the_prespecification_fixes_the_verdict_rule_before_the_run(self):
        text = _normalise(PRESPEC.read_text(encoding="utf-8"))
        assert "The stored-error interval alone is never" in text
        assert "refitting interval governs" in text

    def test_the_results_document_reports_the_measured_width_ratios(self):
        text = RESULTS.read_text(encoding="utf-8")
        for ratio in ("1.97", "1.78", "3.76", "2.09", "1.96", "6.32"):
            assert ratio in text, (
                f"the width ratio {ratio} is in the artefact but not in the "
                "results document")

    def test_the_results_document_states_the_consequence_for_wegler(self):
        text = _normalise(RESULTS.read_text(encoding="utf-8"))
        assert "Wegler pilot's positive findings were declared on the narrower" \
            in text

    def test_the_results_document_does_not_claim_equivalence(self):
        text = _normalise(RESULTS.read_text(encoding="utf-8")).lower()
        assert "demonstration that the relationships are absent" in text

    def test_every_verdict_in_the_artefact_appears_in_the_document(self, payload):
        text = RESULTS.read_text(encoding="utf-8")
        assert str(payload["n_contrasts_supported"]) in text
        for block in payload["results"].values():
            for cell in block.values():
                assert cell["verdict"] in text

    def test_the_fold_errors_in_the_document_match_the_artefact(self, payload):
        text = RESULTS.read_text(encoding="utf-8")
        for block in payload["results"].values():
            for cell in block.values():
                for side in ("arm_score", "comparator_score"):
                    printed = f"{cell[side]['geometric_fold_error']:.3f}"
                    assert printed in text, (
                        f"fold error {printed} is in the artefact and not in "
                        "the results document")
