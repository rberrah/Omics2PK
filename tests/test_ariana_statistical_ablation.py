"""Contract tests for the ARIANA statistical-versus-mechanistic ablation.

This protocol's whole claim to legitimacy is that it adds **no model and no
arm**: it reuses two frozen runners, one fold list, and adds a declared family
with its multiplicity adjustment. These tests exist to make that claim
checkable, and to stop four specific failures:

1. a new arm, panel or fold creeping in, which would make the contrasts
   incomparable to the frozen primary and key secondary;
2. a Family B contrast comparing one omic ratio against two and being read as
   a comparison of representations;
3. a degradation whose adjusted interval excludes zero being filed as "nothing
   detected", or an unadjusted interval being promoted to a claim;
4. the endpoint being described as an AUC, a clearance or an exposure metric.
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
    / "OMICS2PK_ARIANA_STATISTICAL_VERSUS_MECHANISTIC_ABLATION_PRESPECIFICATION_20261110.md"
)
RUNNER = REPO / "scripts" / "run_ariana_statistical_versus_mechanistic_ablation.py"
REPORT = REPO / "reports" / "data" / "ariana_statistical_versus_mechanistic.json"
MECH_REPORT = REPO / "reports" / "data" / "h1_ariana_mechanistic_hybrid.json"
STAT_REPORT = REPO / "reports" / "data" / "h1_ariana_trough.json"

pytest.importorskip("numpy")
pytest.importorskip("pandas")
pytest.importorskip("sklearn")


@pytest.fixture(scope="module")
def run():
    sys.path.insert(0, str(REPO / "src"))
    spec = importlib.util.spec_from_file_location("ariana_stat_ablation", RUNNER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _normalise(text):
    return re.sub(r"\s+", " ", text)


# --------------------------------------------------------------------------
# 1. nothing new is introduced
# --------------------------------------------------------------------------
class TestNothingNewIsIntroduced:
    def test_the_family_has_exactly_eleven_contrasts(self, run):
        assert len(run.FAMILY_A) == 6
        assert len(run.FAMILY_B) == 5
        assert run.FAMILY_SIZE == 11

    def test_every_statistical_arm_comes_from_the_frozen_runner(self, run):
        used = {arm for pair in run.FAMILY_A.values() for arm in pair}
        used |= {new for new, _, _ in run.FAMILY_B.values()}
        assert used <= set(run.H1.ARMS), sorted(used - set(run.H1.ARMS))

    def test_every_mechanistic_arm_comes_from_the_frozen_module(self, run):
        used = {base for _, base, _ in run.FAMILY_B.values()}
        assert used <= set(run.MECH.MECHANISTIC_ARMS)

    def test_every_panel_is_already_declared(self, run):
        used = {panel for _, _, panel in run.FAMILY_B.values()}
        assert used <= set(run.MECH.OMIC_PANELS)

    def test_the_runner_defines_no_arm_of_its_own(self):
        source = RUNNER.read_text(encoding="utf-8")
        # An arm table in this file would mean a new model was defined here.
        assert "ARMS = {" not in source
        assert "CLINICAL = [" not in source

    def test_the_machinery_check_passes_against_the_current_runners(self, run):
        checked = run.check_imported_machinery()
        assert checked["mechanistic_primary_panel"] == run.UNMATCHED_PANEL
        assert len(checked["matched_panel_columns"]) == len(
            checked["statistical_omic_block"]
        )

    def test_a_changed_fold_constant_is_refused(self, run, monkeypatch):
        monkeypatch.setattr(run.H1, "SEED", 1)
        with pytest.raises(run.RunError, match="frozen constants have changed"):
            run.check_imported_machinery()

    def test_a_missing_arm_is_refused(self, run, monkeypatch):
        reduced = {
            name: value
            for name, value in run.H1.ARMS.items()
            if name != "M1a_clinical_genotype"
        }
        monkeypatch.setattr(run.H1, "ARMS", reduced)
        with pytest.raises(run.RunError, match="no longer defines the arms"):
            run.check_imported_machinery()

    def test_a_missing_panel_is_refused(self, run, monkeypatch):
        panels = {
            name: value
            for name, value in run.MECH.OMIC_PANELS.items()
            if name != run.MATCHED_PANEL
        }
        monkeypatch.setattr(run.MECH, "OMIC_PANELS", panels)
        with pytest.raises(run.RunError, match="no longer declares panel"):
            run.check_imported_machinery()


# --------------------------------------------------------------------------
# 2. the matched-panel requirement
# --------------------------------------------------------------------------
class TestPanelMatching:
    def test_the_matched_panel_has_as_many_ratios_as_the_statistical_block(self, run):
        assert len(run.MECH.OMIC_PANELS[run.MATCHED_PANEL]) == len(run.H1.OMIC)

    def test_the_matched_panel_holds_the_same_ratios(self, run):
        assert set(run.MECH.OMIC_PANELS[run.MATCHED_PANEL]) == set(run.H1.OMIC)

    def test_the_unmatched_panel_is_declared_as_unmatched(self, run):
        assert len(run.MECH.OMIC_PANELS[run.UNMATCHED_PANEL]) < len(run.H1.OMIC)

    def test_a_mismatched_matched_panel_is_refused(self, run, monkeypatch):
        panels = dict(run.MECH.OMIC_PANELS)
        panels[run.MATCHED_PANEL] = ("lnMR OH-solanidine/solanidine",)
        monkeypatch.setattr(run.MECH, "OMIC_PANELS", panels)
        with pytest.raises(run.RunError, match="would not be matched"):
            run.check_imported_machinery()

    def test_four_of_the_five_family_b_contrasts_are_matched(self, run):
        matched = [
            key
            for key, (_, _, panel) in run.FAMILY_B.items()
            if panel == run.MATCHED_PANEL
        ]
        assert len(matched) == 4
        assert "B5_the_unmatched_panel" not in matched


# --------------------------------------------------------------------------
# 3. the prohibited predictor
# --------------------------------------------------------------------------
class TestProhibitedPredictor:
    def test_it_appears_in_no_statistical_arm(self, run):
        for arm, features in run.H1.ARMS.items():
            assert run.H1.COL_PROHIBITED not in features, arm

    def test_it_appears_in_no_mechanistic_panel(self, run):
        for panel, columns in run.MECH.OMIC_PANELS.items():
            assert run.H1.COL_PROHIBITED not in columns, panel

    def test_the_runner_asserts_it_before_running(self, run):
        assert run.MECHRUN.assert_no_prohibited_predictor() is None


# --------------------------------------------------------------------------
# 4. the freeze agrees with the document
# --------------------------------------------------------------------------
class TestFreezeAgreesWithTheDocument:
    @pytest.fixture(scope="class")
    def text(self):
        assert PRESPEC.is_file(), f"prespecification absent: {PRESPEC}"
        return PRESPEC.read_text(encoding="utf-8")

    def test_the_check_passes(self, run):
        out = run.check_freeze_against_document()
        assert out["checked"]["Family size"] == 11
        assert out["checked"]["Adjustment"] == "Bonferroni"

    def test_the_family_size_is_in_the_table(self, run, text):
        assert f"| Family size | {run.FAMILY_SIZE} |" in text

    def test_the_seed_is_in_the_table(self, run, text):
        assert f"| Fold seed | `{run.H1.SEED}` |" in text

    def test_every_contrast_is_named(self, run, text):
        for key in list(run.FAMILY_A) + list(run.FAMILY_B):
            assert f"`{key}`" in text, key

    def test_the_notes_named_contrast_is_marked(self, run, text):
        assert run.NOTES_NAMED_CONTRAST == "A5_pgx_on_top_of_metabolomics"
        assert f"| Note's named contrast | `{run.NOTES_NAMED_CONTRAST}` |" in text

    def test_a_changed_family_size_is_refused(self, run, monkeypatch):
        monkeypatch.setattr(run, "FAMILY_SIZE", 12)
        monkeypatch.setitem(
            run._FROZEN_IN_DOCUMENT, "Family size", (lambda: 12, lambda s: int(s))
        )
        with pytest.raises(run.RunError, match="Family size"):
            run.check_freeze_against_document()

    def test_a_document_missing_the_blindness_disclosure_is_refused(self, run, tmp_path):
        reduced = PRESPEC.read_text(encoding="utf-8").replace(
            "Family A below is not a blind prespecification", "Family A is fine"
        )
        path = tmp_path / "prespec.md"
        path.write_text(reduced, encoding="utf-8")
        with pytest.raises(run.RunError, match="not a blind test"):
            run.check_freeze_against_document(path)

    def test_a_document_missing_the_family_b_limit_is_refused(self, run, tmp_path):
        reduced = PRESPEC.read_text(encoding="utf-8").replace(
            "It does not isolate the representation", "It isolates everything"
        )
        path = tmp_path / "prespec.md"
        path.write_text(reduced, encoding="utf-8")
        with pytest.raises(run.RunError, match="what Family B cannot measure"):
            run.check_freeze_against_document(path)


# --------------------------------------------------------------------------
# 5. the document's own content
# --------------------------------------------------------------------------
class TestPrespecificationContent:
    @pytest.fixture(scope="class")
    def text(self):
        return _normalise(PRESPEC.read_text(encoding="utf-8"))

    def test_it_says_no_new_model_is_fitted(self, text):
        assert "fits **no new model**" in text

    def test_it_says_why_a_second_benchmark_was_rejected(self, text):
        assert "rejected" in text
        assert "fifth pass over a cohort" in text

    def test_it_justifies_keeping_the_dose_in_the_baseline(self, text):
        assert (
            "a concentration model without dose is not a weaker model but a "
            "wrong one" in text
        )

    def test_it_refuses_the_literal_alone_arm_with_a_reason(self, text):
        assert "genotype with no dose" in text

    def test_it_denies_the_endpoint_is_an_exposure_metric(self, text):
        assert "**not** an AUC" in text
        assert "**not** a clearance" in text

    def test_it_explains_bonferroni_instead_of_holm(self, text):
        assert "Bonferroni, not Holm" in text
        assert "uniformly more conservative" in text

    def test_it_forbids_promoting_an_unadjusted_interval(self, text):
        assert "No promotion of an unadjusted interval to a claim" in text

    def test_it_states_the_non_detection_rule(self, text):
        assert (
            "no effect was **detected**, never that there is none" in text
        )

    def test_it_forbids_a_clinical_pk_claim(self, text):
        assert "may be presented as a clinical-PK proof" in text

    def test_it_protects_the_frozen_contrasts(self, text):
        assert "No restatement of the frozen primary contrast or key secondary" in text


# --------------------------------------------------------------------------
# 6. the decision rule
# --------------------------------------------------------------------------
class TestDecisionRule:
    def _contrast(self, run, new_mae, base_mae, n=47):
        import numpy as np

        rng = np.random.default_rng(3)
        observed = rng.uniform(20.0, 200.0, size=n)
        new = {
            "absolute": np.full(n, new_mae),
            "log_absolute": np.full(n, new_mae / 100.0),
        }
        base = {
            "absolute": np.full(n, base_mae),
            "log_absolute": np.full(n, base_mae / 100.0),
        }
        return run.one_contrast(new, base, observed, 200, run.FAMILY_SIZE)

    def test_a_large_improvement_carries_a_claim(self, run):
        out = self._contrast(run, 10.0, 50.0)
        assert out["decision"]["carries_a_claim"] is True
        assert out["primary_loss"]["direction"] == "improvement"

    def test_a_degradation_never_carries_a_claim(self, run):
        out = self._contrast(run, 50.0, 10.0)
        assert out["decision"]["carries_a_claim"] is False
        assert out["primary_loss"]["direction"] == "degradation"

    def test_a_tiny_improvement_does_not_reach_the_threshold(self, run):
        out = self._contrast(run, 49.5, 50.0)
        assert out["decision"]["reaches_the_relative_threshold"] is False
        assert out["decision"]["carries_a_claim"] is False

    def test_the_adjusted_interval_is_wider_than_the_unadjusted_one(self, run):
        out = self._contrast(run, 30.0, 50.0)
        unadjusted = out["primary_loss"]["descriptive_95_interval"]
        adjusted = out["multiplicity_adjusted"]["interval"]
        assert adjusted[0] <= unadjusted[0]
        assert adjusted[1] >= unadjusted[1]

    def test_the_adjustment_declares_its_family_size(self, run):
        out = self._contrast(run, 30.0, 50.0)
        assert out["multiplicity_adjusted"]["family_size"] == run.FAMILY_SIZE

    def test_the_rule_is_written_into_the_output(self, run):
        out = self._contrast(run, 30.0, 50.0)
        assert "never a positive finding" in out["decision"]["rule"]
        assert "never that there is none" in out["decision"]["rule"]


# --------------------------------------------------------------------------
# 7. the committed report
# --------------------------------------------------------------------------
class TestCommittedReport:
    @pytest.fixture(scope="class")
    def payload(self):
        if not REPORT.is_file():
            pytest.skip("ablation report not generated")
        return json.loads(REPORT.read_text(encoding="utf-8"))

    def test_it_reports_forty_seven_patients(self, payload):
        assert payload["n_patients"] == 47

    def test_it_declares_the_status_ceiling(self, payload):
        assert "exploratory" in payload["status_ceiling"]
        assert "independent validation" in payload["status_ceiling"]

    def test_it_denies_the_endpoint_is_an_exposure_metric(self, payload):
        assert "NOT an AUC" in payload["target_is"]
        assert "NOT a clearance" in payload["target_is"]

    def test_it_states_that_no_new_model_was_fitted(self, payload):
        assert "imported from the two frozen runners" in payload["no_new_model_was_fitted"]

    def test_it_discloses_that_family_a_was_not_blind(self, payload):
        assert "already computed and visible" in payload[
            "family_a_is_not_a_blind_prespecification"
        ]

    def test_it_discloses_what_family_b_cannot_isolate(self, payload):
        assert "not representation alone" in payload[
            "family_b_cannot_isolate_the_representation"
        ]

    def test_the_folds_are_the_frozen_mechanistic_runs_folds(self, payload):
        if not MECH_REPORT.is_file():
            pytest.skip("mechanistic report not generated")
        mechanistic = json.loads(MECH_REPORT.read_text(encoding="utf-8"))
        assert (
            payload["validation"]["fold_fingerprint_sha256"]
            == mechanistic["validation"]["fold_fingerprint_sha256"]
        )

    def test_the_statistical_losses_match_the_frozen_statistical_run(self, payload):
        if not STAT_REPORT.is_file():
            pytest.skip("statistical report not generated")
        frozen = json.loads(STAT_REPORT.read_text(encoding="utf-8"))["results"]["ridge"]
        for arm, loss in payload["losses"]["statistical_ridge"].items():
            assert loss["mae"] == pytest.approx(frozen[arm]["mae"], rel=1e-9), arm

    def test_every_contrast_carries_both_intervals(self, payload):
        for key, contrast in payload["contrasts"].items():
            assert "descriptive_95_interval" in contrast["primary_loss"], key
            assert "interval" in contrast["multiplicity_adjusted"], key
            assert contrast["multiplicity_adjusted"]["family_size"] == 11, key

    def test_the_three_outcome_lists_partition_the_family(self, payload):
        keys = set(payload["contrasts"])
        claims = set(payload["claims_carried"])
        baseline = set(payload["differences_detected_in_favour_of_the_baseline"])
        nothing = set(payload["nothing_detected"])
        assert claims | baseline | nothing == keys
        assert not (claims & nothing)
        assert not (baseline & nothing)

    def test_a_degradation_is_never_in_the_claims_list(self, payload):
        for key in payload["claims_carried"]:
            assert payload["contrasts"][key]["primary_loss"]["direction"] == "improvement"

    def test_the_notes_question_is_answered_with_its_direction(self, payload):
        answer = payload["the_question_the_note_named"]
        assert answer["contrast"] == "A5_pgx_on_top_of_metabolomics"
        assert answer["direction"] in ("improvement", "degradation")
        if answer["direction"] == "degradation":
            assert "degradation" in answer["answer"]

    def test_the_free_parameter_asymmetry_is_published(self, payload):
        counts = payload["free_parameters_fitted_on_this_cohort"]
        assert counts["M5_mechanistic_pgx_omic"] == 1
        assert counts["M4_mechanistic_pgx"] == 1

    def test_it_refuses_a_clinical_pk_reading(self, payload):
        assert "is a clinical-pharmacokinetic proof" in payload["no_clinical_pk_proof"]
        assert payload["no_clinical_pk_proof"].startswith("No result here")

    def test_it_exports_no_patient_level_data(self, payload):
        assert payload["no_person_identifiers_values_or_predictions_exported"] is True

    def test_family_b_contrasts_declare_whether_the_panels_match(self, payload):
        for key, contrast in payload["contrasts"].items():
            if contrast["family"] == "B":
                assert "panels_are_matched" in contrast, key
                assert "mechanistic_panel" in contrast, key
