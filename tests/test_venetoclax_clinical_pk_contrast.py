"""Tests for the venetoclax clinical-PK contrast.

The design's two load-bearing claims are that the patient is the unit and that
the baseline is the mean rather than a clinical model. Both are easy to lose by
accident and both change what the result means, so both are tested.
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
           / "OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_PRESPECIFICATION_20261110.md")
RESULTS = (ROOT / "docs"
           / "OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_CONTRAST_20261110.md")
ARTEFACT = ROOT / "reports" / "data" / "venetoclax_clinical_pk_contrast.json"


@pytest.fixture(scope="module")
def mod():
    from omics2pk import venetoclax_clinical_pk_contrast
    return venetoclax_clinical_pk_contrast


@pytest.fixture(scope="module")
def payload():
    if not ARTEFACT.is_file():
        pytest.skip("the run has not been executed in this checkout")
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


@pytest.fixture
def cohort():
    """Twelve patients, two visits, with a planted biomarker relationship."""
    rng = np.random.default_rng(11)
    rows = {}
    for index in range(1, 13):
        latent = rng.normal()
        for visit in (1, 2):
            rows[(str(index), visit)] = {
                "dose_mg_per_day": 400.0,
                "concentration_nm": float(
                    10.0 ** (3.0 - 0.5 * latent + 0.1 * rng.normal())),
                "plasma_4b_hc_ng_per_ml": 20.0,
                "plasma_cholesterol_mg_per_ml": 1.0,
                "plasma_4b_hc_cholesterol_ratio": float(
                    10.0 ** (-1.0 + 0.3 * latent)),
            }
    return rows, sorted(rows)


# --------------------------------------------------------------------------
# 1. the endpoint
# --------------------------------------------------------------------------
class TestEndpoint:
    def test_it_divides_by_the_dose(self, mod):
        rows = {("1", 1): {"dose_mg_per_day": 400.0,
                           "concentration_nm": 4000.0},
                ("2", 1): {"dose_mg_per_day": 200.0,
                           "concentration_nm": 2000.0}}
        target = mod.dose_normalised_target(rows, sorted(rows))
        assert target[0] == pytest.approx(target[1]), (
            "two patients with the same concentration-to-dose ratio must get "
            "the same endpoint, or the dose reduction is being read as "
            "pharmacokinetics")

    def test_a_non_positive_concentration_raises(self, mod):
        rows = {("1", 1): {"dose_mg_per_day": 400.0, "concentration_nm": 0.0}}
        with pytest.raises(mod.VenetoclaxContrastError,
                           match="never by a silent drop"):
            mod.dose_normalised_target(rows, sorted(rows))

    def test_a_missing_dose_raises(self, mod):
        rows = {("1", 1): {"dose_mg_per_day": 0.0, "concentration_nm": 100.0}}
        with pytest.raises(mod.VenetoclaxContrastError, match="no dose"):
            mod.dose_normalised_target(rows, sorted(rows))


# --------------------------------------------------------------------------
# 2. the patient is the unit
# --------------------------------------------------------------------------
class TestPatientIsTheUnit:
    def test_both_visits_of_a_patient_are_held_out_together(self, mod, cohort):
        rows, keys = cohort
        predicted, errors, target = mod.leave_one_patient_out(
            rows, keys, ("plasma_4b_hc_cholesterol_ratio",))
        assert predicted.shape == (24,)
        # With both visits held out together and the biomarker identical
        # across a patient's visits, both predictions must coincide.
        for index in range(0, 24, 2):
            assert predicted[index] == pytest.approx(predicted[index + 1]), (
                "a patient's two visits got different predictions although "
                "their predictor is identical, so the fold is leaking")

    def test_the_bootstrap_resamples_patients(self, mod, cohort):
        rows, keys = cohort
        out = mod.patient_level_bootstrap(rows, keys, n_boot=400)
        assert out["resampling_unit"] == "patient"
        assert "not two individuals" in \
            out["why_the_patient_and_not_the_observation"]

    def test_the_bootstrap_refits(self, mod, cohort):
        rows, keys = cohort
        out = mod.patient_level_bootstrap(rows, keys, n_boot=400)
        assert out["refits_the_model_inside_each_draw"] is True

    def test_the_permutation_carries_a_patients_whole_visit_vector(self, mod):
        source = (ROOT / "src" / "omics2pk"
                  / "venetoclax_clinical_pk_contrast.py").read_text(
                      encoding="utf-8")
        flat = re.sub(r"\s+", " ", source)
        assert "Carrying a patient's whole visit vector" in flat
        assert "would make the null easier to beat" in flat

    def test_the_stored_error_interval_is_labelled_as_overstating(
            self, mod, cohort):
        rows, keys = cohort
        block = mod.contrast(rows, keys)
        assert "OVERSTATES independence" in \
            block["stored_error_bootstrap"]["resampling_unit"]

    def test_too_few_patients_is_refused(self, mod):
        rows = {("1", 1): {"dose_mg_per_day": 400.0,
                           "concentration_nm": 100.0,
                           "plasma_4b_hc_cholesterol_ratio": 0.1},
                ("2", 1): {"dose_mg_per_day": 400.0,
                           "concentration_nm": 200.0,
                           "plasma_4b_hc_cholesterol_ratio": 0.2}}
        with pytest.raises(mod.VenetoclaxContrastError,
                           match="at least four patients"):
            mod.leave_one_patient_out(
                rows, sorted(rows), ("plasma_4b_hc_cholesterol_ratio",))


# --------------------------------------------------------------------------
# 3. a planted relationship is detected, so a null result means something
# --------------------------------------------------------------------------
class TestSensitivityOfTheDesign:
    def test_a_strong_planted_relationship_beats_the_mean(self, mod, cohort):
        rows, keys = cohort
        arm, _ = mod.score_arm(rows, keys, "C1_plasma_metabolite")
        base, _ = mod.score_arm(rows, keys, "C0_mean")
        assert arm["mae_log10"] < base["mae_log10"], (
            "the design cannot find a relationship this strong, so a null "
            "result on the real data would say nothing")

    def test_the_baseline_arm_is_the_mean_and_has_no_predictor(self, mod):
        assert mod.ARMS["C0_mean"] == ()

    def test_only_one_predictor_is_used(self, mod):
        assert len(mod.ARMS["C1_plasma_metabolite"]) == 1


# --------------------------------------------------------------------------
# 4. the published anchors gate the extraction
# --------------------------------------------------------------------------
class TestAnchors:
    def _real(self):
        # The dose and concentration columns as published, used to check the
        # anchor logic itself rather than the parser.
        raw = [("1", 1, 400, 2943), ("1", 2, 400, 642), ("2", 1, 300, 1523),
               ("3", 1, 400, 404), ("3", 2, 400, 118), ("4", 1, 400, 1500),
               ("4", 2, 400, 619), ("5", 1, 400, 1670), ("6", 1, 400, 3160),
               ("6", 2, 400, 1760), ("7", 1, 400, 262), ("7", 2, 400, 160),
               ("8", 1, 400, 481), ("8", 2, 400, 1040), ("9", 1, 400, 376),
               ("9", 2, 300, 774), ("10", 1, 400, 438), ("10", 2, 400, 930),
               ("11", 1, 400, 383), ("11", 2, 400, 807), ("12", 1, 400, 2587),
               ("12", 2, 200, 1697)]
        return {(p, v): {"dose_mg_per_day": float(d),
                         "concentration_nm": float(c)}
                for p, v, d, c in raw}

    def test_the_published_figures_reproduce(self, mod):
        out = mod.verify_published_anchors(self._real())
        assert out["fold_variation_recomputed"] == pytest.approx(28.8, abs=0.1)
        assert out["minimum_recomputed"] == pytest.approx(1.02, abs=0.15)
        assert out["maximum_recomputed"] == pytest.approx(29.4, abs=0.15)
        assert out["highest_two_recomputed"] == [["3", 2], ["7", 2]]

    def test_a_shifted_concentration_column_is_caught(self, mod):
        rows = self._real()
        keys = sorted(rows)
        values = [rows[k]["concentration_nm"] for k in keys]
        for key, value in zip(keys, values[1:] + values[:1]):
            rows[key]["concentration_nm"] = value
        with pytest.raises(mod.VenetoclaxContrastError):
            mod.verify_published_anchors(rows)

    def test_it_says_what_the_anchor_does_not_prove(self, mod):
        out = mod.verify_published_anchors(self._real())
        assert "does NOT validate the article's clearance method" in \
            out["what_this_proves"]


# --------------------------------------------------------------------------
# 5. the refusals
# --------------------------------------------------------------------------
class TestRefusals:
    def test_it_says_the_baseline_is_not_clinical(self, mod):
        text = mod.what_this_cannot_claim()["the_baseline_is_not_clinical"]
        assert "by sex group only" in text
        assert "not over clinical knowledge" in text

    def test_it_says_the_endpoint_is_not_an_auc_or_a_clearance(self, mod):
        text = mod.what_this_cannot_claim()["not_an_auc_or_a_clearance"]
        assert "reciprocal" in text

    def test_it_forbids_joining_the_in_vitro_donors(self, mod):
        text = mod.what_this_cannot_claim()[
            "the_in_vitro_donors_are_different_people"]
        assert "no analysis" in text and "may join them" in text

    def test_it_records_that_no_genotype_exists_for_these_patients(self, mod):
        text = mod.what_this_cannot_claim()[
            "not_a_pharmacogenetic_comparison"]
        assert "cannot be formed" in text


# --------------------------------------------------------------------------
# 6. the artefact and the documents
# --------------------------------------------------------------------------
def _normalise(text):
    return re.sub(r"[\s*_`>|]+", " ", text)


class TestArtefactAndDocuments:
    def test_both_documents_exist(self):
        assert PRESPEC.is_file() and RESULTS.is_file()

    def test_the_prespecification_admits_the_answer_was_known(self):
        text = _normalise(PRESPEC.read_text(encoding="utf-8"))
        assert "Not blind" in text
        assert "confirmation with controls" in text

    def test_the_prespecification_declares_the_bias_direction(self):
        text = _normalise(PRESPEC.read_text(encoding="utf-8"))
        assert "biased in favour of the biomarker" in text

    def test_no_patient_level_data_is_exported(self, payload):
        assert payload[
            "no_patient_identifiers_values_or_predictions_exported"] is True

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

    def test_the_per_patient_table_is_outside_the_repository(self, payload):
        where = pathlib.Path(payload[
            "per_patient_table_is_outside_this_repository"])
        assert ROOT not in where.parents and ROOT != where.parent

    def test_the_excluded_observations_are_counted_not_hidden(self, payload):
        assert isinstance(payload[
            "observations_excluded_by_the_declared_rule"], list)
        assert payload["n_observations_all_visits"] < 24, (
            "two concentrations are not determined in the source; if all "
            "24 survived, the missing-data rule did not fire")

    def test_the_primary_has_twelve_patients(self, payload):
        assert payload["n_patients_primary"] == 12

    def test_both_analyses_are_not_supported(self, payload):
        for key in ("primary_visit_1_only", "sensitivity_both_visits"):
            assert payload[key]["verdict"] == "NOT_SUPPORTED"

    def test_the_refitting_interval_governs(self, payload):
        for key in ("primary_visit_1_only", "sensitivity_both_visits"):
            block = payload[key]
            if block["verdict"] == "SUPPORTED":
                assert block["refitting_patient_level_bootstrap"][
                    "interval_excludes_zero"] is True

    def test_the_results_document_carries_the_artefact_figures(self, payload):
        text = RESULTS.read_text(encoding="utf-8")
        for key in ("primary_visit_1_only", "sensitivity_both_visits"):
            block = payload[key]
            for side in ("arm_score", "comparator_score"):
                printed = f"{block[side]['geometric_fold_error']:.3f}"
                assert printed in text, (
                    f"{printed} is in the artefact and not in the document")

    def test_the_results_document_does_not_overclaim_the_negative(self):
        text = _normalise(RESULTS.read_text(encoding="utf-8"))
        assert "a plasma CYP3A marker is useless for venetoclax exposure | not established either" \
            in text or "not established either" in text

    def test_the_results_document_names_the_ariana_boundary(self):
        text = _normalise(RESULTS.read_text(encoding="utf-8"))
        assert "does not generalise on the evidence available" in text

    def test_the_hypothesis_about_marker_class_is_labelled_a_hypothesis(self):
        text = _normalise(RESULTS.read_text(encoding="utf-8"))
        assert "stated as a hypothesis and not a finding" in text
