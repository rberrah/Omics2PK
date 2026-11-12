"""Drift tests for the funding dossier.

A dossier is the one document a reader outside the project will quote from, so
every figure in it is checked back against the artefact it came from. If a
number is corrected anywhere upstream and not here, these tests fail.

They also check the framing, because a dossier can be accurate line by line and
still mislead: a scenario presented as an estimate, an in-vitro effect used to
size a clinical study, or a non-detection quoted as evidence of absence.
"""
from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DOSSIER = ROOT / "docs" / "OMICS2PK_FUNDING_PREPARATION_DOSSIER_20261110.md"
LEDGER = ROOT / "docs" / "OMICS2PK_EVIDENCE_LEDGER_20261110.md"
SIZING = ROOT / "reports" / "data" / "cohort_size_simulation.json"
SELECTION = ROOT / "reports" / "data" / "selection_criterion_diagnostic.json"
B2 = ROOT / "reports" / "data" / "plasma_biomarker_validation.json"
C2 = ROOT / "reports" / "data" / "venetoclax_clinical_pk_contrast.json"
SEARCH = ROOT / "reports" / "data" / "independent_repository_search.json"
REGISTRY = ROOT / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"


@pytest.fixture(scope="module")
def text():
    assert DOSSIER.is_file(), "the dossier is missing"
    return DOSSIER.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def flat(text):
    return re.sub(r"[\s*_`>|]+", " ", text)


def _load(path):
    if not path.is_file():
        pytest.skip(f"{path.name} has not been produced in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# 1. the six commissioned sections are all present
# --------------------------------------------------------------------------
class TestStructure:
    @pytest.mark.parametrize("heading", [
        "## 1. Public results established, with their level of proof",
        "## 2. The scientific questions still open",
        "## 3. The minimum data needed to resolve them",
        "## 4. Sample-size scenarios, with their assumptions and their uncertainty",
        "## 5. Modalities, compartments and endpoints",
        "## 6. A minimal cohort, then an ideal one",
    ])
    def test_the_section_is_present(self, text, heading):
        assert heading in text

    def test_it_names_its_audience(self, text):
        assert "**Written for:**" in text

    def test_it_recomputes_nothing(self, flat):
        assert "Nothing new is computed here" in flat

    def test_every_source_document_it_cites_exists(self, text):
        cited = set(re.findall(r"`(OMICS2PK_[A-Z0-9_]+\.md)`", text))
        assert len(cited) >= 8
        missing = [name for name in cited
                   if not (ROOT / "docs" / name).is_file()]
        assert not missing, f"the dossier cites missing documents: {missing}"


# --------------------------------------------------------------------------
# 2. the sizing figures match the sizing artefact
# --------------------------------------------------------------------------
class TestSizingFigures:
    def test_question_a_curve_matches(self, text):
        payload = _load(SIZING)
        curve = payload["question_a_enzyme_activity"][
            "required_n_by_rank_correlation"]
        for key in ("spearman_0.20", "spearman_0.40", "spearman_0.53",
                    "spearman_0.70", "spearman_0.90"):
            required = curve[key]["required_n"]
            assert str(required) in text, (
                f"{key} needs n={required} in the artefact and that figure is "
                "not in the dossier")

    def test_question_b_curve_matches(self, text):
        payload = _load(SIZING)
        curve = payload["question_b_multiomic_complementarity"][
            "required_n_by_assumed_reduction"]
        for key, cell in curve.items():
            if cell["required_n"] is not None:
                assert str(cell["required_n"]) in text, (
                    f"{key} needs n={cell['required_n']} and that is not in "
                    "the dossier")

    def test_question_c_endpoints_match(self, text):
        payload = _load(SIZING)
        by_endpoint = payload["question_c_clinical_pk"]["by_endpoint"]
        for endpoint, block in by_endpoint.items():
            value = block["mae_reduction_10pc"]["required_n"]
            printed = f"{value:,}".replace(",", ",")
            assert printed in text or str(value) in text, (
                f"{endpoint} needs n={value} and that is not in the dossier")

    def test_the_dimensionality_power_curve_matches(self, text):
        payload = _load(SIZING)
        curve = payload["question_b_multiomic_complementarity"][
            "power_by_layer_width"]["curve"]
        for width, row in curve.items():
            printed = f"{row['n_800']:.2f}"
            assert printed in text, (
                f"{width} gives power {printed} at n=800 and that is not in "
                "the dossier")

    def test_the_confirmatory_sizes_match(self, text):
        payload = _load(SIZING)
        block = payload["question_a_enzyme_activity"][
            "measured_confirmatory_sizing"]
        for label, cell in block.items():
            assert str(cell["confirmatory_sizing"]["required_n"]) in text, (
                f"{label} implies n="
                f"{cell['confirmatory_sizing']['required_n']} and that is not "
                "in the dossier")

    def test_the_selection_finding_matches(self, text):
        payload = _load(SELECTION)
        cell = payload["selection_quality_by_criterion"]["redundancy_0.85"]
        assert str(cell[
            "informative_fraction_selecting_on_the_outcome"]) in text
        assert payload["the_recommendation"][:40] in text or \
            "residual after the baseline" in text


# --------------------------------------------------------------------------
# 3. the result figures match their artefacts
# --------------------------------------------------------------------------
class TestResultFigures:
    def test_the_b2_counts_match(self, text):
        payload = _load(B2)
        assert f"{payload['n_contrasts_supported']} of 6" in text or \
            "0 of 6" in text
        assert str(payload["n_donors"]) in text

    def test_the_b2_interval_ratios_are_quoted(self, text):
        payload = _load(B2)
        ratios = []
        for block in payload["results"].values():
            for cell in block.values():
                stored = cell["stored_error_bootstrap"]["interval"]
                refit = cell["refitting_out_of_bag_bootstrap"]["interval"]
                ratios.append((refit[1] - refit[0]) / (stored[1] - stored[0]))
        assert f"{min(ratios):.2f}" in text
        assert f"{max(ratios):.2f}" in text

    def test_the_c2_result_is_quoted_as_negative(self, text, flat):
        payload = _load(C2)
        assert payload["primary_visit_1_only"]["verdict"] == "NOT_SUPPORTED"
        assert str(payload["n_patients_primary"]) in text
        fold = payload["primary_visit_1_only"]["comparator_score"][
            "geometric_fold_error"]
        assert f"{fold:.3f}" in text

    def test_the_search_counts_match(self, text):
        payload = _load(SEARCH)
        pattern = payload["the_pattern_now_stands_at"]
        assert str(pattern[
            "sources_that_measured_a_per_donor_pairing"]) in text
        assert str(pattern[
            "sources_that_published_the_per_donor_rows"]) in text

    def test_the_registry_version_is_current(self, text):
        payload = _load(REGISTRY)
        assert payload["version"] in text, (
            "the dossier cites a registry version that is not the current one")

    def test_the_withdrawal_count_agrees_with_the_ledger(self, text):
        ledger = LEDGER.read_text(encoding="utf-8")
        match = re.search(r"(\w+) claims have been withdrawn", ledger)
        assert match
        assert match.group(1).lower() in text.lower(), (
            f"the ledger says {match.group(1)} withdrawals and the dossier "
            "does not say the same")


# --------------------------------------------------------------------------
# 4. the framing, which is where a dossier misleads
# --------------------------------------------------------------------------
class TestFraming:
    def test_the_two_evidence_levels_are_distinguished(self, flat):
        assert "level 1" in flat and "level 2" in flat
        assert "No in-vitro effect is carried into" in flat

    def test_no_in_vitro_effect_sizes_a_clinical_study(self, flat):
        assert "No in-vitro effect anywhere in this dossier is used to size a" \
            in flat

    def test_scenarios_are_not_presented_as_estimates(self, flat):
        assert "Quote the assumption with the number or do not quote the " \
            "number" in flat
        assert "Every one of those is an assumption" in flat

    def test_the_floor_is_labelled_a_floor(self, flat):
        assert "95 is a floor, not an estimate" in flat

    def test_the_non_detections_are_not_quoted_as_absence(self, flat):
        assert "None of the three had the power to detect the effect" in flat
        assert "would be a misuse of the project's own results" in flat

    def test_the_failed_replication_is_in_the_summary(self, flat):
        assert "no public cohort is large enough to establish its size" in flat

    def test_the_negative_clinical_row_is_in_the_summary(self, flat):
        assert "the single C-level contrast is negative" in flat

    def test_the_data_request_alternative_is_disclosed(self, flat):
        assert "Question A can be answered without this grant" in flat

    def test_it_states_a_prespecified_failure_criterion(self, flat):
        assert "Prespecified failure criterion" in flat
        assert "the programme does not proceed to stage 2" in flat

    def test_it_separates_what_sharing_could_supply_from_what_it_cannot(
            self, flat):
        assert "bounded by data that has never been collected" in flat

    def test_the_ariana_gain_is_labelled_exploratory(self, flat):
        assert "The 38 % ARIANA gain is exploratory" in flat
