"""Contracts for the catalogue and the hypothesis-dataset matrix.

The matrix is only worth having if it cannot over-claim. The three ways it
could, in order of how much damage they would do:

1. **A functional measurement read as a PK parameter.** An intrinsic clearance
   in microsomes has clearance units and is measured with no absorption, no
   binding and no blood flow. -> `TestPkAndFunctionalAreDisjoint`.
2. **Pairing inferred instead of asserted.** Two layers from different
   publications are not paired. -> `TestPairingIsNeverInferred`.
3. **The legacy H1-H5 tags reused.** The registry's H4 is "when omics become
   harmful"; this programme's H4 is longitudinal. Reusing the tag would claim
   six longitudinal datasets where there are none.
   -> `TestTheTwoVocabulariesStayApart`.

Plus the error this file was written after finding: one `comparator_available`
flag cannot answer five different questions.
-> `TestComparatorIsPerHypothesis`.
"""
from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def cat():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import input_catalogue
    return input_catalogue


@pytest.fixture(scope="module")
def matrix():
    path = ROOT / "reports" / "data" / "hypothesis_dataset_matrix.json"
    if not path.is_file():
        pytest.skip("the matrix artefact has not been generated")
    return json.loads(path.read_text(encoding="utf-8"))


class TestPkAndFunctionalAreDisjoint:
    def test_the_two_vocabularies_do_not_overlap(self, cat):
        assert not (set(cat.PK_PARAMETERS) & set(cat.FUNCTIONAL_ENDPOINTS))

    def test_an_in_vitro_clearance_is_functional_not_pk(self, cat):
        assert cat.endpoint_kind("clint_in_vitro") == "functional"
        assert cat.endpoint_kind("CL") == "pk"

    def test_metabolic_clearance_is_declared_a_component_not_cl(self, cat):
        assert "not CL" in cat.PK_PARAMETERS["CL_metabolic"]

    def test_the_guard_refuses_to_bless_a_real_pk_parameter(self, cat):
        with pytest.raises(Exception):
            cat.assert_not_a_pk_parameter("CL")

    def test_the_guard_accepts_a_functional_measurement(self, cat):
        out = cat.assert_not_a_pk_parameter("clint_in_vitro")
        assert out["is_pk_parameter"] is False
        assert "no absorption" in out["what_it_is"]

    def test_an_undeclared_endpoint_is_refused(self, cat):
        with pytest.raises(Exception):
            cat.endpoint_kind("some_new_thing")

    def test_a_pk_hypothesis_refuses_a_functional_endpoint(self, cat):
        status, reasons = cat.feasibility(
            "H1", endpoint="clint_in_vitro", pairing="same_individuals",
            access="open_public", layers_available=("proteomics",),
            comparator_available=True, n_units=100)
        assert status == "NOT_EVALUABLE"
        assert any("needs a PK parameter" in r for r in reasons)

    def test_a_functional_hypothesis_refuses_a_pk_endpoint(self, cat):
        status, reasons = cat.feasibility(
            "H3", endpoint="AUC", pairing="same_individuals",
            access="open_public", layers_available=("proteomics",),
            comparator_available=True, n_units=100)
        assert status == "NOT_EVALUABLE"
        assert any("belongs to H1, H4 or H5" in r for r in reasons)

    def test_the_matrix_fails_h1_mostly_for_this_reason(self, matrix):
        h1 = [r for r in matrix["matrix"] if r["hypothesis"] == "H1"]
        wrong_kind = [r for r in h1
                      if any("needs a PK parameter" in reason
                             for reason in r["reasons"])]
        assert len(wrong_kind) >= 10, len(wrong_kind)


class TestPairingIsNeverInferred:
    def test_unknown_pairing_blocks_joint_analysis(self, cat):
        status, reasons = cat.feasibility(
            "H3", endpoint="enzyme_activity", pairing="unknown",
            access="open_public", layers_available=("proteomics",),
            comparator_available=True, n_units=100)
        assert status == "NOT_EVALUABLE"
        assert any("not paired" in r for r in reasons)

    def test_different_publications_blocks_joint_analysis(self, cat):
        status, reasons = cat.feasibility(
            "H3", endpoint="enzyme_activity",
            pairing="different_publications", access="open_public",
            layers_available=("proteomics",), comparator_available=True,
            n_units=100)
        assert status == "NOT_EVALUABLE"

    def test_every_matrix_row_carries_a_pairing_provenance(self, matrix):
        for row in matrix["matrix"]:
            assert row["pairing_provenance"], row["dataset"]
            assert len(row["pairing_provenance"]) > 20, row["dataset"]

    def test_the_registry_has_no_pairing_field_to_infer_from(self):
        registry = json.loads(
            (ROOT / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"
             ).read_text(encoding="utf-8"))
        for dataset in registry["datasets"]:
            assert "pairing" not in dataset, dataset["key"]

    def test_h4_refuses_a_cross_sectional_dataset(self, cat):
        status, reasons = cat.feasibility(
            "H4", endpoint="AUC", pairing="same_individuals",
            access="open_public", layers_available=("proteomics",),
            comparator_available=True, n_units=100, longitudinal=False)
        assert status == "NOT_EVALUABLE"
        assert any("repeated occasions" in r for r in reasons)


class TestTheTwoVocabulariesStayApart:
    def test_the_legacy_labels_mean_something_else(self, cat):
        assert "harmful" in cat.LEGACY_OBJECTIVES["H4"].lower()
        assert "longitudinal" not in cat.LEGACY_OBJECTIVES["H4"].lower()
        assert "robustness" in cat.LEGACY_OBJECTIVES["H5"].lower()

    def test_three_legacy_labels_have_no_mapping(self, cat):
        unmapped = [k for k, (target, _) in cat.LEGACY_CROSSWALK.items()
                    if target is None]
        assert set(unmapped) == {"H3", "H4", "H5"}

    def test_asking_without_a_vocabulary_is_refused(self, cat):
        with pytest.raises(Exception):
            cat.assert_hypothesis_vocabulary("H4", programme="whatever")

    def test_a_legacy_lookup_says_it_does_not_map(self, cat):
        out = cat.assert_hypothesis_vocabulary("H4", programme="legacy")
        assert out["maps_to_current"] is None
        assert "would claim" in out["why"]

    def test_the_matrix_does_not_reuse_the_registry_tags(self):
        script = (ROOT / "scripts"
                  / "build_hypothesis_dataset_matrix.py").read_text(
                      encoding="utf-8")
        assert "objectives_supported" not in script
        assert "does not reuse the registry" in script


class TestComparatorIsPerHypothesis:
    """One flag cannot answer five questions, and once it tried."""

    def test_the_declared_comparators_are_dictionaries(self):
        script = (ROOT / "scripts"
                  / "build_hypothesis_dataset_matrix.py").read_text(
                      encoding="utf-8")
        assert '"comparator_available": {"H1"' in script

    def test_ariana_is_not_ready_for_dose_personalisation(self, matrix):
        rows = [r for r in matrix["matrix"]
                if r["dataset"] == "ariana_metoprolol_47"]
        by_hypothesis = {r["hypothesis"]: r["status"] for r in rows}
        assert by_hypothesis["H1"] == "READY"
        assert by_hypothesis["H5"] == "NOT_EVALUABLE"

    def test_two_candidates_out_of_two_is_not_selection(self, cat, matrix):
        assert cat.HYPOTHESES["H2"]["needs_candidate_pool_of_at_least"] >= 3
        status, reasons = cat.feasibility(
            "H2", endpoint="concentration", pairing="same_individuals",
            access="open_public", layers_available=("metabolomics",),
            comparator_available=True, n_units=47, n_candidate_omics=2)
        assert status == "NOT_EVALUABLE"
        assert any("is not selection" in r for r in reasons)

    def test_ariana_is_not_ready_for_selection(self, matrix):
        rows = {r["hypothesis"]: r for r in matrix["matrix"]
                if r["dataset"] == "ariana_metoprolol_47"}
        assert rows["H2"]["status"] == "NOT_EVALUABLE"


class TestMethodsAreNotSources:
    def test_diablo_is_refused_as_a_source(self, cat):
        assert cat.source_kind("DIABLO") == "method"
        with pytest.raises(Exception) as caught:
            cat.assert_is_a_database("DIABLO")
        assert "not a source of biological information" in str(caught.value)

    def test_a_database_is_accepted(self, cat):
        for name in ("STRING", "Reactome", "UniProt", "MaveDB", "PRIDE"):
            assert cat.assert_is_a_database(name) is True

    def test_annotation_tools_are_their_own_kind(self, cat):
        assert cat.source_kind("EnsemblVEP") == "annotation_tool"

    def test_string_confidence_is_disclaimed(self, cat):
        assert "NOT causality" in cat.SOURCES["STRING"][1]


class TestMechanismsAreNotEnzymeOnly:
    def test_transporters_binding_and_absorption_are_present(self, cat):
        for key in ("efflux_transporters", "uptake_transporters",
                    "binding_proteins", "absorption_machinery"):
            assert key in cat.MECHANISTIC_KNOWLEDGE, key

    def test_each_class_names_the_pk_parameters_it_may_inform(self, cat):
        for key, (_, informs) in cat.MECHANISTIC_KNOWLEDGE.items():
            for parameter in informs:
                assert parameter in cat.PK_PARAMETERS, (key, parameter)

    def test_a_pharmacological_target_informs_no_pk_parameter(self, cat):
        assert cat.MECHANISTIC_KNOWLEDGE["targets"][1] == ()

    def test_absorption_informs_f_and_ka_not_clearance(self, cat):
        informs = cat.MECHANISTIC_KNOWLEDGE["absorption_machinery"][1]
        assert set(informs) == {"F", "ka"}


class TestTheMatrixIsHonestAboutWhatIsMissing:
    def test_h4_and_h5_are_empty(self, matrix):
        summary = matrix["summary_by_hypothesis"]
        for hypothesis in ("H4", "H5"):
            assert summary[hypothesis]["READY"] == []
            assert summary[hypothesis]["PARTIAL"] == []
            assert summary[hypothesis]["ACCESS_REQUIRED"] == []

    def test_h1_has_exactly_one_vehicle(self, matrix):
        assert matrix["summary_by_hypothesis"]["H1"]["READY"] == \
            ["ariana_metoprolol_47"]

    def test_every_row_explains_its_status(self, matrix):
        for row in matrix["matrix"]:
            assert row["reasons"], (row["dataset"], row["hypothesis"])

    def test_negative_searches_are_kept(self, matrix):
        geo = [s for s in matrix["searches_run_today"]
               if s["repository"] == "GEO"]
        assert len(geo) >= 3
        assert any(s["on_target"] == 0 for s in geo)

    def test_the_run_declares_that_no_experiment_happened(self, matrix):
        assert matrix["no_experiment_was_run"] is True

    def test_the_imatinib_candidate_records_what_is_unverified(self, matrix):
        candidate = matrix["new_candidates"][
            "pride_pxd066482_imatinib_clearance"]
        assert "NOT yet verified as public" in candidate["what_is_missing"]
        assert candidate["licence"] == "CC0"

    def test_the_document_states_h1_has_one_consumed_vehicle(self):
        text = (ROOT / "docs"
                / "OMICS2PK_CATALOGUE_ET_ENDPOINTS_VERIFIABLES_20261112.md"
                ).read_text(encoding="utf-8")
        assert "qu'un seul véhicule public et il est déjà consommé" in text
        assert "aucun endpoint de transporteur" in text


class TestAxisCoverage:
    """The four gaps the programme is short of data on, counted not asserted.

    "We have no transporter data" is a claim that must be re-derivable from
    the matrix at any time, so coverage is computed. The first version counted
    hypothesis-by-dataset ROWS and multiplied every dataset by five, which
    reported coverage that does not exist.
    """

    def test_the_axes_name_the_pk_parameters_they_inform(self, cat):
        for axis, spec in cat.AXES.items():
            assert spec["informs"], axis
            for parameter in spec["informs"]:
                assert parameter in cat.PK_PARAMETERS, (axis, parameter)

    def test_only_one_axis_informs_metabolic_clearance(self, cat):
        informing = [a for a, s in cat.AXES.items()
                     if s["informs"] == ("CL_metabolic",)]
        assert informing == ["metabolic_clearance"]

    def test_absorption_informs_f_and_ka_and_not_clearance(self, cat):
        assert set(cat.AXES["absorption"]["informs"]) == {"F", "ka"}

    def test_every_mechanism_class_maps_to_an_axis_or_to_none(self, cat):
        for mechanism in cat.MECHANISTIC_KNOWLEDGE:
            axis = cat.axis_of(mechanism)
            assert axis is None or axis in cat.AXES, mechanism

    def test_a_pharmacological_target_belongs_to_no_axis(self, cat):
        assert cat.axis_of("targets") is None

    def test_coverage_counts_datasets_not_rows(self, cat):
        rows = [
            {"dataset": "a", "axes": ("transport",), "endpoint_kind":
             "functional", "status": "NOT_EVALUABLE", "hypothesis": "H1"},
            {"dataset": "a", "axes": ("transport",), "endpoint_kind":
             "functional", "status": "NOT_EVALUABLE", "hypothesis": "H2"},
            {"dataset": "b", "axes": ("transport",), "endpoint_kind": "pk",
             "status": "READY", "hypothesis": "H1"},
        ]
        out = cat.axis_coverage(rows)
        assert out["transport"]["n_datasets"] == 2
        assert out["transport"]["n_with_a_pk_endpoint"] == 1
        assert out["transport"]["n_ready_for_any_hypothesis"] == 1

    def test_the_transport_axis_has_no_pk_endpoint_in_the_matrix(self,
                                                                 matrix):
        for axis in ("transport", "binding_distribution", "absorption",
                     "temporal"):
            block = matrix["axis_coverage"][axis]
            assert block["n_with_a_pk_endpoint"] == 0, axis
            assert block["n_ready_for_any_hypothesis"] == 0, axis

    def test_the_absorption_axis_is_no_longer_empty(self, matrix):
        """The search did fill the measurement side of two gaps."""
        assert matrix["axis_coverage"]["absorption"]["n_datasets"] >= 2
        assert matrix["axis_coverage"]["transport"]["n_datasets"] >= 2

    def test_a_transporter_abundance_does_not_claim_to_inform_f(self, cat):
        assert "does not inform F or ka by itself" in \
            cat.AXES["transport"]["note"]

    def test_the_negative_searches_are_recorded_with_their_reason(self,
                                                                   matrix):
        searches = matrix["searches_2026_11_12"]
        assert len(searches) >= 10
        structural = [s for s in searches
                      if "STRUCTURAL NEGATIVE" in s["verdict"]]
        assert structural, "the PBPK-consumer finding must be recorded"
        tacrolimus = [s for s in searches if s["query"] == "tacrolimus"]
        assert tacrolimus and tacrolimus[0]["on_target"] == 0

    def test_the_longitudinal_design_exists_without_its_endpoint(self,
                                                                 matrix):
        candidate = matrix["new_candidates"][
            "pride_pxd011498_jejunum_before_after_gastrectomy"]
        assert candidate["longitudinal"] is True
        assert "no PK on either occasion" in candidate["what_is_missing"]

    def test_the_document_records_the_two_structural_negatives(self):
        """A negative nobody wrote down gets re-searched by the next run."""
        text = (ROOT / "docs"
                / "OMICS2PK_CATALOGUE_ET_ENDPOINTS_VERIFIABLES_20261112.md"
                ).read_text(encoding="utf-8")
        assert "saturé par les consommateurs" in text
        assert "Aucune omique de tacrolimus dans PRIDE" in text

    def test_the_document_states_the_programme_question_it_raises(self):
        """H1 and H4 lack a DESIGN, which is not an absence of effect."""
        text = (ROOT / "docs"
                / "OMICS2PK_CATALOGUE_ET_ENDPOINTS_VERIFIABLES_20261112.md"
                ).read_text(encoding="utf-8")
        assert "H1 et H4 sont inatteignables" in text
        assert "l'absence du **dessin expérimental**" in text

    def test_the_intestinal_source_is_a_reference_layer_not_a_cohort(self,
                                                                     matrix):
        """An abundance predicting another abundance is axis A, not H3."""
        candidate = matrix["new_candidates"][
            "pride_pxd020987_human_small_intestine"]
        assert "REFERENCE layer" in candidate["what_is_missing"]
        assert "no independent "               "functional endpoint" in candidate["what_is_missing"]
        rows = [r for r in matrix["matrix"]
                if r["dataset"] == "pride_pxd020987_human_small_intestine"]
        assert rows and all(r["status"] == "NOT_EVALUABLE" for r in rows)
