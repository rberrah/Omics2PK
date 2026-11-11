"""Drift tests for the H1-H5 state, the method audit and the refutation protocol.

These three documents decide whether a method gets built. The risks they exist
to control are (a) a figure drifting away from the result it came from, and
(b) the architecture quietly becoming unfalsifiable. Both are tested.

Nothing here executes Fusing; there is nothing to execute yet, by design.
"""
from __future__ import annotations

import json
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
STATE = ROOT / "docs" / "OMICS2PK_H1_H5_SCIENTIFIC_STATE_20261110.md"
AUDIT = (ROOT / "docs"
         / "OMICS2PK_METHOD_AUDIT_AND_HYBRID_FUSING_ASSESSMENT_20261110.md")
PROTOCOL = ROOT / "docs" / "OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md"
CONCEPT = ROOT / "docs" / "HYBRID_FUSING_CONCEPT.md"
PROTOCOL_FROZEN = (ROOT / "docs"
                   / "OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md")
SIZING = ROOT / "reports" / "data" / "cohort_size_simulation.json"
SELECTION = ROOT / "reports" / "data" / "selection_criterion_diagnostic.json"
B2 = ROOT / "reports" / "data" / "plasma_biomarker_validation.json"
REGISTRY = ROOT / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"

SRC = ROOT / "src" / "omics2pk"


def _flat(path):
    return re.sub(r"[\s*_`>|]+", " ", path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def state():
    return STATE.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def audit():
    return AUDIT.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def protocol():
    return PROTOCOL.read_text(encoding="utf-8")


def _load(path):
    if not path.is_file():
        pytest.skip(f"{path.name} absent in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------
# 1. the three documents exist and change no frozen thing
# --------------------------------------------------------------------------
class TestNoFrozenChange:
    def test_all_three_exist(self):
        for path in (STATE, AUDIT, PROTOCOL):
            assert path.is_file(), path.name

    def test_each_declares_it_changes_no_frozen_protocol(self):
        for path in (STATE, AUDIT, PROTOCOL):
            assert "frozen" in _flat(path).lower()
            assert "is modified" in _flat(path) or "no frozen" in \
                _flat(path).lower(), path.name

    def test_the_frozen_protocol_is_still_at_its_version(self):
        text = PROTOCOL_FROZEN.read_text(encoding="utf-8")
        assert "**Version:** 1.4.10" in text, (
            "the five-objective protocol version moved; these documents claim "
            "to read it without changing it")

    def test_the_state_document_cites_that_version(self, state):
        assert "v1.4.10" in state


# --------------------------------------------------------------------------
# 2. the audit's factual claims about this repository
# --------------------------------------------------------------------------
class TestAuditIsTrue:
    def test_the_module_count_is_right(self, audit):
        count = len([p for p in SRC.glob("*.py") if p.name != "__init__.py"])
        assert f"## 1. Audit of Omics2PK — {count} modules" in audit, (
            f"the repository has {count} modules and the audit's heading "
            "says otherwise")

    def test_the_hybrid_module_exists_and_residualises_on_the_baseline(self):
        source = (SRC / "hybrid_mechanistic_selection.py").read_text(
            encoding="utf-8")
        assert "partial_rank_strength" in source
        assert "residualise" in source
        assert "permutation" in source.lower()

    def test_the_hybrid_module_allows_abstention(self):
        source = (SRC / "hybrid_mechanistic_selection.py").read_text(
            encoding="utf-8")
        assert "abstained" in source
        assert "allowed answer" in source

    def test_the_fusion_category_holds_only_what_the_audit_admits(self, audit):
        # The category was empty on 2026-11-10 and now holds constant-weight
        # stacking, inside the diagnostic module. The audit must say so: a
        # function that fits weights and scores them out of sample IS
        # stacking, whatever file it lives in.
        joined = "\n".join(
            p.read_text(encoding="utf-8", errors="replace")
            for p in SRC.glob("*.py")).lower()
        for term in ("ensemble", "model_averag", "meta_learner",
                     "mixture_of_experts"):
            assert term not in joined, (
                f"'{term}' now appears in src/; the audit's inventory is out "
                "of date and the novelty assessment must be revised")
        assert "this is stacking" in audit, (
            "nested_constant_stacking fits constant weights and scores them "
            "out of sample; the audit must not describe it as a ceiling")
        assert "a **mechanistic prior on weights** | **does not exist**" in audit

    #: Stacking is permitted in exactly these modules, both of which exist to
    #: measure it as a comparator. A third would mean the project has started
    #: building fusion rather than testing it, which is what this guard is for.
    STACKING_IS_ALLOWED_IN = {"fusion_headroom_diagnostics.py",
                              "nested_fusion_evaluation.py"}

    def test_stacking_stays_in_the_two_modules_that_measure_it(self):
        offenders = [p.name for p in SRC.glob("*.py")
                     if "stacking" in p.read_text(encoding="utf-8",
                                                  errors="replace").lower()
                     and p.name not in self.STACKING_IS_ALLOWED_IN]
        assert not offenders, (
            f"stacking has spread beyond the modules that measure it: "
            f"{offenders}")

    def test_the_audit_records_that_a_gate_now_exists_and_loses(self, audit):
        """The audit said a gate does not exist. One does now, and it loses.

        Keeping the old sentence would be the same failure the headroom
        document was corrected for: a true statement left standing after it
        stopped being true.
        """
        assert "nested_fusion_evaluation" in audit
        assert "0 of 30" in audit
        assert "a **gate**, w(z_i) | **does not exist**" not in audit

    def test_the_cached_prediction_layer_is_recorded_as_superseded(self):
        audit = AUDIT.read_text(encoding="utf-8")
        assert "Amended 2026-11-11" in audit
        assert "refits the arms inside every outer fold" in audit

    def test_the_audit_names_the_degree_strata_safeguard(self, audit):
        source = (SRC / "mechanistic_graph_selection.py").read_text(
            encoding="utf-8")
        assert "degree_strata" in source
        assert "degree strata" in audit

    def test_it_records_that_axis_a_selects_on_the_raw_outcome(self, audit):
        # A known inconsistency, which must stay disclosed rather than fixed
        # by editing a frozen run.
        assert "do not follow it" in audit
        assert "not fixed by editing frozen runs" in audit


# --------------------------------------------------------------------------
# 3. the audit's figures match their artefacts
# --------------------------------------------------------------------------
class TestAuditFigures:
    def test_the_interval_width_ratios_match_the_b2_artefact(self, audit):
        payload = _load(B2)
        ratios = []
        for block in payload["results"].values():
            for cell in block.values():
                stored = cell["stored_error_bootstrap"]["interval"]
                refit = cell["refitting_out_of_bag_bootstrap"]["interval"]
                ratios.append((refit[1] - refit[0]) / (stored[1] - stored[0]))
        assert f"{min(ratios):.2f}" in audit
        assert f"{max(ratios):.2f}" in audit

    def test_the_selection_finding_matches_its_artefact(self, audit):
        payload = _load(SELECTION)
        cell = payload["selection_quality_by_criterion"]["redundancy_0.85"]
        assert str(cell[
            "informative_fraction_selecting_on_the_outcome"]) in audit
        assert str(cell[
            "informative_fraction_selecting_on_the_residual"]) in audit


# --------------------------------------------------------------------------
# 4. the state document's figures match their artefacts
# --------------------------------------------------------------------------
class TestStateFigures:
    def test_the_dimensionality_curve_matches(self, state):
        payload = _load(SIZING)
        curve = payload["question_b_multiomic_complementarity"][
            "power_by_layer_width"]["curve"]
        for width, row in curve.items():
            assert f"{row['n_800']:.2f}" in state, (
                f"{width} power {row['n_800']:.2f} missing from the state doc")

    def test_the_b2_verdict_matches(self, state):
        payload = _load(B2)
        assert payload["n_contrasts_supported"] == 0
        assert "NOT_SUPPORTED" in state

    def test_the_four_h4_mechanisms_are_all_named(self, state):
        flat = re.sub(r"[\s*_`>|]+", " ", state)
        for phrase in ("Differential shrinkage", "Dimensionality",
                       "Selection on the outcome rewards redundancy",
                       "Interval understatement"):
            assert phrase in flat


# --------------------------------------------------------------------------
# 5. the novelty assessment stays honest
# --------------------------------------------------------------------------
class TestNoveltyHonesty:
    def test_it_states_the_hybrid_already_exists(self, audit):
        flat = re.sub(r"[\s*_`>|]+", " ", audit)
        assert "exists, is tested, and does not beat a univariate filter" in flat

    def test_it_names_the_prior_art_for_each_claimed_element(self, audit):
        for prior in ("MOFA2", "DIABLO", "stacking", "BMA",
                      "mixture of experts", "1991"):
            assert prior in audit, f"{prior} is not credited as prior art"

    def test_it_narrows_novelty_to_one_element(self, audit):
        flat = re.sub(r"[\s*_`>|]+", " ", audit)
        assert "the defensible novelty claim is narrow" in flat.lower()
        assert "Everything else is re-description" in flat

    def test_it_carries_the_metapk_negative_precedent_with_its_number(
            self, audit):
        assert "−0.855" in audit or "-0.855" in audit
        flat = re.sub(r"[\s*_`>|]+", " ", audit)
        assert "concentrated" in flat and "applicable" in flat

    def test_it_states_a_prior_on_the_outcome_before_running(self, audit):
        flat = re.sub(r"[\s*_`>|]+", " ", audit)
        assert "Stated before any run" in flat
        assert "Most likely" in flat and "Least likely" in flat

    def test_it_recommends_the_existing_method_if_the_gate_fails(self, audit):
        flat = re.sub(r"[\s*_`>|]+", " ", audit)
        assert "as the comparator and the solution" in flat or \
            "stacking with a rule filter" in flat


# --------------------------------------------------------------------------
# 6. the protocol can actually refute
# --------------------------------------------------------------------------
class TestRefutability:
    def test_it_states_a_refutation_criterion(self, protocol):
        assert "## 7. The one-sentence refutation criterion" in protocol
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "Contextual fusion is refuted on a cohort when" in flat
        assert "The mechanistic prior is refuted separately" in flat

    def test_the_criterion_rules_out_the_two_corrected_errors(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "Neither is refuted by an interval that merely contains zero"             in flat
        assert "neither is refuted by the constant-weight ceiling, which "             "cannot express a gate" in flat

    def test_there_are_two_co_primary_questions_not_one(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "Two co-primary questions, not one" in flat
        assert "X vs K" in flat and "F vs X" in flat
        assert "the first can succeed while the second fails" in flat

    def test_neither_primary_is_against_uniform_averaging(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "Neither contrast is against U" in flat

    def test_it_distinguishes_the_ceilings_and_names_the_governing_one(
            self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "L ind is the governing ceiling" in flat
        assert "L sel does not bound averaging" in flat
        assert "L const does not bound a gate" in flat

    def test_it_forbids_quoting_the_ceiling_as_a_performance(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "None may be quoted as a performance" in flat
        assert "A permissive ceiling licenses nothing" in flat

    def test_it_closes_the_unfalsifiability_loophole(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "anti-unfalsifiability rule" in flat
        assert "can never lose and therefore can never be shown to win" in flat

    def test_the_degree_permuted_prior_is_a_mandatory_control(self, protocol):
        assert "degree-permuted prior" in protocol
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "how well studied a gene is" in flat

    def test_every_gate_has_a_verdict_including_not_evaluable(self, protocol):
        for verdict in ("REFUTED_CEILING_BELOW_MARGIN", "HEADROOM_PRESENT",
                        "NO_DETECTABLE_HEADROOM_UNDERPOWERED",
                        "REFUTED_CONTEXT_UNINFORMATIVE", "NOT_EVALUABLE"):
            assert verdict in protocol

    def test_the_superseded_verdict_name_is_gone(self, protocol):
        assert "REFUTED_NO_HEADROOM" not in protocol, (
            "the v1.0.0 verdict name implied that a null interval refutes")

    def test_weight_learning_is_required_to_be_out_of_sample(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "leakage even when the base predictions are themselves "             "honestly out-of-fold" in flat

    def test_eight_negative_controls_are_declared(self, protocol):
        section = protocol.split("## 3. Gate 3")[1].split("## 4.")[0]
        rows = [line for line in section.splitlines()
                if line.startswith("| **") and "---" not in line]
        assert len(rows) >= 8, f"only {len(rows)} controls declared"

    def test_implementation_is_deferred_to_step_five(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "The first line of Fusing code belongs in step 5" in flat
        assert "Nothing in steps 1–2 is an implementation" in flat or \
            "Nothing in steps 1-2 is an implementation" in flat


# --------------------------------------------------------------------------
# 7. the data verdicts, and that nothing was overclaimed
# --------------------------------------------------------------------------
class TestDataVerdicts:
    def test_exactly_one_experiment_is_fully_evaluable(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "testable at exactly one place: E1" in flat

    def test_the_three_not_evaluable_experiments_each_give_a_reason(
            self, protocol):
        section = protocol.split("### 5.2")[1].split("### 5.3")[0]
        for label in ("E5", "E6", "E7"):
            assert label in section
        assert section.count("NOT_EVALUABLE") >= 3

    def test_the_pk_fusion_experiment_is_blocked_for_the_right_reason(
            self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "one observation per patient" in flat
        assert "credentialed access" in flat

    def test_the_two_new_sources_are_in_the_registry_as_abstain(self):
        payload = _load(REGISTRY)
        by_key = {e.get("key"): e for e in payload["datasets"]}
        for key in ("lee_venlafaxine_solanidine_20",
                    "kyoto_jmicc_cyp2c19_metabolome_686"):
            assert key in by_key, key
            assert by_key[key]["status"] == "ABSTAIN"
            assert by_key[key]["abstain_reason"] == \
                "ABSTAIN_NO_PER_DONOR_VALUES"

    def test_the_venlafaxine_cohort_is_flagged_as_the_closest_analogue(self):
        payload = _load(REGISTRY)
        entry = next(e for e in payload["datasets"]
                     if e.get("key") == "lee_venlafaxine_solanidine_20")
        assert any("CLOSEST STRUCTURAL ANALOGUE" in c
                   for c in entry["caveats"])

    def test_the_metabolome_endpoint_is_flagged_as_a_genotype_class(self):
        payload = _load(REGISTRY)
        entry = next(e for e in payload["datasets"]
                     if e.get("key") == "kyoto_jmicc_cyp2c19_metabolome_686")
        assert entry["endpoint"]["endpoint_class"] == "biomarker_only"
        assert any("GENOTYPE CLASS, NOT A FUNCTION" in c
                   for c in entry["caveats"])

    def test_the_search_targets_for_next_time_are_recorded(self, protocol):
        assert "### 5.4 The search continues, with these targets" in protocol


# --------------------------------------------------------------------------
# 7b. the protocol and the concept document stay consistent
# --------------------------------------------------------------------------
class TestConceptConsistency:
    def test_the_concept_document_exists(self):
        assert CONCEPT.is_file(), (
            "docs/HYBRID_FUSING_CONCEPT.md states the decision philosophy "
            "this protocol is the experiment for")

    def test_the_protocol_cross_references_it(self, protocol):
        assert "HYBRID_FUSING_CONCEPT.md" in protocol
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "the concept document governs the philosophy" in flat

    def test_both_refuse_to_reinvent_the_failed_weighting(self, protocol):
        concept = CONCEPT.read_text(encoding="utf-8")
        assert "Do not reinvent pi_m under the Fusing name" in concept
        assert "reinvent" in protocol

    def test_both_treat_fusing_as_a_label_not_a_claim(self, audit):
        concept = re.sub(r"[\s*_`>|]+", " ", CONCEPT.read_text(encoding="utf-8"))
        assert "not an originality or performance claim" in concept
        assert "re-description" in audit

    def test_both_require_selection_on_the_residual(self, audit):
        concept = re.sub(r"[\s*_`>|]+", " ", CONCEPT.read_text(encoding="utf-8"))
        assert "baseline residuals" in concept
        assert "residual after the baseline" in audit

    def test_both_permit_not_evaluable_as_a_result(self, protocol):
        concept = CONCEPT.read_text(encoding="utf-8")
        assert "NOT_EVALUABLE" in concept
        assert "NOT_EVALUABLE" in protocol

    def test_both_separate_density_fusion_from_equation_fusion(self, protocol):
        concept = re.sub(r"[\s*_`>|]+", " ", CONCEPT.read_text(encoding="utf-8"))
        assert "not fusing structural ODEs" in concept
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "distinct in kind from distribution fusion" in flat


# --------------------------------------------------------------------------
# 8. no code was imported from the other repositories
# --------------------------------------------------------------------------
class TestNoImport:
    def test_no_module_imports_from_the_forbidden_repositories(self):
        for path in SRC.glob("*.py"):
            text = path.read_text(encoding="utf-8", errors="replace")
            for line in text.splitlines():
                stripped = line.strip()
                if not (stripped.startswith("import ")
                        or stripped.startswith("from ")):
                    continue
                for forbidden in ("pkdt", "metapk", "netpk", "newpk",
                                  "thesis"):
                    assert forbidden not in stripped.lower(), (
                        f"{path.name} imports {forbidden!r}: {stripped}")

    def test_the_audit_declares_the_read_only_relationship(self, audit):
        flat = re.sub(r"[\s*_`>|]+", " ", audit)
        assert "No code was imported from Thesis, NetPK, NewPK or MetaPK" in flat
        assert "nothing from them is copied, adapted or depended upon" in flat

    def test_the_protocol_plans_to_rederive_rather_than_import(self, protocol):
        flat = re.sub(r"[\s*_`>|]+", " ", protocol)
        assert "re-derive independently" in flat
        assert "No code is imported" in flat
