"""Contract tests for the candidate paired-liver-source qualification.

The failure this module guards against is the one a source audit is most prone
to: a verdict that reads as a measurement but is really an assumption. So the
tests check that the rule is applied rather than narrated, that each verdict's
reason is the one the rule returns, that the access verdicts rest on fields
that get re-queried, and that the per-donor verdicts name a file and a
checksum rather than an impression.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "audit_candidate_paired_liver_sources.py"
REPORT = REPO / "reports" / "data" / "candidate_paired_liver_sources.json"
DOCUMENT = (
    REPO / "docs" / "OMICS2PK_CANDIDATE_PAIRED_LIVER_SOURCE_QUALIFICATION_20261110.md"
)


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("candidate_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _normalise(text):
    return re.sub(r"\s+", " ", text)


# --------------------------------------------------------------------------
# 1. the rule is a rule
# --------------------------------------------------------------------------
class TestTheRule:
    def test_the_donor_floor_is_declared_and_justified(self, mod):
        assert mod.DONOR_FLOOR == 10
        offline = mod.audit(check_network=False)
        assert "under eight" in offline["rule"]["why_the_floor"]

    def test_the_admissible_outcomes_match_the_readiness_gate_vocabulary(self, mod):
        for kind in ("drug_specific_clint", "drug_specific_metabolic_activity",
                     "clinical_pk"):
            assert kind in mod.ADMISSIBLE_OUTCOMES

    def test_access_is_checked_before_anything_else(self, mod):
        """A source nobody can open cannot be assessed on its contents."""
        entry = dict(
            reason="ABSTAIN_ACCESS",
            per_donor_layer_and_outcome=True,
            n_donors=500,
            outcome_kind="clinical_pk",
        )
        verdict, reason = mod.qualify(entry)
        assert (verdict, reason) == ("ABSTAIN", "ABSTAIN_ACCESS")

    def test_aggregates_only_is_refused(self, mod):
        entry = dict(
            reason="ABSTAIN_NO_PER_DONOR_VALUES",
            per_donor_layer_and_outcome=False,
            n_donors=128,
            outcome_kind="drug_specific_metabolic_activity",
        )
        assert mod.qualify(entry) == ("ABSTAIN", "ABSTAIN_NO_PER_DONOR_VALUES")

    def test_unknown_per_donor_status_is_refused_not_assumed(self, mod):
        """`None` must not be read as `True`."""
        entry = dict(
            reason="ABSTAIN_NO_PER_DONOR_VALUES",
            per_donor_layer_and_outcome=None,
            n_donors=128,
            outcome_kind="clinical_pk",
        )
        assert mod.qualify(entry)[1] == "ABSTAIN_NO_PER_DONOR_VALUES"

    def test_a_cohort_under_the_floor_is_refused(self, mod):
        entry = dict(
            reason="ABSTAIN_TOO_FEW_DONORS",
            per_donor_layer_and_outcome=True,
            n_donors=4,
            outcome_kind="drug_specific_metabolic_activity",
        )
        assert mod.qualify(entry) == ("ABSTAIN", "ABSTAIN_TOO_FEW_DONORS")

    def test_a_source_at_the_floor_is_not_refused_for_its_size(self, mod):
        entry = dict(
            reason="ABSTAIN_NO_ADMISSIBLE_ENDPOINT",
            per_donor_layer_and_outcome=True,
            n_donors=mod.DONOR_FLOOR,
            outcome_kind=None,
        )
        assert mod.qualify(entry)[1] == "ABSTAIN_NO_ADMISSIBLE_ENDPOINT"

    def test_a_fully_qualifying_source_is_acquired(self, mod):
        entry = dict(
            reason=None,
            per_donor_layer_and_outcome=True,
            n_donors=200,
            outcome_kind="drug_specific_clint",
        )
        assert mod.qualify(entry) == ("ACQUIRED", None)

    def test_a_recorded_reason_that_contradicts_the_rule_stops_the_audit(self, mod):
        """The rule and the record must agree or the audit refuses to report."""
        entry = dict(
            label="a fabricated entry",
            reason="ABSTAIN_TOO_FEW_DONORS",
            per_donor_layer_and_outcome=False,
            n_donors=500,
            outcome_kind="clinical_pk",
        )
        with pytest.raises(mod.AuditError, match="must be resolved"):
            mod.qualify(entry)


# --------------------------------------------------------------------------
# 2. the candidates
# --------------------------------------------------------------------------
class TestTheCandidates:
    def test_there_are_seven(self, mod):
        assert len(mod.CANDIDATES) == 7

    def test_the_three_from_the_commissioning_note_are_present(self, mod):
        for key in ("drozdzik_2018_liver_intestine_9",
                    "almajdoub_2026_kcat_matched_4",
                    "zhang_2015_mppgl_128"):
            assert key in mod.CANDIDATES

    def test_every_candidate_carries_evidence_and_a_reason(self, mod):
        for key, entry in mod.CANDIDATES.items():
            assert entry["evidence"].strip(), key
            assert entry["reason"] in mod.REASONS, key
            assert entry["why_it_matters"].strip(), key
            assert entry["citation"].strip(), key

    def test_every_candidate_declares_whether_its_layers_are_same_donor(self, mod):
        for key, entry in mod.CANDIDATES.items():
            assert isinstance(entry["same_donor_layers"], bool), key

    def test_the_two_paywalled_candidates_claim_nothing_about_contents(self, mod):
        for key in ("drozdzik_2018_liver_intestine_9", "wenzel_2024_jejunum_liver_11"):
            entry = mod.CANDIDATES[key]
            assert entry["per_donor_layer_and_outcome"] is None, (
                f"{key}: the contents could not be checked, so the field must "
                "be None rather than True or False"
            )
            assert entry["usable_as"] == "citation_only"

    def test_the_relayed_claim_about_table_s1_is_marked_unverified(self, mod):
        evidence = mod.CANDIDATES["drozdzik_2018_liver_intestine_9"]["evidence"]
        assert "could NOT be checked" in evidence
        assert "not treated as a fact" in evidence

    def test_the_aggregate_only_verdicts_name_a_file(self, mod):
        for key in ("zhang_2015_mppgl_128", "drozdzik_2023_hcv_58"):
            evidence = mod.CANDIDATES[key]["evidence"]
            assert "SHA-256" in evidence, key
            assert "bytes" in evidence, key

    def test_the_commissioning_notes_donor_count_is_corrected(self, mod):
        entry = mod.CANDIDATES["almajdoub_2026_kcat_matched_4"]
        assert entry["n_donors"] == 4
        assert "three to four is the whole cohort" in entry["evidence"]

    def test_pooled_abundance_disqualifies_before_the_donor_count(self, mod):
        """The pairing was never created, which is not a size problem."""
        entry = mod.CANDIDATES["almajdoub_2026_kcat_matched_4"]
        assert entry["per_donor_layer_and_outcome"] is False
        assert entry["same_donor_layers"] is False
        assert entry["reason"] == "ABSTAIN_NO_PER_DONOR_VALUES"
        assert "pooled" in entry["evidence"].lower()
        assert "LOT NUMBER" in entry["evidence"]
        assert "independent of its size" in entry["evidence"]
        assert "ABSTAIN_TOO_FEW_DONORS" in entry["secondary_disqualification"]

    def test_a_per_donor_outcome_with_a_pooled_layer_is_refused_by_the_rule(self, mod):
        entry = dict(
            label="per-donor outcome, pooled predictor",
            reason="ABSTAIN_NO_PER_DONOR_VALUES",
            per_donor_layer_and_outcome=False,
            n_donors=400,
            outcome_kind="drug_specific_metabolic_activity",
        )
        assert mod.qualify(entry) == ("ABSTAIN", "ABSTAIN_NO_PER_DONOR_VALUES")

    def test_the_rule_text_says_pooled_does_not_count(self, mod):
        offline = mod.audit(check_network=False)
        conditions = " ".join(
            offline["rule"]["a_source_is_usable_only_if_all_four_hold"]
        )
        assert "POOLED layer does not satisfy this" in conditions


# --------------------------------------------------------------------------
# 3. the offline audit
# --------------------------------------------------------------------------
class TestOfflineAudit:
    @pytest.fixture(scope="class")
    def report(self, mod):
        return mod.audit(check_network=False)

    def test_it_fits_nothing(self, report):
        assert report["fits_nothing"] is True

    def test_no_candidate_qualifies_as_a_dataset(self, report):
        assert report["summary"]["n_usable_as_datasets"] == 0

    def test_the_reasons_partition_the_candidates(self, report):
        counts = report["summary"]["by_reason"]
        assert sum(counts.values()) == report["summary"]["n_candidates"] == 7
        assert counts["ABSTAIN_ACCESS"] == 2
        assert counts["ABSTAIN_NO_PER_DONOR_VALUES"] == 5
        assert "ABSTAIN_TOO_FEW_DONORS" not in counts, (
            "the donor floor is no longer the FIRST failure of any candidate; "
            "Al-Majdoub fails earlier, on the pooled abundance"
        )

    def test_the_finding_distinguishes_unshared_from_uncollected(self, report):
        finding = report["summary"]["finding"]
        assert "the measurements exist" in finding.lower()
        assert "individual-level sharing" in finding
        assert "pooled commercial lots" in finding, (
            "the one source whose pairing never existed must be separated "
            "from the four whose pairing exists but is withheld"
        )

    def test_it_says_what_it_does_not_claim(self, report):
        text = " ".join(report["summary"]["what_this_does_not_claim"])
        assert "not an exhaustive census" in text
        assert "normal and legitimate choice" in text

    def test_the_offline_audit_contacts_no_network(self, report):
        for row in report["candidates"].values():
            assert "europe_pmc_reverified" not in row


# --------------------------------------------------------------------------
# 4. the committed artefact
# --------------------------------------------------------------------------
class TestCommittedArtefact:
    @pytest.fixture(scope="class")
    def payload(self):
        if not REPORT.is_file():
            pytest.skip("audit artefact not generated")
        return json.loads(REPORT.read_text(encoding="utf-8"))

    def test_the_live_run_reverified_every_access_field(self, payload):
        for key, row in payload["candidates"].items():
            assert "europe_pmc_reverified" in row, (
                f"{key}: the committed artefact must come from a run that "
                "re-queried Europe PMC"
            )
            assert row["europe_pmc_reverified"], key

    def test_the_paywalled_candidates_are_recorded_as_closed(self, payload):
        for key in ("drozdzik_2018_liver_intestine_9", "wenzel_2024_jejunum_liver_11"):
            assert payload["candidates"][key]["europe_pmc_reverified"][
                "isOpenAccess"
            ] == "N"

    def test_the_open_candidates_are_recorded_as_open(self, payload):
        for key in ("almajdoub_2026_kcat_matched_4", "zhang_2015_mppgl_128",
                    "drozdzik_2023_hcv_58", "drozdzik_2021_pathology_97",
                    "song_2018_hlm_16"):
            row = payload["candidates"][key]
            assert row["europe_pmc_reverified"]["isOpenAccess"] == "Y", key
            assert row["europe_pmc_reverified"]["hasSuppl"] == "Y", key

    def test_the_non_commercial_licence_is_recorded(self, payload):
        assert payload["candidates"]["wenzel_2024_jejunum_liver_11"][
            "licence_as_indexed"
        ] == "cc by-nc-nd"

    def test_none_is_acquired(self, payload):
        for key, row in payload["candidates"].items():
            assert row["verdict"] == "ABSTAIN", key

    def test_it_carries_no_numeric_array(self, payload):
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


# --------------------------------------------------------------------------
# 5. the document
# --------------------------------------------------------------------------
class TestDocument:
    @pytest.fixture(scope="class")
    def text(self):
        assert DOCUMENT.is_file(), f"document absent: {DOCUMENT}"
        return _normalise(DOCUMENT.read_text(encoding="utf-8"))

    def test_it_states_the_sharper_gap(self, text):
        assert "The measurements exist in several cohorts" in text
        # The sentence is inside a blockquote, so whitespace normalisation
        # leaves the "> " marker between its two halves.
        assert "individual-level values are almost never published" in text.replace(
            "> ", ""
        )

    def test_it_records_that_step_two_has_nothing_to_run(self, text):
        assert "has **nothing to run**" in text

    def test_it_refuses_to_compute_the_zhang_concordance(self, text):
        assert "concordance in direction" in text
        assert "not a replication" in text
        assert "deliberately **not** turned into a computation" in text

    def test_it_separates_the_never_created_pairing(self, text):
        assert "the donor-level pairing was never created at all" in text
        assert "rather than created and withheld" in text
        assert "With four hundred donors it would still not support" in text

    def test_its_change_log_records_the_corrected_verdict(self, text):
        assert "One verdict corrected by opening the file" in text
        assert "derived from the same wrong belief" in text

    def test_it_corrects_the_three_framings(self, text):
        assert "2018 rather than 2017" in text
        assert "four donors and three after exclusion" in text
        assert "publishes no individual values at all" in text

    def test_it_names_the_data_request_option_without_taking_it(self, text):
        assert "This audit does not make such a request" in text
        assert "without a request or collaboration" in text

    def test_it_states_that_no_earlier_conclusion_moves(self, text):
        assert "No conclusion from any earlier work moves" in text
        assert "Nothing was fitted here" in text

    def test_it_forbids_digitising_figures(self, text):
        assert "digitising scatter figures, which this project forbids" in text

    def test_it_justifies_the_donor_floor_from_the_design(self, text):
        assert "property of the design, not of the biology" in text
