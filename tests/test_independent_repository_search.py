"""Tests for the repository search.

These do not re-run the live queries; the script does that and fails on drift.
What they check is that the record is internally consistent, that the two
verdicts which turn on a distinction -- access against aggregates -- are not
confused, and that the document states the conclusion it withdraws.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import re
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_independent_repository_search.py"
ARTEFACT = ROOT / "reports" / "data" / "independent_repository_search.json"
DOC = ROOT / "docs" / "OMICS2PK_INDEPENDENT_REPOSITORY_SEARCH_20261110.md"
EARLIER = (ROOT / "docs"
           / "OMICS2PK_CANDIDATE_PAIRED_LIVER_SOURCE_QUALIFICATION_20261110.md")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("repo_search", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["repo_search"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def payload():
    if not ARTEFACT.is_file():
        pytest.skip("the audit has not been run in this checkout")
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


def _normalise(text):
    return re.sub(r"[\s*_`>|]+", " ", text)


# --------------------------------------------------------------------------
# 1. the search is recorded, not asserted
# --------------------------------------------------------------------------
class TestSearchRecord:
    def test_eight_repositories_are_named(self, mod):
        assert len(mod.REPOSITORIES) == 8

    def test_every_repository_records_its_api_and_its_queries(self, mod):
        for name, record in mod.REPOSITORIES.items():
            assert record["api"].startswith("http"), name
            assert record["queries"], f"{name} records no query string"
            assert len(record["outcome"]) > 20, name

    def test_the_productive_index_is_named_as_such(self, mod):
        assert "productive" in \
            mod.REPOSITORIES["BioStudies_ArrayExpress"]["outcome"]


# --------------------------------------------------------------------------
# 2. verdicts are declared and used consistently
# --------------------------------------------------------------------------
class TestVerdicts:
    def test_every_recorded_verdict_is_declared(self, mod):
        for record in list(mod.PRIDE_FINDINGS.values()) + \
                list(mod.ARTICLE_FINDINGS.values()):
            assert record["verdict"] in mod.VERDICTS

    def test_the_wegler_deposition_is_not_counted_as_a_new_dataset(self, mod):
        assert mod.PRIDE_FINDINGS["PXD014131"]["verdict"] == "ALREADY_HELD"

    def test_the_venetoclax_paper_is_not_counted_as_a_second_cohort(self, mod):
        record = mod.ARTICLE_FINDINGS["PMC11645475"]
        assert record["verdict"] == "ALREADY_HELD"
        assert "ONE cohort with two substrates" in record["note"]

    def test_the_proteomics_depositions_are_layer_without_outcome(self, mod):
        for accession in ("PXD020910", "PXD021025", "PXD038776"):
            assert mod.PRIDE_FINDINGS[accession]["verdict"] == \
                "LAYER_WITHOUT_OUTCOME"

    def test_the_171_donor_study_is_not_recorded_as_an_access_failure(self, mod):
        record = mod.ARTICLE_FINDINGS["PMC6222000"]
        assert record["verdict"] == "ABSTAIN_NO_PER_DONOR_VALUES", (
            "the file was obtained; recording it as an access failure would "
            "aim the funding argument at the wrong barrier")
        assert "availability TICKS, not values" in record["note"]
        assert "barrier was never access" in record["note"]

    def test_exactly_two_datasets_are_called_usable(self, mod):
        usable = [k for k, v in mod.ARTICLE_FINDINGS.items()
                  if v["verdict"] == "DATASET_USABLE"]
        assert sorted(usable) == ["PMC9310776", "PMC9926076"]

    def test_the_two_paywalled_sources_are_still_abstain_access(self, mod):
        assert set(mod.STILL_INACCESSIBLE) == {"29205295", "37739780"}


# --------------------------------------------------------------------------
# 3. the access ladder, and the limit on what it proved
# --------------------------------------------------------------------------
class TestAccessLadder:
    def test_the_lesson_records_both_halves(self, payload):
        lesson = payload["the_access_lesson"]
        assert "do not establish" in lesson["finding"]
        assert "did not change the verdict" in lesson["but"]

    def test_the_proof_of_work_challenge_was_not_defeated(self, payload):
        text = payload["the_access_lesson"]["what_was_not_done"]
        assert "not defeated" in text
        assert "requires registration" in text

    def test_the_biostudies_file_sizes_are_recorded(self, payload):
        listings = payload["biostudies_file_listings"]
        assert listings["S-EPMC6222000"]["size_bytes"] == 2547704
        assert listings["S-EPMC9310776"]["size_bytes"] == 11163005

    def test_each_inaccessible_source_records_the_channels_tried(self, payload):
        for pmid, record in payload[
                "still_inaccessible_after_this_search"].items():
            assert len(record["channels_tried"]) >= 4
            assert record["pmcid"] in (None, "")
            assert "PMCID appearing" in record["what_would_change_it"]


# --------------------------------------------------------------------------
# 4. the withdrawn conclusion
# --------------------------------------------------------------------------
class TestWithdrawal:
    def test_the_artefact_records_the_withdrawal(self, payload):
        block = payload["corrects_an_earlier_claim"]
        assert block["status"] == "WITHDRAWN"
        assert "nothing to run" in block["claim"]
        assert "too narrow" in block["why"]

    def test_the_document_states_the_withdrawal_in_its_first_section(self):
        text = _normalise(DOC.read_text(encoding="utf-8"))
        head = text[:2600]
        assert "That conclusion is withdrawn" in head

    def test_the_earlier_document_still_exists_and_is_not_rewritten(self):
        assert EARLIER.is_file(), (
            "the superseded qualification must stay in the record; a withdrawn "
            "claim is evidence about the project, not something to delete")
        assert "nothing to run" in EARLIER.read_text(encoding="utf-8")


# --------------------------------------------------------------------------
# 5. the pattern the funding case argues from
# --------------------------------------------------------------------------
class TestPattern:
    def test_the_counts_are_recorded_and_consistent(self, payload):
        block = payload["the_pattern_now_stands_at"]
        measured = block["sources_that_measured_a_per_donor_pairing"]
        published = block["sources_that_published_the_per_donor_rows"]
        assert measured > published
        assert published == 2
        assert "Wegler" in block["named"] and "Lee" in block["named"]

    def test_the_document_carries_the_same_counts(self, payload):
        text = DOC.read_text(encoding="utf-8")
        block = payload["the_pattern_now_stands_at"]
        assert f"**{block['sources_that_measured_a_per_donor_pairing']}**" in text
        assert f"**{block['sources_that_published_the_per_donor_rows']}**" in text

    def test_the_document_says_the_rows_were_published_by_accident_of_format(self):
        text = _normalise(DOC.read_text(encoding="utf-8"))
        assert "accident of format rather than by policy" in text


# --------------------------------------------------------------------------
# 6. the unresolved items are listed rather than left implied
# --------------------------------------------------------------------------
class TestUnresolved:
    def test_the_clinical_pk_contrast_is_explicitly_not_unblocked(self):
        text = _normalise(DOC.read_text(encoding="utf-8"))
        assert "A clinical-PK contrast is not unblocked" in text
        assert "does not digitise scatter figures" in text

    def test_the_second_dataset_is_recorded_as_having_no_endpoint(self, mod):
        note = mod.ARTICLE_FINDINGS["PMC9310776"]["note"]
        assert "NO pharmacological endpoint" in note

    def test_the_document_declines_to_execute_the_second_dataset(self):
        text = _normalise(DOC.read_text(encoding="utf-8"))
        assert "is not executed here" in text

    def test_the_artefact_contains_no_individual_level_data(self, payload):
        assert payload["no_individual_level_data_in_this_artefact"] is True
        assert payload["this_audit_fits_nothing"] is True
