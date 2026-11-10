"""Contract tests for the Wegler per-donor clearance recovery.

The danger this module exists to prevent is specific and it very nearly
happened: ``pdftotext -layout`` reads Table S3 with the clearance cells
misaligned against the batch labels, so a text-order transcription assigns one
donor's intrinsic clearance to another, silently, and produces a file that
looks entirely plausible. Everything downstream — a protein-to-clearance
model, a replication attempt, a funding claim — would inherit that.

So these tests check, in order of consequence:

1. the geometry definitions cannot silently merge or skip a column;
2. the audit refuses to emit a table that does not reproduce the twenty
   correlation coefficients the article itself published;
3. the per-donor rows cannot be written inside this repository;
4. the committed artefacts say what the documents say they say.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
AUDIT_SCRIPT = REPO / "scripts" / "audit_wegler_table_s3_per_donor_kinetics.py"
COMPARE_SCRIPT = (
    REPO / "scripts" / "compare_wegler_protein_clint_against_published_and_transcript.py"
)
AUDIT_JSON = REPO / "reports" / "data" / "wegler_table_s3_per_donor_kinetics_audit.json"
REPLICATION_JSON = REPO / "reports" / "data" / "wegler_protein_clint_replication.json"
DOCUMENT = (
    REPO / "docs" / "OMICS2PK_WEGLER_PER_DONOR_CLINT_QUALIFICATION_20261110.md"
)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def audit():
    sys.path.insert(0, str(REPO / "src"))
    return _load("wegler_audit", AUDIT_SCRIPT)


@pytest.fixture(scope="module")
def compare():
    sys.path.insert(0, str(REPO / "src"))
    return _load("wegler_compare", COMPARE_SCRIPT)


def _normalise(text):
    return re.sub(r"\s+", " ", text)


# --------------------------------------------------------------------------
# 1. the geometry cannot silently merge or skip a column
# --------------------------------------------------------------------------
class TestGeometryDefinitions:
    def test_there_is_one_band_per_column(self, audit):
        assert len(audit.COLUMN_BANDS) == len(audit.COLUMNS) == 11

    def test_the_bands_are_ordered_and_do_not_overlap(self, audit):
        for (low, high), (next_low, _) in zip(
            audit.COLUMN_BANDS, audit.COLUMN_BANDS[1:]
        ):
            assert low < high
            assert high <= next_low, "two bands overlap, so a column could merge"

    def test_the_bands_leave_no_gap(self, audit):
        for (_, high), (next_low, _) in zip(
            audit.COLUMN_BANDS, audit.COLUMN_BANDS[1:]
        ):
            assert high == next_low, "a gap would silently drop a value"

    def test_there_are_fifteen_donors_and_five_substrates(self, audit):
        assert len(audit.DONORS) == 15
        assert audit.DONORS[0] == "UU1" and audit.DONORS[-1] == "UU15"
        assert len(audit.SUBSTRATES) == 5

    def test_every_substrate_declares_its_cyp(self, audit):
        assert set(audit.SUBSTRATE_CYP) == set(audit.SUBSTRATES)
        assert len(set(audit.SUBSTRATE_CYP.values())) == 5

    def test_the_clearance_columns_are_named_for_their_system(self, audit):
        microsomal = [c for c in audit.COLUMNS if c.startswith("clint_mic")]
        hepatic = [c for c in audit.COLUMNS if c.startswith("clint_hep")]
        assert len(microsomal) == 3, "three microsomal equations"
        assert len(hepatic) == 2, "two hepatocyte equations"

    def test_not_determined_is_a_value_and_not_a_zero(self, audit):
        assert audit._value("ND") is None
        assert audit._value("nd") is None
        assert audit._value("0.0") == 0.0

    def test_an_unparseable_cell_is_refused(self, audit):
        with pytest.raises(audit.AuditError):
            audit._value("~12")


# --------------------------------------------------------------------------
# 2. the statistics owe nothing to a library, so they are checked
# --------------------------------------------------------------------------
class TestCorrelationPrimitives:
    def test_pearson_of_a_perfect_line_is_one(self, audit):
        assert audit._pearson([1, 2, 3, 4], [2, 4, 6, 8]) == pytest.approx(1.0)

    def test_pearson_of_a_perfect_inverse_is_minus_one(self, audit):
        assert audit._pearson([1, 2, 3, 4], [8, 6, 4, 2]) == pytest.approx(-1.0)

    def test_a_constant_series_has_no_correlation(self, audit):
        assert audit._pearson([1, 1, 1, 1], [1, 2, 3, 4]) is None

    def test_spearman_is_one_for_any_monotone_map(self, audit):
        assert audit._spearman([1, 2, 3, 4], [1, 10, 1000, 100000]) == pytest.approx(1.0)

    def test_ties_take_average_ranks(self, audit):
        assert audit._ranks([5, 5, 1]) == [2.5, 2.5, 1.0]

    def test_it_agrees_with_scipy_where_scipy_exists(self, audit):
        stats = pytest.importorskip("scipy.stats")
        xs = [3.0, 1.0, 4.0, 1.5, 5.0, 9.0, 2.0, 6.0]
        ys = [2.0, 7.0, 1.0, 8.0, 2.8, 1.8, 2.8, 4.0]
        assert audit._pearson(xs, ys) == pytest.approx(
            float(stats.pearsonr(xs, ys)[0]), abs=1e-10
        )
        assert audit._spearman(xs, ys) == pytest.approx(
            float(stats.spearmanr(xs, ys)[0]), abs=1e-10
        )

    def test_the_comparison_script_shares_the_same_spearman(self, audit, compare):
        xs = [3.0, 1.0, 4.0, 1.5, 5.0, 9.0, 2.0]
        ys = [2.0, 7.0, 1.0, 8.0, 2.8, 1.8, 2.8]
        assert compare.spearman(xs, ys) == pytest.approx(audit._spearman(xs, ys))


# --------------------------------------------------------------------------
# 3. the refusals
# --------------------------------------------------------------------------
class TestRefusals:
    def test_a_wrong_pdf_checksum_is_refused(self, audit, tmp_path):
        decoy = tmp_path / "not_the_supplement.pdf"
        decoy.write_bytes(b"%PDF-1.4 not the real file")
        with pytest.raises(audit.AuditError, match="checksum mismatch"):
            audit.audit(str(decoy), None)

    def test_writing_the_per_donor_table_into_the_repository_is_refused(self, audit):
        rows = {
            (substrate, donor): dict.fromkeys(audit.COLUMNS, 1.0)
            for substrate in audit.SUBSTRATES
            for donor in audit.DONORS
        }
        with pytest.raises(audit.AuditError, match="must not be written inside"):
            audit.write_local_csv(rows, str(REPO / "reports" / "data"))

    def test_writing_it_outside_the_repository_is_allowed(self, audit, tmp_path):
        rows = {
            (substrate, donor): dict.fromkeys(audit.COLUMNS, 1.0)
            for substrate in audit.SUBSTRATES
            for donor in audit.DONORS
        }
        target = audit.write_local_csv(rows, str(tmp_path / "local"))
        assert target.is_file()
        assert REPO not in target.parents

    def test_all_twenty_anchors_are_required(self, audit):
        assert audit.REQUIRED_ANCHORS == len(audit.PUBLISHED_CORRELATIONS) == 20

    def test_a_mis_assigned_donor_would_fail_the_anchors(self, audit):
        """The point of the verification: rotating one column must break it."""
        rows = _rows_from_committed_csv(audit)
        if rows is None:
            pytest.skip("per-donor csv is local and not present here")
        verification = audit.verify_against_the_article(rows)
        assert verification["all_anchors_agree"] is True

        rotated = dict(rows)
        order = list(audit.DONORS)
        shifted = order[1:] + order[:1]
        for original, replacement in zip(order, shifted):
            row = dict(rows[("Midazolam", original)])
            row["clint_mic_eq7"] = rows[("Midazolam", replacement)]["clint_mic_eq7"]
            rotated[("Midazolam", original)] = row
        broken = audit.verify_against_the_article(rotated)
        assert broken["all_anchors_agree"] is False, (
            "shifting one donor's clearance by one row must be detected, or "
            "the verification is not doing its job"
        )

    def test_the_comparison_refuses_an_unverified_table(self, compare, tmp_path):
        payload = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
        payload["verification"]["all_anchors_agree"] = False
        path = tmp_path / "audit.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(compare.ComparisonError, match="not been verified"):
            compare.run("does-not-matter.csv", str(path))

    def test_the_comparison_refuses_a_different_table(self, compare, tmp_path):
        payload = json.loads(AUDIT_JSON.read_text(encoding="utf-8"))
        payload["per_donor_table_sha256"] = "0" * 64
        path = tmp_path / "audit.json"
        path.write_text(json.dumps(payload), encoding="utf-8")
        with pytest.raises(compare.ComparisonError, match="checksum does not match"):
            compare.run("does-not-matter.csv", str(path))


def _rows_from_committed_csv(audit):
    """Read the per-donor table from the local directory, if it is there."""
    import csv

    candidate = (
        REPO.parent
        / "Wegler2021_PMC8041379_supplement"
        / "extracted"
        / "wegler2021_table_s3_per_donor_kinetics.csv"
    )
    if not candidate.is_file():
        return None
    rows = {}
    with candidate.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows[(row["substrate"], row["donor"])] = {
                name: (None if row[name] == "ND" else float(row[name]))
                for name in audit.COLUMNS
            }
    return rows


# --------------------------------------------------------------------------
# 4. the committed audit artefact
# --------------------------------------------------------------------------
class TestCommittedAudit:
    @pytest.fixture(scope="class")
    def payload(self):
        if not AUDIT_JSON.is_file():
            pytest.skip("audit artefact not generated")
        return json.loads(AUDIT_JSON.read_text(encoding="utf-8"))

    def test_every_anchor_agrees(self, payload):
        verification = payload["verification"]
        assert verification["n_anchors"] == 20
        assert verification["n_anchors_agreeing"] == 20
        assert verification["all_anchors_agree"] is True

    def test_the_published_median_is_reproduced(self, payload):
        median = payload["verification"]["median_midazolam_clint_mic_eq7"]
        assert median["agrees"] is True
        assert median["published_in_results"] == 375.0
        assert median["recomputed"] == pytest.approx(375.0, abs=1.0)
        assert median["n_donors"] == 15

    def test_all_eight_hundred_and_twenty_five_values_are_recovered(self, payload):
        assert payload["n_donors"] == 15
        assert payload["n_substrates"] == 5
        assert len(payload["columns"]) == 11
        assert payload["n_values_recovered"] == 825

    def test_the_table_spans_three_pages(self, payload):
        assert payload["table_s3_pages"] == [5, 6, 7]

    def test_the_per_donor_rows_are_not_in_the_repository(self, payload):
        assert payload["per_donor_rows_are_not_in_this_repository"] is True
        assert re.fullmatch(r"[0-9a-f]{64}", payload["per_donor_table_sha256"])

    def test_the_licence_caveat_carries_the_projects_own_formula(self, payload):
        assert "silence about terms is not permission" in payload["licence_caveat"].lower()

    def test_the_absence_of_per_donor_demographics_is_stated(self, payload):
        text = payload["what_table_s1_does_not_give"]
        assert "aggregates only" in text
        assert "no per-donor" in text.lower()

    def test_it_warns_against_a_text_dump(self, payload):
        assert "mis-assigned" in payload["why_a_text_dump_will_not_do"]

    def test_no_aggregate_hides_a_per_donor_vector(self, payload):
        """Aggregates only: the audit must publish no length-15 vector."""
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

    def test_every_substrate_reports_its_clearance_aggregates(self, payload):
        for substrate, entry in payload["aggregates_only"].items():
            for name in ("clint_mic_eq7", "clint_hep_eq7"):
                column = entry["columns"][name]
                assert column["n"] + column["n_not_determined"] == 15, (
                    f"{substrate}/{name}"
                )

    def test_bufuralol_microsomal_clearance_is_missing_for_two_donors(self, payload):
        column = payload["aggregates_only"]["Bufuralol"]["columns"]["clint_mic_eq7"]
        assert column["n"] == 13
        assert column["n_not_determined"] == 2


# --------------------------------------------------------------------------
# 5. the committed replication artefact
# --------------------------------------------------------------------------
class TestCommittedReplication:
    @pytest.fixture(scope="class")
    def payload(self):
        if not REPLICATION_JSON.is_file():
            pytest.skip("replication artefact not generated")
        return json.loads(REPLICATION_JSON.read_text(encoding="utf-8"))

    def test_it_is_labelled_as_a_replication_attempt(self, payload):
        assert payload["label"] == "INDEPENDENT_REPLICATION_ATTEMPT"

    def test_it_covers_the_five_matched_pairs(self, payload):
        assert payload["summary"]["n_pairs"] == 5
        assert len(payload["per_pair"]) == 5

    def test_every_pair_reports_both_systems_with_an_interval(self, payload):
        for substrate, entry in payload["per_pair"].items():
            for system in ("HLM", "HH"):
                cell = entry["systems"][system]
                assert cell["bootstrap_95_interval"] is not None, f"{substrate}/{system}"
                low, high = cell["bootstrap_95_interval"]
                assert low <= cell["spearman_rho"] <= high, f"{substrate}/{system}"

    def test_the_microsomal_median_is_far_below_the_published_figure(self, payload):
        summary = payload["summary"]
        assert summary["microsomal_median_spearman"] < 0.6
        assert summary["grangeon_published_microsomal_median"] > 0.9

    def test_the_transcript_falls_inside_this_cohorts_interval_everywhere(self, payload):
        inside = payload["summary"][
            "pairs_where_the_hlc_transcript_falls_inside_this_cohorts_interval"
        ]
        assert len(inside) == 5, (
            "the claim in the write-up is 5 of 5; if that changes the document "
            "must change with it"
        )

    def test_the_published_figure_falls_outside_on_some_enzymes(self, payload):
        outside = payload["summary"]["pairs_where_grangeon_falls_outside"]
        assert outside, "at least one enzyme where the published value is excluded"

    def test_it_says_what_it_is_not(self, payload):
        text = " ".join(payload["what_this_is_not"])
        assert "Not a like-for-like replication" in text
        assert "Not a correction of Grangeon's arithmetic" in text
        assert "clinical pharmacokinetics" in text

    def test_it_declares_the_endpoint_difference(self, payload):
        assert "CLEARANCE" in payload["endpoint_difference"]
        assert "ACTIVITY" in payload["endpoint_difference"]

    def test_it_explains_how_a_wide_interval_should_be_read(self, payload):
        assert "rarely EXCLUDES" in payload["how_to_read"]

    def test_it_exports_no_per_donor_data(self, payload):
        assert payload["no_person_identifiers_values_or_predictions_exported"] is True


# --------------------------------------------------------------------------
# 6. the document
# --------------------------------------------------------------------------
class TestDocument:
    @pytest.fixture(scope="class")
    def text(self):
        assert DOCUMENT.is_file(), f"document absent: {DOCUMENT}"
        return _normalise(DOCUMENT.read_text(encoding="utf-8"))

    def test_it_answers_the_question_that_was_asked(self, text):
        assert "They are accessible and they are matchable" in text

    def test_it_states_the_count_of_recovered_values(self, text):
        assert "825 per-donor values recovered" in text

    def test_it_records_that_a_text_dump_mis_assigns_donors(self, text):
        assert "assigns one donor's clearance to another" in text

    def test_it_publishes_the_twenty_anchors(self, text):
        assert "All twenty reproduce" in text
        for published in ("0.9396", "0.3251", "0.1107", "0.4121"):
            assert published in text

    def test_it_states_the_absence_of_per_donor_demographics(self, text):
        assert "No per-donor age, sex, BMI or diagnosis" in text

    def test_it_corrects_the_previous_headline(self, text):
        assert "This is a correction" in text
        assert "0.412" in text and "0.944" in text
        assert "5 of 5" in text

    def test_it_says_what_is_not_being_claimed(self, text):
        assert "Not that Grangeon is wrong" in text
        assert "Not a like-for-like replication" in text

    def test_it_keeps_the_per_donor_rows_out_of_the_repository(self, text):
        assert "silence about terms is not permission" in text.lower()

    def test_it_refuses_a_clinical_pk_reading(self, text):
        assert "never** a clinical-PK proof" in text or (
            "never a clinical-PK proof" in text
        )

    def test_it_states_the_axis_b_design_without_claiming_a_result(self, text):
        assert "the honest prior is that nothing will be detected" in text

    def test_it_forbids_fusion(self, text):
        assert "It must not be fused with anything" in text

    def test_the_superseding_note_is_in_the_older_document(self):
        older = (
            REPO / "docs" / "OMICS2PK_PUBLIC_SOURCE_QUALIFICATION_20261109.md"
        ).read_text(encoding="utf-8")
        assert "SUPERSEDED IN PART, 2026-11-10" in older
        assert "must not be quoted on its own" in _normalise(older)
