"""The evidence ledger must agree with the artefacts it claims to read.

A consolidating document is the one a reviewer reads first and the one most
likely to drift: it restates numbers that live elsewhere, and nothing warns it
when those numbers move. These tests read the ledger's figures back out of the
committed artefacts, so a drift fails the build instead of reaching a reader.

They also check the two structural properties the ledger's credibility rests
on: that it contains no clinical-pharmacokinetic row, and that its register of
withdrawn claims is present and dated.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "docs" / "OMICS2PK_EVIDENCE_LEDGER_20261110.md"
DATA = REPO / "reports" / "data"


def _load(name):
    path = DATA / name
    if not path.is_file():
        pytest.skip(f"artefact not generated: {name}")
    return json.loads(path.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def text():
    assert LEDGER.is_file(), f"ledger absent: {LEDGER}"
    return re.sub(r"\s+", " ", LEDGER.read_text(encoding="utf-8"))


def _signed(value, places=4):
    """Render a number the way the ledger writes it, with a Unicode minus."""
    return f"{value:+.{places}f}".replace("-", "−")


# --------------------------------------------------------------------------
# 1. the ARIANA figures
# --------------------------------------------------------------------------
class TestArianaFigures:
    @pytest.fixture(scope="class")
    def payload(self):
        return _load("ariana_statistical_versus_mechanistic.json")

    def test_the_metabolomics_contrast_matches(self, payload, text):
        contrast = payload["contrasts"]["A3_metabolomics_vs_clinical"]
        assert f"{contrast['primary_loss']['mean_difference']:.2f}".replace(
            "-", "−"
        ) in text
        low, high = contrast["multiplicity_adjusted"]["interval"]
        assert f"{low:.2f}".replace("-", "−") in text
        assert f"{high:.2f}".replace("-", "−") in text

    def test_the_notes_named_contrast_matches(self, payload, text):
        contrast = payload["contrasts"]["A5_pgx_on_top_of_metabolomics"]
        assert f"{contrast['primary_loss']['mean_difference']:.2f}" in text
        assert "+5.68" in text
        assert contrast["primary_loss"]["direction"] == "degradation"

    def test_the_metabolomics_on_top_of_genotype_contrast_matches(self, payload, text):
        contrast = payload["contrasts"]["A6_metabolomics_on_top_of_genotype"]
        assert f"{contrast['primary_loss']['mean_difference']:.2f}".replace(
            "-", "−"
        ) in text

    def test_the_statistical_versus_mechanistic_genotype_row_matches(self, payload, text):
        contrast = payload["contrasts"]["B3_statistical_vs_mechanistic_genotype_only"]
        assert f"{contrast['primary_loss']['mean_difference']:.2f}" in text
        low, high = contrast["multiplicity_adjusted"]["interval"]
        assert f"{low:.2f}" in text and f"{high:.2f}" in text
        assert contrast["multiplicity_adjusted"]["interval_excludes_zero"] is True

    def test_the_matched_panel_row_is_reported_as_undetected(self, payload, text):
        contrast = payload["contrasts"]["B1_statistical_vs_mechanistic_omic_matched"]
        assert contrast["multiplicity_adjusted"]["interval_excludes_zero"] is False
        assert "nothing detected" in text

    def test_the_three_carried_claims_are_the_artefacts_three(self, payload):
        assert len(payload["claims_carried"]) == 3


# --------------------------------------------------------------------------
# 2. the Axis B figures
# --------------------------------------------------------------------------
class TestAxisBFigures:
    @pytest.fixture(scope="class")
    def payload(self):
        return _load("axis_b_wegler_protein_clint_pilot.json")

    def test_the_two_positive_q1_rows_match(self, payload, text):
        for substrate in ("Midazolam", "Omeprazole"):
            contrasts = payload["results"][f"HH::{substrate}::k2::all"]["contrasts"]
            contrast = contrasts["secondary_W1_vs_W0"]
            assert _signed(contrast["mean_difference_log10"]) in text, substrate
            adjusted = contrasts["secondary_W1_vs_W0_bonferroni"]
            assert adjusted["interval_excludes_zero"] is True, substrate
            low, high = adjusted["interval"]
            assert _signed(low) in text, f"{substrate} low"
            assert _signed(high) in text, f"{substrate} high"

    def test_the_microsomal_null_row_matches(self, payload, text):
        contrasts = payload["results"]["HLM::Midazolam::k2::all"]["contrasts"]
        assert _signed(contrasts["secondary_W1_vs_W0"]["mean_difference_log10"]) in text
        adjusted = contrasts["secondary_W1_vs_W0_bonferroni"]
        assert adjusted["interval_excludes_zero"] is False
        assert adjusted["family_size"] == 5

    def test_the_one_significant_degradation_matches(self, payload, text):
        contrasts = payload["results"]["HH::Diclofenac::k2::all"]["contrasts"]
        contrast = contrasts["primary_W2_vs_W1"]
        assert _signed(contrast["mean_difference_log10"]) in text
        assert contrast["direction"] == "degradation"
        adjusted = contrasts["primary_W2_vs_W1_bonferroni"]
        assert adjusted["interval_excludes_zero"] is True
        low, high = adjusted["interval"]
        assert f"{low:+.4f}" in text and f"{high:+.4f}" in text

    def test_no_q2_cell_improves_with_an_interval_excluding_zero(self, payload):
        for key, cell in payload["results"].items():
            contrast = cell["contrasts"]["primary_W2_vs_W1"]
            if contrast["interval_excludes_zero"]:
                assert contrast["direction"] == "degradation", key

    def test_the_paired_system_rows_match(self, payload, text):
        for substrate in ("Omeprazole", "Midazolam", "Diclofenac"):
            entry = payload["paired_system_comparison"][substrate]["W1_vs_W0"]
            assert entry["interval_excludes_zero"] is True, substrate
            assert _signed(entry["mean_difference_log10"]) in text, substrate

    def test_exactly_three_substrates_differ_between_systems(self, payload):
        detected = [
            s for s, e in payload["paired_system_comparison"].items()
            if e["W1_vs_W0"]["interval_excludes_zero"]
        ]
        assert len(detected) == 3


# --------------------------------------------------------------------------
# 3. the Axis A figures
# --------------------------------------------------------------------------
class TestAxisAFigures:
    @pytest.fixture(scope="class")
    def payload(self):
        return _load("axis_a_eatris_masked_modality_reconstruction.json")

    def test_every_cell_is_not_supported(self, payload, text):
        for key, cell in payload["results"].items():
            assert cell["prespecified_decision"]["verdict"] == "NOT_SUPPORTED", key
        assert "`NOT_SUPPORTED` in all four cells" in text

    def test_the_contrast_count_matches(self, payload, text):
        total = sum(len(cell["contrasts"]) for cell in payload["results"].values())
        assert f"**{total}** contrasts" in text, f"the ledger must say {total}"

    def test_exactly_one_contrast_survives_both_controls(self, payload, text):
        survivors = [
            (key, name)
            for key, cell in payload["results"].items()
            for name, contrast in cell["contrasts"].items()
            if contrast.get("interval_excludes_zero")
            and contrast.get("permutation_null", {}).get("better_than_every_replicate")
            is None
            and contrast.get("permutation_null", {}).get("exceeds_every_replicate")
        ]
        assert len(survivors) == 1
        key, name = survivors[0]
        assert "acylcarnitines" in key
        assert name == "all_but_target_A3_vs_A1"
        value = payload["results"][key]["contrasts"][name]["median_delta_r2"]
        assert _signed(value) in text

    def test_the_methylome_artefact_row_matches(self, payload, text):
        contrast = payload["results"]["C11_primary::proteins"]["contrasts"][
            "single_family_B_methylome_vs_A0"
        ]
        assert contrast["interval_excludes_zero"] is True
        assert contrast["permutation_null"]["null_min"] > 0.0, (
            "the ledger's claim is that the whole null is positive"
        )
        assert _signed(contrast["median_delta_r2"]) in text

    def test_the_null_is_positive_in_every_cell(self, payload):
        for key, cell in payload["results"].items():
            median = cell["permutation_control"]["null_median"]
            assert median > 0.0, f"{key}: null median {median}"


# --------------------------------------------------------------------------
# 4. the structural counts
# --------------------------------------------------------------------------
class TestStructuralCounts:
    def test_the_readiness_gate_count_matches(self, text):
        payload = _load("hepatic_multiview_readiness_verdict.json")
        blob = json.dumps(payload)
        assert "fit_allowed` on 0 of 10" in text or "0 of 10" in text
        assert "ABSTAIN" in blob

    def test_the_candidate_audit_count_matches(self, text):
        payload = _load("candidate_paired_liver_sources.json")
        assert payload["summary"]["n_usable_as_datasets"] == 0
        assert payload["summary"]["n_candidates"] == 7
        assert "**0 of 7** usable as datasets" in text
        assert "**2 usable datasets**" in text, (
            "the repository search found two datasets the publisher-first "
            "search missed; the ledger must say so beside the 0 of 7")
        counts = payload["summary"]["by_reason"]
        assert f"{counts['ABSTAIN_ACCESS']} on access" in text
        assert f"{counts['ABSTAIN_NO_PER_DONOR_VALUES']} on the absence" in text

    def test_the_wegler_verification_count_matches(self, text):
        payload = _load("wegler_table_s3_per_donor_kinetics_audit.json")
        verification = payload["verification"]
        assert (
            f"**{verification['n_anchors_agreeing']} of "
            f"{verification['n_anchors']}**" in text
        )
        median = verification["median_midazolam_clint_mic_eq7"]
        assert f"**{median['recomputed']}**" in text
        assert f"**{int(median['published_in_results'])}**" in text

    def test_the_replication_attempt_figures_match(self, text):
        payload = _load("wegler_protein_clint_replication.json")
        summary = payload["summary"]
        assert f"**{summary['microsomal_median_spearman']:.3f}**" in text
        inside = summary[
            "pairs_where_the_hlc_transcript_falls_inside_this_cohorts_interval"
        ]
        assert len(inside) == 5
        assert "**5 of 5**" in text


# --------------------------------------------------------------------------
# 5. the properties the ledger's credibility rests on
# --------------------------------------------------------------------------
class TestLedgerProperties:
    def test_the_five_levels_are_defined_before_the_table(self, text):
        definitions = text.index("The five levels, defined before the table")
        table = text.index("## 1. The ledger")
        assert definitions < table

    def test_it_holds_no_clinical_pk_result(self, text):
        assert "The project holds no C-level result" in text
        assert "No clinical-PK proof" in text

    def test_it_holds_exactly_one_positive_independent_validation(self, text):
        assert "only **positive** V-level row" in text

    def test_the_failed_replication_is_recorded_at_v_level(self, text):
        # An attempted independent validation that did not confirm belongs in
        # the ledger at the level of the attempt, not quietly in a footnote.
        assert "attempted and not confirmed" in text
        assert "0 of 6** contrasts SUPPORTED" in text

    def test_a_correlation_is_not_called_a_performance(self, text):
        assert "a published" in text and "correlation is not a" in text

    def test_the_withdrawal_register_is_present_and_counted(self, text):
        assert "Ten claims have been withdrawn" in text
        section = text.split("register of withdrawn")[1].split("## 3.")[0]
        # Each row of the register ends with a date.
        dates = re.findall(r"2026-11-\d\d", section)
        assert len(dates) >= 8, f"found {len(dates)} dated rows"

    def test_the_withdrawn_variance_fraction_is_in_the_register(self, text):
        assert "squared a **Spearman** coefficient" in text

    def test_the_withdrawn_protein_versus_mrna_headline_is_in_the_register(self, text):
        assert "Interpretation withdrawn" in text
        assert "0.944" in text and "0.412" in text

    def test_it_states_the_interval_caveat_for_the_positive_rows(self, text):
        assert "resamples stored out-of-fold errors without refitting" in text
        assert "narrower than a fully nested interval" in text

    def test_it_confines_the_sample_size_figure(self, text):
        assert "This sizes **that contrast only**" in text

    def test_the_unsizeable_requirements_are_labelled_as_scenarios(self, text):
        # They were "not sized" until the sizing document existed. They are
        # now sized only under declared assumptions, which must never be
        # allowed to read as a measured requirement.
        assert text.count("sized **only as declared scenarios**") == 2
        assert "which are not interchangeable" in text

    def test_it_names_the_data_request_option_without_taking_it(self, text):
        assert "data-sharing request** to the corresponding authors" in text
        assert "leaves the decision where it belongs" in text
        assert "deferred it" in text

    def test_it_distinguishes_itself_from_the_older_matrix(self, text):
        assert "**not** the October 2026" in text
        assert "EVIDENCE_MATRIX.md" in text

    def test_it_recomputes_nothing(self, text):
        assert "Nothing new was computed for this document" in text
        assert "No figure was recomputed for this document" in text

    def test_every_level_appears_in_the_table(self, text):
        for heading in ("Reproduction (R)", "Independent validation (V)",
                        "Functional proof (F)", "Exploratory (E)"):
            assert heading in text
