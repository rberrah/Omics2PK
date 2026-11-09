"""Tests for the 2026-11-09 qualification of three candidate public sources.

These lock the three things a source audit can get wrong in a way that would
mislead a funding case:

* claiming individual-level data that is not published (Grangeon);
* understating a source that is in fact usable (Wegler);
* letting a cohort with no pharmacological endpoint drift into a
  pharmacological claim (EATRIS-Plus).

They also lock the cross-cohort label and the three competing explanations on
the mRNA-versus-protein comparison, because that number is the one most likely
to be quoted without its caveats.
"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

REGISTRY = REPO_ROOT / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"
DOC = REPO_ROOT / "docs" / "OMICS2PK_PUBLIC_SOURCE_QUALIFICATION_20261109.md"
ARTEFACT = REPO_ROOT / "reports" / "data" / (
    "rna_versus_protein_activity_correlation.json")
SCRIPT = REPO_ROOT / "scripts" / (
    "compare_rna_versus_published_protein_activity_correlation.py")


def registry():
    return json.loads(REGISTRY.read_text(encoding="utf-8"))


def dataset(key):
    for entry in registry()["datasets"]:
        if entry["key"] == key:
            return entry
    raise AssertionError(f"no registry entry named {key!r}")


class TestEatrisPlus(unittest.TestCase):
    """A blood-only cohort with no drug must not be able to drift into Axis B or C."""

    def setUp(self):
        self.entry = dataset("eatris_plus_v4")

    def test_it_supports_no_objective(self):
        self.assertEqual(self.entry["objectives_supported"], [])

    def test_it_has_no_endpoint_at_all(self):
        self.assertIsNone(self.entry["endpoint"])
        self.assertIsNone(self.entry["dose_regimen"])

    def test_the_absence_of_a_pharmacological_endpoint_is_stated(self):
        caveats = " ".join(self.entry["caveats"]).lower()
        self.assertIn("no pharmacological endpoint", caveats)
        self.assertIn("no drug", caveats)

    def test_the_compartment_limit_is_stated(self):
        caveats = " ".join(self.entry["caveats"])
        self.assertIn("do not measure hepatic drug-metabolising activity",
                      caveats)

    def test_the_twelve_modalities_are_recorded_with_their_own_n(self):
        layers = self.entry["omics_layers"]
        self.assertEqual(len(layers), 12)
        counts = {layer["layer"]: layer["individuals"] for layer in layers}
        self.assertEqual(counts["transcriptomics_mRNA"], 90)
        self.assertEqual(counts["miRNA_qRT_PCR"], 86)
        self.assertEqual(counts["methylomics_EM_seq"], 117)
        self.assertEqual(counts["proteomics_targeted"], 125)

    def test_the_complete_case_curve_is_recorded(self):
        curve = self.entry["complete_case_by_layer_count"]
        self.assertEqual(curve["all_12_layers"], 55)
        self.assertEqual(curve["11_layers_without_miRNA_qRT_PCR"], 79)
        self.assertEqual(curve["9_layers_without_miRNA_and_without_mRNA"], 114)

    def test_the_metadata_discrepancy_is_registered_not_silently_resolved(self):
        caveats = " ".join(self.entry["caveats"])
        self.assertIn("ISA workbook claims", caveats)
        self.assertIn("over-count", caveats)

    def test_the_union_is_the_analysed_count_not_the_isa_count(self):
        self.assertEqual(self.entry["n_enrolled"], 127)
        self.assertEqual(self.entry["n_analysed"], 125)

    def test_the_downloaded_container_carries_a_checksum(self):
        artifacts = self.entry["local_artifacts"]
        self.assertTrue(artifacts)
        h5mu = next(a for a in artifacts if a["filename"].endswith(".h5mu"))
        self.assertEqual(len(h5mu["sha256"]), 64)
        self.assertNotEqual(set(h5mu["sha256"]), {"0"},
                            "a placeholder checksum must not be committed")
        self.assertEqual(h5mu["bytes"], 141704166)

    def test_it_is_open_without_an_authorisation_request(self):
        self.assertEqual(self.entry["access_tier"], "open_public")
        self.assertEqual(self.entry["access_requirement"]["type"], "none")


class TestWeglerCorrection(unittest.TestCase):
    """The entry understated its own deposit; the correction must stay."""

    def setUp(self):
        self.entry = dataset("wegler_15")

    def test_individual_level_data_is_now_recorded_as_released(self):
        self.assertTrue(self.entry["individual_level_data_released"])

    def test_the_deposit_rather_than_the_supplement_is_named(self):
        self.assertIn("PXD014131", self.entry["original_repository"])
        self.assertIn("proteinGroups.txt", self.entry["access_note"])

    def test_the_partial_flag_is_explained_rather_than_trusted(self):
        self.assertIn("submissionType PARTIAL", self.entry["access_note"])

    def test_the_donor_structure_is_verified_per_preparation(self):
        structure = self.entry["donor_structure_verified"]
        self.assertEqual(structure["homogenate"], 15)
        self.assertEqual(structure["microsomes"], 15)
        self.assertEqual(structure["hepatocytes"], 15)
        self.assertEqual(structure["low_speed_pellet"], 5)
        self.assertEqual(structure["pellet_donors"], [1, 7, 11, 12, 15])

    def test_the_pellet_limitation_is_a_caveat(self):
        self.assertIn("only 5 of the 15 donors",
                      " ".join(self.entry["caveats"]))

    def test_what_is_still_missing_is_named(self):
        self.assertIn("intrinsic-clearance", self.entry["open_requirement"])
        self.assertIn("NOT in", self.entry["access_note"])

    def test_the_licence_is_recorded(self):
        self.assertIn("CC0", self.entry["licence"])


class TestGrangeonIsNotADataset(unittest.TestCase):
    """A published article is not automatically an exploitable base."""

    def setUp(self):
        self.entry = dataset("grangeon_2025_cyp_protein")

    def test_it_abstains_with_a_reason(self):
        self.assertEqual(self.entry["status"], "ABSTAIN")
        self.assertTrue(self.entry["abstain_reason"].strip())
        self.assertEqual(self.entry["objectives_supported"], [])

    def test_individual_level_data_is_recorded_as_not_released(self):
        self.assertFalse(self.entry["individual_level_data_released"])
        self.assertFalse(
            self.entry["omics_layers"][0]["individual_values_published"])
        self.assertFalse(self.entry["endpoint"]["individual_values_published"])

    def test_the_reason_names_the_correlation_table(self):
        reason = self.entry["abstain_reason"]
        self.assertIn("CORRELATION table", reason)
        self.assertIn("digitising scatter figures", reason)

    def test_the_genotype_coverage_is_recorded_honestly(self):
        genotypes = self.entry["genotypes"]
        self.assertTrue(genotypes["individual_values_published"])
        coverage = genotypes["measured_coverage"]
        self.assertIn("10 of 16", coverage)
        self.assertIn("zero variance", coverage)

    def test_it_is_usable_only_as_a_comparator(self):
        self.assertEqual(self.entry["usable_as"], "published_comparator")
        comparator = self.entry["published_comparator_content"]
        self.assertEqual(comparator["hlm_n"], 16)
        self.assertEqual(comparator["hhep_n"], 20)

    def test_the_supplier_provenance_of_the_activity_is_a_caveat(self):
        self.assertIn("provided by the supplier",
                      " ".join(self.entry["caveats"]))


class TestRnaVersusProteinComparison(unittest.TestCase):
    """The headline number must travel with its label and its alternatives."""

    @classmethod
    def setUpClass(cls):
        if not ARTEFACT.is_file():
            raise unittest.SkipTest("the comparison has not been run yet")
        cls.report = json.loads(ARTEFACT.read_text(encoding="utf-8"))
        cls.doc = " ".join(DOC.read_text(encoding="utf-8").split())

    def test_it_is_labelled_cross_cohort(self):
        self.assertEqual(self.report["label"],
                         "DESCRIPTIVE_CROSS_COHORT_COMPARISON")
        self.assertIn("NOT a within-subject comparison",
                      self.report["what_this_is_not"])

    def test_it_refuses_the_pharmacokinetic_reading(self):
        text = self.report["what_this_is_not"]
        self.assertIn("no pharmacokinetic claim", text)
        self.assertIn("functional_in_vitro", text)

    def test_the_published_side_carries_its_own_caveats(self):
        caveats = " ".join(self.report["published_source_caveats"])
        self.assertIn("n = 16", caveats)
        self.assertIn("provided by the commercial supplier", caveats)
        self.assertIn("cannot be recomputed", caveats)

    def test_no_model_was_fitted(self):
        self.assertIn("no regression", self.report["hlc_method"]["no_model"])

    def test_cyp1a2_is_not_evaluable_as_in_every_other_hlc_analysis(self):
        row = self.report["per_endpoint"]["CYP1A2_activity"]
        self.assertEqual(row["status"], "NOT_EVALUABLE")
        self.assertIn("no admissible probe", row["reason"])

    def test_every_measured_endpoint_carries_an_interval(self):
        for endpoint, row in self.report["per_endpoint"].items():
            if row.get("status") != "MEASURED":
                continue
            self.assertLess(row["hlc_ci_low"],
                            row["hlc_transcript_activity_spearman"], endpoint)
            self.assertGreater(row["hlc_ci_high"],
                               row["hlc_transcript_activity_spearman"], endpoint)
            self.assertGreaterEqual(row["hlc_n_donors"], 300, endpoint)

    def test_the_summary_matches_the_per_endpoint_rows(self):
        import statistics

        rows = [r for r in self.report["per_endpoint"].values()
                if r.get("status") == "MEASURED"]
        median = statistics.median(
            r["hlc_transcript_activity_spearman"] for r in rows)
        self.assertAlmostEqual(
            self.report["summary"]["hlc_transcript_rho_median"], median,
            places=10)
        self.assertEqual(self.report["summary"]["endpoints_measured"], len(rows))

    def test_the_strict_count_is_computed_not_asserted(self):
        rows = [r for r in self.report["per_endpoint"].values()
                if r.get("status") == "MEASURED"]
        above = sum(1 for r in rows
                    if r["published_protein"]["hlm_r"] > r["hlc_ci_high"])
        self.assertEqual(
            self.report["summary"][
                "endpoints_where_published_protein_r_exceeds_the_hlc_interval"],
            above)

    def test_the_document_states_three_competing_explanations(self):
        self.assertIn("Three explanations", self.doc)
        for mechanism in ("Biology", "Assay precision", "Shared preparation"):
            self.assertIn(mechanism, self.doc)

    def test_the_document_does_not_claim_biology_alone(self):
        self.assertIn("motivates a measurement rather than settling a question",
                      self.doc)

    def test_the_transcribed_published_values_are_in_the_script(self):
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertIn("0.8168", source)
        self.assertIn("Pharmaceuticals 2025", source)
        self.assertIn("CC BY 4.0", source)


class TestQualificationDocument(unittest.TestCase):
    """Document assertions run on whitespace-normalised text.

    A prose commitment can be line-wrapped at any column by an editor, so
    matching the raw text would make these tests fail on reformatting rather
    than on a change of meaning.
    """

    @classmethod
    def setUpClass(cls):
        raw = DOC.read_text(encoding="utf-8")
        cls.doc = " ".join(raw.split())

    def test_it_gives_a_verdict_per_source(self):
        for verdict in ("USABLE", "COMPARATOR ONLY"):
            self.assertIn(verdict, self.doc)

    def test_it_records_the_eatris_layer_trade_off(self):
        self.assertIn("55", self.doc)
        self.assertIn("79", self.doc)
        self.assertIn("114", self.doc)

    def test_it_forbids_cross_cohort_joining_and_figure_digitisation(self):
        self.assertIn("no analysis joining individuals across cohorts", self.doc)
        self.assertIn("no Grangeon figure digitisation", self.doc)

    def test_it_names_the_data_requirement_rather_than_assuming_it(self):
        self.assertIn("No public dataset has that", self.doc)
        self.assertIn("all different people", self.doc)

    def test_it_keeps_the_registry_amendment_list_explicit(self):
        self.assertIn("What this changes in the registry", self.doc)


class TestReadinessGateExecution(unittest.TestCase):
    """The first execution of a gate that previously had no recorded verdict.

    The first version of the runner invented its own `outcome_kind` vocabulary
    and every source came back `ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME` with zero
    donor counts -- a report that looked plausible and said nothing. These tests
    exist so that failure mode cannot recur silently.
    """

    VERDICT = REPO_ROOT / "reports" / "data" / (
        "hepatic_multiview_readiness_verdict.json")
    RUNNER = REPO_ROOT / "scripts" / (
        "run_hepatic_multiview_readiness_on_registered_sources.py")

    @classmethod
    def setUpClass(cls):
        if not cls.VERDICT.is_file():
            raise unittest.SkipTest("the gate has not been executed yet")
        cls.report = json.loads(cls.VERDICT.read_text(encoding="utf-8"))
        cls.cases = {c["source_id"]: c for c in cls.report["cases"]}

    def test_no_source_is_allowed_to_fit(self):
        self.assertEqual(self.report["verdict"]["fit_allowed"], 0)
        self.assertEqual(self.report["verdict"]["fit_allowed_sources"], [])
        for source, row in self.cases.items():
            self.assertFalse(row["fit_allowed"], source)

    def test_the_cross_cohort_fusion_is_refused(self):
        fusion = self.cases["INVALID_FUSION_hlc_rna_plus_wegler_protein"]
        self.assertFalse(fusion["fit_allowed"])
        self.assertEqual(fusion["status"], "ABSTAIN_NO_SAME_DONOR_INTERSECTION")
        self.assertTrue(self.report["verdict"]["cross_cohort_fusion_refused"])

    def test_the_fusion_case_was_given_real_donor_counts(self):
        """A refusal on empty inputs would prove nothing about the gate."""
        fusion = self.cases["INVALID_FUSION_hlc_rna_plus_wegler_protein"]
        self.assertEqual(fusion["n_rna"], 466)
        self.assertEqual(fusion["n_protein"], 15)
        self.assertGreater(fusion["n_outcome"], 400)
        self.assertEqual(fusion["n_complete_case"], 0)

    def test_the_refusals_are_not_all_the_same_code(self):
        """A uniform code across every source is the signature of a wrong input."""
        codes = {row["status"] for row in self.cases.values()}
        self.assertGreaterEqual(len(codes), 4, sorted(codes))

    def test_each_refusal_cause_is_distinguished(self):
        self.assertEqual(self.cases["hlc_syn4499"]["status"],
                         "ABSTAIN_NO_MATCHED_PROTEIN")
        self.assertEqual(self.cases["wegler_15_microsomes"]["status"],
                         "ABSTAIN_NO_MATCHED_RNA")
        self.assertEqual(self.cases["eatris_plus_v4"]["status"],
                         "ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME")
        self.assertEqual(self.cases["grangeon_2025_cyp_protein"]["status"],
                         "ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED")

    def test_layer_absent_is_not_reported_as_rows_missing(self):
        """ARIANA's rows are held; it is the layer that does not exist."""
        for source in ("ariana_metoprolol_47", "jarmusch_midazolam"):
            row = self.cases[source]
            self.assertTrue(row["inputs"]["source_rows_acquired"], source)
            self.assertEqual(row["status"], "ABSTAIN_NO_MATCHED_RNA", source)

    def test_hlc_is_credited_with_the_transcriptome_it_has(self):
        row = self.cases["hlc_syn4499"]
        self.assertEqual(row["n_rna"], 466)
        self.assertEqual(row["n_protein"], 0)
        self.assertGreater(row["n_outcome"], 400)

    def test_eatris_rna_and_protein_are_recorded_as_same_donor(self):
        row = self.cases["eatris_plus_v4"]
        self.assertTrue(row["inputs"]["rna_donors_supplied"])
        self.assertTrue(row["inputs"]["protein_donors_supplied"])
        self.assertIn("88 individuals carry both", row["note"])

    def test_the_eatris_compartment_limit_travels_with_its_verdict(self):
        row = self.cases["eatris_plus_v4"]
        self.assertIn("BLOOD", row["compartment"])
        self.assertIn("Not hepatic", row["compartment"])

    def test_the_gate_does_not_claim_to_test_compartment(self):
        self.assertIn("does NOT test tissue compartment",
                      self.report["how_to_read"])

    def test_an_abstain_is_called_a_result_not_a_failure(self):
        self.assertIn("recorded scientific result",
                      self.report["how_to_read"])

    def test_the_runner_probes_the_accepted_vocabulary_before_running(self):
        source = self.RUNNER.read_text(encoding="utf-8")
        self.assertIn("def check_outcome_vocabulary", source)
        self.assertIn("check_outcome_vocabulary()", source)
        self.assertIn("would report a uniform refusal that means nothing",
                      source)

    def test_the_accepted_vocabulary_is_the_gate_s_own(self):
        from omics2pk.hepatic_multiview_readiness import (
            assess_same_donor_hepatic_multiview as gate)

        probe = dict(source_id="probe", rna_donors=["a"], protein_donors=["a"],
                     outcome_donors=["a"], source_rows_acquired=True)
        for kind in ("drug_specific_clint", "drug_specific_metabolic_activity",
                     "clinical_pk"):
            self.assertNotEqual(gate(outcome_kind=kind, **probe).status,
                                "ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME", kind)
        self.assertEqual(gate(outcome_kind="invented", **probe).status,
                         "ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME")

    def test_the_runner_refuses_to_write_if_a_fusion_were_allowed(self):
        source = self.RUNNER.read_text(encoding="utf-8")
        self.assertIn("refusing to write a report", source)

    def test_the_headline_states_the_gap_rather_than_a_capability(self):
        headline = self.report["headline"]
        self.assertIn("No registered source permits", headline)
        self.assertIn("EATRIS-Plus is the", headline)


if __name__ == "__main__":
    unittest.main()
