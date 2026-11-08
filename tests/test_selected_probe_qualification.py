"""Contract tests for the selected-probe qualification.

The point of the qualifier is that it decides from the public record and not
from the look of a symbol. These tests pin both halves of that: a non-human
accession is called non-human however biological its symbol looks, and a human
record is left alone however much its symbol resembles a pseudogene.

No test here reaches the network. The classifier is exercised against records
supplied inline, and the committed report is checked as data.
"""
from __future__ import annotations

import importlib.util
import json
import unittest
from pathlib import Path

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "qualify_selected_probes.py"
REPORT = REPO_ROOT / "docs" / "H2A_SELECTED_PROBE_QUALIFICATION.json"
ERRATUM = REPO_ROOT / "docs" / "OMICS2PK_H2A_ERRATUM_20261109.md"


def load():
    spec = importlib.util.spec_from_file_location("qualify_probes", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


Q = load()


def record(organism, molecule="mRNA", definition="x"):
    return {"organism": organism, "molecule": molecule, "definition": definition}


class TestTheVerdictComesFromTheRecord(unittest.TestCase):
    def test_a_macaque_accession_is_non_human_however_the_symbol_looks(self):
        verdict, reason = Q.classify(
            "CYP1A2", ["AB185338"], {"AB185338": record("Macaca fuscata")})
        self.assertEqual(verdict, "non_human")
        self.assertIn("Macaca fuscata", reason)

    def test_a_human_accession_is_human_however_the_symbol_looks(self):
        # A pseudogene-looking symbol must not be condemned by its name.
        verdict, _ = Q.classify(
            "CYP2B7P1", ["M29873"], {"M29873": record("Homo sapiens")})
        self.assertEqual(verdict, "human_transcript")

    def test_a_human_genomic_record_is_distinguished_from_a_transcript(self):
        verdict, reason = Q.classify(
            "CYP2A7P1", ["NG_000008"],
            {"NG_000008": record("Homo sapiens", molecule="DNA")})
        self.assertEqual(verdict, "human_genomic")
        self.assertIn("not a transcript", reason)

    def test_accessions_spanning_organisms_are_flagged_not_resolved(self):
        verdict, reason = Q.classify(
            "CYP2E1", ["NM_000773", "AY635465"],
            {"NM_000773": record("Homo sapiens"),
             "AY635465": record("Macaca mulatta")})
        self.assertEqual(verdict, "mixed_organism")
        self.assertIn("Homo sapiens", reason)

    def test_a_technical_control_is_recognised_without_a_lookup(self):
        verdict, reason = Q.classify("NegCont_210", [], {})
        self.assertEqual(verdict, "technical_control")
        self.assertIn("control feature", reason)

    def test_control_patterns_cover_the_platform_s_families(self):
        for name in ("NegCont_1", "PosCont_9", "Corner", "DarkCorner",
                     "SpikeIn_3", "GE_BrightCorner", "E1A_r60_1", "ETG01"):
            self.assertEqual(Q.classify(name, [], {})[0], "technical_control",
                             msg=name)

    def test_a_missing_record_is_unverified_not_a_guess(self):
        verdict, reason = Q.classify("MYSTERY1", ["ZZ999999"], {})
        self.assertEqual(verdict, "unverified")
        self.assertIn("no public record", reason)

    def test_a_record_without_an_organism_is_unverified(self):
        verdict, _ = Q.classify("MYSTERY1", ["ZZ999999"],
                                {"ZZ999999": record(None)})
        self.assertEqual(verdict, "unverified")


class TestSelectionHarvesting(unittest.TestCase):
    def test_every_arm_and_endpoint_that_chose_a_symbol_is_recorded(self):
        result = {"endpoints": {
            "A_activity": {"first_fold_selection": {
                "univariate_spearman": {"1": ["GENEA"], "5": ["GENEA", "GENEB"]},
                "netpk_target_seeded": {"1": ["GENEA"]}}},
            "B_activity": {"first_fold_selection": {
                "xgboost_gain": {"1": ["GENEB"]}}},
        }}
        chosen = Q.selected_symbols(result)
        self.assertEqual(set(chosen), {"GENEA", "GENEB"})
        self.assertEqual(chosen["GENEA"]["arms"],
                         {"univariate_spearman", "netpk_target_seeded"})
        self.assertEqual(chosen["GENEB"]["endpoints"],
                         {"A_activity", "B_activity"})

    def test_an_endpoint_without_selections_is_skipped_quietly(self):
        self.assertEqual(Q.selected_symbols({"endpoints": {"A": {}}}), {})


class TestTheCommittedReport(unittest.TestCase):
    def report(self):
        return json.loads(REPORT.read_text(encoding="utf-8"))

    def test_the_report_exists_and_states_its_source_of_truth(self):
        self.assertTrue(REPORT.is_file())
        self.assertIn("never inferred from the", self.report()["source_of_truth"])

    def test_the_target_gene_of_one_endpoint_is_non_human(self):
        # The defect the erratum withdraws. Pinned so it cannot be forgotten.
        entry = self.report()["probes"]["CYP1A2"]
        self.assertEqual(entry["verdict"], "non_human")
        self.assertTrue(entry["organism"].startswith("Macaca"))
        self.assertEqual(len(entry["accessions"]), 4)

    def test_the_other_negative_endpoint_s_target_probe_is_human(self):
        entry = self.report()["probes"]["CYP2E1"]
        self.assertEqual(entry["verdict"], "mixed_organism")
        self.assertIn("NM_000773", entry["accessions"])

    def test_the_statistical_arms_selected_inadmissible_probes(self):
        by_arm = self.report()["by_arm"]
        for arm in ("univariate_spearman", "xgboost_gain", "lasso_entry_order"):
            self.assertGreater(len(by_arm[arm]["inadmissible"]), 0, arm)

    def test_unverified_probes_are_counted_rather_than_assumed_fine(self):
        counts = self.report()["counts_by_verdict"]
        self.assertGreater(counts.get("unverified", 0), 0)
        self.assertGreater(counts.get("human_transcript", 0), 0)

    def test_the_report_refuses_both_over_readings(self):
        text = self.report()["how_to_read"]
        self.assertIn("neither can be reported as a distinct biomarker", text)
        self.assertIn("NOT dismissed because its symbol looks like a", text)


class TestTheErratumWithdrawsWhatItMust(unittest.TestCase):
    def text(self):
        return ERRATUM.read_text(encoding="utf-8")

    def test_the_erratum_exists(self):
        self.assertTrue(ERRATUM.is_file())

    def test_it_withdraws_the_cyp1a2_conclusion_explicitly(self):
        flat = " ".join(self.text().split())
        self.assertIn("is withdrawn", flat)
        self.assertIn("never measured human CYP1A2 messenger RNA", flat)
        self.assertIn("It was never tested", flat)

    def test_it_keeps_the_endpoint_that_survives_the_check(self):
        flat = " ".join(self.text().split())
        self.assertIn("CYP2E1 survives the same check", flat)
        self.assertIn("Not two.", flat)

    def test_it_requires_an_organism_check_before_the_next_experiment(self):
        flat = " ".join(self.text().split())
        self.assertIn("must verify the organism of the probe", flat)
        self.assertIn("Required before H2d", flat)

    def test_it_says_each_restriction_cuts_both_ways(self):
        flat = " ".join(self.text().split())
        self.assertIn("moves the comparison in both directions at once", flat)
        self.assertIn("All three are reported", flat)

    def test_it_recomputes_no_frozen_number(self):
        flat = " ".join(self.text().split())
        self.assertIn("no h2a number is recomputed", flat.lower())


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
