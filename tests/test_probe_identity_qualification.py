"""Contract tests for graded probe identity and the qualified pools.

The verdict has to be graded, because an earlier binary version called two
unequal things identified: a real concordance between two annotations, and a
lone annotation with nothing to agree with it. And two different accessions are
not a conflict when they designate transcripts of the same gene -- comparing
accession strings manufactured several times more conflicts than exist.
"""
from __future__ import annotations

import json
import unittest
from pathlib import Path

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

from omics2pk.probe_identity_qualification import (
    ADMISSIBLE_BY_DEFAULT,
    VERDICTS,
    ProbeIdentityError,
    build_probe_pool,
    qualify_probe,
    summarise,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
POOLS = REPO_ROOT / "registry" / "QUALIFIED_PROBE_POOLS.json"

HUMAN = {
    "1565": ("CYP2D6", frozenset({"CPD6", "CYP2D"})),
    "1571": ("CYP2E1", frozenset()),
    "1576": ("CYP3A4", frozenset({"CYP3A3"})),
}
ACCESSIONS = {
    "NM_000106": ["1565"],
    "NM_000773": ["1571"],
    "NM_017460": ["1576"],
    "D00003": ["1576"],
}


def grade(feature="f1", **kwargs):
    base = {"symbol": "", "archive_accession": "", "platform_accession": "",
            "platform_entrez": "",
            "accession_to_human_genes": ACCESSIONS,
            "human_gene_symbols": HUMAN}
    base.update(kwargs)
    return qualify_probe(feature, **base)


class TestTheEvidenceIsGraded(unittest.TestCase):
    def test_two_sources_on_the_same_gene_are_concordant(self):
        identity = grade(symbol="CYP2D6", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertEqual(identity.verdict, "concordant_human")
        self.assertEqual(identity.human_gene, "1565")
        self.assertTrue(identity.symbol_matches_gene)
        self.assertTrue(identity.admissible)

    def test_different_accessions_of_one_gene_are_not_a_conflict(self):
        # The correction that mattered quantitatively: the archive cites one
        # transcript, the platform another, and both resolve to gene 1576.
        identity = grade(symbol="CYP3A4", archive_accession="D00003",
                         platform_accession="NM_017460", platform_entrez="1576")
        self.assertEqual(identity.verdict, "concordant_human")
        self.assertEqual(identity.human_gene, "1576")

    def test_one_source_alone_is_not_a_concordance(self):
        identity = grade(symbol="CYP2E1", archive_accession="NM_000773")
        self.assertEqual(identity.verdict, "single_source_human")
        self.assertIn("no second annotation", identity.reason)
        self.assertTrue(identity.admissible)

    def test_a_platform_only_human_assignment_is_also_single_source(self):
        identity = grade(symbol="CYP2D6", platform_entrez="1565")
        self.assertEqual(identity.verdict, "single_source_human")
        self.assertIn("only the platform", identity.reason)

    def test_two_sources_naming_different_genes_conflict(self):
        identity = grade(symbol="CYP2D6", archive_accession="NM_000106",
                         platform_accession="NM_000773", platform_entrez="1571")
        self.assertEqual(identity.verdict, "conflicting")
        self.assertFalse(identity.admissible)

    def test_a_non_human_gene_assignment_against_a_human_accession_conflicts(self):
        identity = grade(symbol="CYP2D6", archive_accession="NM_000106",
                         platform_entrez="479559")
        self.assertEqual(identity.verdict, "conflicting")
        self.assertIn("not a human gene", identity.reason)

    def test_a_non_human_assignment_with_nothing_else_conflicts(self):
        identity = grade(symbol="CYP1A2", archive_accession="AF123056",
                         platform_accession="DN745607",
                         platform_entrez="479559")
        self.assertEqual(identity.verdict, "conflicting")

    def test_nothing_resolvable_is_undetermined_not_wrong(self):
        identity = grade(symbol="MYSTERY1", archive_accession="ZZ999999")
        self.assertEqual(identity.verdict, "undetermined")
        self.assertFalse(identity.admissible)

    def test_a_control_needs_no_annotation(self):
        self.assertEqual(grade(symbol="NegCont_7").verdict, "technical_control")

    def test_a_symbol_that_disagrees_with_the_gene_is_flagged_not_rejected(self):
        identity = grade(symbol="SOMETHINGELSE", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertEqual(identity.verdict, "concordant_human")
        self.assertFalse(identity.symbol_matches_gene)

    def test_a_listed_synonym_counts_as_a_match(self):
        identity = grade(symbol="CPD6", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertTrue(identity.symbol_matches_gene)

    def test_an_accession_version_suffix_does_not_defeat_the_lookup(self):
        identity = grade(symbol="CYP2D6", archive_accession="NM_000106.6",
                         platform_accession="NM_000106.5",
                         platform_entrez="1565")
        self.assertEqual(identity.verdict, "concordant_human")

    def test_the_verdicts_are_ordered_strongest_first(self):
        self.assertEqual(VERDICTS[0], "concordant_human")
        self.assertEqual(ADMISSIBLE_BY_DEFAULT,
                         ("concordant_human", "single_source_human"))


class TestThePoolIsBuiltByProbe(unittest.TestCase):
    def identities(self):
        return [
            grade("f1", symbol="CYP2E1", archive_accession="NM_000773",
                  platform_accession="NM_000773", platform_entrez="1571"),
            grade("f2", symbol="CYP2E1", archive_accession="AY635465",
                  platform_accession="CO608737", platform_entrez="479559"),
            grade("f3", symbol="CYP2D6", archive_accession="NM_000106",
                  platform_accession="NM_000106", platform_entrez="1565"),
            grade("f4", symbol="NegCont_1"),
        ]

    def test_a_mixed_gene_keeps_its_admissible_probe(self):
        pool = build_probe_pool(self.identities())
        self.assertIn("CYP2E1", pool["by_symbol"])
        self.assertEqual(pool["by_symbol"]["CYP2E1"], ["f1"])
        self.assertIn("CYP2E1", pool["symbols_keeping_only_some_probes"])

    def test_a_control_never_enters_a_pool(self):
        pool = build_probe_pool(self.identities())
        self.assertNotIn("NegCont_1", pool["by_symbol"])
        self.assertEqual(pool["probes_dropped_by_verdict"]["technical_control"], 1)

    def test_admitting_a_control_is_refused(self):
        with self.assertRaisesRegex(ProbeIdentityError, "never admissible"):
            build_probe_pool(self.identities(),
                             admit=("concordant_human", "technical_control"))

    def test_an_unknown_verdict_is_refused(self):
        with self.assertRaisesRegex(ProbeIdentityError, "unknown verdict"):
            build_probe_pool(self.identities(), admit=("probably_fine",))

    def test_the_strict_pool_admits_no_more_than_the_default_one(self):
        strict = build_probe_pool(self.identities(), admit=("concordant_human",))
        default = build_probe_pool(self.identities())
        self.assertLessEqual(strict["probes"], default["probes"])

    def test_restricting_to_a_gene_set_intersects_rather_than_replaces(self):
        pool = build_probe_pool(self.identities(), restrict_to_genes=["CYP2D6"])
        self.assertEqual(set(pool["by_symbol"]), {"CYP2D6"})
        self.assertTrue(pool["restricted_to_gene_set"])

    def test_the_pool_records_why_it_groups_after_grading(self):
        pool = build_probe_pool(self.identities())
        self.assertIn("probe by probe, then grouped", pool["note"])


class TestTheSummaryRefusesToOverclaim(unittest.TestCase):
    def many(self):
        return [
            grade("f1", symbol="CYP2D6", archive_accession="NM_000106",
                  platform_accession="NM_000106", platform_entrez="1565"),
            grade("f2", symbol="CYP2E1", archive_accession="NM_000773"),
            grade("f3", symbol="MYSTERY1"),
            grade("f4", symbol="NegCont_2"),
        ]

    def test_it_says_what_identified_does_not_mean(self):
        text = summarise(self.many())["what_identified_does_not_mean"]
        self.assertIn("hybridisation specificity", text)
        self.assertIn("identified by annotation", text)
        self.assertIn("not a wrong measurement", text)

    def test_the_admissible_count_is_the_two_strongest_levels(self):
        overview = summarise(self.many())
        self.assertEqual(
            overview["admissible_by_default"],
            overview["counts_by_verdict"]["concordant_human"]
            + overview["counts_by_verdict"]["single_source_human"])


class TestTheCommittedPools(unittest.TestCase):
    def pools(self):
        return json.loads(POOLS.read_text(encoding="utf-8"))

    def test_the_pool_file_exists_and_covers_the_whole_array(self):
        self.assertTrue(POOLS.is_file())
        self.assertEqual(self.pools()["overview"]["probes"], 40638)

    def test_every_source_carries_a_checksum(self):
        provenance = self.pools()["provenance"]
        for key in ("platform_record", "human_gene_reference",
                    "accession_to_gene_map"):
            self.assertEqual(len(provenance[key]["sha256"]), 64, key)

    def test_the_rule_version_is_recorded(self):
        self.assertEqual(self.pools()["rule_version"], "probe-identity-1.0.0")

    def test_the_complete_admissible_list_is_committed_not_just_a_summary(self):
        features = self.pools()["admissible_features"]
        for name in ("A_qualified_human", "A_strict_concordant_only",
                     "B_qualified_human_on_network"):
            self.assertGreater(len(features[name]), 1000, name)

    def test_the_strata_are_kept_rather_than_deleted(self):
        strata = self.pools()["strata_kept_separate"]
        self.assertGreater(strata["undetermined_probes"], 0)
        self.assertIn("not a wrong measurement", strata["why"])
        self.assertIn("states which strata it admits", strata["why"])

    def test_the_conflict_count_is_the_corrected_one(self):
        # Comparing accession strings produced about seven thousand apparent
        # conflicts; resolving accessions to genes leaves far fewer.
        counts = self.pools()["overview"]["counts_by_verdict"]
        self.assertLess(counts["conflicting"], 3000)
        self.assertGreater(counts["concordant_human"], 20_000)

    def test_pool_b_is_a_subset_of_pool_a(self):
        features = self.pools()["admissible_features"]
        self.assertTrue(set(features["B_qualified_human_on_network"])
                        <= set(features["A_qualified_human"]))

    def test_the_strict_pool_is_a_subset_of_the_default_pool(self):
        features = self.pools()["admissible_features"]
        self.assertTrue(set(features["A_strict_concordant_only"])
                        <= set(features["A_qualified_human"]))

    def test_the_four_probes_of_the_withdrawn_target_are_all_inadmissible(self):
        pools = self.pools()
        conflicting = set(pools["conflicting_features"])
        for feature in ("10033668809", "10033668817", "10033668906",
                        "10033668921"):
            self.assertIn(feature, conflicting, msg=feature)
            for name in pools["admissible_features"]:
                self.assertNotIn(feature, pools["admissible_features"][name])


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
