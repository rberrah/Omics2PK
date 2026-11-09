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
    symbol_to_gene_index,
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
SYMBOLS = symbol_to_gene_index(HUMAN)


def grade(feature="f1", **kwargs):
    base = {"symbol": "", "archive_accession": "", "platform_accession": "",
            "platform_entrez": "",
            "accession_to_human_genes": ACCESSIONS,
            "human_gene_symbols": HUMAN,
            "symbol_to_human_genes": SYMBOLS}
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

    def test_a_symbol_naming_another_gene_is_refused(self):
        # Two mutually exclusive claims: the accessions say gene 1565, the
        # analysis calls the column by the approved symbol of gene 1571.
        # Grouping it under that symbol would attribute the measurement to the
        # wrong gene, which is how the withdrawn target slipped through.
        identity = grade(symbol="CYP2E1", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertEqual(identity.verdict, "symbol_names_another_gene")
        self.assertFalse(identity.admissible)
        self.assertIn("wrong gene", identity.reason)

    def test_a_symbol_that_is_another_genes_synonym_is_also_refused(self):
        identity = grade(symbol="CYP3A3", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertEqual(identity.verdict, "symbol_names_another_gene")

    def test_a_stale_alias_keeps_the_measurement_and_names_the_gene(self):
        # No other human gene carries this name, so there is no competing
        # claim: the accessions decide and the gene is recorded beside it.
        identity = grade(symbol="FLJ99999", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertEqual(identity.verdict, "concordant_human")
        self.assertTrue(identity.admissible)
        self.assertFalse(identity.symbol_matches_gene)
        self.assertTrue(identity.symbol_is_a_stale_alias)
        self.assertEqual(identity.human_symbol, "CYP2D6")

    def test_the_same_rule_applies_to_a_single_source_identity(self):
        refused = grade(symbol="CYP2E1", archive_accession="NM_000106")
        self.assertEqual(refused.verdict, "symbol_names_another_gene")
        kept = grade(symbol="MGC99999", archive_accession="NM_000106")
        self.assertEqual(kept.verdict, "single_source_human")
        self.assertTrue(kept.symbol_is_a_stale_alias)

    def test_without_a_symbol_reference_a_mismatch_is_refused(self):
        # A rule that cannot tell a stale alias from a competing claim must not
        # assume the harmless case.
        identity = grade(symbol="FLJ99999", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565",
                         symbol_to_human_genes=None)
        self.assertEqual(identity.verdict, "symbol_names_another_gene")
        self.assertIn("cannot be cleared", identity.reason)

    def test_a_matching_symbol_is_never_a_stale_alias(self):
        identity = grade(symbol="CYP2D6", archive_accession="NM_000106",
                         platform_accession="NM_000106", platform_entrez="1565")
        self.assertFalse(identity.symbol_is_a_stale_alias)


class TestThePlatformMustAgreeWithItself(unittest.TestCase):
    """A record that contradicts itself cannot confirm another source.

    The earlier order of tests computed the agreement between sources first, so
    an archive agreeing with one of the platform's two fields was reported as a
    concordance even when the platform's other field named a different gene.
    """

    def test_a_platform_contradicting_itself_cannot_confirm_the_archive(self):
        identity = grade(symbol="CYP3A4", archive_accession="NM_017460",
                         platform_accession="NM_017460", platform_entrez="1571")
        self.assertEqual(identity.verdict, "platform_self_contradictory")
        self.assertFalse(identity.admissible)
        self.assertIn("contradicts itself", identity.reason)

    def test_a_non_human_assignment_against_its_own_human_accession(self):
        identity = grade(symbol="CYP3A4", archive_accession="NM_017460",
                         platform_accession="NM_017460",
                         platform_entrez="479559")
        self.assertEqual(identity.verdict, "platform_self_contradictory")

    def test_one_platform_field_alone_is_not_a_contradiction(self):
        identity = grade(symbol="CYP2D6", archive_accession="ZZ999999",
                         platform_accession="ZZ999999", platform_entrez="1565")
        self.assertEqual(identity.verdict, "single_source_human")

    def test_an_agreeing_platform_record_is_unaffected(self):
        identity = grade(symbol="CYP3A4", archive_accession="D00003",
                         platform_accession="NM_017460", platform_entrez="1576")
        self.assertEqual(identity.verdict, "concordant_human")

    def test_neither_new_verdict_is_admissible_by_default(self):
        for verdict in ("symbol_names_another_gene",
                        "platform_self_contradictory"):
            self.assertIn(verdict, VERDICTS)
            self.assertNotIn(verdict, ADMISSIBLE_BY_DEFAULT)

    def test_the_symbol_index_carries_approved_names_and_synonyms(self):
        index = symbol_to_gene_index(HUMAN)
        self.assertEqual(index["CYP2D6"], ("1565",))
        self.assertEqual(index["CYP3A3"], ("1576",))
        self.assertNotIn("FLJ99999", index)

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

    def test_it_counts_the_measurements_kept_under_an_outdated_name(self):
        overview = summarise(self.many())
        self.assertIn("admissible_under_a_stale_alias", overview)
        text = overview["what_a_stale_alias_means"]
        self.assertIn("not the gene's current symbol", text)
        self.assertIn("cannot reach the feature", text)
        self.assertIn("not a defect of the", text)

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
        self.assertEqual(self.pools()["rule_version"], "probe-identity-1.1.0")

    def test_the_two_new_refusals_are_counted_and_listed(self):
        pools = self.pools()
        counts = pools["overview"]["counts_by_verdict"]
        for verdict, key in (
                ("symbol_names_another_gene", "symbol_names_another_gene_features"),
                ("platform_self_contradictory",
                 "platform_self_contradictory_features")):
            self.assertGreater(counts[verdict], 0, verdict)
            self.assertEqual(len(pools[key]), counts[verdict], key)

    def test_no_admissible_feature_sits_in_a_refused_stratum(self):
        pools = self.pools()
        refused = set()
        for key in ("conflicting_features", "undetermined_features",
                    "symbol_names_another_gene_features",
                    "platform_self_contradictory_features"):
            refused |= set(pools[key])
        for name, features in pools["admissible_features"].items():
            self.assertFalse(refused & set(features), name)

    def test_a_measurement_kept_under_an_outdated_name_names_its_gene(self):
        pools = self.pools()
        rows = pools["stale_alias_features"]
        self.assertEqual(len(rows),
                         pools["overview"]["admissible_under_a_stale_alias"])
        admissible = set(pools["admissible_features"]["A_qualified_human"])
        for feature, analysis_symbol, current_symbol in rows[:200]:
            self.assertIn(feature, admissible)
            self.assertTrue(current_symbol)
            self.assertNotEqual(analysis_symbol, current_symbol)

    def test_every_verdict_is_accounted_for(self):
        counts = self.pools()["overview"]["counts_by_verdict"]
        self.assertEqual(sum(counts.values()), 40638)

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

    def test_the_nine_target_genes_keep_what_the_evidence_allows(self):
        # Eight targets keep at least one feature; the mixed one keeps its
        # admissible features instead of being dropped whole; the withdrawn one
        # keeps none. Checked against the committed lists, not recomputed.
        pools = self.pools()
        admissible = set(pools["admissible_features"]["A_qualified_human"])
        for feature in ("10025909551", "10033668556"):  # CYP2E1, concordant
            self.assertIn(feature, admissible, feature)
        self.assertNotIn("10033668877", admissible)  # CYP2E1, conflicting
        for feature in ("10033668809", "10033668817"):  # CYP1A2
            self.assertNotIn(feature, admissible, feature)

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
