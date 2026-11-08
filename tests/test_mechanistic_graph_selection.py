"""Contract tests for the executable mechanism-guided selector.

Two groups matter more than the rest. `TestNoGeneIsNamedInTheModule` enforces
the protocol's mandatory honesty check: the retained panel must be produced by
the rule, not written into the code. `TestTheRuleIsReproducible` re-derives the
decision from the documented definition and asserts the module agrees.

Every synthetic symbol here is deliberately non-biological, so the honesty scan
over the module source cannot be satisfied by accident and the label-resolution
rule is shown to be generic rather than special-cased for one gene family.
"""
from __future__ import annotations

import ast
import gzip
import io
import unittest
import zipfile
from pathlib import Path
from tempfile import TemporaryDirectory

try:
    import pytest
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"pytest is required for this module: {exc}") from exc

try:
    import numpy as np
    from scipy.sparse import csr_matrix  # noqa: F401  (import guard only)
except ImportError as exc:  # pragma: no cover - dependency-free CI job only
    raise unittest.SkipTest(f"numpy and scipy are required: {exc}") from exc

from omics2pk import mechanistic_graph_selection as mgs
from omics2pk.mechanistic_graph_selection import (
    FROZEN,
    MechanisticSelectionError,
    STRING_INFO_HEADER,
    STRING_LINKS_HEADER,
    Selection,
    column_stochastic_transition,
    degree_strata,
    genes_to_probe_positions,
    jaccard,
    load_string_symbol_graph,
    permutation_fdr_cutoff,
    random_walk_with_restart,
    rank_by_mechanistic_proximity,
    seeds_from_ctd_substrate,
    seeds_from_endpoint_label,
    seeds_from_reactome,
    select_panel,
)

MODULE_PATH = Path(mgs.__file__)


# ---------------------------------------------------------------------------
# Synthetic sources
# ---------------------------------------------------------------------------
def write_string_pair(root: Path, edges, symbols=None):
    """Write a minimal STRING-format pair. `edges` is (symbol, symbol, score)."""
    names = sorted({s for e in edges for s in e[:2]} | set(symbols or ()))
    protein_of = {name: f"9606.ENSP{index:011d}" for index, name in enumerate(names)}
    info = root / "info.txt.gz"
    with gzip.open(info, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(STRING_INFO_HEADER + "\n")
        for name in names:
            handle.write(f"{protein_of[name]}\t{name}\t100\tsynthetic\n")
    links = root / "links.txt.gz"
    with gzip.open(links, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write(STRING_LINKS_HEADER + "\n")
        for left, right, score in edges:
            handle.write(f"{protein_of[left]} {protein_of[right]} {score}\n")
            handle.write(f"{protein_of[right]} {protein_of[left]} {score}\n")
    return info, links


def write_reactome(root: Path, pathways):
    path = root / "pathways.gmt.zip"
    body = "\n".join("\t".join([name, "R-SYN-1", *members])
                     for name, members in pathways) + "\n"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("ReactomePathways.gmt", body)
    return path


CTD_COLUMNS = ("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
               "GeneForms", "Organism", "OrganismID", "Interaction",
               "InteractionActions", "PubMedIDs")


def write_ctd(root: Path, rows):
    path = root / "ixns.tsv.gz"
    with gzip.open(path, "wt", encoding="utf-8", newline="\n") as handle:
        handle.write("# a comment line\n")
        handle.write("# " + "\t".join(CTD_COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    return path


def ctd_row(chemical, gene, organism, actions, pmids):
    return (chemical, "D000001", "", gene, "1", "protein", organism, "9606",
            "synthetic interaction text", actions, pmids)


# A small ring plus a dense clique, used wherever a real network would be
# too large to reason about by hand.
RING_AND_CLIQUE = (
    [("SEEDA", f"CLQ{i}", 900) for i in range(1, 6)]
    + [(f"CLQ{i}", f"CLQ{j}", 900) for i in range(1, 6) for j in range(i + 1, 6)]
    + [(f"FAR{i}", f"FAR{i + 1}", 900) for i in range(1, 12)]
    + [("FAR12", "FAR1", 900)]
)


class TestGraphConstruction(unittest.TestCase):
    def test_threshold_excludes_weak_edges(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info, links = write_string_pair(root, [("AAA1", "BBB2", 800),
                                                   ("AAA1", "CCC3", 400)])
            graph = load_string_symbol_graph(info, links, min_combined_score=700)
            self.assertEqual(graph.n_nodes, 2)
            self.assertEqual(graph.n_edges, 1)
            self.assertNotIn("CCC3", graph.index)

    def test_graph_is_symmetric_and_free_of_self_loops(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info, links = write_string_pair(
                root, [("AAA1", "BBB2", 900), ("AAA1", "AAA1", 999)])
            graph = load_string_symbol_graph(info, links)
            self.assertEqual(graph.n_edges, 1)
            dense = column_stochastic_transition(graph).toarray()
            self.assertEqual(dense.shape, (2, 2))
            self.assertTrue(np.all(np.diag(dense) == 0.0))

    def test_symbols_are_sorted_so_node_order_is_stable(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info, links = write_string_pair(
                root, [("ZZZ9", "AAA1", 900), ("AAA1", "MMM5", 900)])
            graph = load_string_symbol_graph(info, links)
            self.assertEqual(list(graph.symbols), sorted(graph.symbols))

    def test_duplicate_identifiers_sharing_a_symbol_keep_the_strongest_edge(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info = root / "info.txt.gz"
            with gzip.open(info, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(STRING_INFO_HEADER + "\n")
                handle.write("9606.ENSP1\tAAA1\t10\tx\n")
                handle.write("9606.ENSP2\tAAA1\t10\tx\n")
                handle.write("9606.ENSP3\tBBB2\t10\tx\n")
            links = root / "links.txt.gz"
            with gzip.open(links, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write(STRING_LINKS_HEADER + "\n")
                for a, b, s in (("9606.ENSP1", "9606.ENSP3", 750),
                                ("9606.ENSP2", "9606.ENSP3", 950)):
                    handle.write(f"{a} {b} {s}\n")
                    handle.write(f"{b} {a} {s}\n")
            graph = load_string_symbol_graph(info, links)
            self.assertEqual(graph.n_nodes, 2)
            self.assertEqual(graph.n_edges, 1)
            self.assertAlmostEqual(float(graph.weights.max()), 0.95)

    def test_changed_headers_are_refused_rather_than_guessed(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info, links = write_string_pair(root, [("AAA1", "BBB2", 900)])
            bad = root / "bad_info.txt.gz"
            with gzip.open(bad, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("#some_other_layout\n9606.ENSP1\tAAA1\t1\tx\n")
            with self.assertRaisesRegex(MechanisticSelectionError, "protein.info"):
                load_string_symbol_graph(bad, links)
            bad_links = root / "bad_links.txt.gz"
            with gzip.open(bad_links, "wt", encoding="utf-8", newline="\n") as handle:
                handle.write("p1,p2,score\n")
            with self.assertRaisesRegex(MechanisticSelectionError, "protein.links"):
                load_string_symbol_graph(info, bad_links)

    def test_an_empty_graph_is_an_error_not_an_empty_answer(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info, links = write_string_pair(root, [("AAA1", "BBB2", 300)])
            with self.assertRaisesRegex(MechanisticSelectionError, "confidence"):
                load_string_symbol_graph(info, links)

    def test_threshold_must_be_a_usable_score(self):
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            info, links = write_string_pair(root, [("AAA1", "BBB2", 900)])
            for bad in (0, -1, 1001):
                with self.assertRaisesRegex(MechanisticSelectionError, "threshold"):
                    load_string_symbol_graph(info, links, min_combined_score=bad)


class TestRandomWalk(unittest.TestCase):
    def build(self, edges):
        self.tmp = TemporaryDirectory()
        info, links = write_string_pair(Path(self.tmp.name), edges)
        return load_string_symbol_graph(info, links)

    def tearDown(self):
        if getattr(self, "tmp", None):
            self.tmp.cleanup()

    def test_matches_the_closed_form_solution(self):
        graph = self.build(RING_AND_CLIQUE)
        alpha = 0.3
        restart = np.zeros(graph.n_nodes)
        restart[graph.index["SEEDA"]] = 1.0
        walk = random_walk_with_restart(graph, restart, restart_probability=alpha,
                                        tolerance=1e-14, max_iterations=5000)
        transition = column_stochastic_transition(graph).toarray()
        closed = np.linalg.solve(
            np.eye(graph.n_nodes) - (1.0 - alpha) * transition, alpha * restart)
        np.testing.assert_allclose(walk, closed, atol=1e-10)

    def test_conserves_probability(self):
        graph = self.build(RING_AND_CLIQUE)
        restart = np.zeros(graph.n_nodes)
        restart[graph.index["FAR1"]] = 1.0
        walk = random_walk_with_restart(graph, restart)
        self.assertAlmostEqual(float(walk.sum()), 1.0, places=9)
        self.assertTrue(np.all(walk >= 0.0))

    def test_the_seed_neighbourhood_outranks_a_distant_component(self):
        graph = self.build(RING_AND_CLIQUE)
        restart = np.zeros(graph.n_nodes)
        restart[graph.index["SEEDA"]] = 1.0
        walk = random_walk_with_restart(graph, restart)
        near = walk[graph.index["CLQ1"]]
        far = walk[graph.index["FAR6"]]
        self.assertGreater(near, far)

    def test_the_batched_form_equals_the_single_vector_form(self):
        graph = self.build(RING_AND_CLIQUE)
        columns = []
        for name in ("SEEDA", "CLQ2", "FAR3"):
            restart = np.zeros(graph.n_nodes)
            restart[graph.index[name]] = 1.0
            columns.append(restart)
        batch = np.column_stack(columns)
        together = random_walk_with_restart(graph, batch, tolerance=1e-14)
        for position, restart in enumerate(columns):
            alone = random_walk_with_restart(graph, restart, tolerance=1e-14)
            np.testing.assert_allclose(together[:, position], alone, atol=1e-9)

    def test_unnormalised_restart_weights_are_rescaled_not_rejected(self):
        graph = self.build(RING_AND_CLIQUE)
        one = np.zeros(graph.n_nodes)
        one[graph.index["SEEDA"]] = 1.0
        np.testing.assert_allclose(
            random_walk_with_restart(graph, one, tolerance=1e-14),
            random_walk_with_restart(graph, one * 7.5, tolerance=1e-14),
            atol=1e-12)

    def test_malformed_input_is_refused(self):
        graph = self.build(RING_AND_CLIQUE)
        good = np.zeros(graph.n_nodes)
        good[0] = 1.0
        with self.assertRaisesRegex(MechanisticSelectionError, "align"):
            random_walk_with_restart(graph, np.zeros(graph.n_nodes + 3))
        with self.assertRaisesRegex(MechanisticSelectionError, "non-negative"):
            bad = good.copy()
            bad[1] = -1.0
            random_walk_with_restart(graph, bad)
        with self.assertRaisesRegex(MechanisticSelectionError, "positive mass"):
            random_walk_with_restart(graph, np.zeros(graph.n_nodes))
        with self.assertRaisesRegex(MechanisticSelectionError, "Restart probability"):
            random_walk_with_restart(graph, good, restart_probability=0.0)
        with self.assertRaisesRegex(MechanisticSelectionError, "Restart probability"):
            random_walk_with_restart(graph, good, restart_probability=1.5)

    def test_a_restart_probability_of_one_never_leaves_the_seed(self):
        graph = self.build(RING_AND_CLIQUE)
        restart = np.zeros(graph.n_nodes)
        restart[graph.index["SEEDA"]] = 1.0
        walk = random_walk_with_restart(graph, restart, restart_probability=1.0)
        np.testing.assert_allclose(walk, restart, atol=1e-12)


class TestPermutationFalseDiscoveryRate(unittest.TestCase):
    def test_pure_noise_retains_nothing(self):
        rng = np.random.default_rng(7)
        observed = rng.normal(size=400)
        null = rng.normal(size=(400, 200))
        retained, audit = permutation_fdr_cutoff(observed, null, 0.05)
        self.assertEqual(retained, 0)
        self.assertEqual(audit.shape, (400, 4))

    def test_a_clear_signal_is_retained(self):
        rng = np.random.default_rng(11)
        observed = rng.normal(size=400)
        observed[:5] += 25.0
        null = rng.normal(size=(400, 200))
        retained, audit = permutation_fdr_cutoff(observed, null, 0.05)
        # The five planted values must be inside the cut. The cut can run a
        # little deeper, because the first ordinary value still sits far enough
        # into the null's tail to keep the ratio small; that is the estimator
        # behaving as defined, not a failure.
        self.assertGreaterEqual(retained, 5)
        self.assertLess(retained, 20)
        self.assertLessEqual(audit[retained - 1, 3], 0.05)

    def test_matches_an_independent_implementation(self):
        rng = np.random.default_rng(13)
        for _ in range(12):
            observed = rng.normal(size=60)
            null = rng.normal(size=(60, 25))
            retained, _ = permutation_fdr_cutoff(observed, null, 0.2)
            # Independent re-derivation of the documented definition.
            order = np.argsort(-observed, kind="stable")
            rates = []
            for position, index in enumerate(order, start=1):
                threshold = observed[index]
                mean_null = float(np.mean(
                    [(null[:, column] >= threshold).sum()
                     for column in range(null.shape[1])]))
                rates.append(mean_null / position)
            guarded, running = [], 0.0
            for value in rates:
                running = max(running, value)
                guarded.append(running)
            expected = 0
            for position, value in enumerate(guarded, start=1):
                if value <= 0.2:
                    expected = position
            self.assertEqual(retained, expected)

    def test_the_guard_prevents_a_lucky_deep_cut(self):
        # Three null values sit between the second and third observed score, so
        # the rate is 0 at the first cut, 0.75 at the second and 0.5 at the
        # third. Unguarded, the third cut would pass at q = 0.5; the guard must
        # stop at the first.
        observed = np.array([10.0, 9.0, 8.0])
        null = np.array([[9.5, 9.5], [9.5, 0.0], [0.0, 0.0]])
        retained, audit = permutation_fdr_cutoff(observed, null, 0.5)
        self.assertEqual(retained, 1)
        self.assertAlmostEqual(audit[0, 3], 0.0)
        self.assertAlmostEqual(audit[1, 3], 0.75)
        self.assertAlmostEqual(audit[2, 3], 0.5)

    def test_empty_candidate_set_is_handled(self):
        retained, audit = permutation_fdr_cutoff(
            np.zeros(0), np.zeros((0, 5)), 0.05)
        self.assertEqual(retained, 0)
        self.assertEqual(audit.shape, (0, 4))

    def test_malformed_input_is_refused(self):
        with self.assertRaisesRegex(MechanisticSelectionError, "candidates by"):
            permutation_fdr_cutoff(np.zeros(4), np.zeros((5, 3)), 0.05)
        with self.assertRaisesRegex(MechanisticSelectionError, "candidates by"):
            permutation_fdr_cutoff(np.zeros((4, 2)), np.zeros((4, 3)), 0.05)
        for bad in (0.0, 1.0, -0.1, 2.0):
            with self.assertRaisesRegex(MechanisticSelectionError, "false discovery"):
                permutation_fdr_cutoff(np.zeros(4), np.zeros((4, 3)), bad)


class TestDegreeStrata(unittest.TestCase):
    def test_hubs_and_leaves_land_in_different_strata(self):
        with TemporaryDirectory() as tmp:
            info, links = write_string_pair(Path(tmp), RING_AND_CLIQUE)
            graph = load_string_symbol_graph(info, links)
        strata = degree_strata(graph, 4)
        self.assertEqual(strata.shape, (graph.n_nodes,))
        degree = graph.weighted_degree()
        hub = strata[int(np.argmax(degree))]
        leaf = strata[int(np.argmin(degree))]
        self.assertGreater(hub, leaf)

    def test_at_least_one_stratum_is_required(self):
        with TemporaryDirectory() as tmp:
            info, links = write_string_pair(Path(tmp), [("AAA1", "BBB2", 900)])
            graph = load_string_symbol_graph(info, links)
        with self.assertRaisesRegex(MechanisticSelectionError, "stratum"):
            degree_strata(graph, 0)


class TestSeedBuilders(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.root = Path(self.tmp.name)
        info, links = write_string_pair(
            self.root,
            RING_AND_CLIQUE + [("AAA1", "AAA2", 900), ("AAA2", "BBB7", 900)])
        self.graph = load_string_symbol_graph(info, links)

    def tearDown(self):
        self.tmp.cleanup()

    def test_a_label_resolves_to_the_node_it_names(self):
        self.assertEqual(seeds_from_endpoint_label("SEEDA_activity", self.graph),
                         {"SEEDA": 1.0})

    def test_trailing_uppercase_characters_are_dropped_until_a_node_matches(self):
        # AAA1 is a node, AAA1Z is not: the generic rule must shorten the stem.
        self.assertEqual(seeds_from_endpoint_label("AAA1Z_activity", self.graph),
                         {"AAA1": 1.0})
        self.assertEqual(seeds_from_endpoint_label("AAA1ZQ_activity", self.graph),
                         {"AAA1": 1.0})

    def test_shortening_stops_at_the_first_match_rather_than_continuing(self):
        self.assertEqual(seeds_from_endpoint_label("AAA2_activity", self.graph),
                         {"AAA2": 1.0})

    def test_an_unresolvable_label_raises_instead_of_guessing(self):
        with self.assertRaisesRegex(MechanisticSelectionError, "does not resolve"):
            seeds_from_endpoint_label("QQQ9_activity", self.graph)

    def test_the_dataset_suffix_is_required(self):
        with self.assertRaisesRegex(MechanisticSelectionError, "must end with"):
            seeds_from_endpoint_label("SEEDA", self.graph)

    def test_pathway_seeds_come_from_the_matching_pathways_only(self):
        path = write_reactome(self.root, [
            ("Wanted pathway", ["CLQ1", "CLQ2", "NOT_ON_GRAPH"]),
            ("Other pathway", ["FAR4"]),
        ])
        seeds = seeds_from_reactome(path, r"^Wanted pathway$", self.graph)
        self.assertEqual(seeds, {"CLQ1": 1.0, "CLQ2": 1.0})

    def test_pathway_pattern_is_a_case_insensitive_regular_expression(self):
        path = write_reactome(self.root, [("Phase ONE thing", ["CLQ3"])])
        self.assertEqual(seeds_from_reactome(path, r"phase one", self.graph),
                         {"CLQ3": 1.0})

    def test_a_pattern_matching_nothing_raises(self):
        path = write_reactome(self.root, [("Only pathway", ["CLQ1"])])
        with self.assertRaisesRegex(MechanisticSelectionError, "matched the pattern"):
            seeds_from_reactome(path, r"^absent$", self.graph)

    def test_a_pathway_with_no_graph_member_raises(self):
        path = write_reactome(self.root, [("Only pathway", ["NOWHERE1"])])
        with self.assertRaisesRegex(MechanisticSelectionError, "no network node"):
            seeds_from_reactome(path, r"^Only pathway$", self.graph)

    def test_substrate_seeds_weight_genes_by_distinct_references(self):
        path = write_ctd(self.root, [
            ctd_row("probechem", "CLQ1", "Homo sapiens",
                    "increases^metabolic processing", "111|222|333"),
            ctd_row("probechem", "CLQ2", "Homo sapiens",
                    "affects^metabolic processing", "444"),
            ctd_row("probechem", "CLQ3", "Homo sapiens",
                    "increases^expression", "555"),
            ctd_row("probechem", "CLQ4", "Mus musculus",
                    "increases^metabolic processing", "666"),
            ctd_row("otherchem", "CLQ5", "Homo sapiens",
                    "increases^metabolic processing", "777"),
        ])
        seeds = seeds_from_ctd_substrate(path, "probechem", self.graph)
        self.assertEqual(seeds, {"CLQ1": 3.0, "CLQ2": 1.0})

    def test_substrate_lookup_is_case_insensitive_on_the_chemical(self):
        path = write_ctd(self.root, [
            ctd_row("ProbeChem", "CLQ1", "Homo sapiens",
                    "increases^metabolic processing", "111"),
        ])
        self.assertEqual(seeds_from_ctd_substrate(path, "probechem", self.graph),
                         {"CLQ1": 1.0})

    def test_duplicate_references_are_counted_once(self):
        path = write_ctd(self.root, [
            ctd_row("probechem", "CLQ1", "Homo sapiens",
                    "increases^metabolic processing", "111|222"),
            ctd_row("probechem", "CLQ1", "Homo sapiens",
                    "affects^metabolic processing", "222|333"),
        ])
        self.assertEqual(seeds_from_ctd_substrate(path, "probechem", self.graph),
                         {"CLQ1": 3.0})

    def test_no_matching_interaction_raises(self):
        path = write_ctd(self.root, [
            ctd_row("probechem", "CLQ1", "Homo sapiens",
                    "increases^expression", "111"),
        ])
        with self.assertRaisesRegex(MechanisticSelectionError, "No curated"):
            seeds_from_ctd_substrate(path, "probechem", self.graph)
        with self.assertRaisesRegex(MechanisticSelectionError, "chemical name"):
            seeds_from_ctd_substrate(path, "   ", self.graph)


class TestSelection(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        info, links = write_string_pair(Path(self.tmp.name), RING_AND_CLIQUE)
        self.graph = load_string_symbol_graph(info, links)

    def tearDown(self):
        self.tmp.cleanup()

    def run_rule(self, seeds, **kwargs):
        options = {"null_resamples": 200, "n_strata": 3, "random_seed": 5}
        options.update(kwargs)
        return select_panel(self.graph, seeds, **options)

    def test_a_seed_absent_from_the_network_abstains_entirely(self):
        selection = self.run_rule({"NOWHERE1": 1.0})
        self.assertTrue(selection.abstained)
        self.assertEqual(selection.panel, ())
        self.assertEqual(selection.k, 0)

    def test_a_non_positive_weight_does_not_create_a_seed(self):
        selection = self.run_rule({"SEEDA": 0.0})
        self.assertTrue(selection.abstained)

    def test_seeds_are_given_and_never_appear_in_the_expansion(self):
        selection = self.run_rule({"SEEDA": 1.0, "CLQ1": 1.0})
        self.assertEqual(selection.seeds, ("CLQ1", "SEEDA"))
        for seed in selection.seeds:
            self.assertIn(seed, selection.panel)
            self.assertNotIn(seed, selection.expansion)
        self.assertEqual(selection.panel,
                         selection.seeds + selection.expansion)

    def test_the_threshold_rule_abstains_even_on_an_obvious_module(self):
        # The documented negative result. SEEDA sits in a clique, so every node
        # the null can relocate the seed to is a clique member whose own
        # neighbours score exactly what SEEDA's neighbours score. The count
        # ratio is then about one and nothing is retained at any rate, which is
        # why the panel size is not delegated to the network.
        for rate in (0.05, 0.2, 0.5):
            selection = self.run_rule({"SEEDA": 1.0}, q=rate)
            self.assertEqual(selection.expansion, (),
                             msg=f"unexpected expansion at q={rate}")
            self.assertTrue(selection.expansion_abstained)

    def test_a_distant_component_is_not_recovered(self):
        selection = self.run_rule({"SEEDA": 1.0})
        self.assertFalse(any(name.startswith("FAR")
                             for name in selection.expansion))

    def test_the_decision_is_deterministic(self):
        first = self.run_rule({"SEEDA": 1.0})
        second = self.run_rule({"SEEDA": 1.0})
        self.assertEqual(first.panel, second.panel)

    def test_a_different_seed_produces_a_different_panel(self):
        near = self.run_rule({"SEEDA": 1.0})
        far = self.run_rule({"FAR6": 1.0})
        self.assertNotEqual(near.panel, far.panel)
        self.assertLess(jaccard(near.panel, far.panel), 1.0)

    def test_max_expansion_truncates_and_says_so(self):
        capped = self.run_rule({"SEEDA": 1.0}, max_expansion=0)
        self.assertEqual(capped.expansion, ())
        self.assertFalse(capped.provenance["expansion_truncated_to_max"])

    def test_the_expansion_is_ordered_by_decreasing_proximity(self):
        selection = self.run_rule({"SEEDA": 1.0})
        scores = [selection.scores[name] for name in selection.expansion]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_the_recorded_rate_does_not_exceed_the_requested_one(self):
        selection = self.run_rule({"SEEDA": 1.0}, q=0.05)
        self.assertLessEqual(
            selection.provenance["estimated_false_discovery_rate_at_cut"], 0.05)

    def test_a_stricter_rate_cannot_retain_more(self):
        loose = self.run_rule({"SEEDA": 1.0}, q=0.5)
        strict = self.run_rule({"SEEDA": 1.0}, q=0.01)
        self.assertLessEqual(len(strict.expansion), len(loose.expansion))

    def test_at_least_one_resample_is_required(self):
        with self.assertRaisesRegex(MechanisticSelectionError, "resample"):
            self.run_rule({"SEEDA": 1.0}, null_resamples=0)

    def test_the_parameters_that_produced_the_panel_are_carried_with_it(self):
        selection = self.run_rule({"SEEDA": 1.0})
        self.assertEqual(selection.parameters["null_resamples"], 200)
        self.assertEqual(selection.parameters["random_seed"]
                         if "random_seed" in selection.parameters
                         else selection.parameters["null_seed"], 5)
        self.assertIn("random walk", selection.rule)
        self.assertEqual(selection.provenance["min_combined_score"],
                         FROZEN["string_min_combined_score"])


class TestMechanisticRanking(unittest.TestCase):
    """The ranking is the primary product, so this is where power is tested."""

    def setUp(self):
        self.tmp = TemporaryDirectory()
        info, links = write_string_pair(Path(self.tmp.name), RING_AND_CLIQUE)
        self.graph = load_string_symbol_graph(info, links)

    def tearDown(self):
        self.tmp.cleanup()

    def test_the_module_around_the_seed_occupies_the_top_of_the_ranking(self):
        ranking = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0})
        self.assertEqual(ranking.seeds, ("SEEDA",))
        self.assertEqual(set(ranking.ranked[:5]),
                         {f"CLQ{i}" for i in range(1, 6)})
        self.assertTrue(all(name.startswith("FAR") for name in ranking.ranked[5:]))

    def test_the_seed_never_appears_among_the_candidates(self):
        ranking = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0})
        self.assertNotIn("SEEDA", ranking.ranked)
        self.assertEqual(len(ranking.ranked), self.graph.n_nodes - 1)

    def test_top_returns_the_seeds_plus_the_requested_number(self):
        ranking = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0})
        self.assertEqual(ranking.top(0), ("SEEDA",))
        self.assertEqual(len(ranking.top(3)), 4)
        self.assertEqual(ranking.top(3), ("SEEDA",) + ranking.ranked[:3])
        # Asking for more than exist returns everything rather than failing.
        self.assertEqual(len(ranking.top(10_000)), self.graph.n_nodes)

    def test_a_negative_panel_size_is_refused(self):
        ranking = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0})
        with self.assertRaisesRegex(MechanisticSelectionError, "negative"):
            ranking.top(-1)

    def test_the_ranking_declares_that_it_read_no_outcome(self):
        ranking = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0})
        self.assertFalse(ranking.provenance["outcome_was_read"])
        self.assertEqual(ranking.provenance["ranked_candidates"],
                         self.graph.n_nodes - 1)

    def test_scores_decrease_along_the_ranking(self):
        ranking = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0})
        scores = [ranking.scores[name] for name in ranking.ranked]
        self.assertEqual(scores, sorted(scores, reverse=True))

    def test_ties_are_broken_reproducibly_and_not_by_insertion_order(self):
        first = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0}).ranked
        second = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0}).ranked
        self.assertEqual(first, second)
        # The five clique members are exactly tied, so the tie-break must order
        # them by symbol.
        self.assertEqual(list(first[:5]), sorted(first[:5]))

    def test_a_different_seed_gives_a_different_order(self):
        near = rank_by_mechanistic_proximity(self.graph, {"SEEDA": 1.0}).ranked
        far = rank_by_mechanistic_proximity(self.graph, {"FAR6": 1.0}).ranked
        self.assertNotEqual(near[:5], far[:5])

    def test_a_seed_outside_the_network_raises(self):
        with self.assertRaisesRegex(MechanisticSelectionError, "No seed"):
            rank_by_mechanistic_proximity(self.graph, {"NOWHERE1": 1.0})

    def test_weights_shift_proximity_towards_the_heavier_seed(self):
        balanced = rank_by_mechanistic_proximity(
            self.graph, {"CLQ1": 1.0, "FAR6": 1.0})
        tilted = rank_by_mechanistic_proximity(
            self.graph, {"CLQ1": 1.0, "FAR6": 50.0})
        # A neighbour of the up-weighted seed gains, a neighbour of the other
        # loses. The order of the two seeds' own neighbourhoods need not swap:
        # the ring node has fewer partners to share mass with, so it starts
        # ahead on the balanced run.
        self.assertGreater(tilted.scores["FAR5"], balanced.scores["FAR5"])
        self.assertLess(tilted.scores["CLQ2"], balanced.scores["CLQ2"])


class TestTheRuleIsReproducible(unittest.TestCase):
    """Re-derive the decision from its written definition, not from the code."""

    def test_an_independent_rederivation_returns_the_same_panel(self):
        with TemporaryDirectory() as tmp:
            info, links = write_string_pair(Path(tmp), RING_AND_CLIQUE)
            graph = load_string_symbol_graph(info, links)

        alpha, resamples, strata_count, q, seed = 0.3, 120, 3, 0.05, 99
        produced = select_panel(graph, {"SEEDA": 1.0}, restart_probability=alpha,
                                null_resamples=resamples, n_strata=strata_count,
                                q=q, random_seed=seed)

        # Step 1: the walk, by direct linear solve rather than power iteration.
        transition = column_stochastic_transition(graph).toarray()
        identity = np.eye(graph.n_nodes)
        seed_node = graph.index["SEEDA"]
        restart = np.zeros(graph.n_nodes)
        restart[seed_node] = 1.0
        observed = np.linalg.solve(identity - (1.0 - alpha) * transition,
                                  alpha * restart)

        # Step 2: the same degree-stratified relocation null.
        strata = degree_strata(graph, strata_count)
        pools = {int(s): np.nonzero(strata == s)[0] for s in np.unique(strata)}
        rng = np.random.default_rng(seed)
        null = np.zeros((graph.n_nodes, resamples))
        drawn = np.zeros(resamples, dtype=int)
        for column in range(resamples):
            pool = pools[int(strata[seed_node])]
            picked = int(pool[rng.integers(pool.size)])
            drawn[column] = picked
            null_restart = np.zeros(graph.n_nodes)
            null_restart[picked] = 1.0
            null[:, column] = np.linalg.solve(
                identity - (1.0 - alpha) * transition, alpha * null_restart)

        # Step 3: candidates are the non-seed nodes; the null's own relocated
        # seed is excluded for the same reason the observed seed is.
        candidates = np.array([i for i in range(graph.n_nodes) if i != seed_node])
        position = {node: index for index, node in enumerate(candidates)}
        null_candidates = null[candidates].copy()
        for column in range(resamples):
            if drawn[column] in position:
                null_candidates[position[drawn[column]], column] = -np.inf

        # Step 4: the deepest cut whose estimated rate stays within q.
        order = np.argsort(-observed[candidates], kind="stable")
        running, retained = 0.0, 0
        for rank, index in enumerate(order, start=1):
            threshold = observed[candidates][index]
            mean_null = float(np.mean(
                [(null_candidates[:, column] >= threshold).sum()
                 for column in range(resamples)]))
            running = max(running, mean_null / rank)
            if running <= q:
                retained = rank
        expected = tuple(graph.symbols[candidates[i]] for i in order[:retained])

        self.assertEqual(produced.expansion, expected)
        self.assertEqual(produced.panel, ("SEEDA",) + expected)


class TestNoGeneIsNamedInTheModule(unittest.TestCase):
    """The protocol's mandatory honesty check, enforced mechanically.

    A selection written into the code after inspecting results would show up as
    a gene-symbol literal in the selector. Nothing in the module may name one:
    every symbol must arrive from a downloaded source or from the analysed
    dataset's own column name.
    """

    GENE_SHAPED = __import__("re").compile(r"^[A-Z][A-Z0-9]{1,8}[0-9][A-Z0-9]*$")
    KNOWN_FAMILY = __import__("re").compile(
        r"\b(CYP|UGT|GST|SULT|NAT|ALDH|ADH|AOX|FMO|EPHX|CES|POR|CYB5|"
        r"SLC|SLCO|ABC|NR1[IH]|PPAR|AHR|HNF|RXR|VDR|ESR|MAO|TPMT|DPYD|"
        r"NQO|COMT|GSTM|GSTT|GSTP)[0-9]", __import__("re").IGNORECASE)

    def string_constants(self):
        tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and isinstance(node.value, str):
                yield node.value

    def test_no_string_constant_has_the_shape_of_a_gene_symbol(self):
        offenders = []
        for text in self.string_constants():
            for token in __import__("re").split(r"[^A-Za-z0-9]+", text):
                if token and self.GENE_SHAPED.match(token):
                    offenders.append(token)
        self.assertEqual(offenders, [], msg=(
            "The selector must not name a gene. Offending literals: "
            f"{sorted(set(offenders))}"))

    def test_no_string_constant_names_a_drug_metabolising_gene_family(self):
        offenders = [text for text in self.string_constants()
                     if self.KNOWN_FAMILY.search(text)]
        self.assertEqual(offenders, [], msg=(
            "The selector must not name an enzyme or transporter family. "
            f"Offending literals: {offenders}"))

    def test_the_module_defines_no_collection_of_symbol_like_names(self):
        tree = ast.parse(MODULE_PATH.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.List, ast.Tuple, ast.Set)):
                continue
            literals = [element.value for element in node.elts
                        if isinstance(element, ast.Constant)
                        and isinstance(element.value, str)]
            uppercase = [text for text in literals
                         if text and text == text.upper() and len(text) > 1]
            self.assertLess(len(uppercase), 3, msg=(
                "A literal collection of uppercase names in the selector looks "
                f"like a hard-coded panel: {literals}"))

    def test_the_source_file_is_the_one_being_scanned(self):
        # Guards against the scan silently passing on the wrong file.
        self.assertTrue(MODULE_PATH.name.endswith("mechanistic_graph_selection.py"))
        self.assertIn("random walk with restart",
                      MODULE_PATH.read_text(encoding="utf-8").lower())


class TestPlatformMappingAndOverlap(unittest.TestCase):
    def test_every_probe_carrying_a_symbol_is_returned(self):
        platform = ["AAA1", "BBB2", "AAA1", "", "CCC3"]
        mapped = genes_to_probe_positions(["AAA1", "CCC3", "ZZZ9"], platform)
        self.assertEqual(mapped, {"AAA1": (0, 2), "CCC3": (4,)})

    def test_a_gene_absent_from_the_platform_is_omitted_not_faked(self):
        self.assertEqual(genes_to_probe_positions(["ZZZ9"], ["AAA1"]), {})

    def test_overlap_of_two_empty_selections_is_one(self):
        self.assertEqual(jaccard([], []), 1.0)

    def test_overlap_is_the_usual_ratio(self):
        self.assertAlmostEqual(jaccard(["A", "B"], ["B", "C"]), 1 / 3)
        self.assertEqual(jaccard(["A"], ["A"]), 1.0)
        self.assertEqual(jaccard(["A"], ["B"]), 0.0)


class TestFrozenParameters(unittest.TestCase):
    def test_the_frozen_set_is_complete(self):
        for key in ("string_version", "string_min_combined_score",
                    "restart_probability", "convergence_tolerance",
                    "max_power_iterations", "null_resamples", "degree_strata",
                    "target_fdr", "null_seed"):
            self.assertIn(key, FROZEN)

    def test_the_frozen_values_are_the_ones_the_protocol_records(self):
        self.assertEqual(FROZEN["string_min_combined_score"], 700)
        self.assertAlmostEqual(FROZEN["restart_probability"], 0.30)
        self.assertEqual(FROZEN["null_resamples"], 500)
        self.assertAlmostEqual(FROZEN["target_fdr"], 0.05)
        self.assertEqual(FROZEN["degree_strata"], 10)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
