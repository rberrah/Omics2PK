"""Executable mechanism-guided variable selection: knowledge graph to decision.

Why this module exists
----------------------
Until now the project's mechanistic selection was a *prioritisation* step: a
network ranking was computed, but the final retain/exclude decision was written
into the code as a literal gene list. A result obtained that way demonstrates
the value of a human reading biology, not the value of an algorithm. The
validation protocol therefore forbids claiming H2 on an autonomy basis until an
executable rule reproduces the retained set mechanically.

This module is that rule. It takes

  * a public protein-protein functional network (STRING),
  * a seed set obtained automatically from a public resource or from the
    dataset's own endpoint label, and
  * frozen numeric parameters,

and returns an ordering of candidate genes, together with the evidence that
produced it. It never sees an outcome variable, so it cannot be tuned against
the thing it is meant to predict. That is a deliberate restriction with a
visible cost: an outcome-blind ranking cannot recognise that a gene it ranks
highly carries no measurable signal, and the H2/H4 analyses must report that
cost rather than hide it.

What the network can and cannot decide
--------------------------------------
Building this module established a negative result that shapes the whole H2
design, and it was established from the network alone, before any endpoint was
predicted. A propagation score cannot be given a calibrated significance
threshold: every null that preserves the degree sequence also preserves how
many nodes reach a given proximity, so a count-based or score-based criterion
retains almost nothing at any parameter setting. The argument and the measured
behaviour are recorded in `permutation_fdr_cutoff`.

The consequence is that this module supplies a *ranking*, not a panel size. The
number of genes to keep is chosen by the analysis inside training folds only,
exactly as a data-driven competitor chooses its own. That keeps the comparison
fair -- every arm is allowed to use the outcome in the same way, and the only
difference between arms is where the ordering came from -- and it keeps the
test folds untouched by any selection step.

Two properties are enforced by tests rather than by promise:

1. **No gene is named in this file.** `tests/test_mechanistic_graph_selection.py`
   scans the module source for gene-symbol literals and fails if any appears.
   Every symbol that reaches a panel comes from a downloaded source file or
   from the analysed dataset's own column name.
2. **The decision is reproducible from the documented rule.** The test
   re-implements the retain rule independently and asserts the same output.

Data handling
-------------
The source files are downloaded separately and are not redistributed here. Only
gene symbols, scores and counts cross this boundary; no participant-level
value, and no source archive, enters the repository.

Frozen parameters live in `FROZEN`. Changing one is an amendment, not a tweak.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import gzip
import re
import zipfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
from scipy.sparse import csr_matrix

# --------------------------------------------------------------------------
# Frozen analysis parameters. An amendment must be recorded in the registry.
# --------------------------------------------------------------------------
FROZEN = {
    "string_version": "v12.0",
    "string_min_combined_score": 700,
    "restart_probability": 0.30,
    "convergence_tolerance": 1e-11,
    "max_power_iterations": 400,
    "null_resamples": 500,
    "degree_strata": 10,
    "target_fdr": 0.05,
    "null_seed": 20261108,
}

STRING_INFO_HEADER = "#string_protein_id\tpreferred_name\tprotein_size\tannotation"
STRING_LINKS_HEADER = "protein1 protein2 combined_score"
STRING_MAX_SCORE = 1000.0

# The dataset's own column-naming convention, not a biological fact.
ENDPOINT_LABEL_SUFFIX = "_activity"
_TRAILING_UPPERCASE = re.compile(r"[A-Z]$")


class MechanisticSelectionError(ValueError):
    """Raised when a source file, a seed set or a parameter is unusable."""


# --------------------------------------------------------------------------
# Graph
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class SymbolGraph:
    """An undirected, weighted, symbol-level functional network.

    `indptr`/`indices`/`weights` are a CSR adjacency over `symbols`. Weights are
    STRING combined scores rescaled to (0, 1]. The graph is symmetric by
    construction and carries no self-loops.
    """

    symbols: tuple[str, ...]
    index: Mapping[str, int]
    indptr: np.ndarray
    indices: np.ndarray
    weights: np.ndarray
    provenance: Mapping[str, object]

    @property
    def n_nodes(self) -> int:
        return len(self.symbols)

    @property
    def n_edges(self) -> int:
        return int(self.indices.size // 2)

    def weighted_degree(self) -> np.ndarray:
        out = np.zeros(self.n_nodes, dtype=np.float64)
        np.add.at(out, np.repeat(np.arange(self.n_nodes), np.diff(self.indptr)),
                  self.weights)
        return out


def _open_text(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def load_string_symbol_graph(info_path, links_path, *,
                             min_combined_score: int | None = None) -> SymbolGraph:
    """Build a symbol-level graph from the two STRING human download files.

    Several STRING protein identifiers can share a preferred gene symbol. Those
    are collapsed onto one node and the strongest edge is kept, which is the
    conservative choice: it never invents a link that STRING does not assert.
    """
    threshold = (FROZEN["string_min_combined_score"] if min_combined_score is None
                 else int(min_combined_score))
    if not 0 < threshold <= STRING_MAX_SCORE:
        raise MechanisticSelectionError(
            "STRING combined score threshold must lie in (0, 1000]")

    protein_to_symbol: dict[str, str] = {}
    with _open_text(info_path) as handle:
        header = next(handle).rstrip("\n")
        if header != STRING_INFO_HEADER:
            raise MechanisticSelectionError(
                "STRING protein.info header changed; the loader must be re-verified")
        for line in handle:
            parts = line.split("\t", 2)
            if len(parts) < 2:
                continue
            symbol = parts[1].strip()
            if symbol:
                protein_to_symbol[parts[0]] = symbol

    if not protein_to_symbol:
        raise MechanisticSelectionError("STRING protein.info carried no symbols")

    best: dict[tuple[int, int], float] = {}
    symbol_to_node: dict[str, int] = {}
    symbols: list[str] = []
    kept_lines = 0
    with _open_text(links_path) as handle:
        header = next(handle).rstrip("\n")
        if header != STRING_LINKS_HEADER:
            raise MechanisticSelectionError(
                "STRING protein.links header changed; the loader must be re-verified")
        for line in handle:
            left, right, raw = line.split(" ")
            score = int(raw)
            if score < threshold:
                continue
            a = protein_to_symbol.get(left)
            b = protein_to_symbol.get(right)
            if a is None or b is None or a == b:
                continue
            kept_lines += 1
            ia = symbol_to_node.get(a)
            if ia is None:
                ia = symbol_to_node[a] = len(symbols)
                symbols.append(a)
            ib = symbol_to_node.get(b)
            if ib is None:
                ib = symbol_to_node[b] = len(symbols)
                symbols.append(b)
            key = (ia, ib) if ia < ib else (ib, ia)
            weight = score / STRING_MAX_SCORE
            if weight > best.get(key, 0.0):
                best[key] = weight

    if not best:
        raise MechanisticSelectionError(
            "No STRING edge survived the confidence threshold")

    order = np.argsort(np.asarray(symbols, dtype=object), kind="stable")
    remap = np.empty(len(symbols), dtype=np.int64)
    remap[order] = np.arange(len(symbols))
    sorted_symbols = tuple(symbols[i] for i in order)

    rows = np.empty(2 * len(best), dtype=np.int64)
    cols = np.empty(2 * len(best), dtype=np.int64)
    vals = np.empty(2 * len(best), dtype=np.float64)
    for position, ((ia, ib), weight) in enumerate(best.items()):
        ra, rb = int(remap[ia]), int(remap[ib])
        rows[2 * position], cols[2 * position], vals[2 * position] = ra, rb, weight
        rows[2 * position + 1] = rb
        cols[2 * position + 1] = ra
        vals[2 * position + 1] = weight

    n = len(sorted_symbols)
    sort_key = np.lexsort((cols, rows))
    rows, cols, vals = rows[sort_key], cols[sort_key], vals[sort_key]
    counts = np.bincount(rows, minlength=n)
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])

    provenance = {
        "string_version": FROZEN["string_version"],
        "min_combined_score": threshold,
        "proteins_with_symbol": len(protein_to_symbol),
        "link_lines_above_threshold": kept_lines,
        "nodes": n,
        "undirected_edges": len(best),
        "symbol_collapse": "strongest edge kept when identifiers share a symbol",
    }
    return SymbolGraph(symbols=sorted_symbols,
                       index={s: i for i, s in enumerate(sorted_symbols)},
                       indptr=indptr, indices=cols, weights=vals,
                       provenance=provenance)


# --------------------------------------------------------------------------
# Random walk with restart
# --------------------------------------------------------------------------
def column_stochastic_transition(graph: SymbolGraph) -> csr_matrix:
    """Return W with W[i, j] = weight(i, j) / weighted degree of j.

    The walker leaving node j picks a neighbour in proportion to edge
    confidence, so each column sums to one and total probability is conserved.
    """
    degree = graph.weighted_degree()
    safe = np.where(degree > 0.0, degree, 1.0)
    data = graph.weights / safe[graph.indices]
    return csr_matrix((data, graph.indices, graph.indptr),
                      shape=(graph.n_nodes, graph.n_nodes))


def random_walk_with_restart(graph: SymbolGraph, restart: np.ndarray, *,
                             restart_probability: float | None = None,
                             tolerance: float | None = None,
                             max_iterations: int | None = None,
                             transition: csr_matrix | None = None) -> np.ndarray:
    """Stationary distribution of a random walk that restarts on `restart`.

    `restart` is a non-negative vector over graph nodes, or a matrix whose
    columns are independent restart distributions evaluated in one pass. Power
    iteration on p <- (1 - a) W p + a r, with W column-stochastic.
    """
    alpha = (FROZEN["restart_probability"] if restart_probability is None
             else float(restart_probability))
    if not 0.0 < alpha <= 1.0:
        raise MechanisticSelectionError("Restart probability must lie in (0, 1]")
    tol = FROZEN["convergence_tolerance"] if tolerance is None else float(tolerance)
    iterations = (FROZEN["max_power_iterations"] if max_iterations is None
                  else int(max_iterations))

    restart = np.asarray(restart, dtype=np.float64)
    vector_input = restart.ndim == 1
    if vector_input:
        restart = restart[:, None]
    if restart.ndim != 2 or restart.shape[0] != graph.n_nodes:
        raise MechanisticSelectionError("Restart vector must align with graph nodes")
    if np.any(restart < 0.0):
        raise MechanisticSelectionError("Restart vector must be non-negative")
    totals = restart.sum(axis=0)
    if not np.all(np.isfinite(totals)) or np.any(totals <= 0.0):
        raise MechanisticSelectionError("Restart vector must carry positive mass")
    restart = restart / totals

    matrix = column_stochastic_transition(graph) if transition is None else transition
    p = restart.copy()
    for _ in range(iterations):
        nxt = (1.0 - alpha) * (matrix @ p) + alpha * restart
        # Any mass lost on a node the walker cannot leave returns to the seeds.
        deficit = 1.0 - nxt.sum(axis=0)
        positive = deficit > 0.0
        if np.any(positive):
            nxt = nxt + restart * np.where(positive, deficit, 0.0)
        if np.abs(nxt - p).sum(axis=0).max() < tol:
            p = nxt
            break
        p = nxt
    return p[:, 0] if vector_input else p


# --------------------------------------------------------------------------
# Seed construction. Each builder reads a public source or the dataset label.
# --------------------------------------------------------------------------
def seeds_from_endpoint_label(label: str, graph: SymbolGraph) -> dict[str, float]:
    """Resolve a dataset endpoint column name to a seed, using the graph as the
    dictionary.

    The rule, frozen: drop the dataset's documented suffix; if the remainder is
    not a node of the network, drop trailing single uppercase characters one at
    a time until it is. No gene list is consulted, and no family is special-
    cased. Where a dataset distinguishes two assays of one enzyme by appending
    a letter to its symbol, the shortening step recovers the symbol because the
    lengthened form is absent from the network.
    """
    text = str(label).strip()
    if not text.endswith(ENDPOINT_LABEL_SUFFIX):
        raise MechanisticSelectionError(
            f"Endpoint label must end with {ENDPOINT_LABEL_SUFFIX!r}")
    stem = text[: -len(ENDPOINT_LABEL_SUFFIX)]
    while stem:
        if stem in graph.index:
            return {stem: 1.0}
        if not _TRAILING_UPPERCASE.search(stem):
            break
        stem = stem[:-1]
    raise MechanisticSelectionError(
        "Endpoint label does not resolve to any network node")


def seeds_from_reactome(gmt_zip_path, pathway_pattern: str, graph: SymbolGraph,
                        *, member: str | None = None) -> dict[str, float]:
    """Seed from every Reactome pathway whose name matches `pathway_pattern`.

    The pattern is a case-insensitive regular expression over pathway names. All
    members of all matching pathways become seeds with equal weight, so the seed
    set carries no endpoint-specific information whatsoever.
    """
    path = Path(gmt_zip_path)
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        chosen = member or (names[0] if len(names) == 1 else None)
        if chosen is None or chosen not in names:
            raise MechanisticSelectionError(
                "Reactome archive member is ambiguous; name it explicitly")
        raw = archive.read(chosen).decode("utf-8", "replace")

    pattern = re.compile(pathway_pattern, re.IGNORECASE)
    seeds: dict[str, float] = {}
    matched: list[str] = []
    for row in raw.splitlines():
        if not row.strip():
            continue
        fields = row.split("\t")
        if len(fields) < 3 or not pattern.search(fields[0]):
            continue
        matched.append(fields[0])
        for symbol in fields[2:]:
            symbol = symbol.strip()
            if symbol in graph.index:
                seeds[symbol] = 1.0
    if not matched:
        raise MechanisticSelectionError(
            "No Reactome pathway name matched the pattern")
    if not seeds:
        raise MechanisticSelectionError(
            "Matching Reactome pathways contain no network node")
    return seeds


def seeds_from_ctd_substrate(ctd_path, chemical_name: str, graph: SymbolGraph, *,
                             organism: str = "Homo sapiens",
                             action_token: str = "metabolic processing"
                             ) -> dict[str, float]:
    """Seed from a probe substrate, via curated chemical-gene interactions.

    A gene becomes a seed when the source records an interaction between the
    named chemical and that gene carrying `action_token`. The weight is the
    number of distinct supporting references, so a single-report association
    does not outweigh a repeatedly replicated one.

    This is a literature lookup and must be described as one: the module does
    not discover that an enzyme metabolises a drug, it reads that a public
    database asserts it. What the rest of this module adds is the expansion
    around that assertion.
    """
    wanted = str(chemical_name).strip().lower()
    if not wanted:
        raise MechanisticSelectionError("A chemical name is required")
    token = str(action_token).strip().lower()
    references: dict[str, set[str]] = {}
    columns: list[str] | None = None
    with _open_text(ctd_path) as handle:
        for line in handle:
            if line.startswith("#"):
                if "ChemicalName" in line:
                    columns = line.lstrip("# ").rstrip("\n").split("\t")
                continue
            if columns is None:
                raise MechanisticSelectionError(
                    "Chemical-gene source carried no column header")
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(columns):
                continue
            record = dict(zip(columns, fields))
            if record["ChemicalName"].strip().lower() != wanted:
                continue
            if record["Organism"].strip() != organism:
                continue
            if token not in record["InteractionActions"].lower():
                continue
            symbol = record["GeneSymbol"].strip()
            if symbol not in graph.index:
                continue
            pmids = {p for p in record["PubMedIDs"].split("|") if p}
            references.setdefault(symbol, set()).update(pmids or {record["Interaction"]})
    if not references:
        raise MechanisticSelectionError(
            "No curated interaction matched the chemical, organism and action")
    return {symbol: float(len(pmids)) for symbol, pmids in references.items()}


# --------------------------------------------------------------------------
# The decision
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class MechanisticRanking:
    """An ordered list of candidate genes, produced without any outcome value.

    This is the module's primary product. It says what the network considers
    close to the seeds, in order, and nothing about how many to keep. Choosing
    how many is a separate decision that belongs to the analysis, because -- see
    `permutation_fdr_cutoff` -- the network cannot make it.
    """

    seeds: tuple[str, ...]
    ranked: tuple[str, ...]
    scores: Mapping[str, float]
    provenance: Mapping[str, object]

    def top(self, k: int) -> tuple[str, ...]:
        """The seeds plus the `k` closest non-seed genes."""
        if int(k) < 0:
            raise MechanisticSelectionError("A panel size cannot be negative")
        return self.seeds + self.ranked[: int(k)]


def rank_by_mechanistic_proximity(graph: SymbolGraph, seeds: Mapping[str, float],
                                  *, restart_probability: float | None = None
                                  ) -> MechanisticRanking:
    """Order every non-seed gene by network proximity to the seed set.

    No outcome variable is read, so the ordering is identical in every
    cross-validation fold and cannot have been tuned against the quantity it is
    later used to predict. The cost of that property is stated in the module
    docstring: an outcome-blind ordering cannot know that a gene it ranks highly
    carries no measurable signal.
    """
    alpha = (FROZEN["restart_probability"] if restart_probability is None
             else float(restart_probability))
    present = {s: float(w) for s, w in seeds.items()
               if s in graph.index and float(w) > 0.0}
    if not present:
        raise MechanisticSelectionError("No seed is a node of the network")
    seed_names = tuple(sorted(present))
    seed_nodes = np.array([graph.index[s] for s in seed_names], dtype=np.int64)
    restart = np.zeros(graph.n_nodes, dtype=np.float64)
    restart[seed_nodes] = [present[s] for s in seed_names]
    observed = random_walk_with_restart(graph, restart, restart_probability=alpha)
    candidate = np.ones(graph.n_nodes, dtype=bool)
    candidate[seed_nodes] = False
    nodes = np.nonzero(candidate)[0]
    # Ties are broken by symbol, which is arbitrary but reproducible; the
    # alternative, breaking them by anything data-derived, would leak.
    nodes = nodes[np.lexsort((np.array([graph.symbols[i] for i in nodes],
                                       dtype=object), -observed[nodes]))]
    provenance = dict(graph.provenance)
    provenance.update({"restart_probability": alpha,
                       "ranked_candidates": int(nodes.size),
                       "outcome_was_read": False})
    return MechanisticRanking(
        seeds=seed_names,
        ranked=tuple(graph.symbols[i] for i in nodes),
        scores={graph.symbols[i]: float(observed[i])
                for i in np.concatenate([seed_nodes, nodes])},
        provenance=provenance,
    )


def permutation_fdr_cutoff(observed: np.ndarray, null: np.ndarray, q: float
                           ) -> tuple[int, np.ndarray]:
    """Largest retained set whose estimated false discovery rate is at most `q`.

    `observed` holds the candidate scores; `null` is candidates by resamples.
    At a score threshold t the estimated rate is the average number of null
    candidates reaching t, divided by the number of observed candidates
    reaching it -- the permutation estimator of Tusher, Tibshirani and Chu
    (2001), transposed to network proximity.

    Why not an empirical p-value with Benjamini-Hochberg: a p-value from B
    resamples cannot fall below 1/(B+1), while a Benjamini-Hochberg threshold
    over sixteen thousand hypotheses is of order 1e-6. That combination can
    never reject anything, whatever the data. The count ratio has no such
    resolution floor because its precision comes from the number of nodes, not
    from the number of resamples. This was found by running the first version
    of the rule, which abstained on every endpoint; it is recorded here so the
    change is not mistaken for tuning.

    A second and deeper limitation survives that repair, and the analysis must
    report it rather than work around it. Any null that preserves the degree
    sequence also preserves, approximately, *how many* nodes reach a given
    proximity to a seed: relocating a seed to a comparably connected node gives
    its new neighbours about the scores the real seed's neighbours have. The
    count ratio is therefore close to one wherever the walk is informative, and
    this criterion retains almost nothing however the parameters are set. That
    is not conservatism to be tuned away; it says that network proximity
    carries information about *which* genes are partners and almost none about
    whether being a partner is unusual. Significance of a propagation score
    against a topological null is not identifiable, so the decision of how many
    genes to keep cannot be delegated to the network. This function is kept as
    the reported diagnostic of that fact, not as the selection rule;
    `rank_by_mechanistic_proximity` is the primary product.

    Returns the number of retained candidates and a per-threshold audit table
    with columns: threshold, observed count, mean null count, estimated rate.
    """
    observed = np.asarray(observed, dtype=np.float64)
    null = np.asarray(null, dtype=np.float64)
    if observed.ndim != 1 or null.ndim != 2 or null.shape[0] != observed.size:
        raise MechanisticSelectionError(
            "Null scores must be candidates by resamples and align with observed")
    if not 0.0 < float(q) < 1.0:
        raise MechanisticSelectionError("The false discovery rate must lie in (0, 1)")
    if observed.size == 0:
        return 0, np.zeros((0, 4), dtype=np.float64)

    order = np.argsort(-observed, kind="stable")
    thresholds = observed[order]
    hits = np.arange(1, observed.size + 1, dtype=np.float64)
    # Mean number of null candidates at or above each threshold, per resample.
    flat = np.sort(null.reshape(-1))
    reached = flat.size - np.searchsorted(flat, thresholds, side="left")
    mean_null = reached / float(null.shape[1])
    rate = np.where(hits > 0.0, mean_null / hits, 1.0)
    # Conservative monotonisation: a cut inherits the worst rate seen at any
    # shallower cut, so a locally favourable fluctuation cannot open the gate.
    guarded = np.maximum.accumulate(rate)
    passing = np.nonzero(guarded <= float(q))[0]
    retained = int(passing[-1] + 1) if passing.size else 0
    audit = np.column_stack([thresholds, hits, mean_null, rate])
    return retained, audit


def degree_strata(graph: SymbolGraph, n_strata: int | None = None) -> np.ndarray:
    """Assign each node to a stratum of comparable connectivity.

    Strata are quantiles of log10(1 + weighted degree). A null seed drawn from a
    node's own stratum has comparable reach, so a gene cannot appear significant
    merely for being a hub.
    """
    count = FROZEN["degree_strata"] if n_strata is None else int(n_strata)
    if count < 1:
        raise MechanisticSelectionError("At least one degree stratum is required")
    scale = np.log10(1.0 + graph.weighted_degree())
    edges = np.quantile(scale, np.linspace(0.0, 1.0, count + 1)[1:-1])
    return np.searchsorted(edges, scale, side="right").astype(np.int64)


@dataclass
class Selection:
    """What the rule decided, and the evidence it decided on."""

    seeds: tuple[str, ...]
    expansion: tuple[str, ...]
    panel: tuple[str, ...]
    rate_at_cut: Mapping[str, float]
    scores: Mapping[str, float]
    abstained: bool
    expansion_abstained: bool
    rule: str
    parameters: Mapping[str, object]
    provenance: Mapping[str, object] = field(default_factory=dict)

    @property
    def k(self) -> int:
        return len(self.panel)


def select_panel(graph: SymbolGraph, seeds: Mapping[str, float], *,
                 restart_probability: float | None = None,
                 null_resamples: int | None = None,
                 n_strata: int | None = None,
                 q: float | None = None,
                 random_seed: int | None = None,
                 max_expansion: int | None = None) -> Selection:
    """Rank the network around `seeds` and decide, automatically, what to retain.

    The rule, frozen before any run:

    1. Run a random walk with restart from the seed distribution.
    2. Build a null by relocating every seed to a random node of the same
       connectivity stratum and re-running the walk, `null_resamples` times.
    3. At every candidate score threshold, estimate the false discovery rate
       as the mean number of null nodes reaching it divided by the number of
       observed nodes reaching it.
    4. Retain the deepest cut whose estimated rate stays at or below `q`.

    Seeds are *given*, not selected, so they are excluded from the testing
    family and always present in the panel. The algorithm's own output is the
    expansion. If nothing survives, the expansion is empty and that is recorded
    as an abstention rather than silently filled.
    """
    alpha = (FROZEN["restart_probability"] if restart_probability is None
             else float(restart_probability))
    resamples = (FROZEN["null_resamples"] if null_resamples is None
                 else int(null_resamples))
    rate = FROZEN["target_fdr"] if q is None else float(q)
    base_seed = FROZEN["null_seed"] if random_seed is None else int(random_seed)
    if resamples < 1:
        raise MechanisticSelectionError("At least one null resample is required")

    present = {s: float(w) for s, w in seeds.items()
               if s in graph.index and float(w) > 0.0}
    if not present:
        return Selection(seeds=(), expansion=(), panel=(), rate_at_cut={}, scores={},
                         abstained=True, expansion_abstained=True,
                         rule="no seed is a node of the network",
                         parameters=dict(FROZEN), provenance=dict(graph.provenance))

    seed_names = tuple(sorted(present))
    seed_nodes = np.array([graph.index[s] for s in seed_names], dtype=np.int64)
    seed_weights = np.array([present[s] for s in seed_names], dtype=np.float64)

    restart = np.zeros(graph.n_nodes, dtype=np.float64)
    restart[seed_nodes] = seed_weights
    transition = column_stochastic_transition(graph)
    observed = random_walk_with_restart(graph, restart,
                                        restart_probability=alpha,
                                        transition=transition)

    strata = degree_strata(graph, n_strata)
    by_stratum = {int(s): np.nonzero(strata == s)[0] for s in np.unique(strata)}
    rng = np.random.default_rng(base_seed)
    null_restarts = np.zeros((graph.n_nodes, resamples), dtype=np.float64)
    drawn_nodes = np.empty((seed_nodes.size, resamples), dtype=np.int64)
    for draw in range(resamples):
        for position, (node, weight) in enumerate(zip(seed_nodes, seed_weights)):
            pool = by_stratum[int(strata[node])]
            picked = int(pool[rng.integers(pool.size)])
            drawn_nodes[position, draw] = picked
            null_restarts[picked, draw] += weight
    null_scores = random_walk_with_restart(graph, null_restarts,
                                           restart_probability=alpha,
                                           transition=transition)

    candidate = np.ones(graph.n_nodes, dtype=bool)
    candidate[seed_nodes] = False
    candidate_nodes = np.nonzero(candidate)[0]
    position_of = np.full(graph.n_nodes, -1, dtype=np.int64)
    position_of[candidate_nodes] = np.arange(candidate_nodes.size)

    # The observed seeds are given rather than discovered, so the null's own
    # relocated seeds are excluded too. Without this the null would be charged
    # for nodes that the observed analysis never had to find.
    null_candidates = null_scores[candidate_nodes].copy()
    for draw in range(resamples):
        rows = position_of[drawn_nodes[:, draw]]
        null_candidates[rows[rows >= 0], draw] = -np.inf

    retained_count, audit = permutation_fdr_cutoff(
        observed[candidate_nodes], null_candidates, rate)
    ranked = candidate_nodes[np.argsort(-observed[candidate_nodes], kind="stable")]
    retained_nodes = ranked[:retained_count]
    truncated = False
    if max_expansion is not None and retained_nodes.size > int(max_expansion):
        retained_nodes = retained_nodes[: int(max_expansion)]
        truncated = True

    expansion = tuple(graph.symbols[i] for i in retained_nodes)
    achieved = float(audit[retained_count - 1, 3]) if retained_count else 0.0
    provenance = dict(graph.provenance)
    provenance.update({
        "candidates_tested": int(candidate_nodes.size),
        "estimated_false_discovery_rate_at_cut": achieved,
        "expected_false_leads_at_cut": (
            float(audit[retained_count - 1, 2]) if retained_count else 0.0),
        "expansion_truncated_to_max": truncated,
    })
    return Selection(
        seeds=seed_names,
        expansion=expansion,
        panel=seed_names + expansion,
        rate_at_cut={graph.symbols[n]: float(audit[i, 3])
                     for i, n in enumerate(ranked[:retained_count])},
        scores={graph.symbols[n]: float(observed[n])
                for n in np.concatenate([seed_nodes, retained_nodes])},
        abstained=False,
        expansion_abstained=bool(expansion == ()),
        rule=("random walk with restart on the public functional network; "
              "degree-stratified seed-relocation null; permutation false "
              "discovery rate on node counts; seeds given, expansion selected"),
        parameters={
            "restart_probability": alpha,
            "null_resamples": resamples,
            "degree_strata": (FROZEN["degree_strata"] if n_strata is None
                              else int(n_strata)),
            "target_fdr": rate,
            "null_seed": base_seed,
            "max_expansion": max_expansion,
        },
        provenance=provenance,
    )


# --------------------------------------------------------------------------
# Mapping a panel onto an expression platform
# --------------------------------------------------------------------------
def genes_to_probe_positions(symbols: Iterable[str],
                             platform_symbols: Sequence[str]
                             ) -> dict[str, tuple[int, ...]]:
    """Map gene symbols to every probe position annotated with that symbol.

    Several probes can carry one symbol. Choosing among them needs data and is
    therefore left to the caller, which must do it inside a training fold.
    """
    platform = np.asarray(platform_symbols, dtype=str)
    out: dict[str, tuple[int, ...]] = {}
    for symbol in symbols:
        positions = np.nonzero(platform == symbol)[0]
        if positions.size:
            out[symbol] = tuple(int(p) for p in positions)
    return out


def jaccard(left: Iterable[str], right: Iterable[str]) -> float:
    """Overlap of two selections; 1.0 when both are empty."""
    a, b = set(left), set(right)
    if not a and not b:
        return 1.0
    union = a | b
    return len(a & b) / len(union) if union else 1.0
