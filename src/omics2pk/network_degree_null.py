"""Conditional top-SHAP ↔ pharmacology STRING bridge null on an observed graph.

Degree-stratified *node-label* permutations hold the observed graph, number of
genes in each role, and the number of SHAP labels at each EXACT node degree
constant. This is a graph-internal descriptive check, NOT an omics-universe
background test, independent pathway enrichment, or evidence of causal PK.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, product
from math import prod, comb, isfinite


class NetworkNullError(ValueError):
    pass


def degree_conditioned_bridge_null(genes, shap_genes, drug_prior_genes, edges,
                                   *, max_exact_assignments=100_000):
    """Exact conditioned-label reference distribution, abstaining if too large.

    No new biological genes or STRING connections are invented. Edges are
    STRING *functional* associations with score in [0.4, 1.0].
    """
    genes = tuple(genes)
    shap = set(shap_genes)
    prior = set(drug_prior_genes)
    if (len(genes) < 4 or len(set(genes)) != len(genes)
            or not all(isinstance(g, str) and g for g in genes)
            or not shap or not prior or shap & prior
            or shap | prior != set(genes)):
        raise NetworkNullError("Distinct exhaustive measured-source roles required")
    if (not isinstance(max_exact_assignments, int)
            or isinstance(max_exact_assignments, bool)
            or max_exact_assignments < 2):
        raise NetworkNullError("Positive exact enumeration limit required")

    normalized = []
    seen = set()
    degrees = Counter({name: 0 for name in genes})
    for edge in edges:
        if not isinstance(edge, dict):
            raise NetworkNullError("Original source edge mapping required")
        a, b = edge.get("gene_A"), edge.get("gene_B")
        w = edge.get("confidence_functional_association")
        if (a not in degrees or b not in degrees or a == b
                or isinstance(w, bool) or not isinstance(w, (float, int))
                or not isfinite(w) or not 0.4 <= w <= 1.0):
            raise NetworkNullError("Unqualified STRING edge")
        pair = tuple(sorted((a, b)))
        if pair in seen:
            raise NetworkNullError("Duplicated STRING graph edge")
        seen.add(pair)
        degrees[a] += 1
        degrees[b] += 1
        normalized.append((a, b, float(w)))

    strata = defaultdict(list)
    for gene in genes:
        strata[degrees[gene]].append(gene)
    values = [tuple(sorted(group)) for _, group in sorted(strata.items())]
    counts = [sum(gene in shap for gene in group) for group in values]
    n_assignments = prod(comb(len(group), take) for group, take in zip(values, counts))
    observed_bridges = sum((a in shap) != (b in shap) for a, b, _ in normalized)
    observed_weight = sum(w for a, b, w in normalized if (a in shap) != (b in shap))
    result = {
        "status": "EXACT_DEGREE_CONDITIONED_FIXED_17_GENE_GRAPH_DESCRIPTIVE",
        "n_observed_nodes": len(genes),
        "n_observed_string_functional_edges": len(normalized),
        "n_selected_shap_nodes": len(shap),
        "n_prespecified_prior_nodes": len(prior),
        "observed_shap_to_prior_edges": observed_bridges,
        "observed_shap_to_prior_edge_weight_sum": round(observed_weight, 6),
        "n_conditioned_role_assignments": n_assignments,
        "null_preserves_exact_node_degree_histogram_by_role": True,
        "null_uses_only_original_selected_genes_not_all_measured_omics": True,
        "no_external_background_degree_or_expression_matching_performed": True,
        "preselected_shap_and_drug_priors_invalidate_confirmatory_p_value": True,
        "direct_transport_or_pk_causality_demonstrated": False,
    }
    if n_assignments < 2:
        return {**result, "status": "ABSTAIN_NO_EXCHANGEABLE_DEGREE_MATCHED_LABELS"}
    if n_assignments > max_exact_assignments:
        return {**result, "status": "ABSTAIN_EXACT_ROLE_ENUMERATION_TOO_LARGE"}
    frequencies = Counter()
    weight_sums = []
    for choice in product(*(combinations(group, take) for group, take in zip(values, counts))):
        alternative = set(gene for subset in choice for gene in subset)
        bridges = sum((a in alternative) != (b in alternative) for a, b, _ in normalized)
        weights = sum(w for a, b, w in normalized
                      if (a in alternative) != (b in alternative))
        frequencies[bridges] += 1
        weight_sums.append(weights)
    if sum(frequencies.values()) != n_assignments:
        raise NetworkNullError("Incomplete exact null enumeration")
    return {
        **result,
        "expected_number_of_cross_role_edges_under_conditional_relabeling":
            round(sum(k * v for k, v in frequencies.items()) / n_assignments, 6),
        "fraction_of_conditional_assignments_with_at_least_observed_bridges":
            round(sum(v for k, v in frequencies.items() if k >= observed_bridges)
                  / n_assignments, 6),
        "expected_cross_role_string_score_sum_under_conditional_relabeling":
            round(sum(weight_sums) / n_assignments, 6),
        "fraction_of_conditional_assignments_with_at_least_observed_weight_sum":
            round(sum(w >= observed_weight - 1e-12 for w in weight_sums)
                  / n_assignments, 6),
        "conditional_cross_role_edge_count_histogram":
            {str(k): frequencies[k] for k in sorted(frequencies)},
        "interpretation": "DESCRIPTIVE selected-network topology only; not an independent pathway test",
    }
