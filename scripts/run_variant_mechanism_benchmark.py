#!/usr/bin/env python
"""First executable demonstrator: can the mechanism tell us which assay to trust?

Runs on public deep-mutational-scanning data in which the **same variant
library** was measured by two assays, one reporting protein abundance or
stability and one reporting activity. That pairing is what makes the question
answerable at all, and it is rare.

Data, all CC0, downloaded from the MaveDB API:

=============================  =====  =============================================
`urn:mavedb:00000055-a-1`      3100   NUDT15 protein stability, VAMP-seq
`urn:mavedb:00000055-b-1`      3100   NUDT15 **activity**, thiopurine cytotoxicity
`urn:mavedb:00000095-b-1`      6821   CYP2C9 abundance, VAMP-seq in HEK293T
`urn:mavedb:00000095-a-1`      6524   CYP2C9 **activity**, Click-seq in yeast
=============================  =====  =============================================

**NUDT15 is the primary gene.** Both of its assays were run on human cells
from one library. The CYP2C9 activity assay is in humanised yeast while its
abundance assay is in HEK293T, so a discordance there can be an assay-system
artefact, and every CYP2C9 result is reported as a sensitivity analysis with
that confound named.

Every split is grouped by protein **position**, and repeated with contiguous
**region** blocks, because several substitutions at one residue are not
independent observations.

Usage:
    python -I scripts/run_variant_mechanism_benchmark.py \
        --mavedb-dir <dir-with-the-four-csv-files> \
        --out reports/data/variant_mechanism_benchmark.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8")

from omics2pk import structural_annotation as SA  # noqa: E402
from omics2pk import variant_mechanism_omic_selection as VMS  # noqa: E402

#: Verified against the UniProt REST API on 2026-11-11. Embedded rather than
#: fetched so the benchmark runs offline and the annotation it used is
#: auditable; `--verify-uniprot` re-fetches and fails on any disagreement.
#: Note how THIN this is: CYP2C9 has exactly one annotated binding residue.
#: The mechanistic layer available from public annotation is sparse and
#: asymmetric between proteins, and that is a finding, not a configuration
#: detail.
ANNOTATION = {
    "NUDT15": {
        "uniprot": "Q9NV35",
        "length": 164,
        "functional_sites": (47, 49, 63, 67),
        "motifs": ((48, 69),),
        "motif_names": ("Nudix box",),
        "activity_urn": "urn:mavedb:00000055-b-1",
        "abundance_urn": "urn:mavedb:00000055-a-1",
        "role": "primary",
        "assay_confound": None,
    },
    # Pre-registered in OMICS2PK_BRAF_PREREGISTRATION_20261111.md and
    # declared here before any BRAF score was downloaded. UniProt P15056
    # annotates one active site (576, proton acceptor), eight point binding
    # sites in the cysteine-rich region (235-280, zinc), one binding range
    # (463-471, the ATP glycine-rich loop) and binding site 483.
    #
    # Mapping, to keep the earlier two genes' convention: point sites go to
    # functional_sites, a SHORT functional range goes to motifs. The protein
    # kinase domain 457-717 is a domain, not a motif, and is recorded for the
    # reader rather than used as a feature -- a near-constant indicator over
    # an assayed region that lies inside it would carry no information while
    # looking like knowledge.
    "BRAF": {
        "uniprot": "P15056",
        "length": 766,
        "functional_sites": (235, 248, 251, 261, 264, 269, 272, 280, 483, 576),
        "motifs": ((463, 471),),
        "motif_names": ("ATP glycine-rich loop",),
        "domains_not_used_as_features": ((155, 227), (457, 717)),
        "domain_names": ("RBD", "protein kinase"),
        "activity_urn": "urn:mavedb:00001315-d-1",
        "abundance_urn": "urn:mavedb:00001315-c-1",
        "role": "pre_registered_independent_test",
        "assay_confound": (
            "both assays are the C-terminal construct at baseline, so they "
            "share a construct and a condition. The pooled standardised "
            "score sets mix the N- and C-terminal constructs and are a "
            "declared sensitivity only."),
    },
    "CYP2C9": {
        "uniprot": "P11712",
        "length": 490,
        "functional_sites": (435,),
        "motifs": (),
        "motif_names": (),
        "activity_urn": "urn:mavedb:00000095-a-1",
        "abundance_urn": "urn:mavedb:00000095-b-1",
        "role": "sensitivity",
        "assay_confound": (
            "the abundance assay is VAMP-seq in HEK293T and the activity "
            "assay is Click-seq in humanised yeast. A discordance between "
            "them may be a difference of system rather than of biology, and "
            "this design cannot separate the two."),
    },
}

#: Two genes, two groupings, two conditional contrasts. Declared before the run.
PRIMARY_FAMILY_SIZE = 8


class RunError(RuntimeError):
    pass


def load_scores(path):
    """`hgvs_pro` -> score, keeping only variants this module can parse."""
    scores = {}
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if "hgvs_pro" not in (reader.fieldnames or []):
            raise RunError(f"{path.name} has no hgvs_pro column; got "
                           f"{reader.fieldnames}")
        for row in reader:
            parsed = VMS.parse_protein_variant(row.get("hgvs_pro"))
            raw = row.get("score")
            if parsed is None or raw in (None, "", "NA"):
                continue
            try:
                value = float(raw)
            except ValueError:
                continue
            if not np.isfinite(value):
                continue
            scores[row["hgvs_pro"].strip()] = (parsed, value)
    return scores


def spearman(left, right):
    left, right = np.asarray(left, float), np.asarray(right, float)
    def rank(values):
        order = values.argsort()
        out = np.empty(values.size)
        out[order] = np.arange(1, values.size + 1)
        return out
    a, b = rank(left) - (left.size + 1) / 2.0, rank(right) - (right.size + 1) / 2.0
    denominator = float(np.sqrt((a ** 2).sum() * (b ** 2).sum()))
    return float((a * b).sum() / denominator) if denominator else float("nan")


def verify_uniprot(gene, annotation):
    """Re-fetch the annotation and fail loudly on any disagreement."""
    import urllib.request
    url = (f"https://rest.uniprot.org/uniprotkb/{annotation['uniprot']}.json"
           "?fields=accession,sequence,ft_act_site,ft_binding,ft_domain,"
           "ft_motif,ft_site")
    with urllib.request.urlopen(url, timeout=60) as handle:
        payload = json.load(handle)
    length = len(payload["sequence"]["value"])
    sites, motifs = set(), set()
    for feature in payload.get("features", []):
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]
        if feature["type"] in ("Binding site", "Active site", "Site"):
            sites.update(range(start, end + 1))
        elif feature["type"] == "Motif":
            motifs.add((start, end))
    problems = []
    if length != annotation["length"]:
        problems.append(f"length {length} != {annotation['length']}")
    if sites != set(annotation["functional_sites"]):
        problems.append(f"sites {sorted(sites)} != "
                        f"{sorted(annotation['functional_sites'])}")
    if motifs != set(annotation["motifs"]):
        problems.append(f"motifs {sorted(motifs)} != "
                        f"{sorted(annotation['motifs'])}")
    if problems:
        raise RunError(f"{gene}: embedded annotation is stale: "
                       + "; ".join(problems))
    return {"verified": True, "length": length,
            "functional_sites": sorted(sites), "motifs": sorted(motifs)}


#: The structural columns that involve no threshold at all: two continuous
#: log-distances and the flag for a position the crystal does not resolve.
CONTINUOUS_COLUMNS = [0, 1, 5]

#: The cofactor half only, dropping every substrate column. For BRAF the
#: substrate is a protein chain, so `to_substrate` is an interface distance
#: rather than a pocket distance; this arm shows what the ruler does without
#: it. Columns: log distance to cofactor, cofactor contact, second shell,
#: unresolved. `in_second_shell` is EXCLUDED: it is computed from
#: min(cofactor, substrate), so keeping it would smuggle the substrate
#: distance into an arm whose point is to be without it.
COFACTOR_ONLY_COLUMNS = [0, 2, 5]


def design_effect(row_wise, clustered):
    """How much wider the clustered interval is, and whether a verdict moved.

    Reported rather than asserted. A factor near 1 would mean the rows were
    effectively independent after all and the objection did not bite; a large
    factor means the published row-wise interval was too narrow, in the
    direction that manufactures a verdict.
    """
    width_row = row_wise["interval"][1] - row_wise["interval"][0]
    width_cluster = clustered["interval"][1] - clustered["interval"][0]
    return {
        "row_wise_interval": row_wise["interval"],
        "clustered_interval": clustered["interval"],
        "width_ratio": float(width_cluster / width_row)
            if width_row > 0 else float("nan"),
        "row_wise_verdict": row_wise["verdict"],
        "clustered_verdict": clustered["verdict"],
        "verdict_changed": bool(row_wise["verdict"] != clustered["verdict"]),
        "n_rows": int(row_wise["n_resampling_units"]),
        "n_clusters": int(clustered["n_resampling_units"]),
    }


def analyse(gene, directory, *, verify=False, pdb_dir=None):
    annotation = ANNOTATION[gene]
    activity_path = (directory
                     / f"mavedb_{annotation['activity_urn'].split(':')[-1]}"
                       "_scores.csv")
    abundance_path = (directory
                      / f"mavedb_{annotation['abundance_urn'].split(':')[-1]}"
                        "_scores.csv")
    for path in (activity_path, abundance_path):
        if not path.is_file():
            raise RunError(f"missing score file: {path}")

    activity = load_scores(activity_path)
    abundance = load_scores(abundance_path)
    shared = sorted(set(activity) & set(abundance))
    if len(shared) < 200:
        return {"verdict": "NOT_EVALUABLE",
                "reason": f"only {len(shared)} variants measured by both assays"}

    parsed = [activity[v][0] for v in shared]
    y = np.array([activity[v][1] for v in shared], dtype=float)
    x = np.array([abundance[v][1] for v in shared], dtype=float)
    positions = np.array([p[1] for p in parsed], dtype=int)
    is_nonsense = np.array([p[2] is None for p in parsed])
    is_synonymous = np.array([p[2] == "=" for p in parsed])

    keep = ~is_synonymous
    parsed = [p for p, k in zip(parsed, keep) if k]
    y, x, positions = y[keep], x[keep], positions[keep]
    is_nonsense = is_nonsense[keep]

    substitution, mechanism = VMS.build_design(parsed, annotation)
    position_folds = VMS.position_groups(positions, VMS.FROZEN["position_folds"],
                                         VMS.FROZEN["seed"])
    region_folds = VMS.region_groups(positions, VMS.FROZEN["region_blocks"],
                                     annotation["length"])
    # The folds being grouped does not make the INTERVAL grouped, and only the
    # interval decides a verdict. Up to nineteen substitutions share one
    # residue here, so a row bootstrap counts correlated observations as
    # independent and narrows every interval in the direction that
    # manufactures a verdict. Primary unit: the position, where the
    # dependence is mechanical. Sensitivity: ten-residue blocks, which also
    # absorb local structural context.
    boot_units = {
        "protein_position": positions,
        "contiguous_block_of_10_residues": VMS.contiguous_blocks(positions, 10),
    }

    out = {
        "gene": gene,
        "role": annotation["role"],
        "uniprot": annotation["uniprot"],
        "assay_confound": annotation["assay_confound"],
        "n_activity": len(activity),
        "n_abundance": len(abundance),
        "n_shared_parsed": int(y.size),
        "n_distinct_positions": int(np.unique(positions).size),
        "n_nonsense": int(is_nonsense.sum()),
        "n_synonymous_dropped": int((~keep).sum()),
        # Written as text, not a list. This project forbids a numeric
        # list of ten or more entries in a committed artefact, and
        # BRAF has ten annotated sites. The rule exists to stop
        # per-row data leaking into the repository; declared
        # annotation is not that, but a guard that holds only when
        # convenient is not a guard.
        "annotated_functional_sites": ", ".join(
            str(int(site)) for site in annotation["functional_sites"]),
        "n_annotated_functional_sites": len(
            annotation["functional_sites"]),
        "annotated_motifs": [list(m) for m in annotation["motifs"]],
        "spearman_abundance_activity": spearman(x, y),
        "spearman_excluding_nonsense": spearman(x[~is_nonsense],
                                                y[~is_nonsense]),
        "activity_score_file_sha256": hashlib.sha256(
            activity_path.read_bytes()).hexdigest(),
        "abundance_score_file_sha256": hashlib.sha256(
            abundance_path.read_bytes()).hexdigest(),
    }
    if verify:
        out["uniprot_verification"] = verify_uniprot(gene, annotation)

    # The discordant quadrant: abundance preserved, activity lost.
    hi_x = float(np.quantile(x, 0.75))
    lo_y = float(np.quantile(y, 0.25))
    discordant = (x >= hi_x) & (y <= lo_y)
    out["discordant_share_abundance_high_activity_low"] = float(
        discordant.mean())
    out["discordant_positions_most_frequent"] = [
        [int(p), int(c)] for p, c in sorted(
            zip(*np.unique(positions[discordant], return_counts=True)),
            key=lambda row: -row[1])[:8]]

    # Does the mechanistic hypothesis hold WHERE the annotation exists? This
    # tests the commissioning note's own table row directly: "missense near a
    # catalytic site -> protein present but inactive". It also separates a
    # weak method from a weak knowledge layer, which is the difference that
    # decides what to build next.
    sites = annotation["functional_sites"]
    if sites:
        distance = np.array([min(abs(int(p) - s) for s in sites)
                             for p in positions], dtype=float)
        near = distance <= 3.0
        far = ~near
        seq_intervals = {
            "variant_iid_too_narrow": SA.discordance_ratio_interval(
                near, discordant, n_boot=2000, seed=VMS.FROZEN["seed"] + 5),
        }
        for unit, labels in boot_units.items():
            seq_intervals[unit] = SA.discordance_ratio_interval(
                near, discordant, groups=labels, n_boot=2000,
                seed=VMS.FROZEN["seed"] + 5)
        out["discordance_near_annotated_sites"] = {
            "n_variants_within_3_residues": int(near.sum()),
            "n_variants_beyond": int(far.sum()),
            "discordance_rate_near": float(discordant[near].mean())
                if near.any() else float("nan"),
            "discordance_rate_far": float(discordant[far].mean()),
            "enrichment_ratio": (float(discordant[near].mean()
                                       / discordant[far].mean())
                                 if near.any() and discordant[far].mean() > 0
                                 else float("nan")),
            # Headline: the position-clustered interval. The row-wise one is
            # kept visible so the design effect can be read off, not hidden.
            "enrichment_interval": seq_intervals["protein_position"],
            "enrichment_interval_by_resampling_unit": seq_intervals,
            "resampling_unit": "protein_position",
            "share_of_discordant_variants_the_annotation_reaches": float(
                near[discordant].mean()) if discordant.any() else float("nan"),
            "how_to_read": (
                "an enrichment above 1 supports the hypothesis that a variant "
                "near annotated machinery keeps abundance and loses activity. "
                "The LAST figure is the one that decides what to build next: "
                "it is the share of discordant variants the annotation can "
                "even see. A high enrichment with a low share means the rule "
                "is right and the knowledge layer is too sparse to use it."),
        }

    # ---- the structural layer: a 3D ruler where the sequence ruler failed
    # The earlier diagnosis said annotation COVERAGE was the binding
    # constraint. That is refutable: if it was right, 3D proximity to the
    # catalytic machinery should reach the discordant residues that sequence
    # proximity misses. If it does not, the diagnosis was wrong.
    structural = None
    if pdb_dir is not None and gene in SA.CATALYTIC_CENTRES:
        centre = SA.CATALYTIC_CENTRES[gene]
        path = pathlib.Path(pdb_dir) / f"{centre['pdb']}.pdb"
        if path.is_file():
            text = path.read_text(encoding="utf-8", errors="replace")
            # A chain must be named for a complex: in 6PP9 the substrate MEK1
            # is chain B and carries its own AMP-PNP and Mg, so unrestricted
            # ligands would measure the distance to the OTHER kinase's
            # nucleotide and never say so.
            residues, ligands = SA.parse_pdb(text, chain=centre.get("chain"))
            guard = SA.verify_catalytic_geometry(residues, ligands, centre)
            substrate_points = (
                SA.polymer_chain_points(text, centre["substrate_chain"])
                if centre.get("substrate_chain") else None)
            distances = SA.residue_distances(
                residues, ligands, centre, substrate_points=substrate_points)
            structural = SA.build_structural_block(positions, distances)
            reach = SA.coverage(positions, discordant, distances)
            out["structural_layer"] = {
                "pdb": centre["pdb"],
                "resolution_angstrom": centre["resolution_angstrom"],
                "title": centre["title"],
                "chain": centre.get("chain"),
                "cofactor_ligands": list(centre["cofactor_ligands"]),
                "substrate_ligands": list(centre["substrate_ligands"]),
                "substrate_chain": centre.get("substrate_chain"),
                "substrate_is_a_protein": centre.get("substrate_is_a_protein"),
                "numbering_guard": guard,
                "structure_sha256": hashlib.sha256(
                    path.read_bytes()).hexdigest(),
                "coverage": reach,
                "enrichment_interval": SA.enrichment_interval(
                    positions, discordant, distances,
                    groups=boot_units["protein_position"]),
                "enrichment_interval_by_resampling_unit": {
                    "variant_iid_too_narrow": SA.enrichment_interval(
                        positions, discordant, distances),
                    **{unit: SA.enrichment_interval(
                        positions, discordant, distances, groups=labels)
                       for unit, labels in boot_units.items()},
                },
                "resampling_unit": "protein_position",
                # The 8 A cutoff was declared before this run, but it was
                # declared AFTER the discordant residues of these two proteins
                # had been looked at. A single threshold cannot answer that
                # objection; a sweep can show whether the threshold is
                # carrying the result or merely sitting inside a broad plateau.
                "enrichment_by_cutoff": {
                    f"{cut:g}": dict(
                        SA.coverage(positions, discordant, distances,
                                    cutoff=cut),
                        interval_position_clustered=SA.enrichment_interval(
                            positions, discordant, distances, cutoff=cut,
                            groups=boot_units["protein_position"]))
                    for cut in (4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0)},
                "distance_for_the_residues_the_sequence_ruler_missed": {
                    str(int(p)): {
                        "to_cofactor": distances[int(p)]["to_cofactor"],
                        "to_substrate": distances[int(p)]["to_substrate"],
                        "sequence_distance_to_nearest_annotation": min(
                            abs(int(p) - site)
                            for site in annotation["functional_sites"]),
                    }
                    for p, _ in out["discordant_positions_most_frequent"]
                    if int(p) in distances},
                "what_this_is_not": SA.what_this_module_does_not_do(),
            }

    runs, contrasts = {}, {}
    for label, folds in (("grouped_by_position", position_folds),
                         ("grouped_by_contiguous_region", region_folds)):
        real = VMS.evaluate(y, x, substitution, mechanism, folds, positions)
        no_abundance = VMS.evaluate(y, x, substitution, mechanism, folds,
                                    positions, permute_abundance=True)
        no_mechanism = VMS.evaluate(y, x, substitution, mechanism, folds,
                                    positions, permute_mechanism=True)
        # Specificity control: the same information predicting ABUNDANCE. If
        # the variant and mechanism blocks do that just as well, they are not
        # telling us anything specific to activity.
        on_abundance = VMS.evaluate(x, y, substitution, mechanism, folds,
                                    positions)
        runs[label] = {
            # scored_rows is a per-variant mask: it belongs in memory for the
            # pairing assertion, never in a committed artefact.
            "real": {k: v for k, v in real.items()
                     if k not in ("losses", "scored_rows")},
            "permuted_abundance_control": no_abundance["mean_loss"],
            "permuted_mechanism_control": no_mechanism["mean_loss"],
            "predicting_abundance_instead_specificity_control": {
                "mean_loss": on_abundance["mean_loss"],
                "relative_gain_of_mechanism_over_intercept": float(
                    1.0 - on_abundance["mean_loss"]["M3_mechanism_only"]
                    / on_abundance["mean_loss"]["M0_intercept"]),
                "what_must_happen": (
                    "the mechanism block should help LESS for abundance than "
                    "for activity. If it helps equally, it is describing the "
                    "variant in general and not its functional consequence."),
            },
        }
        cluster = boot_units["protein_position"]
        if structural is not None:
            richer = np.column_stack([mechanism, structural])
            with_structure = VMS.evaluate(y, x, substitution, richer, folds,
                                          positions)
            threshold_free = VMS.evaluate(
                y, x, substitution,
                np.column_stack([mechanism, structural[:, CONTINUOUS_COLUMNS]]),
                folds, positions)
            cofactor_only = VMS.evaluate(
                y, x, substitution,
                np.column_stack([mechanism,
                                 structural[:, COFACTOR_ONLY_COLUMNS]]),
                folds, positions)
            runs[label]["with_structural_layer"] = {
                "mean_loss": with_structure["mean_loss"],
                "mean_abundance_blend_weight":
                    with_structure["mean_abundance_blend_weight"],
                "M5_vs_M4": VMS.contrast(
                    with_structure, "M5_conditional_interaction",
                    "M4_fixed_union", family_size=PRIMARY_FAMILY_SIZE,
                    groups=cluster),
                "M3_vs_M0": VMS.contrast(
                    with_structure, "M3_mechanism_only", "M0_intercept",
                    family_size=PRIMARY_FAMILY_SIZE, groups=cluster),
                "M4_vs_M1": VMS.contrast(
                    with_structure, "M4_fixed_union", "M1_abundance_only",
                    family_size=PRIMARY_FAMILY_SIZE, groups=cluster),
                # THE DIRECT CONTRAST. M4 beating M1 says the union beats the
                # best single measurement; it does not say what the six
                # structural variables contributed, because the union already
                # held the sequence mechanism block. This pairs the SAME arm,
                # on the SAME rows and the SAME folds, with and without them.
                "ablation_of_the_structural_block": {
                    arm: VMS.ablation(with_structure, real, arm,
                                      label="structural",
                                      family_size=PRIMARY_FAMILY_SIZE,
                                      groups=cluster)
                    for arm in ("M3_mechanism_only", "M4_fixed_union",
                                "M5_conditional_interaction")},
                "ablation_resampled_by_10_residue_blocks": {
                    arm: VMS.ablation(
                        with_structure, real, arm, label="structural",
                        family_size=PRIMARY_FAMILY_SIZE,
                        groups=boot_units["contiguous_block_of_10_residues"])
                    for arm in ("M3_mechanism_only", "M4_fixed_union")},
                # The 5 A and 8 A thresholds were declared before this run but
                # AFTER these proteins' discordant residues had been looked
                # at, so a reader is entitled to ask whether they are carrying
                # the result. This arm removes all three threshold indicators
                # and keeps only the two continuous log-distances and the
                # unresolved flag. If the gain survives with no threshold in
                # the model, no threshold can be carrying it.
                "ablation_of_a_threshold_free_block": {
                    arm: VMS.ablation(threshold_free, real, arm,
                                      label="structural_continuous_only",
                                      family_size=PRIMARY_FAMILY_SIZE,
                                      groups=cluster)
                    for arm in ("M3_mechanism_only", "M4_fixed_union")},
                "threshold_free_block_columns": [
                    SA.STRUCTURAL_NAMES[i] for i in CONTINUOUS_COLUMNS],
                # Declared sensitivity for a PROTEIN substrate: an interface
                # distance is not a pocket distance, so the ruler must also
                # be reported with the substrate column removed entirely.
                "ablation_of_a_cofactor_only_block": {
                    arm: VMS.ablation(cofactor_only, real, arm,
                                      label="structural_cofactor_only",
                                      family_size=PRIMARY_FAMILY_SIZE,
                                      groups=cluster)
                    for arm in ("M3_mechanism_only", "M4_fixed_union")},
                "cofactor_only_block_columns": [
                    SA.STRUCTURAL_NAMES[i] for i in COFACTOR_ONLY_COLUMNS],
                "relative_change_in_M3_from_adding_structure": float(
                    1.0 - with_structure["mean_loss"]["M3_mechanism_only"]
                    / real["mean_loss"]["M3_mechanism_only"]),
                "relative_change_in_M4_from_adding_structure": float(
                    1.0 - with_structure["mean_loss"]["M4_fixed_union"]
                    / real["mean_loss"]["M4_fixed_union"]),
                "what_must_happen": (
                    "if annotation coverage was the binding constraint, the "
                    "mechanism arm M3 must improve when the ruler becomes "
                    "three-dimensional. If it does not, the diagnosis was "
                    "wrong."),
            }

        contrasts[label] = {
            "M5_vs_M4": VMS.contrast(real, "M5_conditional_interaction",
                                     "M4_fixed_union",
                                     family_size=PRIMARY_FAMILY_SIZE,
                                     groups=cluster),
            "M6_vs_M4": VMS.contrast(real, "M6_trust_gate", "M4_fixed_union",
                                     family_size=PRIMARY_FAMILY_SIZE,
                                     groups=cluster),
            "M4_vs_M1": VMS.contrast(real, "M4_fixed_union",
                                     "M1_abundance_only",
                                     family_size=PRIMARY_FAMILY_SIZE,
                                     groups=cluster),
            "M1_vs_M0": VMS.contrast(real, "M1_abundance_only",
                                     "M0_intercept",
                                     family_size=PRIMARY_FAMILY_SIZE,
                                     groups=cluster),
            "M3_vs_M0": VMS.contrast(real, "M3_mechanism_only", "M0_intercept",
                                     family_size=PRIMARY_FAMILY_SIZE,
                                     groups=cluster),
            "M5_vs_M4_permuted_mechanism": VMS.contrast(
                no_mechanism, "M5_conditional_interaction", "M4_fixed_union",
                family_size=PRIMARY_FAMILY_SIZE, groups=cluster),
            "resampling_unit": "protein_position (cluster bootstrap)",
            # How much the row bootstrap was understating the width, measured
            # on the headline contrast rather than asserted.
            "design_effect_on_M4_vs_M1": design_effect(
                VMS.contrast(real, "M4_fixed_union", "M1_abundance_only",
                             family_size=PRIMARY_FAMILY_SIZE),
                VMS.contrast(real, "M4_fixed_union", "M1_abundance_only",
                             family_size=PRIMARY_FAMILY_SIZE,
                             groups=cluster)),
        }
    out["runs"] = runs
    out["contrasts"] = contrasts

    # Does the mechanism block predict WHERE abundance fails? This is the
    # selection question stated directly, without any blending.
    from omics2pk.fusion_headroom_diagnostics import ridge_predict
    standard_x = (x - x.mean()) / (x.std() or 1.0)
    standard_y = (y - y.mean()) / (y.std() or 1.0)
    gap = np.abs(standard_x - standard_y)
    without = np.column_stack([substitution, mechanism])
    predicted = np.full(gap.size, np.nan)
    for fold in np.unique(position_folds):
        held = np.flatnonzero(position_folds == fold)
        kept = np.flatnonzero(position_folds != fold)
        if kept.size < 10:
            continue
        predicted[held] = ridge_predict(without[kept], gap[kept],
                                        without[held])
    ok = np.isfinite(predicted)
    rng = np.random.default_rng(VMS.FROZEN["seed"])
    draws = []
    index = np.flatnonzero(ok)
    for _ in range(600):
        picked = rng.choice(index, size=index.size, replace=True)
        draws.append(spearman(predicted[picked], gap[picked]))
    out["can_we_anticipate_where_abundance_fails"] = {
        "target": ("|standardised abundance - standardised activity|, the "
                   "per-variant disagreement between the two assays"),
        "spearman_out_of_sample": spearman(predicted[ok], gap[ok]),
        "spearman_interval": [float(np.percentile(draws, 2.5)),
                              float(np.percentile(draws, 97.5))],
        "n": int(ok.sum()),
        "grouping": "protein position",
        "how_to_read": (
            "this is the user's question asked directly: can variant and "
            "mechanism information predict, for an unseen residue, that "
            "abundance will disagree with activity? It is a correlation, not "
            "a performance, and it carries no margin verdict."),
    }
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mavedb-dir", required=True)
    parser.add_argument("--genes", nargs="*", default=list(ANNOTATION))
    parser.add_argument("--verify-uniprot", action="store_true")
    parser.add_argument("--pdb-dir", default=None,
                        help="directory holding the declared PDB structures; "
                             "omitted means the structural layer is skipped")
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "variant_mechanism_benchmark.json"))
    args = parser.parse_args(argv)

    directory = pathlib.Path(args.mavedb_dir).resolve()
    if ROOT == directory or ROOT in directory.parents:
        raise RunError("the downloaded score files must stay outside the "
                       "repository")

    genes = {}
    for gene in args.genes:
        if gene not in ANNOTATION:
            raise RunError(f"unknown gene {gene}; known: {sorted(ANNOTATION)}")
        print(f"  {gene} [{ANNOTATION[gene]['role']}] ...", flush=True)
        genes[gene] = analyse(gene, directory, verify=args.verify_uniprot,
                              pdb_dir=args.pdb_dir)

    payload = {
        "run": "variant_mechanism_omic_selection",
        "generated": "2026-11-11",
        "question": ("Given a variant, can its presumed mechanism tell us "
                     "which biological measurement carries the information "
                     "about its functional consequence?"),
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in VMS.FROZEN.items()},
        "primary_family_size": PRIMARY_FAMILY_SIZE,
        "primary_gene": "NUDT15",
        "structural_layer_declared": {
            gene: {"pdb": c["pdb"], "uniprot": c["uniprot"],
                   "cofactor": list(c["cofactor_ligands"]),
                   "substrate": list(c["substrate_ligands"])}
            for gene, c in SA.CATALYTIC_CENTRES.items()},
        "positive_control": VMS.self_test(),
        "what_this_is_not": VMS.what_this_module_does_not_do(),
        "genes": genes,
        "no_individual_level_data_in_this_artefact": True,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              default=float) + "\n", encoding="utf-8")

    for gene, block in genes.items():
        if block.get("verdict") == "NOT_EVALUABLE":
            print(f"\n{gene}: NOT_EVALUABLE — {block['reason']}")
            continue
        print(f"\n===== {gene} [{block['role']}] n={block['n_shared_parsed']} "
              f"over {block['n_distinct_positions']} positions =====")
        print(f"  Spearman(abundance, activity) "
              f"{block['spearman_abundance_activity']:+.4f}  "
              f"(excluding nonsense {block['spearman_excluding_nonsense']:+.4f})")
        print(f"  abundance high / activity low: "
              f"{block['discordant_share_abundance_high_activity_low']:.1%}")
        probe = block["can_we_anticipate_where_abundance_fails"]
        print(f"  predicting WHERE the assays disagree: rho "
              f"{probe['spearman_out_of_sample']:+.4f} (n={probe['n']})")
        for label, cell in block["contrasts"].items():
            print(f"  --- {label}")
            loss = block["runs"][label]["real"]["mean_loss"]
            print("      " + "  ".join(
                f"{k.split('_')[0]}={v:.3f}" for k, v in loss.items()))
            for key, value in cell.items():
                if not (isinstance(value, dict) and "relative_gain" in value):
                    continue
                print(f"      {key:<32} {value['relative_gain']:+7.2%} "
                      f"[{value['interval'][0]:+.4f},"
                      f"{value['interval'][1]:+.4f}] {value['verdict']}")
            effect = cell.get("design_effect_on_M4_vs_M1")
            if effect:
                print(f"      {'(row bootstrap would give)':<32} "
                      f"[{effect['row_wise_interval'][0]:+.4f},"
                      f"{effect['row_wise_interval'][1]:+.4f}] "
                      f"{effect['row_wise_verdict']}  "
                      f"width x{effect['width_ratio']:.2f} from "
                      f"{effect['n_rows']} rows to "
                      f"{effect['n_clusters']} positions")
            layer = block["runs"][label].get("with_structural_layer")
            if layer:
                print("      --- with the structural block")
                for key in ("M3_vs_M0", "M4_vs_M1", "M5_vs_M4"):
                    value = layer[key]
                    print(f"      {key:<32} {value['relative_gain']:+7.2%} "
                          f"[{value['interval'][0]:+.4f},"
                          f"{value['interval'][1]:+.4f}] {value['verdict']}")
                for arm, value in layer[
                        "ablation_of_the_structural_block"].items():
                    print(f"      ABLATE {arm:<25} "
                          f"{value['relative_gain']:+7.2%} "
                          f"[{value['interval'][0]:+.4f},"
                          f"{value['interval'][1]:+.4f}] {value['verdict']}")
    print(f"\nartefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
