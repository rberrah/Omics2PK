#!/usr/bin/env python
"""Grade every candidate probe, then build the qualified pools H2d will use.

Replaces `qualify_all_candidate_probes.py`, whose single `identified_human`
verdict covered both a real concordance between two annotations and a lone
annotation with nothing to agree with it, and which treated any difference
between two accession strings as a conflict even when both designate transcripts
of the same gene.

Outputs, all written with the provenance needed to reconstruct every decision:

* the graded verdict for each of the 40,638 features, with its evidence chain;
* **pool A**, the qualified human pool, built probe by probe then grouped;
* **pool B**, pool A intersected with the public network's genes;
* the undetermined stratum, kept separate rather than deleted;
* the frozen provenance block: source files, their SHA-256, the rule version,
  and the complete list of admissible feature identifiers.

    python scripts/build_qualified_probe_pools.py \
        --hlc "<archive dir>" --platform <GPL table> \
        --gene-info <Homo_sapiens.gene_info.gz> \
        --accession-map <human_accession_to_gene.tsv> \
        --knowledge "<knowledge dir>" \
        --out registry/QUALIFIED_PROBE_POOLS.json \
        --full-out "<local path for the per-probe table>"
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from datetime import date
from pathlib import Path
from zipfile import ZipFile

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk.probe_identity_qualification import (  # noqa: E402
    ADMISSIBLE_BY_DEFAULT,
    VERDICTS,
    ProbeIdentityError,
    build_probe_pool,
    qualify_probe,
    summarise,
    symbol_to_gene_index,
)

RULE_VERSION = "probe-identity-1.1.0"
EXPECTED_PROBES = 40638
CHUNK = 1 << 20


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            block = handle.read(CHUNK)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def read_human_genes(path: Path):
    by_id = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).lstrip("#").rstrip("\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            synonyms = frozenset(
                s for s in fields[index["Synonyms"]].split("|") if s != "-")
            by_id[fields[index["GeneID"]]] = (fields[index["Symbol"]], synonyms)
    if len(by_id) < 50_000:
        raise ProbeIdentityError(
            f"the human gene reference looks truncated: {len(by_id)} records")
    return by_id


def read_accession_map(path: Path):
    mapping: dict[str, list[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            accession, _, gene = line.rstrip("\n").partition("\t")
            if accession and gene:
                genes = mapping.setdefault(accession, [])
                if gene not in genes:
                    genes.append(gene)
    if len(mapping) < 100_000:
        raise ProbeIdentityError(
            f"the accession map looks truncated: {len(mapping)} accessions")
    return mapping


def read_archive(hlc_dir: Path):
    import pandas as pd

    with ZipFile(hlc_dir / "curatedExpressionLiver.zip") as archive:
        with archive.open("features.txt") as handle:
            features = pd.read_csv(handle, sep="\t", dtype=str)
    if len(features) != EXPECTED_PROBES:
        raise ProbeIdentityError(f"the archive annotation changed: {len(features)}")
    out = {}
    for row in features.itertuples(index=False):
        out[str(row.feature_id)] = (
            row.genesymbol if isinstance(row.genesymbol, str) else "",
            row.accessioncode if isinstance(row.accessioncode, str) else "")
    return out


def read_platform(path: Path):
    out, inside = {}, False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                header = next(handle).rstrip("\n").split("\t")
                index = {name: position for position, name in enumerate(header)}
                inside = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not inside:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(index):
                fields += [""] * (len(index) - len(fields))
            out[fields[index["ID"]]] = (fields[index["EntrezGeneID"]].strip(),
                                        fields[index["GB_ACC"]].strip())
    if not out:
        raise ProbeIdentityError("the platform record carried no table")
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--gene-info", required=True)
    parser.add_argument("--accession-map", required=True)
    parser.add_argument("--knowledge", default=None,
                        help="Knowledge directory; needed for pool B.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--full-out", default=None)
    args = parser.parse_args(argv)

    platform_path = Path(args.platform)
    gene_info_path = Path(args.gene_info)
    accession_path = Path(args.accession_map)

    human_genes = read_human_genes(gene_info_path)
    # The symbol index decides whether a symbol that disagrees with the
    # resolved gene is a competing claim or merely a name the reference no
    # longer knows. Built from the same gene reference as the verdicts.
    symbol_index = symbol_to_gene_index(human_genes)
    accession_to_genes = read_accession_map(accession_path)
    archive = read_archive(Path(args.hlc))
    platform = read_platform(platform_path)

    identities = []
    for feature_id, (symbol, accession) in archive.items():
        entrez, platform_accession = platform.get(feature_id, ("", ""))
        identities.append(qualify_probe(
            feature_id, symbol=symbol, archive_accession=accession,
            platform_accession=platform_accession, platform_entrez=entrez,
            accession_to_human_genes=accession_to_genes,
            human_gene_symbols=human_genes,
            symbol_to_human_genes=symbol_index))

    overview = summarise(identities)
    pool_a = build_probe_pool(identities, admit=ADMISSIBLE_BY_DEFAULT)
    pool_concordant = build_probe_pool(identities, admit=("concordant_human",))

    pools = {"A_qualified_human": pool_a,
             "A_strict_concordant_only": pool_concordant}
    network_genes = None
    if args.knowledge:
        from omics2pk.mechanistic_graph_selection import load_string_symbol_graph

        knowledge = Path(args.knowledge)
        graph = load_string_symbol_graph(
            knowledge / "9606.protein.info.v12.0.txt.gz",
            knowledge / "9606.protein.links.v12.0.txt.gz")
        network_genes = sorted(graph.index)
        pools["B_qualified_human_on_network"] = build_probe_pool(
            identities, admit=ADMISSIBLE_BY_DEFAULT,
            restrict_to_genes=network_genes)

    undetermined = sorted(i.feature_id for i in identities
                          if i.verdict == "undetermined")
    conflicting = sorted(i.feature_id for i in identities
                         if i.verdict == "conflicting")
    wrong_gene = sorted(i.feature_id for i in identities
                        if i.verdict == "symbol_names_another_gene")
    incoherent = sorted(i.feature_id for i in identities
                        if i.verdict == "platform_self_contradictory")
    stale_alias = sorted(
        (i.feature_id, i.symbol, i.human_symbol)
        for i in identities if i.symbol_is_a_stale_alias)

    report = {
        "report": "Graded probe identity and the qualified pools",
        "built_on": date.today().isoformat(),
        "rule_version": RULE_VERSION,
        "verdict_order_strongest_first": list(VERDICTS),
        "admissible_by_default": list(ADMISSIBLE_BY_DEFAULT),
        "provenance": {
            "cohort_archive_member": "features.txt in curatedExpressionLiver.zip",
            "platform_record": {
                "file": platform_path.name, "sha256": sha256_of(platform_path)},
            "human_gene_reference": {
                "file": gene_info_path.name, "sha256": sha256_of(gene_info_path)},
            "accession_to_gene_map": {
                "file": accession_path.name,
                "sha256": sha256_of(accession_path),
                "derived_from": ("NCBI gene2accession.gz, streamed and filtered "
                                 "to tax_id 9606, keeping RNA, protein and "
                                 "genomic accessions without version suffix"),
                "accessions": len(accession_to_genes)},
            "network_genes": (len(network_genes) if network_genes else None),
        },
        "overview": overview,
        "pools": {name: {key: value for key, value in pool.items()
                         if key != "by_symbol"}
                  for name, pool in pools.items()},
        "strata_kept_separate": {
            "undetermined_probes": len(undetermined),
            "conflicting_probes": len(conflicting),
            "symbol_names_another_gene_probes": len(wrong_gene),
            "platform_self_contradictory_probes": len(incoherent),
            "admissible_under_a_stale_alias": len(stale_alias),
            "why": ("Neither stratum is deleted. An undetermined probe is not a "
                    "wrong measurement, only one whose biological meaning the "
                    "available annotation does not pin down, and dropping the "
                    "stratum everywhere would quietly remove real biology. Each "
                    "sensitivity analysis states which strata it admits."),
        },
        "admissible_features": {
            name: sorted(f for members in pool["by_symbol"].values()
                         for f in members)
            for name, pool in pools.items()
        },
        "undetermined_features": undetermined,
        "conflicting_features": conflicting,
        "symbol_names_another_gene_features": wrong_gene,
        "platform_self_contradictory_features": incoherent,
        "stale_alias_features": [list(row) for row in stale_alias],
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    print(f"wrote {out}")
    for verdict, count in overview["counts_by_verdict"].items():
        print(f"  {verdict:24s} {count:6d}")
    for name, pool in pools.items():
        print(f"  pool {name:30s} {pool['symbols']:6d} symbols, "
              f"{pool['probes']:6d} probes, "
              f"{len(pool['symbols_keeping_only_some_probes']):5d} partial")

    if args.full_out:
        path = Path(args.full_out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(
            {i.feature_id: {
                "symbol": i.symbol, "verdict": i.verdict, "reason": i.reason,
                "archive_accession": i.archive_accession,
                "archive_genes": list(i.archive_genes),
                "platform_accession": i.platform_accession,
                "platform_accession_genes": list(i.platform_accession_genes),
                "platform_entrez": i.platform_entrez,
                "human_gene": i.human_gene, "human_symbol": i.human_symbol,
                "symbol_matches_gene": i.symbol_matches_gene}
             for i in identities}, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8", newline="\n")
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
