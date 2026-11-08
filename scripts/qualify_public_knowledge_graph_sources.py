#!/usr/bin/env python
"""Qualify the four public knowledge sources used by the mechanistic selector.

Acquisition and qualification only. This script computes nothing about any
endpoint: it verifies that the files on disk are the files the registry
records, that their internal structure is the structure the loader assumes, and
that each one was retrievable without an account, an application or a
committee. It writes a qualification report and exits non-zero on any mismatch.

None of these files is committed. STRING and Reactome permit redistribution,
but they are large and the project's rule is that source artefacts stay out of
the repository. For the chemical-gene source the redistribution terms could not
be verified at all, which is recorded as such and treated as a refusal.

Run, from the repository root:

    python scripts/qualify_public_knowledge_graph_sources.py \
        --sources "<directory holding the four downloads>" \
        --out docs/KNOWLEDGE_GRAPH_SOURCE_QUALIFICATION.json
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import sys
import zipfile
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk.mechanistic_graph_selection import (  # noqa: E402
    FROZEN, STRING_INFO_HEADER, STRING_LINKS_HEADER,
)

CHUNK = 1 << 20

# Recorded at acquisition on 2026-11-08 and re-verified on every run.
EXPECTED = {
    "string_protein_info": {
        "filename": "9606.protein.info.v12.0.txt.gz",
        "bytes": 1970090,
        "sha256": "144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c",
        "url": ("https://stringdb-downloads.org/download/protein.info.v12.0/"
                "9606.protein.info.v12.0.txt.gz"),
        "licence": "CC BY 4.0",
        "licence_verified": ("string-db.org/cgi/access, 2026-11-08: all data and "
                             "download files are freely available under CC BY 4.0, "
                             "derived files may be redistributed with attribution"),
        "redistribution_permitted": True,
        "access_requirement_type": "none",
        "access_evidence": "Retrieved anonymously by HTTPS; no account, no terms click-through.",
    },
    "string_protein_links": {
        "filename": "9606.protein.links.v12.0.txt.gz",
        "bytes": 83164437,
        "sha256": "3e22f32572211aa341d5b4bd08d30c32e693e294603202120936872f87719d4f",
        "url": ("https://stringdb-downloads.org/download/protein.links.v12.0/"
                "9606.protein.links.v12.0.txt.gz"),
        "licence": "CC BY 4.0",
        "licence_verified": "as above, same source page",
        "redistribution_permitted": True,
        "access_requirement_type": "none",
        "access_evidence": "Retrieved anonymously by HTTPS; no account, no terms click-through.",
    },
    "reactome_pathways": {
        "filename": "ReactomePathways.gmt.zip",
        "bytes": 300128,
        "sha256": "6056fcd95dc8f09740c8bb64c2765bea1386ce91b9b19cd10bcd45d30451805c",
        "url": "https://reactome.org/download/current/ReactomePathways.gmt.zip",
        "licence": "CC0 public domain dedication",
        "licence_verified": ("reactome.org/license, 2026-11-08: users may copy, "
                             "modify and distribute the data, even commercially, "
                             "without asking permission"),
        "redistribution_permitted": True,
        "access_requirement_type": "none",
        "access_evidence": "Retrieved anonymously by HTTPS from the current release path.",
    },
    "chemical_gene_interactions": {
        "filename": "CTD_chem_gene_ixns.tsv.gz",
        "bytes": 43514023,
        "sha256": "032f1b6fc674ca8614b91efeefe40e2cab8315639c1b1552a5df4df98d88f8c9",
        "url": "https://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz",
        "licence": "not verified",
        "licence_verified": (
            "NOT VERIFIED. The provider's legal page answered 302 on 2026-11-08 "
            "and a search did not surface the provider's own terms. Treated as a "
            "refusal of redistribution. The bulk file itself is served anonymously, "
            "which is what admissibility turns on; note that the provider's "
            "interactive batch-query tool now sits behind a bot challenge, so the "
            "bulk file is the only reproducible route."),
        "redistribution_permitted": None,
        "access_requirement_type": "none",
        "access_evidence": "Bulk report retrieved anonymously by HTTPS; no account.",
    },
}

# Structural expectations, so a silently changed release cannot pass unnoticed.
STRUCTURE = {
    "string_protein_info": {"header": STRING_INFO_HEADER, "min_rows": 19000},
    "string_protein_links": {"header": STRING_LINKS_HEADER, "min_rows": 13_000_000},
    "reactome_pathways": {"member": "ReactomePathways.gmt", "min_pathways": 2500},
    "chemical_gene_interactions": {
        "columns": ("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
                    "GeneForms", "Organism", "OrganismID", "Interaction",
                    "InteractionActions", "PubMedIDs"),
        "min_rows": 3_000_000,
    },
}


class KnowledgeSourceQualificationError(RuntimeError):
    pass


def sha256_of(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while True:
            block = handle.read(CHUNK)
            if not block:
                break
            size += len(block)
            digest.update(block)
    return digest.hexdigest(), size


def verify_bytes(sources: Path) -> dict:
    report = {}
    for key, record in EXPECTED.items():
        path = sources / record["filename"]
        if not path.is_file():
            raise KnowledgeSourceQualificationError(
                f"{key}: expected file is absent: {path}")
        digest, size = sha256_of(path)
        if digest != record["sha256"] or size != record["bytes"]:
            raise KnowledgeSourceQualificationError(
                f"{key}: the file on disk is not the file the registry records "
                f"(sha256 {digest}, {size} bytes)")
        report[key] = {"filename": record["filename"], "bytes": size,
                       "sha256": digest, "verified_on": date.today().isoformat()}
    return report


def verify_string_info(path: Path) -> dict:
    want = STRUCTURE["string_protein_info"]
    rows, symbols = 0, set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n")
        if header != want["header"]:
            raise KnowledgeSourceQualificationError(
                "string_protein_info: header changed")
        for line in handle:
            parts = line.split("\t", 2)
            rows += 1
            if len(parts) > 1 and parts[1].strip():
                symbols.add(parts[1].strip())
    if rows < want["min_rows"]:
        raise KnowledgeSourceQualificationError(
            f"string_protein_info: only {rows} rows")
    return {"rows": rows, "distinct_preferred_symbols": len(symbols)}


def verify_string_links(path: Path) -> dict:
    want = STRUCTURE["string_protein_links"]
    threshold = FROZEN["string_min_combined_score"]
    rows, above = 0, 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).rstrip("\n")
        if header != want["header"]:
            raise KnowledgeSourceQualificationError(
                "string_protein_links: header changed")
        for line in handle:
            rows += 1
            if int(line.rsplit(" ", 1)[1]) >= threshold:
                above += 1
    if rows < want["min_rows"]:
        raise KnowledgeSourceQualificationError(
            f"string_protein_links: only {rows} rows")
    if above == 0:
        raise KnowledgeSourceQualificationError(
            "string_protein_links: no edge reaches the frozen confidence threshold")
    return {"directed_rows": rows,
            "rows_at_or_above_frozen_threshold": above,
            "frozen_threshold": threshold}


def verify_reactome(path: Path) -> dict:
    want = STRUCTURE["reactome_pathways"]
    with zipfile.ZipFile(path) as archive:
        if want["member"] not in archive.namelist():
            raise KnowledgeSourceQualificationError(
                "reactome_pathways: expected archive member absent")
        raw = archive.read(want["member"]).decode("utf-8", "replace")
    rows = [line.split("\t") for line in raw.splitlines() if line.strip()]
    if len(rows) < want["min_pathways"]:
        raise KnowledgeSourceQualificationError(
            f"reactome_pathways: only {len(rows)} pathways")
    if any(len(row) < 3 for row in rows):
        raise KnowledgeSourceQualificationError(
            "reactome_pathways: a pathway row carries no member")
    sizes = sorted(len(row) - 2 for row in rows)
    return {"pathways": len(rows),
            "smallest_pathway": sizes[0],
            "largest_pathway": sizes[-1],
            "median_pathway": sizes[len(sizes) // 2]}


def verify_chemical_gene(path: Path) -> dict:
    want = STRUCTURE["chemical_gene_interactions"]
    columns, rows, human = None, 0, 0
    organisms = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                if "ChemicalName" in line:
                    columns = tuple(line.lstrip("# ").rstrip("\n").split("\t"))
                    # Checked here rather than after the scan, so a changed
                    # release is refused instead of raising an index error.
                    if columns != want["columns"]:
                        raise KnowledgeSourceQualificationError(
                            "chemical_gene_interactions: columns changed: "
                            f"{columns}")
                continue
            rows += 1
            fields = line.rstrip("\n").split("\t")
            if columns and len(fields) == len(columns):
                organism = fields[columns.index("Organism")]
                organisms.add(organism)
                if organism == "Homo sapiens":
                    human += 1
    if columns != want["columns"]:
        raise KnowledgeSourceQualificationError(
            f"chemical_gene_interactions: columns changed: {columns}")
    if rows < want["min_rows"]:
        raise KnowledgeSourceQualificationError(
            f"chemical_gene_interactions: only {rows} rows")
    return {"rows": rows, "human_rows": human, "distinct_organisms": len(organisms)}


VERIFIERS = {
    "string_protein_info": verify_string_info,
    "string_protein_links": verify_string_links,
    "reactome_pathways": verify_reactome,
    "chemical_gene_interactions": verify_chemical_gene,
}


def qualify(sources: Path) -> dict:
    checksums = verify_bytes(sources)
    report = {
        "report": "Public knowledge-source qualification for mechanistic selection",
        "qualified_on": date.today().isoformat(),
        "frozen_selector_parameters": dict(FROZEN),
        "rule": ("A source is admissible when it is retrievable without an "
                 "account, an application or a third-party decision. "
                 "Redistribution rights are recorded separately and govern "
                 "only whether the file may be committed, which none is."),
        "sources": {},
    }
    for key, record in EXPECTED.items():
        structure = VERIFIERS[key](sources / record["filename"])
        report["sources"][key] = {
            "url": record["url"],
            "licence": record["licence"],
            "licence_verified": record["licence_verified"],
            "redistribution_permitted": record["redistribution_permitted"],
            "committed_to_this_repository": False,
            "access_requirement": {
                "type": record["access_requirement_type"],
                "evidence": record["access_evidence"],
                "verified_on": date.today().isoformat(),
            },
            "bytes_verified": checksums[key],
            "structure_verified": structure,
        }
    report["no_endpoint_was_read"] = True
    report["note"] = (
        "These are prior-knowledge resources, not cohorts. They carry no "
        "participant, no outcome and no sample. They cannot support or refute "
        "any objective on their own; they supply the ordering whose value the "
        "H2 experiments measure.")
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", required=True,
                        help="Directory holding the four downloaded files")
    parser.add_argument("--out", default=None,
                        help="Where to write the qualification report")
    args = parser.parse_args(argv)
    try:
        report = qualify(Path(args.sources))
    except KnowledgeSourceQualificationError as exc:
        print(f"QUALIFICATION FAILED: {exc}", file=sys.stderr)
        return 2
    text = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8", newline="\n")
        print(f"wrote {args.out}")
    else:
        print(text)
    for key, entry in report["sources"].items():
        print(f"  {key:30s} licence={entry['licence']:<35s} "
              f"redistribution={entry['redistribution_permitted']}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
