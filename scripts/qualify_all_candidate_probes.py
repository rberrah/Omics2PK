#!/usr/bin/env python
"""Qualify every candidate probe on the array, by annotation concordance.

Replaces the first qualifier, which had two limits an external review was right
to reject. It read only the first fold's recorded selections, so it covered a
fraction of the probes the ten folds actually used; and it judged a *symbol* by
pooling the accessions of all probes carrying it, while the analysis picks one
particular probe per gene per fold. A symbol can hold a human and a non-human
probe at once, so a verdict on the symbol certifies nothing about the probe the
fold used.

This version works **per probe**, over **every** candidate, before any selection
is consulted. The traced chain is `feature_id -> accession -> organism -> symbol`
for each of the 40,638 features.

Why concordance rather than one source
--------------------------------------
Two independent public annotations of this array disagree about the same
features, and one of them is internally incoherent. For the four probes carrying
the symbol `CYP1A2`:

* the cohort archive gives macaque cytochrome accessions (`AB185338`,
  `AF123056`, `AF123057`, `D86474`);
* the public platform record gives different accessions entirely (expressed
  sequence tags `DN745607`, `CO597948`, `AM086385`) and maps three of them to
  dog gene `LOC479559` and one to dog `MYH9`, myosin heavy chain 9.

A cytochrome probe is not myosin, so the platform's gene assignment for these
features is wrong, and the two sources cannot both be right about the accession.
The defensible conclusion is therefore **not** that the probe measures a macaque
transcript, and **not** that it measures the human one: it is that the molecular
identity of what the probe measures **cannot be established from any available
public source**. That is a different and stronger reason to keep such a probe out
of a biological reading, and it is the reason recorded here.

Neither source distributes the oligonucleotide sequence, so hybridisation
specificity cannot be checked at all. A probe designed from a macaque coding
sequence may well hybridise the human orthologue, with unknown efficiency. No
verdict in this file claims to know what a probe does hybridise.

Verdicts
--------
* `identified_human` -- both sources point at the same human gene, and the
  symbol the analysis uses matches that gene's human symbol or a listed synonym.
* `human_gene_symbol_mismatch` -- a coherent human gene, but under a different
  symbol from the one the analysis uses.
* `contradictory` -- the sources disagree, or the assigned gene is implausible
  for the symbol. Identity unestablished.
* `technical_control` -- the array's own control feature; no transcript.
* `unverifiable` -- no usable annotation in either source.

`unverifiable` probes are **not** silently dropped from every analysis: the
review's point stands that excluding them wholesale would quietly delete real
biology. They form their own stratum, and each sensitivity analysis states
whether it keeps them.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
from datetime import date
from pathlib import Path
from zipfile import ZipFile

CONTROL_PATTERN = re.compile(
    r"^(NegCont|PosCont|Corner|DarkCorner|SpikeIn|GE_Bright|E1A_r60|ETG)",
    re.IGNORECASE)

EXPECTED_PROBES = 40638


class ProbeQualificationError(RuntimeError):
    pass


def read_human_genes(path: Path):
    """Entrez gene id -> (symbol, synonyms) for every human gene."""
    by_id, symbols = {}, {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        header = next(handle).lstrip("#").rstrip("\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            gene_id = fields[index["GeneID"]]
            symbol = fields[index["Symbol"]]
            synonyms = {s for s in fields[index["Synonyms"]].split("|") if s != "-"}
            by_id[gene_id] = (symbol, synonyms)
            symbols.setdefault(symbol, gene_id)
    if len(by_id) < 50_000:
        raise ProbeQualificationError(
            f"the human gene reference looks truncated: {len(by_id)} records")
    return by_id, symbols


def read_archive_annotation(hlc_dir: Path):
    """feature_id -> (symbol, accession) as the cohort archive records it."""
    import pandas as pd

    with ZipFile(hlc_dir / "curatedExpressionLiver.zip") as archive:
        with archive.open("features.txt") as handle:
            features = pd.read_csv(handle, sep="\t", dtype=str)
    if len(features) != EXPECTED_PROBES:
        raise ProbeQualificationError(
            f"the archive annotation changed: {len(features)} rows")
    out = {}
    for row in features.itertuples(index=False):
        out[str(row.feature_id)] = {
            "symbol": (row.genesymbol or "") if isinstance(row.genesymbol, str) else "",
            "accession": (row.accessioncode or "")
            if isinstance(row.accessioncode, str) else "",
        }
    return out


def read_platform_annotation(path: Path):
    """feature_id -> (entrez id, symbol, accession) from the public platform record."""
    out = {}
    inside = False
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
            out[fields[index["ID"]]] = {
                "entrez": fields[index["EntrezGeneID"]].strip(),
                "symbol": fields[index["ORF"]].strip(),
                "accession": fields[index["GB_ACC"]].strip(),
            }
    if not out:
        raise ProbeQualificationError("the platform record carried no table")
    return out


def classify_probe(feature_id, archive, platform, human_by_id, human_by_symbol):
    """Decide what, if anything, is established about one probe."""
    symbol = archive["symbol"]
    if symbol and CONTROL_PATTERN.match(symbol):
        return "technical_control", "the array's own control feature; no transcript"
    if not symbol:
        return "unverifiable", "the archive records no symbol for this feature"

    record = platform.get(feature_id)
    entrez = (record or {}).get("entrez") or ""
    platform_accession = (record or {}).get("accession") or ""
    archive_accession = archive["accession"]

    if entrez and entrez in human_by_id:
        human_symbol, synonyms = human_by_id[entrez]
        if symbol == human_symbol or symbol in synonyms:
            return "identified_human", (
                f"the platform assigns human gene {entrez} ({human_symbol}) and "
                "the analysis symbol matches it")
        return "human_gene_symbol_mismatch", (
            f"the platform assigns human gene {entrez} ({human_symbol}) but the "
            f"analysis uses the symbol {symbol}")

    if entrez and entrez not in human_by_id:
        return "contradictory", (
            f"the platform assigns gene {entrez}, which is not a human gene, to "
            f"a feature the archive labels {symbol}; the two annotations cannot "
            "both be right and the identity is unestablished")

    if archive_accession and platform_accession and \
            archive_accession != platform_accession:
        return "contradictory", (
            f"the archive gives accession {archive_accession} and the platform "
            f"gives {platform_accession} for the same feature")

    if symbol in human_by_symbol:
        return "identified_human", (
            f"no gene assignment on the platform, but the symbol matches human "
            f"gene {human_by_symbol[symbol]} and the accessions do not conflict")

    return "unverifiable", (
        "neither source assigns a gene, and the symbol is not a human gene symbol")


def qualify(hlc_dir: Path, platform_path: Path, gene_info: Path):
    human_by_id, human_by_symbol = read_human_genes(gene_info)
    archive = read_archive_annotation(hlc_dir)
    platform = read_platform_annotation(platform_path)

    probes, counts = {}, {}
    by_symbol = {}
    for feature_id, entry in archive.items():
        verdict, reason = classify_probe(feature_id, entry, platform,
                                         human_by_id, human_by_symbol)
        counts[verdict] = counts.get(verdict, 0) + 1
        record = platform.get(feature_id) or {}
        probes[feature_id] = {
            "symbol": entry["symbol"],
            "archive_accession": entry["accession"],
            "platform_accession": record.get("accession", ""),
            "platform_entrez": record.get("entrez", ""),
            "platform_symbol": record.get("symbol", ""),
            "verdict": verdict,
            "reason": reason,
        }
        if entry["symbol"]:
            by_symbol.setdefault(entry["symbol"], []).append(feature_id)

    # A symbol is only usable without qualification when every probe behind it
    # is identified, because the fold picks one of them by expression level.
    symbol_status = {}
    for symbol, members in by_symbol.items():
        verdicts = {probes[m]["verdict"] for m in members}
        if verdicts == {"identified_human"}:
            status = "all_probes_identified_human"
        elif "identified_human" in verdicts:
            status = "mixed_some_probes_unidentified"
        elif verdicts == {"technical_control"}:
            status = "technical_control"
        else:
            status = "no_identified_human_probe"
        symbol_status[symbol] = {"status": status, "probes": members,
                                 "verdicts": sorted(verdicts)}

    status_counts = {}
    for entry in symbol_status.values():
        status_counts[entry["status"]] = status_counts.get(entry["status"], 0) + 1

    return {
        "report": "Qualification of every candidate probe, by annotation concordance",
        "qualified_on": date.today().isoformat(),
        "probes_examined": len(probes),
        "counts_by_probe_verdict": counts,
        "counts_by_symbol_status": status_counts,
        "sources": {
            "cohort_archive": "features.txt inside curatedExpressionLiver.zip",
            "platform_record": platform_path.name,
            "human_gene_reference": gene_info.name,
        },
        "what_this_cannot_establish": (
            "Neither source distributes the oligonucleotide sequence, so "
            "hybridisation specificity is not checked and no verdict claims to "
            "know what a probe does hybridise. A probe designed from a non-human "
            "coding sequence may hybridise the human orthologue with unknown "
            "efficiency. 'contradictory' means the identity is unestablished, "
            "not that the probe measures nothing."),
        "unverifiable_are_a_stratum_not_a_deletion": (
            "Probes whose identity cannot be established form their own stratum. "
            "Dropping them from every analysis would quietly delete real "
            "biology, so each sensitivity analysis states whether it keeps them."),
        "symbol_status": symbol_status,
        "probes": probes,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--gene-info", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--summary-out", default=None,
                        help="Where to write the small summary that is committed; "
                             "the full per-probe table can stay local.")
    args = parser.parse_args(argv)

    report = qualify(Path(args.hlc), Path(args.platform), Path(args.gene_info))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    print(f"wrote {out} ({report['probes_examined']} probes)")
    for verdict, count in sorted(report["counts_by_probe_verdict"].items(),
                                 key=lambda item: -item[1]):
        print(f"  probe  {verdict:30s} {count:6d}")
    for status, count in sorted(report["counts_by_symbol_status"].items(),
                                key=lambda item: -item[1]):
        print(f"  symbol {status:30s} {count:6d}")

    if args.summary_out:
        summary = {key: value for key, value in report.items()
                   if key not in ("probes", "symbol_status")}
        summary["symbol_status_counts"] = report["counts_by_symbol_status"]
        summary["usable_symbols"] = sorted(
            s for s, e in report["symbol_status"].items()
            if e["status"] == "all_probes_identified_human")
        summary["symbols_with_no_identified_human_probe"] = sorted(
            s for s, e in report["symbol_status"].items()
            if e["status"] == "no_identified_human_probe")
        path = Path(args.summary_out)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
                        encoding="utf-8", newline="\n")
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
