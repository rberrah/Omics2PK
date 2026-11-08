#!/usr/bin/env python
"""Qualify the probes the H2a arms actually selected, against public records.

Why this exists. A statistical selector ranked `NegCont_210`, one of the array's
own negative controls, into its top five on one endpoint, and its top choices on
two others were symbols like `CYP2D17` and `CYP2A24`. Neither observation can be
interpreted from the expression table alone: a pseudogene-looking symbol may
still mark real biology at that locus, and the opposite assumption would be just
as unfounded. So the probe's accession is looked up at the public sequence
archive and the answer is recorded, rather than guessed.

Three outcomes matter and are kept separate:

  * `technical_control` -- the array's own control features, which carry no
    transcript at all;
  * `non_human` -- probes whose accession belongs to another organism, so on
    human liver RNA any signal is cross-hybridisation with the human orthologue
    or noise, and cannot be independent human biology;
  * `human_transcript` / `human_genomic` -- a real human record, which must NOT
    be dismissed merely because its symbol looks like a pseudogene.

Only probes a selection actually reached are looked up, so the number of
requests stays small and every one is auditable. Results are cached in the
output file; a rerun with the cache present performs no network access.

    python scripts/qualify_selected_probes.py \
        --result reports/data/h2a_hlc_selection.json \
        --hlc "<archive directory>" \
        --out docs/H2A_SELECTED_PROBE_QUALIFICATION.json
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path
from zipfile import ZipFile

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
TOOL = "omics2pk"
REQUEST_PAUSE_SECONDS = 0.4

CONTROL_PATTERN = re.compile(
    r"^(NegCont|PosCont|Corner|DarkCorner|SpikeIn|GE_Bright|E1A_r60|ETG)",
    re.IGNORECASE)

MECHANISTIC_ARMS = ("netpk_target_seeded", "netpk_wrong_mechanism",
                    "netpk_pathway_seeded", "netpk_substrate_seeded")
STATISTICAL_ARMS = ("univariate_spearman", "lasso_entry_order", "xgboost_gain")


class ProbeQualificationError(RuntimeError):
    pass


def selected_symbols(result):
    """Every symbol any arm selected, with the arms and endpoints that chose it."""
    chosen = {}
    for endpoint, block in result["endpoints"].items():
        for arm, by_size in (block.get("first_fold_selection") or {}).items():
            for size, symbols in by_size.items():
                for symbol in symbols:
                    entry = chosen.setdefault(
                        str(symbol), {"arms": set(), "endpoints": set()})
                    entry["arms"].add(arm)
                    entry["endpoints"].add(endpoint)
    return chosen


def accession_of(symbols, hlc_dir: Path):
    import pandas as pd

    path = hlc_dir / "curatedExpressionLiver.zip"
    with ZipFile(path) as archive:
        with archive.open("features.txt") as handle:
            features = pd.read_csv(handle, sep="\t", dtype=str)
    features["genesymbol"] = features["genesymbol"].fillna("")
    out = {}
    for symbol in symbols:
        rows = features[features["genesymbol"] == symbol]
        out[symbol] = sorted({str(a) for a in rows["accessioncode"].fillna("")
                              if str(a)})
    return out


def fetch_record(accession):
    """Read the organism and definition for one accession, from the archive."""
    query = urllib.parse.urlencode({"db": "nuccore", "id": accession,
                                    "rettype": "gb", "retmode": "text",
                                    "tool": TOOL})
    request = urllib.request.Request(f"{EUTILS}?{query}",
                                     headers={"User-Agent": TOOL})
    with urllib.request.urlopen(request, timeout=60) as response:
        text = response.read().decode("utf-8", "replace")
    organism, definition, molecule = None, None, None
    for line in text.splitlines():
        if line.startswith("LOCUS") and molecule is None:
            molecule = "mRNA" if " mRNA " in line else (
                "DNA" if " DNA " in line else None)
        elif line.startswith("DEFINITION") and definition is None:
            definition = line[len("DEFINITION"):].strip()
        elif line.strip().startswith("ORGANISM") and organism is None:
            organism = line.strip()[len("ORGANISM"):].strip()
        if organism and definition and molecule:
            break
    if organism is None and definition is None:
        raise ProbeQualificationError(f"{accession}: no record returned")
    return {"organism": organism, "definition": definition,
            "molecule": molecule}


def classify(symbol, accessions, records):
    if CONTROL_PATTERN.match(symbol):
        return "technical_control", "the array's own control feature"
    organisms = {records[a]["organism"] for a in accessions
                 if a in records and records[a]["organism"]}
    if not organisms:
        return "unverified", "no public record was retrievable for its accession"
    if organisms == {"Homo sapiens"}:
        molecules = {records[a]["molecule"] for a in accessions if a in records}
        if molecules == {"DNA"}:
            return "human_genomic", "a human genomic record, not a transcript"
        return "human_transcript", "a human sequence record"
    if "Homo sapiens" in organisms:
        return "mixed_organism", f"accessions span {sorted(organisms)}"
    return "non_human", f"the accession belongs to {sorted(organisms)}"


def qualify(result, hlc_dir: Path, cache=None):
    chosen = selected_symbols(result)
    accessions = accession_of(sorted(chosen), hlc_dir)
    records = dict(cache or {})
    wanted = sorted({a for group in accessions.values() for a in group})
    for accession in wanted:
        if accession in records or CONTROL_PATTERN.match(accession):
            continue
        try:
            records[accession] = fetch_record(accession)
        except Exception as exc:  # noqa: BLE001 - a lookup failure is data
            records[accession] = {"organism": None, "definition": None,
                                  "molecule": None, "error": str(exc)}
        time.sleep(REQUEST_PAUSE_SECONDS)

    report = {
        "report": "Qualification of the probes the H2a arms selected",
        "qualified_on": date.today().isoformat(),
        "source_of_truth": ("NCBI Nucleotide, read by accession for every probe "
                            "a selection reached. Organism and molecule type are "
                            "taken from the record, never inferred from the "
                            "symbol."),
        "probes": {},
        "accession_records": records,
    }
    counts = {}
    for symbol, entry in sorted(chosen.items()):
        group = accessions.get(symbol, [])
        verdict, reason = classify(symbol, group, records)
        counts[verdict] = counts.get(verdict, 0) + 1
        report["probes"][symbol] = {
            "verdict": verdict,
            "reason": reason,
            "accessions": group,
            "definition": next((records[a]["definition"] for a in group
                                if a in records and records[a]["definition"]), None),
            "organism": next((records[a]["organism"] for a in group
                              if a in records and records[a]["organism"]), None),
            "selected_by_mechanistic_arms": sorted(
                entry["arms"] & set(MECHANISTIC_ARMS)),
            "selected_by_statistical_arms": sorted(
                entry["arms"] & set(STATISTICAL_ARMS)),
            "endpoints": sorted(entry["endpoints"]),
        }
    report["counts_by_verdict"] = counts

    inadmissible = {"technical_control", "non_human"}
    by_arm = {}
    for symbol, entry in report["probes"].items():
        for arm in (entry["selected_by_mechanistic_arms"]
                    + entry["selected_by_statistical_arms"]):
            bucket = by_arm.setdefault(arm, {"total": 0, "inadmissible": [],
                                             "unverified": []})
            bucket["total"] += 1
            if entry["verdict"] in inadmissible:
                bucket["inadmissible"].append(
                    {"symbol": symbol, "verdict": entry["verdict"],
                     "endpoints": entry["endpoints"]})
            elif entry["verdict"] == "unverified":
                bucket["unverified"].append(symbol)
    report["by_arm"] = by_arm
    report["how_to_read"] = (
        "A non-human probe cannot carry independent human biology: on human "
        "liver RNA its signal is cross-hybridisation with the human orthologue, "
        "or noise. A technical control carries no transcript at all. Neither "
        "makes the arm that chose it wrong about prediction -- averaging two "
        "cross-hybridising probes of one transcript is a legitimate way to "
        "reduce noise -- but neither can be reported as a distinct biomarker. A "
        "human record is NOT dismissed because its symbol looks like a "
        "pseudogene; that question is left to the record, and where the record "
        "is silent the verdict is 'unverified' rather than a guess."
    )
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--result", required=True)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    result = json.loads(Path(args.result).read_text(encoding="utf-8"))
    out = Path(args.out)
    cache = None
    if out.is_file():
        cache = json.loads(out.read_text(encoding="utf-8")).get("accession_records")
    report = qualify(result, Path(args.hlc), cache=cache)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8", newline="\n")
    print(f"wrote {out}")
    for verdict, count in sorted(report["counts_by_verdict"].items()):
        print(f"  {verdict:20s} {count:4d}")
    for arm, bucket in sorted(report["by_arm"].items()):
        if bucket["inadmissible"]:
            names = ", ".join(f"{i['symbol']}({i['verdict']})"
                              for i in bucket["inadmissible"])
            print(f"  {arm}: {len(bucket['inadmissible'])} of "
                  f"{bucket['total']} inadmissible -> {names}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
