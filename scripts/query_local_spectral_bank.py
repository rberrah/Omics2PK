#!/usr/bin/env python3
"""Run an offline, source-controlled spectral-candidate search against local references.

Usage (licensed reference JSON, no patients):
 python -m scripts.query_local_spectral_bank --db bank.sqlite --seed-literature \
   --reference-json /private/licensed_reference.json \
   --query-json /private/unknown_ei_spectrum.json --mz-tolerance-da 0.3 \
   --ri-tolerance 100

This script NEVER fetches reference libraries by itself.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from omics2pk.spectral_knowledge_bank import (
    REFERENCE_CATALOGUE, SpectralKnowledgeError, add_reference,
    make_db, search_spectral_candidates,
)
from omics2pk.spectral_biochemical_seed import seed_bibliography


def _local_json(path):
    p=Path(path)
    if not p.is_file() or p.stat().st_size>2_000_000:
        raise SpectralKnowledgeError("Local original reference/query JSON missing or over 2MB")
    value=json.loads(p.read_text(encoding="utf-8"))
    if not isinstance(value,dict):
        raise SpectralKnowledgeError("Original spectral JSON must be a mapping")
    return value


def main(argv=None):
    ap=argparse.ArgumentParser(description="Research-only spectral evidence bank, no PK claims")
    ap.add_argument("--db",help="Local SQLite path; never put licensed original spectra in GitHub")
    ap.add_argument("--seed-literature",action="store_true")
    ap.add_argument("--reference-json",help="Single approved library reference spectral record")
    ap.add_argument("--query-json",help="Locally measured unknown spectrum with ionization metadata")
    ap.add_argument("--mz-tolerance-da",type=float,default=.5)
    ap.add_argument("--ri-tolerance",type=float)
    ap.add_argument("--precursor-tolerance-da",type=float)
    ap.add_argument("--min-matched-peaks",type=int,default=6)
    ap.add_argument("--minimum-cosine",type=float,default=.7)
    ap.add_argument("--list-sources",action="store_true")
    args=ap.parse_args(argv)
    if args.list_sources:
        return {"sources":REFERENCE_CATALOGUE}
    if not args.db:
        ap.error("A local SQLite --db file is required for experiments")
    db=make_db(args.db)
    try:
        seeded=seed_bibliography(db) if args.seed_literature else None
        if args.reference_json:
            add_reference(db,_local_json(args.reference_json))
        if args.query_json:
            result=search_spectral_candidates(
                db,_local_json(args.query_json),
                mass_tolerance_da=args.mz_tolerance_da,
                retention_tolerance=args.ri_tolerance,
                precursor_tolerance_da=args.precursor_tolerance_da,
                min_matched_peaks=args.min_matched_peaks,
                minimum_cosine=args.minimum_cosine)
        else:
            result={"status":"REFERENCE_REGISTRY_UPDATED_NO_UNKNOWN_SPECTRUM_QUERIED",
                    "candidates":[]}
        return {"literature_seed":seeded,
                "n_registered_reference_spectra":db.execute("SELECT COUNT(*) FROM spectral_reference").fetchone()[0],
                "results":result,
                "patient_pk_or_biological_function_validated":False}
    finally:
        db.close()


if __name__=="__main__":
    try:
        print(json.dumps(main(),indent=2,sort_keys=True))
    except (SpectralKnowledgeError,OSError,ValueError,TypeError,json.JSONDecodeError) as exc:
        print(json.dumps({"status":"ABSTAIN_SOURCE_UNVERIFIED_OR_SPECTRAL_DATA_INSUFFICIENT",
                          "reason":type(exc).__name__,"molecular_identity_claimed":False}))
        raise SystemExit(2)
