"""Download ORIGINAL Synapse HLC files ONLY with user's own approved access.

Use personal local Synapse account/profile or SYNAPSE_AUTH_TOKEN environment,
never submit token or original source data to a chat, public repo or CI.

Run locally:
  python -m pip install synapseclient
  python -m scripts.download_authorized_synapse_hlc
  python -m scripts.download_authorized_synapse_hlc --include-genotype

Each file goes to private/hlc_synapse (already gitignored); controlled data
remain governed by Synapse file conditions, NOT redistributable by default.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

from omics2pk.hlc_authorized_local_source import (
    SOURCES,LocalHLCSourceError,audit_local_hlc_sources
)


def main(argv=None):
    parser=argparse.ArgumentParser(description=
        "Download only legally authorized HLC files and inspect aggregate source structure")
    parser.add_argument("--out",default="private/hlc_synapse",
                        help="Private gitignored source directory, not the repository root")
    parser.add_argument("--include-genotype",action="store_true")
    args=parser.parse_args(argv)
    destination=Path(args.out).expanduser().resolve()
    try:
        import synapseclient
    except ImportError:
        raise SystemExit("Install Synapse Python SDK locally: python -m pip install synapseclient")
    # Reject root/public docs location. User explicitly opts into PRIVATE files.
    if "private" not in destination.parts:
        raise SystemExit("Source download must be inside a private directory excluded from Git")
    destination.mkdir(parents=True,exist_ok=True)
    try:
        syn=synapseclient.login(silent=True)
    except Exception:
        raise SystemExit("Local Synapse login unavailable. Sign in or configure an authorized user token locally; never upload tokens.")
    target={k:SOURCES[k] for k in (
        ("phenotype","expression","genotype") if args.include_genotype
        else ("phenotype","expression")
    )}
    paths={}
    errors={}
    for role,source_id in target.items():
        try:
            # Synapse's client enforces current account access policies.
            entity=syn.get(source_id,downloadLocation=str(destination))
            path=getattr(entity,"path",None)
            if not path:raise LocalHLCSourceError("Authorized source did not return a local file")
            paths[role]=path
        except Exception as err:
            # Exception might contain personal details or token-bearing URLs;
            # never print it. Only emit broad access/transport failure.
            errors[role]={"source":source_id,
                          "status":"DOWNLOAD_DENIED_OR_UNAVAILABLE",
                          "error_class":type(err).__name__,
                          "next_action":"In Synapse, check Download or Request Access on this exact entity"}
    result={
        "status":"PRIVATE_HLC_AUTHORIZED_DOWNLOAD_ATTEMPT",
        "requested_source_ids":target,
        "n_files_authorized_and_downloaded":len(paths),
        "blocked_sources":errors,
        "no_local_paths_original_person_values_or_credentials_emitted":True,
    }
    if paths:
        result["source_structure_qc"]=audit_local_hlc_sources(paths)
    print(json.dumps(result,indent=2,sort_keys=True))
    return 0 if not errors else 2


if __name__=="__main__":
    raise SystemExit(main())
