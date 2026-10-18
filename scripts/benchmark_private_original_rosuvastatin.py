"""Private rosuvastatin E3 benchmark: no individual rows leave the local machine.

Requires authentic one-human-per-row original combined CSV, explicit features,
an externally pinned original SHA-256 and four separately checked source facts.
This is a fixed, exploratory n=40 LODO benchmark, not a clinical dose tool.
"""
from __future__ import annotations

import argparse
import json
from omics2pk.rosuvastatin_heldout import (
    ATTESTATIONS, BenchmarkGateError, evaluate_original_combined,
)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-combined", required=True, help="Local original-author CSV; never sent to GitHub")
    p.add_argument("--original-sha256", required=True, help="SHA256 pinned after independent source authentication")
    p.add_argument("--person-id-column", required=True)
    p.add_argument("--auc-column", required=True)
    p.add_argument("--cmax-column", required=True)
    p.add_argument("--omic-column", action="append", required=True, help="Predeclared quantitative predose column; repeat")
    p.add_argument("--clinical-column", action="append", default=[], help="Predeclared observed numeric baseline; repeat")
    for name in sorted(ATTESTATIONS):
        p.add_argument("--attest-" + name.replace("_", "-"), action="store_true",
                       help="Manually verified from original study source; NOT automatic")
    args = p.parse_args(argv)
    attest = {name: bool(getattr(args, "attest_" + name)) for name in ATTESTATIONS}
    try:
        result = evaluate_original_combined(
            args.original_combined,
            person_id_col=args.person_id_column,
            auc_col=args.auc_column,
            cmax_col=args.cmax_column,
            omic_columns=args.omic_column,
            clinical_columns=args.clinical_column,
            expected_sha256=args.original_sha256,
            attestations=attest,
        )
    except (BenchmarkGateError, OSError, ValueError) as exc:
        p.exit(2, "ABSTAIN: " + str(exc) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
