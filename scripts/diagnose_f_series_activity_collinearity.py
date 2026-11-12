#!/usr/bin/env python
"""Post-hoc: are individual ridge coefficients on activity columns identified?

**This is a diagnostic, added after the F-series run and labelled as such.** It
changes no arm, no loss, no threshold and no contrast, and it produces no new
performance number. It exists because the run's biological-relevance block
reported that the fitted coefficient sign agrees with CollecTRI's recorded sign
for only 15 of 49 regulators -- 30.6 %, two-sided binomial p = 0.009, which is
*below* chance rather than above it.

That number has two incompatible readings:

1. the resource's signs are wrong, or
2. the coefficient signs are not identified, because the activity columns are
   near-duplicates of each other and a ridge distributes weight among
   collinear columns with signs that carry no individual meaning.

Reading 2 is a measurable claim and this script measures it: the pairwise
correlation of the activity columns, the condition number of the standardised
design, and the effective rank. If the columns are strongly collinear, the
30.6 % says nothing about the biology and must not be reported as though it
did.

The regulon overlap is reported alongside, because it is the mechanism that
would produce the collinearity: CollecTRI regulons here hold a median of 195
admissible members drawn from the same 21,580-gene universe, so two regulons
can share most of their members.

Run from the repository root with the same sources as the main runner:

    python scripts/diagnose_f_series_activity_collinearity.py \
        --hlc "<HLC directory>" --network "<collectri.tsv>" \
        --identity "<probe_identity_full_table.json>" \
        --out reports/data/f_series_activity_collinearity.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from datetime import date
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk import functional_omics as FO  # noqa: E402


def load_runner():
    path = REPO_ROOT / "scripts" / "run_f_series_hlc_functional_omics.py"
    spec = importlib.util.spec_from_file_location("_f_series_runner", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--network", required=True)
    parser.add_argument("--identity", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--fold", type=int, default=0)
    args = parser.parse_args(argv)

    runner = load_runner()
    import pandas as pd

    identity = json.loads(Path(args.identity).read_text(encoding="utf-8"))
    resolved = {str(k): str((v or {}).get("human_symbol") or "")
                for k, v in identity.items()}
    data = runner.load_hlc(Path(args.hlc))
    data["resolved"] = resolved
    frame, _ = runner.load_network(Path(args.network))
    universes = runner.build_universes(data)

    signed = FO.resolve_signed_network(frame)
    unsigned = FO.resolve_signed_network(frame, drop_cytochrome_targets=False)
    regulons = signed.regulons()
    enzyme = FO.ENDPOINT_ENZYME[FO.PRIMARY_ENDPOINT]
    candidates = unsigned.regulators_of(enzyme)
    evaluable, _ = FO.evaluable_regulators(
        regulons, candidates, set(universes["scoring"]))

    folds = runner.H2A.repeated_folds(
        len(data["donors"]), n_splits=FO.FROZEN["outer_folds"],
        n_repeats=FO.FROZEN["outer_repeats"], seed=FO.FROZEN["seed"])
    train, _test = folds[args.fold]

    chosen = runner.choose_probes(data["values"], universes["scoring"], train)
    genes, block = runner.assemble(
        data["values"], chosen, np.arange(len(data["donors"])), train)
    standardised = FO.standardise_on_training(block, train)
    weights = FO.weight_matrix(regulons, evaluable, genes)
    activities = FO.ulm_activities(standardised, weights)

    # Correlation of the activity columns on the training donors of this fold.
    scores = activities[train]
    centred = scores - scores.mean(axis=0, keepdims=True)
    spread = centred.std(axis=0)
    spread = np.where(spread > 0.0, spread, 1.0)
    z = centred / spread
    correlation = (z.T @ z) / z.shape[0]
    off = correlation[~np.eye(len(evaluable), dtype=bool)]

    singular = np.linalg.svd(z / np.sqrt(z.shape[0]), compute_uv=False)
    normalised = singular ** 2 / (singular ** 2).sum()
    entropy = float(-np.sum(normalised * np.log(normalised + 1e-300)))

    # Regulon overlap: the mechanism that would produce the collinearity.
    members = {r: set(g for g in regulons[r] if g in universes["scoring"])
               for r in evaluable}
    jaccard = []
    for i, a in enumerate(evaluable):
        for b in evaluable[i + 1:]:
            union = members[a] | members[b]
            if union:
                jaccard.append(len(members[a] & members[b]) / len(union))
    jaccard = np.asarray(jaccard, dtype=np.float64)

    report = {
        "report": "post-hoc: identifiability of ridge coefficients on F-series activity columns",
        "label": "POST_HOC_DIAGNOSTIC",
        "changes_no_result": (
            "No arm, loss, threshold or contrast is affected. This script "
            "produces no performance number. It exists to decide whether the "
            "run's coefficient-sign agreement can be read as a biological "
            "statement."),
        "built_on": date.today().isoformat(),
        "endpoint": FO.PRIMARY_ENDPOINT,
        "fold": int(args.fold),
        "n_training_donors": int(train.size),
        "n_activity_columns": len(evaluable),
        "activity_column_correlation": {
            "max_absolute_offdiagonal": float(np.abs(off).max()),
            "median_absolute_offdiagonal": float(np.median(np.abs(off))),
            "fraction_above_0_9": float((np.abs(off) > 0.9).mean()),
            "fraction_above_0_7": float((np.abs(off) > 0.7).mean()),
        },
        "design": {
            "condition_number": float(singular[0] / singular[-1]),
            "effective_rank_shannon": round(float(np.exp(entropy)), 2),
            "variance_in_first_component": float(normalised[0]),
            "components_for_95_percent_variance": int(
                np.searchsorted(np.cumsum(normalised), 0.95) + 1),
        },
        "regulon_overlap_jaccard": {
            "pairs": int(jaccard.size),
            "median": float(np.median(jaccard)),
            "max": float(jaccard.max()),
            "fraction_above_0_5": float((jaccard > 0.5).mean()),
        },
        "how_to_read": (
            "A high condition number with an effective rank far below the "
            "column count means the individual coefficients are not identified: "
            "ridge distributes weight among near-duplicate columns and the sign "
            "of any one of them is arbitrary. In that case the run's "
            "coefficient-sign agreement is a statement about the design, not "
            "about the resource's biology, and must not be reported as evidence "
            "that the signs are wrong."),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
