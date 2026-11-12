#!/usr/bin/env python
"""EXPLORATORY: is the CYP2D6 result a network, or one transcript in a panel?

    python -I scripts/run_cyp2d6_network_ablation.py \
        --hlc-expression "<dir>/curatedExpressionLiver.zip" \
        --hlc-phenotype "<dir>/curatedPhenotypesynapse.zip"

**This is not part of any frozen protocol and changes none of them.** It is run
separately, written to its own artefact and labelled EXPLORATORY throughout,
because it was designed after seeing protocol 1b's result.

Why it exists
-------------
Protocol 1b found that the 25-probe mechanistic pool beats 0 of 1000
size-matched random pools on CYP2D6 activity, and that the panel selected
inside the folds contains POR, the obligate electron donor of every microsomal
CYP. That is biologically coherent. It is **not** evidence that the network
adds anything, and the per-fold record says why: of the five folds, two chose
**k = 1**, which is the CYP2D6 probe alone. POR appears in **2 folds of 5**,
CYP2E1 in 3.

So the honest question is which of four things is true.

======================================  ====================================
Arm                                     What it would mean if it wins
======================================  ====================================
`CYP2D6_alone`                          the target transcript explains
                                        essentially the whole gain, and
                                        "pharmacological reasoning" reduces
                                        to "measure the gene you study"
`CYP2D6_plus_POR`                       the electron donor adds information
                                        beyond the target transcript
`pool_without_CYP2D6`                   the network carries signal that is
                                        independent of the target transcript
`full_pool`                             the reference, protocol 1b's arm
======================================  ====================================

`baseline_only` is run too: age and sex with no transcript at all, so every
number has a floor to be measured against.

What this cannot settle
-----------------------
It is one cohort, one endpoint, and it was specified after the result it
examines. A win for `CYP2D6_alone` would bound the claim; a win for the full
pool would not establish the network, because both arms were chosen by looking
at the same donors. **Automated expert reasoning must be able to find a
mechanism beyond the target gene without assuming that the mechanism improves
prediction**, and these four arms say which of those happened here.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omics2pk import supervised_screening_selection as SSS  # noqa: E402

ENDPOINT = "CYP2D6_activity"
TARGET_GENE = "CYP2D6"
DONOR_GENE = "POR"
FAMILY_SIZE = 6


class RunError(RuntimeError):
    pass


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc-expression", required=True)
    parser.add_argument("--hlc-phenotype", required=True)
    parser.add_argument("--out", default=str(
        ROOT / "reports" / "data" / "cyp2d6_network_ablation.json"))
    args = parser.parse_args(argv)

    sys.path.insert(0, str(ROOT / "scripts"))
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "screen_runner", ROOT / "scripts"
        / "run_supervised_screening_selection.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    from omics2pk import hlc_real_cyp_validation as HLC

    loaded = HLC.load_original_hlc(args.hlc_expression, args.hlc_phenotype)
    pools, skipped, _ = runner.hlc_fullarray_pools(
        args.hlc_expression, args.hlc_phenotype, (ENDPOINT,))
    if ENDPOINT not in pools:
        raise RunError(f"{ENDPOINT} absent: {skipped.get(ENDPOINT)}")
    pool = pools[ENDPOINT]
    design, target = pool["design"], pool["target"]
    groups, baseline = pool["groups"], pool["baseline_columns"]
    full = np.asarray(pool["knowledge_order"], dtype=int)[:pool["pool_size"]]

    def genes_of(columns):
        return sorted({runner.gene_of_probe(loaded, int(c), pool["names"])
                       for c in columns})

    target_probes = np.asarray(
        [int(c) for c in full
         if runner.gene_of_probe(loaded, int(c), pool["names"]) == TARGET_GENE],
        dtype=int)
    donor_probes = np.asarray(
        [int(c) for c in full
         if runner.gene_of_probe(loaded, int(c), pool["names"]) == DONOR_GENE],
        dtype=int)
    if target_probes.size == 0:
        raise RunError(f"no {TARGET_GENE} probe survived the completeness "
                       "filter, so this ablation has no subject")
    without_target = np.asarray([int(c) for c in full
                                 if int(c) not in set(target_probes.tolist())],
                                dtype=int)

    arms = {
        "baseline_only": np.zeros(0, dtype=int),
        "CYP2D6_alone": target_probes,
        "CYP2D6_plus_POR": np.asarray(
            sorted(set(target_probes.tolist()) | set(donor_probes.tolist())),
            dtype=int),
        "pool_without_CYP2D6": without_target,
        "full_pool": full,
    }
    results, losses = {}, {}
    for name, columns in arms.items():
        out = SSS.pool_arm_loss(
            target, design, groups, baseline_columns=baseline,
            pool_for_fold=lambda rows, fold, c=columns: c)
        losses[name] = out["losses"]
        results[name] = {
            "n_probes": int(columns.size),
            "genes": genes_of(columns),
            "mean_loss": out["mean_loss"],
            "k_chosen_per_fold": out["k_chosen_per_fold"],
            "n_folds_that_abstained": out["n_folds_that_abstained"],
            "selected_genes_per_fold": [genes_of(sel)
                                        for sel in out["selected_per_fold"]],
        }
        print(f"  {name:<22} n={columns.size:>2}  loss {out['mean_loss']:.4f}  "
              f"k={out['k_chosen_per_fold']}", flush=True)

    def paired(better, worse):
        shell = {"losses": losses,
                 "scored_rows": np.ones(target.size, dtype=bool)}
        # Every arm scored the same rows: the folds and the cohort are fixed
        # and only the candidate pool differs.
        sizes = {name: v.size for name, v in losses.items()}
        if len(set(sizes.values())) != 1:
            return {"verdict": "NOT_EVALUABLE",
                    "reason": f"arms scored different row counts: {sizes}"}
        shell["scored_rows"] = np.ones(losses[better].size, dtype=bool)
        return SSS.contrast(shell, better, worse, groups=None,
                            family_size=FAMILY_SIZE)

    contrasts = {
        # The question that bounds the claim.
        "full_pool_vs_CYP2D6_alone": paired("full_pool", "CYP2D6_alone"),
        "CYP2D6_plus_POR_vs_CYP2D6_alone": paired("CYP2D6_plus_POR",
                                                  "CYP2D6_alone"),
        "full_pool_vs_CYP2D6_plus_POR": paired("full_pool",
                                               "CYP2D6_plus_POR"),
        # Does the network carry anything without the target transcript?
        "pool_without_CYP2D6_vs_baseline_only": paired("pool_without_CYP2D6",
                                                       "baseline_only"),
        "CYP2D6_alone_vs_baseline_only": paired("CYP2D6_alone",
                                                "baseline_only"),
        "full_pool_vs_baseline_only": paired("full_pool", "baseline_only"),
    }
    non_inferiority = {
        # If the single transcript is non-inferior to the whole pool, the
        # claim is "measure the gene you study", not "automate a network".
        "CYP2D6_alone_non_inferior_to_full_pool": SSS.non_inferiority(
            {"losses": losses,
             "scored_rows": np.ones(losses["full_pool"].size, dtype=bool)},
            "CYP2D6_alone", "full_pool", family_size=FAMILY_SIZE),
    }

    payload = {
        "run": "cyp2d6_network_ablation",
        "generated": "2026-11-11",
        "status": "EXPLORATORY",
        "why_exploratory": (
            "specified after seeing protocol 1b's CYP2D6 result, on the same "
            "cohort and the same donors. It bounds the claim; it cannot "
            "confirm it."),
        "not_part_of_any_frozen_protocol": True,
        "endpoint": ENDPOINT,
        "question": ("is the CYP2D6 result a network, or one transcript "
                     "travelling in a panel?"),
        "from_protocol_1b": {
            "pool_beats_0_of_1000_matched_pools": True,
            "folds_that_chose_k_equal_1": 2,
            "folds_containing_POR": 2,
            "folds_containing_CYP2D6": 5,
            "why_this_ablation_follows": (
                "two folds of five selected the CYP2D6 probe alone, so the "
                "panel is a single probe in 40 % of the folds and POR appears "
                "in 2 of 5. Biological coherence is not evidence of added "
                "information."),
        },
        "arms": results,
        "contrasts": contrasts,
        "non_inferiority": non_inferiority,
        "what_this_cannot_settle": (
            "one cohort, one endpoint, specified after the result it examines. "
            "A win for CYP2D6_alone bounds the claim; a win for the full pool "
            "does not establish the network, because both arms were chosen by "
            "looking at the same donors."),
        "no_individual_level_data_in_this_artefact": True,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              default=float) + "\n", encoding="utf-8")
    print(f"\nartefact: {out}")
    print("\n=== contrasts ===")
    for name, value in contrasts.items():
        if value.get("verdict") == "NOT_EVALUABLE":
            print(f"  {name:<44} NOT_EVALUABLE {value['reason']}")
            continue
        print(f"  {name:<44} {value['relative_gain']:+7.2%} "
              f"[{value['interval'][0]:+.4f},{value['interval'][1]:+.4f}] "
              f"{value['verdict']}")
    for name, value in non_inferiority.items():
        print(f"  {name:<44} {value['relative_gap']:+7.2%} {value['verdict']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
