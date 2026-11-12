#!/usr/bin/env python
"""The frozen Latham protocol, run once.

    python -I scripts/run_latham_imatinib_protocol.py

Protocol: `docs/OMICS2PK_LATHAM_IMATINIB_PREREGISTRATION_20261112.md` v1.2.0.
Nothing here chooses a comparator, a margin or a model: all three were frozen
before the values were read, and the extraction bug that produced a false
claim about the source was fixed and recorded before this ran.

Why a bespoke loop rather than `pool_arm_loss`
----------------------------------------------
The protocol requires everything learned to be learned **inside the fold**,
including the centring and scaling, and at thirteen donors with singleton PGx
categories that is not a formality. A fixed design matrix cannot express it,
so the design is rebuilt for every fold from the training rows alone.

No per-donor value enters the artefact. PMC permits text mining and otherwise
fair use, so the numbers may be analysed and not redistributed.
"""
from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omics2pk import latham_imatinib_cohort as LATHAM  # noqa: E402
from omics2pk.fusion_headroom_diagnostics import (  # noqa: E402
    paired_gap, ridge_predict)
from omics2pk.variant_mechanism_omic_selection import _ladder  # noqa: E402

TABLE = pathlib.Path(r"C:\Users\abdel\Documents\Omics"
                     r"\Latham2025_imatinib_PHH"
                     r"\latham2025_table4_phh_per_donor.csv")

#: Frozen in the protocol, repeated here so a reader of the artefact sees them.
FROZEN = {
    "margin": 0.05,
    "family_size": 4,
    "ridge_penalty": 1.0,
    "seed": 20261112,
    "k_grid": (0, 1, 2, 3),
    "validation": "leave-one-donor-out, outer and inner",
    "loss": "absolute error",
    "resampling_unit": ("the donor, which IS the row here, so the cluster and "
                        "row bootstraps coincide"),
}

#: The declared endpoints. CYP3A activity is midazolam; testosterone is a
#: sensitivity and the two probes are not two independent validations.
ENDPOINTS = {
    "primary_clearance": ("clearance", "H1-in-vitro"),
    "h3_activity_cyp2c8": ("activity_cyp2c8", "H3"),
    "h3_activity_cyp3a_midazolam": ("activity_cyp3a_midazolam", "H3"),
    "secondary_metabolite": ("metabolite", "secondary"),
    "sensitivity_activity_cyp3a_testosterone": (
        "activity_cyp3a_testosterone", "sensitivity"),
}


def _rank(values):
    order = np.asarray(values, dtype=float).argsort(kind="mergesort")
    out = np.empty(order.size, dtype=float)
    out[order] = np.arange(1, order.size + 1, dtype=float)
    return out


def _abs_spearman(columns, target):
    ranked = np.column_stack([_rank(columns[:, j])
                              for j in range(columns.shape[1])])
    rt = _rank(target)
    centred = ranked - ranked.mean(axis=0)
    rt = rt - rt.mean()
    scale = np.sqrt((centred ** 2).sum(axis=0)) * np.sqrt((rt ** 2).sum())
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.abs(centred.T @ rt) / scale
    return np.nan_to_num(out, nan=0.0)


def run_arm(loaded, rows, endpoint_role, *, baseline_blocks, pool_block):
    """Leave-one-donor-out loss, with the design rebuilt inside every fold."""
    parsed = loaded["parsed"]
    target = np.array([parsed[endpoint_role][i] for i in rows], dtype=float)
    target = np.log(target)
    n = len(rows)
    scored = np.full(n, np.nan)
    budgets, chosen = [], []

    def design_for(train_positions):
        train_rows = [rows[p] for p in train_positions]
        blocks = LATHAM.encode_blocks(loaded, rows, train_rows=train_rows)
        baseline = np.column_stack([blocks[b] for b in baseline_blocks])
        pool = blocks[pool_block] if pool_block else np.zeros((n, 0))
        return baseline, pool

    for held in range(n):
        train = [p for p in range(n) if p != held]
        baseline, pool = design_for(train)
        # inner leave-one-out on the training donors picks k
        inner_loss = {}
        for k in FROZEN["k_grid"]:
            if k > pool.shape[1]:
                continue
            losses = []
            for inner_held in train:
                inner_train = [p for p in train if p != inner_held]
                b_in, p_in = design_for(inner_train)
                if k and p_in.shape[1]:
                    strength = _abs_spearman(p_in[inner_train],
                                             target[inner_train])
                    keep = strength.argsort()[::-1][:k]
                    columns = np.column_stack([b_in, p_in[:, keep]])
                else:
                    columns = b_in
                prediction = ridge_predict(
                    columns[inner_train], target[inner_train],
                    columns[[inner_held]], penalty=FROZEN["ridge_penalty"])
                losses.append(abs(float(prediction[0]) - target[inner_held]))
            if losses:
                inner_loss[k] = float(np.mean(losses))
        best_k = min(sorted(inner_loss), key=lambda k: (inner_loss[k], k))
        if best_k and pool.shape[1]:
            strength = _abs_spearman(pool[train], target[train])
            keep = strength.argsort()[::-1][:best_k]
            columns = np.column_stack([baseline, pool[:, keep]])
        else:
            keep = np.zeros(0, dtype=int)
            columns = baseline
        scored[held] = ridge_predict(columns[train], target[train],
                                     columns[[held]],
                                     penalty=FROZEN["ridge_penalty"])[0]
        budgets.append(int(best_k))
        chosen.append([int(c) for c in keep])
    losses = np.abs(scored - target)
    return {"losses": losses, "mean_loss": float(losses.mean()),
            "k_chosen_per_fold": budgets,
            "n_folds_that_abstained": int(sum(1 for k in budgets if k == 0)),
            "protein_columns_chosen_per_fold": chosen}


def contrast(richer, plainer, *, label, n):
    gap = paired_gap(plainer["losses"], richer["losses"],
                     family_size=FROZEN["family_size"], n_boot=2000,
                     alpha=0.05, seed=FROZEN["seed"])
    reference = float(np.mean(plainer["losses"]))
    threshold = FROZEN["margin"] * reference
    return {"contrast": label, "n": n,
            "verdict": _ladder(gap["interval"], threshold),
            "relative_gain": float(gap["gap"] / reference),
            "minimum_interesting_gain": float(threshold),
            "reference_mean_loss": reference,
            "richer_mean_loss": float(np.mean(richer["losses"])), **gap}


def main(argv=None):
    if not TABLE.is_file():
        raise SystemExit(f"per-donor table absent: {TABLE}")
    loaded = LATHAM.load_table4(TABLE)
    rows, donors = LATHAM.primary_cohort(loaded)
    print(f"primary cohort n={len(rows)}: {', '.join(donors)}", flush=True)

    results, influence = {}, {}
    for name, (role, hypothesis) in ENDPOINTS.items():
        baseline_blocks = (("clinical", "pgx") if hypothesis == "H1-in-vitro"
                           or hypothesis == "secondary" else ("pgx",))
        plain = run_arm(loaded, rows, role, baseline_blocks=baseline_blocks,
                        pool_block=None)
        rich = run_arm(loaded, rows, role, baseline_blocks=baseline_blocks,
                       pool_block="proteomics")
        out = contrast(rich, plain, label=f"{name}_proteomics_added",
                       n=len(rows))
        out["hypothesis"] = hypothesis
        out["baseline_blocks"] = list(baseline_blocks)
        # Written as text: this project forbids a numeric list of ten or more
        # entries in a committed artefact, and there are thirteen folds. The
        # per-k tally is added because it is what a reader actually wants.
        out["k_chosen_per_fold"] = ", ".join(
            str(k) for k in rich["k_chosen_per_fold"])
        out["k_tally"] = {str(k): rich["k_chosen_per_fold"].count(k)
                          for k in sorted(set(rich["k_chosen_per_fold"]))}
        out["n_folds_that_abstained"] = rich["n_folds_that_abstained"]
        results[name] = out
        print(f"  {name:<44} {out['relative_gain']:+7.2%} "
              f"[{out['interval'][0]:+.4f},{out['interval'][1]:+.4f}] "
              f"{out['verdict']}  k={rich['k_chosen_per_fold']}", flush=True)

        # mandatory influence analysis, one re-run per donor dropped
        verdicts = {}
        for drop in range(len(rows)):
            subset = [r for j, r in enumerate(rows) if j != drop]
            p2 = run_arm(loaded, subset, role,
                         baseline_blocks=baseline_blocks, pool_block=None)
            r2 = run_arm(loaded, subset, role,
                         baseline_blocks=baseline_blocks,
                         pool_block="proteomics")
            verdicts[donors[drop]] = contrast(
                r2, p2, label=name, n=len(subset))["verdict"]
        influence[name] = {
            "verdict_with_every_donor": out["verdict"],
            "verdict_per_donor_dropped": verdicts,
            "n_exclusions_that_change_the_verdict": sum(
                1 for v in verdicts.values() if v != out["verdict"]),
            "reading": ("if any exclusion changes the verdict, the contrast "
                        "is reported INCONCLUSIVE rather than positive with a "
                        "caveat, as declared in the protocol"),
        }
        changed = influence[name]["n_exclusions_that_change_the_verdict"]
        print(f"       influence: {changed} of {len(rows)} exclusions change "
              f"the verdict", flush=True)

    payload = {
        "run": "latham_imatinib_protocol",
        "generated": "2026-11-12",
        "protocol": ("OMICS2PK_LATHAM_IMATINIB_PREREGISTRATION_20261112.md "
                     "v1.2.0"),
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in FROZEN.items()},
        "cohort": {"n": len(rows), "donors": donors,
                   "excluded_by_the_data": loaded["completeness"][
                       "proteomics"]["donors_missing"],
                   "completeness_per_block": {
                       k: v["n_complete"]
                       for k, v in loaded["completeness"].items()}},
        "endpoint_is_not_a_pk_parameter": (
            "hepatocyte intrinsic clearance, measured with no absorption, no "
            "plasma binding, no hepatic blood flow and no renal excretion. "
            "Nothing here is extrapolated to CL, F, V or an AUC in patients."),
        "h2_status": ("NOT_EVALUABLE: three quantified CYP proteins against a "
                      "floor of five independent candidates, set before this "
                      "dataset existed"),
        "contrasts": results,
        "influence_analysis": influence,
        "what_this_is_not": LATHAM.what_this_module_does_not_do(),
        "no_individual_level_data_in_this_artefact": True,
    }
    out = ROOT / "reports" / "data" / "latham_imatinib_protocol.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              default=float) + "\n", encoding="utf-8")
    print(f"\nartefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
