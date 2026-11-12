#!/usr/bin/env python
"""Protocol 3's one prediction that can fail on its own: HSP90 inhibition.

    python -I scripts/run_araf_hsp90_proxy.py \
        --mavedb-dir "<dir>"

§3.2 of `OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md` predicted that
**abundance under HSP90 inhibition becomes a better proxy for activity than
abundance under vehicle**. It is the only ARAF prediction that could fail while
every other one held, which is why it is worth running: the BRAF refutation was
attributed to a gain-of-function mechanism in which abundance is an inverted
proxy, and activating kinase variants that are destabilised depend on HSP90. If
this fails while the rest holds, that explanation is wrong and the BRAF failure
goes back to being unexplained.

All three score sets are the **same C-terminal construct and the same variant
library**, so the two abundance arms differ only in the condition measured:
DMSO vehicle against HSP90 inhibitor. Folds are grouped by protein position and
intervals resample positions, as everywhere else in this project.

Also reported, because §3.2 predicted it too: whether this library is a
gain-of-function screen, read off the mean activity of missense against
synonymous variants. That is a property of the assay, not of any model.
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
from omics2pk import variant_mechanism_omic_selection as VMS  # noqa: E402

URNS = {
    "activity": "00001314-c-1",
    "abundance_dmso": "00001314-a-1",
    "abundance_hsp90i": "00001314-b-1",
}


class RunError(RuntimeError):
    pass


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mavedb-dir", required=True)
    parser.add_argument("--out", default=str(
        ROOT / "reports" / "data" / "araf_hsp90_proxy.json"))
    args = parser.parse_args(argv)

    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "variant_bench", ROOT / "scripts" / "run_variant_mechanism_benchmark.py")
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)

    directory = pathlib.Path(args.mavedb_dir)
    tables = {}
    for role, urn in URNS.items():
        path = directory / f"mavedb_{urn}_scores.csv"
        if not path.is_file():
            raise RunError(f"missing declared score file: {path}")
        tables[role] = bench.load_scores(path)

    shared = sorted(set(tables["activity"]) & set(tables["abundance_dmso"])
                    & set(tables["abundance_hsp90i"]))
    if len(shared) < 200:
        return _not_evaluable(args.out, len(shared))

    parsed = [tables["activity"][v][0] for v in shared]
    activity = np.array([tables["activity"][v][1] for v in shared], dtype=float)
    dmso = np.array([tables["abundance_dmso"][v][1] for v in shared],
                    dtype=float)
    hsp90i = np.array([tables["abundance_hsp90i"][v][1] for v in shared],
                      dtype=float)
    positions = np.array([p[1] for p in parsed], dtype=int)
    kind = np.array(["nonsense" if p[2] is None
                     else "synonymous" if p[2] == "=" else "missense"
                     for p in parsed])

    # ---- is this a gain-of-function library? A property of the assay.
    polarity = {
        role: {
            label: float(np.mean(values[kind == label]))
            for label in ("nonsense", "synonymous", "missense")
            if (kind == label).any()}
        for role, values in (("activity", activity), ("abundance_dmso", dmso),
                             ("abundance_hsp90i", hsp90i))}
    polarity["missense_above_synonymous_in_activity"] = bool(
        polarity["activity"]["missense"] > polarity["activity"]["synonymous"])
    polarity["how_to_read"] = (
        "nonsense lowest in every assay means the scales run the same way. "
        "Missense ABOVE synonymous in activity means the average substitution "
        "makes the kinase more active than wild type, which is the signature "
        "of a gain-of-function screen on an autoinhibited kinase.")

    # ---- keep missense and nonsense, drop synonymous, as elsewhere
    keep = kind != "synonymous"
    activity, dmso, hsp90i = activity[keep], dmso[keep], hsp90i[keep]
    positions = positions[keep]

    def rank(values):
        order = values.argsort(kind="mergesort")
        out = np.empty_like(order, dtype=float)
        out[order] = np.arange(1, values.size + 1, dtype=float)
        return out

    def spearman(a, b):
        ra, rb = rank(a), rank(b)
        ra, rb = ra - ra.mean(), rb - rb.mean()
        return float((ra @ rb) / np.sqrt((ra @ ra) * (rb @ rb)))

    correlations = {
        "dmso_vs_activity": spearman(dmso, activity),
        "hsp90i_vs_activity": spearman(hsp90i, activity),
        "dmso_vs_hsp90i": spearman(dmso, hsp90i),
    }

    # ---- the declared comparison: which abundance predicts activity better?
    folds = VMS.position_groups(positions, VMS.FROZEN["position_folds"],
                                VMS.FROZEN["seed"])
    arms = {}
    for name, predictor in (("abundance_dmso", dmso),
                            ("abundance_hsp90i", hsp90i),
                            ("both_abundances", None)):
        design = (np.column_stack([dmso, hsp90i]) if predictor is None
                  else predictor.reshape(-1, 1))
        out = SSS.pool_arm_loss(
            activity, design, positions, baseline_columns=(),
            pool_for_fold=lambda rows, fold, n=design.shape[1]: np.arange(n),
            k_grid=tuple(k for k in (0, 1, 2) if k <= design.shape[1]))
        arms[name] = out
        print(f"  {name:<20} loss {out['mean_loss']:.4f}  "
              f"k={out['k_chosen_per_fold']}", flush=True)

    def paired(better, worse):
        shell = {"losses": {k: v["losses"] for k, v in arms.items()},
                 "scored_rows": np.asarray(arms[better]["scored_rows"],
                                           dtype=bool)}
        if arms[better]["losses"].size != arms[worse]["losses"].size:
            return {"verdict": "NOT_EVALUABLE",
                    "reason": "the two arms scored different rows"}
        return SSS.contrast(shell, better, worse,
                            groups=positions, family_size=3)

    contrasts = {
        # THE declared prediction.
        "hsp90i_better_than_dmso": paired("abundance_hsp90i",
                                          "abundance_dmso"),
        "both_better_than_dmso": paired("both_abundances", "abundance_dmso"),
        "both_better_than_hsp90i": paired("both_abundances",
                                          "abundance_hsp90i"),
    }

    payload = {
        "run": "araf_hsp90_proxy",
        "generated": "2026-11-11",
        "protocol": ("OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md "
                     "§3.2, the one prediction that can fail alone"),
        "declared_urns": URNS,
        "same_construct_and_library": (
            "all three sets are the C-terminal construct on one variant "
            "library, so the two abundance arms differ only in the condition "
            "measured: DMSO vehicle against HSP90 inhibitor."),
        "n_shared_all_three": len(shared),
        "n_analysed_without_synonymous": int(activity.size),
        "n_positions": int(np.unique(positions).size),
        "assay_polarity": polarity,
        "spearman": correlations,
        "mean_loss": {k: v["mean_loss"] for k, v in arms.items()},
        "k_chosen_per_fold": {k: v["k_chosen_per_fold"]
                              for k, v in arms.items()},
        "contrasts": contrasts,
        "prediction": (
            "§3.2 predicted that abundance under HSP90 inhibition becomes a "
            "BETTER proxy for activity than abundance under vehicle. If this "
            "fails while the other ARAF predictions hold, the HSP90 "
            "explanation of the BRAF refutation is wrong and that failure is "
            "unexplained again."),
        "no_individual_level_data_in_this_artefact": True,
    }
    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              default=float) + "\n", encoding="utf-8")
    print(f"\nartefact: {out}")
    print("\nassay polarity (mean score):")
    for role in ("activity", "abundance_dmso", "abundance_hsp90i"):
        print("   %-18s %s" % (role, {k: round(v, 3)
                                      for k, v in polarity[role].items()}))
    print("   missense above synonymous in activity:",
          polarity["missense_above_synonymous_in_activity"])
    print("\nSpearman:", {k: round(v, 4) for k, v in correlations.items()})
    print("\ncontrasts:")
    for name, value in contrasts.items():
        if value.get("verdict") == "NOT_EVALUABLE":
            print(f"   {name:<28} NOT_EVALUABLE")
            continue
        print(f"   {name:<28} {value['relative_gain']:+7.2%} "
              f"[{value['interval'][0]:+.4f},{value['interval'][1]:+.4f}] "
              f"{value['verdict']}")
    return 0


def _not_evaluable(path, n):
    payload = {"run": "araf_hsp90_proxy", "verdict": "NOT_EVALUABLE",
               "reason": f"only {n} variants are present in all three sets"}
    out = pathlib.Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"NOT_EVALUABLE: only {n} shared variants")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
