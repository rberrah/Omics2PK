#!/usr/bin/env python
"""How well does mRNA explain CYP activity, against a published protein benchmark?

**Descriptive, cross-cohort, and labelled as such.** This script computes one
prespecified statistic -- the Spearman correlation between a cytochrome's own
admissible transcript and its measured microsomal activity, in the Human Liver
Cohort -- for every endpoint, with no model, no selection and no tuning. It then
places it beside the *published* protein-to-activity Spearman correlations of
Grangeon 2025 (Pharmaceuticals 18(12):1789, CC BY 4.0), which were measured in
**different donors**.

Why this comparison is worth making, and what it cannot be:

* Grangeon measured protein abundance and probe-substrate activity in the same
  16 genotyped human liver microsome donors and published a Spearman r per
  enzyme. Those are whole-cohort descriptive correlations.
* The Human Liver Cohort measured transcript and the same probe-substrate
  activities in the same ~390 donors. The same descriptive correlation is
  computable here.
* **The donors are not the same people.** This is therefore a comparison of two
  published descriptive statistics from two cohorts, not a within-subject
  contrast of two layers. It cannot license a causal or per-patient claim, and
  no individual from one cohort is ever joined to the other.
* Grangeon's n is 16 and 20. A Spearman r at n = 16 has a wide interval, so the
  comparison is reported with a bootstrap interval for the side this project
  controls and with the published point estimate for the side it does not.

The one thing it does establish, if the gap is large and consistent: how much
of the enzyme-activity variance an mRNA measurement leaves on the table that a
protein measurement captures. That is the quantity a funding case for paired
proteomics has to state.

Run from the repository root:

    python scripts/compare_rna_versus_published_protein_activity_correlation.py \
        --hlc "<HLC directory>" \
        --out reports/data/rna_versus_protein_activity_correlation.json
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path
from zipfile import ZipFile

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

POOL_PATH = REPO_ROOT / "registry" / "QUALIFIED_PROBE_POOLS.json"
POOL_NAME = "A_qualified_human"
N_BOOTSTRAP = 2000
SEED = 20261109

ENDPOINT_ENZYME = {
    "CYP1A2_activity": "CYP1A2",
    "CYP2A6_activity": "CYP2A6",
    "CYP2B6_activity": "CYP2B6",
    "CYP2C8_activity": "CYP2C8",
    "CYP2C9_activity": "CYP2C9",
    "CYP2C19_activity": "CYP2C19",
    "CYP2D6_activity": "CYP2D6",
    "CYP2E1_activity": "CYP2E1",
    "CYP3A4M_activity": "CYP3A4",
    "CYP3A4T_activity": "CYP3A4",
}

# Grangeon 2025, Pharmaceuticals 18(12):1789, Table 3B (genotyped human liver
# microsomes) and Table 3A (individual human hepatocytes). Spearman r between
# absolute-quantified protein content and probe-substrate activity, in the SAME
# donors. Transcribed from the published tables; CC BY 4.0. These are the
# published numbers and are not recomputed here.
GRANGEON_PROTEIN_ACTIVITY = {
    "CYP1A2_activity": {"probe": "phenacetin", "hlm_r": 0.6814, "hlm_p": "0.0046",
                        "hlm_n": 16, "hhep_r": 0.7880, "hhep_p": "<0.0001",
                        "hhep_n": 20},
    "CYP2A6_activity": {"probe": "coumarin", "hlm_r": 0.9471, "hlm_p": "<0.0001",
                        "hlm_n": 16, "hhep_r": 0.7985, "hhep_p": "<0.0001",
                        "hhep_n": 20},
    "CYP2B6_activity": {"probe": "bupropion", "hlm_r": 0.9794, "hlm_p": "<0.0001",
                        "hlm_n": 16, "hhep_r": 0.7424, "hhep_p": "0.0002",
                        "hhep_n": 20},
    "CYP2C8_activity": {"probe": "amodiaquine", "hlm_r": 0.9176,
                        "hlm_p": "<0.0001", "hlm_n": 16, "hhep_r": 0.7574,
                        "hhep_p": "0.0001", "hhep_n": 20},
    "CYP2C9_activity": {"probe": "diclofenac", "hlm_r": 0.6500, "hlm_p": "0.0078",
                        "hlm_n": 16, "hhep_r": 0.3805, "hhep_p": "0.098",
                        "hhep_n": 20},
    "CYP2C19_activity": {"probe": "mephenytoin", "hlm_r": 0.9500,
                         "hlm_p": "<0.0001", "hlm_n": 16, "hhep_r": 0.7003,
                         "hhep_p": "0.0006", "hhep_n": 20},
    "CYP2D6_activity": {"probe": "dextromethorphan", "hlm_r": 0.8291,
                        "hlm_p": "0.0002", "hlm_n": 16, "hhep_r": 0.8015,
                        "hhep_p": "<0.0001", "hhep_n": 20},
    "CYP2E1_activity": {"probe": "chlorzoxazone", "hlm_r": 0.9441,
                        "hlm_p": "<0.0001", "hlm_n": 16, "hhep_r": 0.2915,
                        "hhep_p": "0.2124", "hhep_n": 20},
    "CYP3A4M_activity": {"probe": "midazolam", "hlm_r": 0.8168, "hlm_p": "0.0002",
                         "hlm_n": 16, "hhep_r": 0.5353, "hhep_p": "0.015",
                         "hhep_n": 20},
    "CYP3A4T_activity": {"probe": "testosterone", "hlm_r": 0.9559,
                         "hlm_p": "<0.0001", "hlm_n": 16, "hhep_r": 0.7684,
                         "hhep_p": "<0.0001", "hhep_n": 20},
}

GRANGEON_CITATION = (
    "Grangeon A, et al. Protein Levels of 16 Cytochrome P450s and 2 Carboxyl "
    "Esterases Using Absolute Quantitative Proteomics. Pharmaceuticals 2025, "
    "18(12), 1789. PMID 41471278. CC BY 4.0. Table 3A/3B.")

EXPECTED_PROBES = 40638
EXPECTED_DONORS = 466
AGE_ROW = "AGE_(YRS)"


class CompareError(RuntimeError):
    pass


def spearman(x, y):
    """Spearman's rho via Pearson on ranks, average ranks for ties."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    if x.size < 3:
        return float("nan")
    rx, ry = _rank(x), _rank(y)
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    denominator = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / denominator) if denominator > 0 else float("nan")


def _rank(values):
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    ranks[order] = np.arange(1, values.size + 1, dtype=np.float64)
    # average ranks within ties
    sorted_values = values[order]
    start = 0
    for index in range(1, values.size + 1):
        if index == values.size or sorted_values[index] != sorted_values[start]:
            if index - start > 1:
                ranks[order[start:index]] = ranks[order[start:index]].mean()
            start = index
    return ranks


def bootstrap_interval(x, y, *, seed=SEED, draws=N_BOOTSTRAP):
    rng = np.random.default_rng(seed)
    n = x.size
    out = np.empty(draws, dtype=np.float64)
    for index in range(draws):
        pick = rng.integers(0, n, size=n)
        out[index] = spearman(x[pick], y[pick])
    low, high = np.nanpercentile(out, [2.5, 97.5])
    return float(low), float(high)


def load_hlc(hlc_dir: Path):
    import pandas as pd

    expression_zip = hlc_dir / "curatedExpressionLiver.zip"
    phenotype_zip = hlc_dir / "curatedPhenotypesynapse.zip"
    for path in (expression_zip, phenotype_zip):
        if not path.is_file():
            raise CompareError(f"expected archive absent: {path}")
    with ZipFile(expression_zip) as archive:
        with archive.open("features.txt") as handle:
            features = pd.read_csv(handle, sep="\t", dtype=str)
        with archive.open("expression.txt") as handle:
            expression = pd.read_csv(handle, sep="\t", index_col=0)
    with ZipFile(phenotype_zip) as archive:
        with archive.open("phenotype.txt") as handle:
            phenotype = pd.read_csv(handle, sep="\t", index_col=0, dtype=str)
    if expression.shape != (EXPECTED_PROBES, EXPECTED_DONORS):
        raise CompareError(f"expression dimensions changed: {expression.shape}")
    donors = expression.columns.astype(str).to_numpy()
    phenotype = phenotype.loc[:, donors]
    return {"values": expression.to_numpy(dtype=np.float64),
            "symbols": features["genesymbol"].fillna("").to_numpy(dtype=str),
            "feature_ids": features["feature_id"].astype(str).to_numpy(dtype=str),
            "donors": donors, "phenotype": phenotype}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    import pandas as pd

    pool = json.loads(POOL_PATH.read_text(encoding="utf-8"))
    admissible = set(map(str, pool["admissible_features"][POOL_NAME]))
    data = load_hlc(Path(args.hlc))

    by_symbol: dict[str, list[int]] = {}
    for position, identifier in enumerate(data["feature_ids"]):
        symbol = str(data["symbols"][position])
        if symbol and str(identifier) in admissible:
            by_symbol.setdefault(symbol, []).append(position)

    rows = {}
    for endpoint, enzyme in ENDPOINT_ENZYME.items():
        activity = pd.to_numeric(data["phenotype"].loc[endpoint],
                                 errors="coerce").to_numpy(float)
        probes = by_symbol.get(enzyme, [])
        if not probes:
            rows[endpoint] = {
                "enzyme": enzyme, "status": "NOT_EVALUABLE",
                "reason": ("no admissible probe for this enzyme; every probe is "
                           "graded conflicting by the frozen identity rule"),
                "published_protein": GRANGEON_PROTEIN_ACTIVITY[endpoint]}
            continue
        # One statistic, no selection: the probe with the highest mean, the same
        # rule every other runner in this repository uses.
        means = [float(np.nanmean(data["values"][p])) for p in probes]
        chosen = probes[int(np.argmax(means))]
        transcript = data["values"][chosen]
        keep = np.isfinite(transcript) & np.isfinite(activity)
        x, y = transcript[keep], activity[keep]
        rho = spearman(x, y)
        low, high = bootstrap_interval(x, y)
        published = GRANGEON_PROTEIN_ACTIVITY[endpoint]
        rows[endpoint] = {
            "enzyme": enzyme,
            "status": "MEASURED",
            "hlc_transcript_activity_spearman": rho,
            "hlc_ci_low": low,
            "hlc_ci_high": high,
            "hlc_n_donors": int(keep.sum()),
            "hlc_probes_admissible": len(probes),
            "hlc_probe_used": str(data["feature_ids"][chosen]),
            "published_protein": published,
            "gap_protein_minus_transcript_hlm": published["hlm_r"] - rho,
            "published_hlm_r_above_hlc_interval": bool(
                published["hlm_r"] > high),
        }

    measured = [r for r in rows.values() if r.get("status") == "MEASURED"]
    transcript_rhos = [r["hlc_transcript_activity_spearman"] for r in measured]
    protein_rhos = [r["published_protein"]["hlm_r"] for r in measured]
    above = sum(1 for r in measured if r["published_hlm_r_above_hlc_interval"])

    report = {
        "report": ("mRNA versus published protein as an explanation of CYP "
                   "activity"),
        "label": "DESCRIPTIVE_CROSS_COHORT_COMPARISON",
        "built_on": date.today().isoformat(),
        "what_this_is": (
            "Two whole-cohort Spearman correlations placed side by side: "
            "transcript-to-activity measured here in the Human Liver Cohort, and "
            "protein-to-activity as published by Grangeon 2025 in genotyped human "
            "liver microsomes."),
        "what_this_is_not": (
            "NOT a within-subject comparison. The donors are different people in "
            "different cohorts, and no individual from one is joined to the "
            "other. It licenses no causal claim, no per-patient claim and no "
            "pharmacokinetic claim. The endpoint class on both sides is "
            "functional_in_vitro."),
        "published_source": GRANGEON_CITATION,
        "published_source_caveats": [
            "n = 16 (microsomes) and n = 20 (hepatocytes), so each published r "
            "carries a wide interval that the article does not report.",
            "Activity values were provided by the commercial supplier, not "
            "measured by the authors.",
            "No per-donor numeric value for protein or activity is published, so "
            "the published correlations cannot be recomputed or re-analysed.",
        ],
        "hlc_method": {
            "statistic": "Spearman rho, average ranks for ties",
            "probe_rule": ("one probe per enzyme, the admissible probe with the "
                           "highest mean; pool A_qualified_human"),
            "interval": f"{N_BOOTSTRAP}-draw donor bootstrap, seed {SEED}",
            "no_model": ("no regression, no selection, no tuning, no "
                         "cross-validation: a descriptive whole-cohort "
                         "correlation, which is what the published number is"),
        },
        "summary": {
            "endpoints_measured": len(measured),
            "endpoints_not_evaluable": len(rows) - len(measured),
            "hlc_transcript_rho_median": float(np.median(transcript_rhos)),
            "hlc_transcript_rho_min": float(np.min(transcript_rhos)),
            "hlc_transcript_rho_max": float(np.max(transcript_rhos)),
            "published_protein_rho_median": float(np.median(protein_rhos)),
            "published_protein_rho_min": float(np.min(protein_rhos)),
            "published_protein_rho_max": float(np.max(protein_rhos)),
            "endpoints_where_published_protein_r_exceeds_the_hlc_interval": above,
            "how_to_read": (
                "The last count is the strict comparison: on how many enzymes "
                "the published protein correlation lies entirely above the "
                "bootstrap interval of this cohort's transcript correlation. It "
                "is still a cross-cohort statement."),
        },
        "per_endpoint": rows,
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print(f"written {out}\n")
    print(f"{'endpoint':20s} {'n':>4s} {'mRNA rho':>9s} {'95% CI':>18s} "
          f"{'protein r':>9s} {'gap':>7s} above?")
    for endpoint, row in rows.items():
        if row.get("status") != "MEASURED":
            print(f"{endpoint:20s}  NOT_EVALUABLE  "
                  f"(published protein r = {row['published_protein']['hlm_r']})")
            continue
        print(f"{endpoint:20s} {row['hlc_n_donors']:4d} "
              f"{row['hlc_transcript_activity_spearman']:9.3f} "
              f"[{row['hlc_ci_low']:6.3f},{row['hlc_ci_high']:6.3f}] "
              f"{row['published_protein']['hlm_r']:9.4f} "
              f"{row['gap_protein_minus_transcript_hlm']:7.3f} "
              f"{row['published_hlm_r_above_hlc_interval']}")
    s = report["summary"]
    print(f"\nmedian mRNA rho {s['hlc_transcript_rho_median']:.3f} "
          f"(range {s['hlc_transcript_rho_min']:.3f} to "
          f"{s['hlc_transcript_rho_max']:.3f})")
    print(f"median published protein r {s['published_protein_rho_median']:.3f} "
          f"(range {s['published_protein_rho_min']:.3f} to "
          f"{s['published_protein_rho_max']:.3f})")
    print(f"published protein r above the HLC interval on "
          f"{s['endpoints_where_published_protein_r_exceeds_the_hlc_interval']} "
          f"of {s['endpoints_measured']} enzymes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
