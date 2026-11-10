#!/usr/bin/env python
"""Size the three open questions, each at its own level of evidence.

Protocol: ``docs/OMICS2PK_COHORT_SIZING_PRESPECIFICATION_20261110.md``.

Question A is sized from **transferable published effects** plus this
project's own two cohorts. Questions B and C are sized from **declared
scenarios**, because no public source gives this project an out-of-sample
incremental effect for either, and no in-vitro effect is carried into them.

Usage:
    python -I scripts/run_cohort_size_simulation.py \
        --per-donor <outside-repo>/lee_jackson_cyp3a_panel_per_donor.csv \
        --out reports/data/cohort_size_simulation.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.stdout.reconfigure(encoding="utf-8")

from omics2pk import cohort_size_simulation as CSS  # noqa: E402
from omics2pk import plasma_biomarker_validation as PBV  # noqa: E402

PRESPEC = ROOT / "docs" / "OMICS2PK_COHORT_SIZING_PRESPECIFICATION_20261110.md"

# --------------------------------------------------------------------------
# Level 1 inputs: published rank correlations between an enzyme's abundance
# and its own probe reaction, in individual human preparations. These enter
# as SCALAR parameters. No individual from any of these sources is joined to
# any individual from any other.
# --------------------------------------------------------------------------
PUBLISHED_ABUNDANCE_ACTIVITY = {
    "source": ("Pharmaceuticals 2025;18:1789, PMCID PMC12735975, CC BY, "
               "Table 3: Spearman correlation between CYP450 activity and "
               "CYP450 content in individual preparations."),
    "hepatocytes_n20": {
        "CYP1A2_phenacetin": 0.788, "CYP2A6_coumarin": 0.7985,
        "CYP2B6_bupropion": 0.7424, "CYP2C8_amodiaquine": 0.7574,
        "CYP2C9_diclofenac": 0.3805, "CYP2C19_mephenytoin": 0.7003,
        "CYP2D6_dextromethorphan": 0.8015, "CYP2E1_chlorzoxazone": 0.2915,
        "CYP3A4_midazolam": 0.5353, "CYP3A4_testosterone": 0.7684,
    },
    "genotyped_microsomes_n16": {
        "CYP1A2_phenacetin": 0.6814, "CYP2A6_coumarin": 0.9471,
        "CYP2B6_bupropion": 0.9794, "CYP2C8_amodiaquine": 0.9176,
        "CYP2C9_diclofenac": 0.65, "CYP2C19_mephenytoin": 0.95,
        "CYP2D6_dextromethorphan": 0.8291, "CYP2E1_chlorzoxazone": 0.9441,
        "CYP3A4_midazolam": 0.8168, "CYP3A4_testosterone": 0.9559,
        "CYP4A11_lauric_acid": 0.7176,
    },
    "intestinal_mucosa_n15": {
        "CYP2C9_diclofenac": -0.3071, "CYP2C19_mephenytoin": 0.7036,
        "CYP2D6_dextromethorphan": 0.5771, "CYP3A4_midazolam": 0.1321,
        "CYP3A4_testosterone": 0.2321, "CYP2J2_astemizole": 0.1679,
        "CES2_irinotecan": -0.0643,
    },
    "why_these_are_transferable": (
        "Same kind of quantity, same kind of endpoint and the same kind of "
        "preparation as question A: an absolutely quantified enzyme abundance "
        "against its own probe reaction in individual human material. They "
        "are IN-SAMPLE rank correlations, so they are used only to set the "
        "generating correlation of a simulation whose own out-of-sample "
        "performance is then measured -- never read as performance."),
    "why_they_are_not_a_result": (
        "they are published aggregates with no per-donor row, so nothing can "
        "be recomputed from them and they cannot validate anything"),
}

#: The rank-correlation grid question A is sized over, chosen to span the
#: published range rather than to flatter it.
RHO_GRID = (0.20, 0.40, 0.53, 0.70, 0.80, 0.90)

#: Question A reference scenario. The analytical dispersion is on the log10
#: scale and 0.15 is about a 40 % coefficient of variation, which is the order
#: of magnitude targeted proteomics reports for absolute CYP quantification.
A_REFERENCE = {
    "rho_pearson_log": CSS.spearman_to_pearson(0.53),
    "analytical_cv_log": 0.15,
    "missing_rate": 0.0,
}
A_VARIATIONS = {
    "rho_pearson_log": [CSS.spearman_to_pearson(r) for r in (0.20, 0.70, 0.90)],
    "analytical_cv_log": [0.05, 0.30, 0.50],
    "missing_rate": [0.05, 0.20],
}

#: Questions B and C reference scenario. Every number here is a DECLARED
#: ASSUMPTION; none is an estimate from data.
B_REFERENCE = {
    "baseline_r2": 0.30,
    "target_mae_reduction": 0.10,
    "n_features": 200,
    "n_informative": 10,
    "redundancy": 0.30,
    "analytical_cv": 0.30,
    "missing_rate": 0.05,
    "selected": 10,
}
B_VARIATIONS = {
    "target_mae_reduction": [0.05, 0.20],
    "redundancy": [0.0, 0.60, 0.85],
    "n_features": [20, 2000],
    "analytical_cv": [0.10, 0.60],
    "baseline_r2": [0.10, 0.50],
    "missing_rate": [0.20],
    "selected": [2, 30],
}

#: The dimensionality probe runs at a reduced replicate count, declared here
#: rather than chosen quietly: a 20,000-feature block is roughly fifty times
#: the cost per replicate of the 200-feature reference.
DIMENSIONALITY_REPLICATES = 100

#: Question C: the three endpoints are not interchangeable, and the difference
#: between them is encoded as how much of the outcome a clinical and
#: pharmacogenetic baseline already explains, and how noisy the measurement is.
C_ENDPOINTS = {
    "single_trough_concentration": {
        "baseline_r2": 0.20, "analytical_cv": 0.40,
        "why": ("one timed sample carries assay error, sampling-time error and "
                "adherence error on top of between-subject variability, so a "
                "clinical plus pharmacogenetic baseline explains least here. "
                "ARIANA's endpoint is of this kind."),
    },
    "auc_from_dense_sampling": {
        "baseline_r2": 0.40, "analytical_cv": 0.20,
        "why": ("averaging over a profile removes much of the sampling-time "
                "error, so the same covariates explain more and the "
                "measurement is less noisy -- at the cost of a far more "
                "expensive protocol."),
    },
    "clearance_from_a_population_model": {
        "baseline_r2": 0.50, "analytical_cv": 0.15,
        "why": ("a model-derived clearance is the least noisy and the most "
                "explained by covariates, because the model has already "
                "absorbed covariate effects -- which means an omic layer must "
                "beat a baseline that is partly the same model, and the "
                "apparent gain depends on the structural model chosen."),
    },
}


class RunError(RuntimeError):
    pass


def read_per_donor(path):
    rows, order = {}, []
    with open(path, newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle):
            donor = record.pop("donor")
            parsed = {}
            for key, value in record.items():
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    parsed[key] = value
            parsed["cyp3a_total_abundance"] = (parsed["cyp3a4_abundance"]
                                               + parsed["cyp3a5_abundance"])
            rows[donor] = parsed
            order.append(donor)
    return rows, order


def measured_effects(rows, donors):
    """The paired effect and dispersion this project has actually measured.

    Used for the CONFIRMATORY sizing, which is a different question from
    discovery sizing and is reported separately.
    """
    out = {}
    for label, endpoint in (("midazolam_1_hydroxylation",
                             PBV.FROZEN["primary_endpoint"]),
                            ("ibrutinib_m37",
                             PBV.FROZEN["secondary_endpoint"])):
        _, arm_errors = PBV.score_arm(rows, donors, endpoint,
                                      PBV.ARMS["V1_hepatic_protein"])
        _, base_errors = PBV.score_arm(rows, donors, endpoint,
                                       PBV.ARMS["V0_mean"])
        difference = np.asarray(arm_errors) - np.asarray(base_errors)
        effect = float(difference.mean())
        dispersion = float(difference.std(ddof=1))
        sized = CSS.validation_cohort_size(effect, dispersion)
        out[label] = {
            "contrast": "V1_hepatic_protein vs V0_mean, leave-one-donor-out",
            "mean_paired_difference_log10": effect,
            "dispersion_of_paired_difference_log10": dispersion,
            "n_donors_observed": int(difference.size),
            "confirmatory_sizing": sized,
        }
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-donor", required=True)
    parser.add_argument("--replicates", type=int,
                        default=CSS.FROZEN["replicates"])
    parser.add_argument("--dim-replicates", type=int,
                        default=DIMENSIONALITY_REPLICATES,
                        help="replicates for the layer-width power curve only")
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "cohort_size_simulation.json"))
    args = parser.parse_args(argv)

    if not PRESPEC.is_file():
        raise RunError(f"the frozen protocol is missing: {PRESPEC}")

    per_donor = pathlib.Path(args.per_donor).resolve()
    if ROOT == per_donor.parent or ROOT in per_donor.parents:
        raise RunError("the per-donor table must stay outside this repository")
    rows, donors = read_per_donor(per_donor)

    kwargs = {"replicates": args.replicates}

    # ---- question A: level 1, transferable -------------------------------
    rho_curve = {}
    for rho in RHO_GRID:
        scenario = dict(A_REFERENCE)
        scenario["rho_pearson_log"] = CSS.spearman_to_pearson(rho)
        result = CSS.required_sample_size(
            CSS.simulate_enzyme_activity, scenario, **kwargs)
        rho_curve[f"spearman_{rho:.2f}"] = {
            "pearson_on_log_assumed": round(
                CSS.spearman_to_pearson(rho), 4),
            "required_n": result["required_n"],
            "verdict": result["verdict"],
        }
    a_sensitivity = CSS.one_at_a_time_sensitivity(
        CSS.simulate_enzyme_activity, A_REFERENCE, A_VARIATIONS, **kwargs)

    # ---- question B: level 2, declared scenarios ------------------------
    b_curve = {}
    for reduction in (0.05, 0.10, 0.20, 0.30):
        scenario = dict(B_REFERENCE)
        scenario["target_mae_reduction"] = reduction
        result = CSS.required_sample_size(
            CSS.simulate_incremental_layer, scenario, **kwargs)
        b_curve[f"mae_reduction_{int(reduction * 100)}pc"] = {
            "required_n": result["required_n"], "verdict": result["verdict"]}
    b_sensitivity = CSS.one_at_a_time_sensitivity(
        CSS.simulate_incremental_layer, B_REFERENCE, B_VARIATIONS, **kwargs)

    # Dimensionality is reported as a power curve at fixed n rather than as a
    # required n, for one reason that is declared and not aesthetic: a
    # 20,000-feature block costs about 1.6 s per replicate at n = 800, so a
    # bisection over n would cost more than the rest of this run together. A
    # curve at four fixed sizes shows the same effect and says what it cost.
    dimensionality = {}
    for features in (200, 2000, 20000):
        row = {}
        for n in (100, 200, 400, 800):
            scenario = dict(B_REFERENCE)
            scenario["n_features"] = features
            row[f"n_{n}"] = round(CSS.power_at(
                n, CSS.simulate_incremental_layer, scenario,
                replicates=args.dim_replicates), 3)
        dimensionality[f"features_{features}"] = row

    # ---- question C: level 2, three endpoints that are not the same -----
    c_results = {}
    for endpoint, settings in C_ENDPOINTS.items():
        per_endpoint = {"why_this_endpoint_differs": settings["why"]}
        for reduction in (0.10, 0.20):
            scenario = dict(B_REFERENCE)
            scenario["baseline_r2"] = settings["baseline_r2"]
            scenario["analytical_cv"] = settings["analytical_cv"]
            scenario["target_mae_reduction"] = reduction
            result = CSS.required_sample_size(
                CSS.simulate_incremental_layer, scenario, **kwargs)
            per_endpoint[f"mae_reduction_{int(reduction * 100)}pc"] = {
                "required_n": result["required_n"],
                "verdict": result["verdict"]}
        c_results[endpoint] = per_endpoint

    payload = {
        "run": "cohort_size_simulation",
        "generated": "2026-11-10",
        "protocol": PRESPEC.name,
        "protocol_sha256": hashlib.sha256(PRESPEC.read_bytes()).hexdigest(),
        "frozen": {k: (list(v) if isinstance(v, tuple) else v)
                   for k, v in CSS.FROZEN.items()},
        "replicates_used": args.replicates,
        "validation_scheme": CSS.VALIDATION_NOTE,
        "question_a_enzyme_activity": {
            "evidence_level": "transferable_literature",
            "published_inputs": PUBLISHED_ABUNDANCE_ACTIVITY,
            "required_n_by_rank_correlation": rho_curve,
            "sensitivity": a_sensitivity,
            "reference_scenario": A_REFERENCE,
            "measured_confirmatory_sizing": measured_effects(rows, donors),
        },
        "question_b_multiomic_complementarity": {
            "evidence_level": "declared_scenario",
            "why_not_transferable": (
                "no public source gives this project an OUT-OF-SAMPLE "
                "incremental effect for an omic layer beyond clinical "
                "covariates and pharmacogenetics on a pharmacokinetic "
                "endpoint. ARIANA is the only internal estimate and it is one "
                "drug, 47 patients, one trough concentration and two "
                "metabolite ratios -- it cannot stand for 'an omic layer'. No "
                "in-vitro effect is carried in."),
            "reference_scenario": B_REFERENCE,
            "required_n_by_assumed_reduction": b_curve,
            "sensitivity": b_sensitivity,
            "power_by_layer_width": {
                "replicates": args.dim_replicates,
                "why_a_power_curve_and_not_a_required_n": (
                    "a 20,000-feature block costs about fifty times the "
                    "reference per replicate, so a bisection over n would "
                    "cost more than the rest of this run together. The curve "
                    "at four fixed sizes shows the same effect."),
                "curve": dimensionality,
            },
        },
        "question_c_clinical_pk": {
            "evidence_level": "declared_scenario",
            "endpoints_are_not_interchangeable": (
                "a single timed concentration, an AUC and a model-derived "
                "clearance differ in how much a clinical and pharmacogenetic "
                "baseline already explains and in how noisy the measurement "
                "is. Sizing one does not size another."),
            "refusal": (
                "no effect measured on in-vitro intrinsic clearance or on an "
                "in-vitro activity rate is transferred to a clinical "
                "pharmacokinetic endpoint here. The in-vitro work sizes the "
                "in-vitro question only."),
            "by_endpoint": c_results,
        },
        "how_to_read_any_required_n": CSS.how_to_read_any_required_n(),
        "no_individual_level_data_in_this_artefact": True,
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    print("QUESTION A -- required n by published rank correlation")
    for key, row in rho_curve.items():
        print(f"  {key:<18} -> n = {row['required_n']}")
    print("  dominant assumptions, widest required-n first:")
    for key in a_sensitivity["ordered_by_required_n_descending"][:5]:
        print(f"    {key:<34} n = "
              f"{a_sensitivity['variations'][key]['required_n']}")
    print("\nQUESTION B -- required n by ASSUMED MAE reduction (scenarios)")
    for key, row in b_curve.items():
        print(f"  {key:<24} -> n = {row['required_n']}")
    print("  dominant assumptions, widest required-n first:")
    for key in b_sensitivity["ordered_by_required_n_descending"][:6]:
        print(f"    {key:<34} n = "
              f"{b_sensitivity['variations'][key]['required_n']}")
    print("\nQUESTION C -- required n by endpoint (scenarios)")
    for endpoint, block in c_results.items():
        sizes = {k: v["required_n"] for k, v in block.items()
                 if k.startswith("mae_reduction")}
        print(f"  {endpoint:<36} {sizes}")
    print("\nCONFIRMATORY sizing from measured effects")
    for label, block in payload["question_a_enzyme_activity"][
            "measured_confirmatory_sizing"].items():
        print(f"  {label:<28} effect "
              f"{block['mean_paired_difference_log10']:+.4f} "
              f"sd {block['dispersion_of_paired_difference_log10']:.4f} "
              f"-> n = {block['confirmatory_sizing']['required_n']}")
    print(f"\nartefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
