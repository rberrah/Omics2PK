"""Research-only original mavoglurant adaptive GH E3 NLME benchmark.

Fits are executed on the training subjects only. For held-out subjects,
the log prediction is the unconditional population prior expectation
(no EBEs), while the prospective joint marginal predictive density
integrates over unknown random effects of each test subject.
"""
from __future__ import annotations

from hashlib import sha256
from json import dumps
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen

from omics2pk.adaptive_iv_nlme import (
    MarginalLikelihoodError, adaptive_subject_loglik, adaptive_total_loglik,
)
from omics2pk.iv_pk_structural import subject_disjoint_folds
from scripts.benchmark_mavoglurant_iv_nlme import _data, score_population_heldout
from scripts.benchmark_mavoglurant_iv_structural import (
    _fit, _predict, _validate, qualified_original_profiles,
)
from scripts.validate_mavoglurant_120_real_occasion_nca import SOURCE, EXPECTED_SHA


def fit_adaptive(train, pooled=None, nodes=11, maxiter=180):
    import numpy as np
    from scipy.optimize import minimize
    data = _data(train)
    pooled = np.asarray(pooled if pooled is not None else _fit("2C", train), dtype=float)
    if pooled.shape != (4,) or not np.all(np.isfinite(pooled)) or np.any(pooled <= 0):
        raise ValueError("Invalid training-only 2C initializer")
    # Only training outcomes enter this residual initializer.
    residuals = np.log(_predict("2C", pooled, data[0])) - data[4]
    sigma0 = float(np.clip(np.sqrt(np.mean(residuals**2)), 0.12, 1.5))
    lower = np.array([0.1, 1.0, 0.05, 1.0, 0.01, 0.03])
    upper = np.array([300, 5000, 500, 5000, 1.5, 2.5])
    bounds = list(zip(np.log(lower), np.log(upper)))

    def loss(logparams):
        try:
            nll = -adaptive_total_loglik(np.exp(logparams), data, nodes)
            return nll if isfinite(nll) else 1e30
        except (MarginalLikelihoodError, ValueError, OverflowError, FloatingPointError):
            return 1e30

    candidates = []
    for omega in (0.18, 0.50):
        start = np.clip(np.r_[pooled, omega, sigma0], lower*1.001, upper/1.001)
        result = minimize(loss, np.log(start), bounds=bounds, method="L-BFGS-B",
                          options={"maxiter": maxiter, "ftol": 1e-7})
        if np.isfinite(result.fun):
            candidates.append((float(result.fun), np.exp(result.x), bool(result.success)))
    if not candidates:
        raise ValueError("No finite adaptive marginal likelihood solution")
    candidates.sort(key=lambda item: item[0])
    _, theta, converged = candidates[0]
    if not converged:
        raise ValueError("Adaptive ML optimizer not converged; abstain")
    return theta


def evaluate(profiles, folds=5, nodes=11, check_nodes=(7, 15, 21)):
    import numpy as np
    people = {p["subject"] for p in profiles}
    splits = subject_disjoint_folds(people, folds)
    result_rows = []
    for train_ids, test_ids in splits:
        train = [p for p in profiles if p["subject"] in train_ids]
        heldout = [p for p in profiles if p["subject"] in test_ids]
        if {p["subject"] for p in train} != train_ids or {p["subject"] for p in heldout} != test_ids:
            raise ValueError("Missing distinct biological person in partition")
        fixed = _fit("2C", train)
        theta = fit_adaptive(train, fixed, nodes)
        fixed_score = _validate("2C", fixed, heldout)
        mixed_score = score_population_heldout(theta, heldout, nodes)
        # Re-score the SAME fitted parameters at multiple quadrature orders:
        # this checks integration convergence rather than conflating changes
        # in optimizer solutions with changes in approximation order.
        test_data = _data(heldout)
        n_humans = len(test_data[1])
        density = {
            n: float(-adaptive_total_loglik(theta, test_data, nodes=n)/n_humans)
            for n in sorted(set((*check_nodes, nodes)))
        }
        result_rows.append({
            "n_people": n_humans, "n_occasions": len(heldout),
            "n_concentrations": fixed_score["n_measured_heldout_concentrations"],
            "fixed_mse": fixed_score["person_mse_log"],
            "mixed_mse": mixed_score["person_mse_log"],
            "joint_negative_log_density_per_person": density,
        })
    total = sum(row["n_people"] for row in result_rows)
    total_obs = sum(row["n_concentrations"] for row in result_rows)
    if total != len(people) or total_obs != sum(len(p["points"]) for p in profiles):
        raise ValueError("Subject-disjoint original source coverage differs")
    def mse(name):
        return sum(row[name] * row["n_people"] for row in result_rows)/total
    results = {
        "status": "RESEARCH_2C_ADAPTIVE_NLME_ORIGINAL_MEASURED_PK_NO_OMICS",
        "n_humans": total, "n_occasions": len(profiles),
        "n_concentrations": total_obs, "n_folds": folds,
        "models": {
            "2C_pooled": {"RMSLE": round(float(np.sqrt(mse("fixed_mse"))), 6)},
            "2C_adaptive_CL_IIV": {"RMSLE": round(float(np.sqrt(mse("mixed_mse"))), 6)},
        },
        "heldout_joint_NLL_per_human_by_nodes": {
            str(n): round(sum(row["joint_negative_log_density_per_person"][n]*row["n_people"]
                              for row in result_rows)/total, 6)
            for n in sorted(set((*check_nodes, nodes)))
        },
        "population_point_prediction": "prior-only/no test-subject EBE or conditioning",
        "marginal_density": "per-person adaptive posterior-proposal GH over shared CL ETA",
        "not_SAEM_or_FOCEi": True, "no_IOV_BLQ_or_omics": True,
    }
    vals = list(results["heldout_joint_NLL_per_human_by_nodes"].values())
    results["quadrature_range_NLL_per_human"] = round(max(vals)-min(vals), 6)
    return results


def run_original():
    import pyreadr
    with urlopen(Request(SOURCE, headers={"User-Agent": "Omics2PK-adaptive-nlme-research/1.0"}), timeout=45) as resp:
        payload = resp.read(3_000_001)
    if sha256(payload).hexdigest() != EXPECTED_SHA:
        raise ValueError("Unrecognized original PK source, ABSTAIN")
    with TemporaryDirectory() as folder:
        original = Path(folder) / "mavoglurant.rda"
        original.write_bytes(payload)
        rows = pyreadr.read_r(str(original))["mavoglurant"].to_dict("records")
    result = evaluate(qualified_original_profiles(rows))
    if (result["n_humans"], result["n_occasions"], result["n_concentrations"]) != (120,196,2403):
        raise ValueError("Changed original source cardinalities")
    result["sha256_verified"] = EXPECTED_SHA
    return result


if __name__ == "__main__":
    print(dumps(run_original(), indent=2, sort_keys=True))
