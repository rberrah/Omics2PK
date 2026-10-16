"""Original-source gated marginal-likelihood 2C IV NLME proof of concept.

One subject-level lognormal inter-individual random effect on clearance.
Both observed occasions share ETA. Predictions for unseen people integrate
the population prior without EBEs or held-out concentration conditioning.
Research-only: no IOV, BLQ likelihood, SAEM, clinical or omics covariates.
"""
from __future__ import annotations

from hashlib import sha256
import json
from math import log, pi, sqrt
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen

from omics2pk.iv_pk_structural import subject_disjoint_folds
from scripts.benchmark_mavoglurant_iv_structural import (
    _arrays, _fit, _predict, _validate, qualified_original_profiles,
)
from scripts.validate_mavoglurant_120_real_occasion_nca import SOURCE, EXPECTED_SHA


def _data(profiles):
    """Group measurements by biological subject across occasions."""
    import numpy as np
    if not profiles:
        raise ValueError("Require measured PK profiles")
    if any(not str(p["subject"]).strip() for p in profiles):
        raise ValueError("Unidentified patient is not eligible for NLME")
    pack = _arrays(profiles)
    if len(pack[1]) == 0 or np.any(~np.isfinite(pack[1])) or np.any(pack[1] <= 0):
        raise ValueError("Only actual positive measured concentrations; BLQ needs a separate likelihood")
    if np.any(~np.isfinite(pack[0])) or np.any(pack[0] <= 0):
        raise ValueError("Invalid original time origin")
    people, indices = np.unique(np.asarray(pack[5], dtype=str), return_inverse=True)
    counts = np.bincount(indices, minlength=len(people))
    if np.any(counts < 3):
        raise ValueError("Too few measured concentrations for a subject")
    return pack, people, indices, counts, np.log(pack[1])


def _grid_prediction(parameters, pack, gh_nodes=9):
    """Log C(t; CL*exp(ETA)) at Gaussian prior quadrature nodes, not EBEs."""
    import numpy as np
    if gh_nodes < 3 or gh_nodes % 2 == 0:
        raise ValueError("Quadrature order must be odd and >=3")
    theta = np.asarray(parameters, dtype=float)
    if theta.shape != (6,) or np.any(~np.isfinite(theta)) or np.any(theta <= 0):
        raise ValueError("Require positive CL,Vc,Q,Vp,omega_CL,sigma_log")
    cl, vc, q, vp, omega, _ = theta
    x, w = np.polynomial.hermite.hermgauss(gh_nodes)
    weights = w / sqrt(pi)
    log_pred = []
    for eta in sqrt(2) * omega * x:
        pred = _predict("2C", [cl * np.exp(eta), vc, q, vp], pack)
        if np.any(~np.isfinite(pred)) or np.any(pred <= 0):
            raise ValueError("Nonpositive or nonfinite conditional concentration")
        log_pred.append(np.log(pred))
    return np.asarray(log_pred), np.log(weights)


def _marginal_loglike(parameters, data, gh_nodes=9):
    """Integrate product of all conditional residual densities per subject."""
    import numpy as np
    from scipy.special import logsumexp
    pack, people, indices, counts, log_dv = data
    log_predictions, log_weights = _grid_prediction(parameters, pack, gh_nodes)
    sigma = float(parameters[5])
    per_node = []
    constant = -counts * (log(sigma) + 0.5 * log(2 * pi))
    for node in range(len(log_weights)):
        squared = (log_dv - log_predictions[node]) ** 2
        rss = np.bincount(indices, weights=squared, minlength=len(people))
        per_node.append(constant - 0.5 * rss / (sigma ** 2) + log_weights[node])
    return float(np.sum(logsumexp(np.asarray(per_node), axis=0)))


def fit_iiv_cl(train_profiles, pooled_parameters=None, gh_nodes=9, maxiter=160):
    """Marginal ML for typical 2C PK, omega_CL and log residual sigma."""
    import numpy as np
    from scipy.optimize import minimize
    data = _data(train_profiles)
    pooled = np.asarray(pooled_parameters if pooled_parameters is not None
                        else _fit("2C", train_profiles), dtype=float)
    if pooled.shape != (4,) or np.any(~np.isfinite(pooled)) or np.any(pooled <= 0):
        raise ValueError("Invalid training-only pooled 2C initializer")
    pred = _predict("2C", pooled, data[0])
    sigma0 = float(np.clip(np.sqrt(np.mean((np.log(pred) - data[4]) ** 2)), 0.12, 1.5))
    lower = np.array([0.1, 1.0, 0.05, 1.0, 0.01, 0.03])
    upper = np.array([300, 5000, 500, 5000, 1.5, 2.5])
    bounds = list(zip(np.log(lower), np.log(upper)))
    def objective(log_theta):
        try:
            val = -_marginal_loglike(np.exp(log_theta), data, gh_nodes)
            return float(val) if np.isfinite(val) else 1e30
        except (FloatingPointError, ValueError, OverflowError):
            return 1e30
    candidates = []
    for omega0 in (0.15, 0.5):
        start = np.clip(np.r_[pooled, omega0, sigma0], lower * 1.001, upper / 1.001)
        fit = minimize(objective, np.log(start), method="L-BFGS-B",
                       bounds=bounds, options={"maxiter": maxiter, "ftol": 1e-7})
        if np.isfinite(fit.fun):
            candidates.append((float(fit.fun), np.exp(fit.x), bool(fit.success)))
    if not candidates:
        raise ValueError("Marginal NLME optimization did not yield a finite result")
    candidates.sort(key=lambda row: row[0])
    nll, parameters, converged = candidates[0]
    if not converged:
        raise ValueError("NLME optimizer did not converge; abstain")
    return parameters


def score_population_heldout(parameters, heldout_profiles, gh_nodes=9):
    """Population-prior geometric prediction only, without individual EBEs."""
    import numpy as np
    pack, people, indices, counts, log_dv = _data(heldout_profiles)
    grid, log_weights = _grid_prediction(parameters, pack, gh_nodes)
    prediction = np.sum(grid * np.exp(log_weights)[:, None], axis=0)
    residuals = prediction - log_dv
    person_mse = np.bincount(indices, weights=residuals ** 2) / counts
    person_bias = np.bincount(indices, weights=residuals) / counts
    return {
        "n_independent_heldout_humans": len(people),
        "n_heldout_occasions": len(heldout_profiles),
        "n_measured_heldout_concentrations": len(log_dv),
        "person_mse_log": float(np.mean(person_mse)),
        "person_mean_bias_log": float(np.mean(person_bias)),
        "heldout_joint_negative_log_likelihood_per_person":
            -_marginal_loglike(parameters, (pack, people, indices, counts, log_dv), gh_nodes) / len(people),
    }


def heldout_iiv_cl_comparison(profiles, n_folds=5, gh_nodes=9, maxiter=160):
    """Prespecified fixed-2C versus 2C IIV-CL, whole-patient heldouts."""
    import numpy as np
    people = {p["subject"] for p in profiles}
    partitions = subject_disjoint_folds(people, n_folds=n_folds)
    scores = {"2C_pooled": [], "2C_NLME_IIV_CL": []}
    for train_people, test_people in partitions:
        train = [p for p in profiles if p["subject"] in train_people]
        test = [p for p in profiles if p["subject"] in test_people]
        if {p["subject"] for p in train} != train_people or {p["subject"] for p in test} != test_people:
            raise ValueError("Missing biological subject in fold")
        pooled = _fit("2C", train)
        parameters = fit_iiv_cl(train, pooled, gh_nodes, maxiter)
        scores["2C_pooled"].append(_validate("2C", pooled, test))
        scores["2C_NLME_IIV_CL"].append(score_population_heldout(parameters, test, gh_nodes))
    n_people = len(people)
    n_obs = sum(len(p["points"]) for p in profiles)
    if sum(s["n_independent_heldout_humans"] for s in scores["2C_pooled"]) != n_people:
        raise ValueError("Whole-person validation coverage failed")
    if sum(s["n_measured_heldout_concentrations"] for s in scores["2C_NLME_IIV_CL"]) != n_obs:
        raise ValueError("Original observation coverage failed")
    report = {
        "status": "MEASURED_PK_ONLY_2C_NLME_IIV_CL_INTERNAL_PATIENT_HELDOUT_NO_OMICS",
        "n_biological_humans": n_people, "n_qualified_occasions": len(profiles),
        "n_measured_original_concentrations": n_obs,
        "folds": n_folds, "gh_nodes": gh_nodes,
        "random_effects": "shared subject ETA on log(CL), including both occasions",
        "residual_model": "Gaussian residual on log DV, actual positive MDV0 measurements only",
        "likelihood": "per-subject Gauss-Hermite-integrated marginal maximum likelihood",
        "heldout": "unconditional population-prior geometric forecast; no EBEs or test-subject concentration fitting",
        "limitations": "no IOV, BLQ likelihood, SAEM, uncertainty, covariates, omics or independent external validation",
        "no_omics_in_source": True, "not_clinically_validated": True,
    }
    for model, pieces in scores.items():
        count = sum(x["n_independent_heldout_humans"] for x in pieces)
        mse = sum(x["person_mse_log"] * x["n_independent_heldout_humans"] for x in pieces) / count
        bias = sum(x["person_mean_bias_log"] * x["n_independent_heldout_humans"] for x in pieces) / count
        report[model] = {
            "heldout_root_mean_person_mse_log": round(float(np.sqrt(mse)), 6),
            "heldout_mean_person_log_bias": round(float(bias), 6),
            "heldout_people": count,
        }
        if model == "2C_NLME_IIV_CL":
            report[model]["joint_negative_log_likelihood_per_person"] = round(
                sum(x["heldout_joint_negative_log_likelihood_per_person"] * x["n_independent_heldout_humans"]
                    for x in pieces) / count, 6)
    report["NLME_minus_pooled_RMSLE"] = round(
        report["2C_NLME_IIV_CL"]["heldout_root_mean_person_mse_log"]
        - report["2C_pooled"]["heldout_root_mean_person_mse_log"], 6)
    return report


def run_original():
    import pyreadr
    with urlopen(Request(SOURCE, headers={"User-Agent": "Omics2PK-NLME-research/1.0"}), timeout=45) as resp:
        payload = resp.read(3_000_001)
    if sha256(payload).hexdigest() != EXPECTED_SHA:
        raise ValueError("Original RDA source SHA mismatch; abstain")
    with TemporaryDirectory() as directory:
        path = Path(directory) / "mavoglurant.rda"
        path.write_bytes(payload)
        rows = pyreadr.read_r(str(path))["mavoglurant"].to_dict("records")
    report = heldout_iiv_cl_comparison(qualified_original_profiles(rows))
    if (report["n_biological_humans"], report["n_qualified_occasions"],
        report["n_measured_original_concentrations"]) != (120, 196, 2403):
        raise ValueError("Original subject/occasion/DV count mismatch; abstain")
    report["source_sha256_verified"] = EXPECTED_SHA
    return report


if __name__ == "__main__":
    print(json.dumps(run_original(), sort_keys=True, indent=2))
