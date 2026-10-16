"""Vectorized adaptive Gauss-Hermite integration for a narrow 2C IV IIV-CL model.

Only a one-dimensional *subject-shared* CL effect is supported. Per-subject
posterior proposals are centered using the complete TRAIN subject likelihood.
Do not apply subject-specific centering to held-out point prediction: this
module is solely for marginal training and prospective joint score density.
"""
from __future__ import annotations

from math import pi, sqrt, log

import numpy as np
from scipy.special import logsumexp


class MarginalLikelihoodError(ValueError):
    pass


def _positive_theta(theta):
    theta = np.asarray(theta, dtype=float)
    if theta.shape != (6,) or not np.all(np.isfinite(theta)) or np.any(theta <= 0):
        raise MarginalLikelihoodError("Need six strictly positive finite parameters")
    return theta


def log_concentration_per_record(theta, pack, indices, eta_by_subject):
    """Finite infusion 2C analytic PK with donor-level CL multiplier."""
    cl, vc, q, vp, _, _ = _positive_theta(theta)
    t, _, dose, duration, _, _ = pack
    eta = np.asarray(eta_by_subject, dtype=float)
    if not np.all(np.isfinite(eta)) or len(indices) != len(t):
        raise MarginalLikelihoodError("Invalid random effects or subject record mapping")
    if np.any(indices < 0) or np.any(indices >= len(eta)):
        raise MarginalLikelihoodError("Subject ID outside random-effect vector")
    if np.any(t <= 0) or np.any(dose <= 0) or np.any(duration <= 0):
        raise MarginalLikelihoodError("Invalid observed IV event")
    k10 = (cl / vc) * np.exp(eta[indices])
    k12, k21 = q / vc, q / vp
    trace = k10 + k12 + k21
    discriminant = np.maximum(trace * trace - 4 * k10 * k21, 0)
    diff = np.sqrt(discriminant)
    alpha = (trace + diff) / 2
    beta = (trace - diff) / 2
    if np.any(beta <= 0) or np.any(diff <= 1e-12):
        raise MarginalLikelihoodError("Non-identifiable biexponential rates")
    fraction_fast = (alpha - k21) / diff
    fraction_slow = (k21 - beta) / diff
    before = np.minimum(t, duration)
    after = np.maximum(t - duration, 0)
    fast = (-np.expm1(-alpha * before) / alpha) * np.exp(-alpha * after)
    slow = (-np.expm1(-beta * before) / beta) * np.exp(-beta * after)
    pred = 1000.0 * dose / (duration * vc) * (
        fraction_fast * fast + fraction_slow * slow
    )
    if np.any(~np.isfinite(pred)) or np.any(pred <= 0):
        raise MarginalLikelihoodError("Nonpositive marginal IV predictions")
    return np.log(pred)


def _joint(theta, data, eta):
    """Normalized log joint likelihood for each *biological person*."""
    _, people, indices, counts, log_dv = data
    pack = data[0]
    _, _, _, _, omega, sigma = _positive_theta(theta)
    eta = np.asarray(eta, dtype=float)
    if eta.shape != (len(people),):
        raise MarginalLikelihoodError("Each human must have exactly one ETA")
    predicted = log_concentration_per_record(theta, pack, indices, eta)
    rss = np.bincount(indices, weights=(predicted - log_dv) ** 2,
                      minlength=len(people))
    normal = sqrt(2 * pi)
    return (
        -0.5 * rss / sigma ** 2
        - counts * log(sigma * normal)
        - 0.5 * (eta / omega) ** 2
        - log(omega * normal)
    )


def _posterior_proposal(theta, data, iterations=12):
    """Stable vectorized posterior mode/curvature approximation.

    Multi-start grid avoids false local maxima; small-step Newton updates
    include a monotone backtracking check. All observations of each person's
    two occasions contribute to the same posterior.
    """
    _, people, _, _, _ = data
    omega = float(_positive_theta(theta)[4])
    grid = np.asarray([-3, -1.5, 0, 1.5, 3.0]) * omega
    scores = np.stack([_joint(theta, data, np.full(len(people), e)) for e in grid])
    modes = grid[np.argmax(scores, axis=0)].copy()
    difference = max(0.012 * omega, 0.0001)
    for _ in range(iterations):
        mid = _joint(theta, data, modes)
        plus = _joint(theta, data, modes + difference)
        minus = _joint(theta, data, modes - difference)
        gradient = (plus - minus) / (2 * difference)
        negative_curvature = (2 * mid - plus - minus) / difference ** 2
        negative_curvature = np.maximum(negative_curvature, 0.05 / omega ** 2)
        step = np.clip(gradient / negative_curvature, -omega, omega)
        candidate = np.clip(modes + step, -6 * omega, 6 * omega)
        step_scale = 1.0
        for _backtracking in range(9):
            proposed = modes + step_scale * (candidate - modes)
            good = _joint(theta, data, proposed) >= mid - 1e-9
            modes = np.where(good, proposed, modes)
            if np.all(good):
                break
            step_scale *= 0.5
    mid = _joint(theta, data, modes)
    plus = _joint(theta, data, modes + difference)
    minus = _joint(theta, data, modes - difference)
    curvature = (2 * mid - plus - minus) / difference ** 2
    if np.any(~np.isfinite(curvature)) or np.any(curvature <= 0):
        raise MarginalLikelihoodError("Nonpositive posterior curvature; abstain")
    scales = 1.0 / np.sqrt(curvature)
    return modes, scales


def adaptive_subject_loglik(theta, data, nodes=9):
    """Gauss-Hermite quadrature around each individual's posterior mode.

    Integrates out ETA from all occasions *jointly* and reports one natural-log
    probability density per person (not subject-independent observations).
    """
    if isinstance(nodes, bool) or not isinstance(nodes, int) or nodes < 3 or nodes % 2 == 0:
        raise MarginalLikelihoodError("Require odd Gauss-Hermite nodes >=3")
    theta = _positive_theta(theta)
    centers, scales = _posterior_proposal(theta, data)
    x, w = np.polynomial.hermite.hermgauss(nodes)
    pieces = []
    for xx, ww in zip(x, w):
        eta = centers + sqrt(2) * scales * xx
        logjoint = _joint(theta, data, eta)
        # q(eta) = N(centers, scales²). GH integrates wrt q;
        # retain the true normalized prior within logjoint.
        pieces.append(
            logjoint + xx ** 2 + np.log(scales * sqrt(2 * pi))
            + log(ww / sqrt(pi))
        )
    value = logsumexp(np.stack(pieces), axis=0)
    if np.any(~np.isfinite(value)):
        raise MarginalLikelihoodError("Nonfinite adaptive marginal density")
    return value


def adaptive_total_loglik(theta, data, nodes=9):
    return float(np.sum(adaptive_subject_loglik(theta, data, nodes)))
