"""Analytic one-/two-compartment IV finite-infusion reference.

This is a fixed-effects structural reference, NOT a mixed-effects PopPK model.
Units: administered mg, infusion hours, V in L, CL/Q in L/h,
result in plasma ng/mL. Single dose and linear elimination only.
No post-hoc individual fit and no omic effects are estimated.
"""
from __future__ import annotations

from hashlib import sha256
from math import exp, expm1, isfinite, sqrt


class StructuralPKError(ValueError):
    pass


def _positive(**values):
    for name, value in values.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise StructuralPKError(f"{name} must be a finite positive number")
        if not isfinite(value) or value <= 0:
            raise StructuralPKError(f"{name} must be a finite positive number")


def _infusion_term(k: float, t: float, duration_h: float) -> float:
    """Integral of exp(-k t) during IV infusion; stable at small k."""
    within = min(t, duration_h)
    after = max(0.0, t - duration_h)
    return (-expm1(-k * within) / k) * exp(-k * after)


def iv_1c_ng_ml(time_h: float, dose_mg: float, duration_h: float, cl_l_h: float, v_l: float) -> float:
    _positive(dose_mg=dose_mg, duration_h=duration_h, cl_l_h=cl_l_h, v_l=v_l)
    if not isinstance(time_h, (int, float)) or isinstance(time_h, bool) or not isfinite(time_h) or time_h < 0:
        raise StructuralPKError("Time must be finite and nonnegative relative to original IV dose")
    k = cl_l_h / v_l
    return 1000.0 * dose_mg / duration_h / v_l * _infusion_term(k, time_h, duration_h)


def iv_2c_ng_ml(
    time_h: float, dose_mg: float, duration_h: float,
    cl_l_h: float, vc_l: float, q_l_h: float, vp_l: float
) -> float:
    _positive(dose_mg=dose_mg, duration_h=duration_h, cl_l_h=cl_l_h,
              vc_l=vc_l, q_l_h=q_l_h, vp_l=vp_l)
    if not isinstance(time_h, (int, float)) or isinstance(time_h, bool) or not isfinite(time_h) or time_h < 0:
        raise StructuralPKError("Time must be finite and nonnegative relative to original IV dose")
    k10, k12, k21 = cl_l_h / vc_l, q_l_h / vc_l, q_l_h / vp_l
    s = k10 + k12 + k21
    discriminant = s * s - 4 * k10 * k21
    alpha = (s + sqrt(max(discriminant, 0.0))) / 2
    beta = (s - sqrt(max(discriminant, 0.0))) / 2
    if beta <= 0 or alpha <= beta:
        raise StructuralPKError("Non-identifiable biexponential rates")
    frac_alpha = (alpha - k21) / (alpha - beta)
    frac_beta = (k21 - beta) / (alpha - beta)
    return (1000.0 * dose_mg / duration_h / vc_l *
            (frac_alpha * _infusion_term(alpha, time_h, duration_h) +
             frac_beta * _infusion_term(beta, time_h, duration_h)))


def subject_disjoint_folds(subjects, n_folds: int = 5, seed: str = "omics2pk-mavo-v1") -> list[tuple[set[str], set[str]]]:
    """Deterministic disjoint human folds, not occasion or time-point splits."""
    if not isinstance(n_folds, int) or isinstance(n_folds, bool) or n_folds < 2:
        raise StructuralPKError("Require at least two whole-subject folds")
    names = [str(x) for x in subjects]
    if any(not x.strip() for x in names) or len(names) != len(set(names)):
        raise StructuralPKError("Biological subject IDs must be distinct and nonempty")
    if len(names) < n_folds:
        raise StructuralPKError("Fewer subjects than validation folds")
    ordered = sorted(names, key=lambda x: sha256((seed + "|" + x).encode("utf-8")).hexdigest())
    fold_sets = [set(ordered[i::n_folds]) for i in range(n_folds)]
    full = set(ordered)
    result = [(full - test, test) for test in fold_sets]
    if any(train & test or not train or not test for train, test in result):
        raise StructuralPKError("Invalid subject-level partition")
    if set().union(*(test for _, test in result)) != full:
        raise StructuralPKError("Subject-level coverage incomplete")
    return result
