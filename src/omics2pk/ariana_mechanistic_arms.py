"""Mechanistic and hybrid arms for the ARIANA metoprolol cohort (protocol §3).

What this module is
-------------------
The fixed-coefficient mechanistic model (arms **M4** and **M5**), the
covariate PopPK arms (**P1**, **P2**, **P3**) and the hybrid arm (**HY**) of
protocol §3, implemented so that every coefficient is either

* **transferred from a named publication** and never estimated on ARIANA, or
* **estimated on the training patients of an outer fold only**, and counted.

Every arm therefore declares its number of free parameters fitted on the
evaluated cohort, which is the requirement protocol §3 imposes on all arms.

The target
----------
The endpoint is the **observed pre-dose metoprolol concentration** in ng/mL.
It is never called an AUC, a clearance or an exposure. ``lnMR
OH-metoprolol/metoprolol`` shares its denominator with that concentration and
appears in **no** design here, in any arm; the runner enforces it.

MAP-Bayesian estimation is excluded by construction: with one concentration
per subject, conditioning the individual random effect on that concentration
and then predicting the same concentration is circular. No arm in this module
uses the evaluated patient's concentration. No NLME likelihood is claimed:
between-subject and residual variability are not separately identifiable from
one observation per subject.

Provenance of the structure
---------------------------
The model structure and every fixed constant are **transferred verbatim** from
the documented NetPK-metoprolol module of the project's prior ARIANA bundle
(``MODULES_NETPK_METOPROLOL.md``, ``R/netpk_metoprolol.R``,
``R/poppk_covariate.R``, 2026-11-07), which records the publication behind each
value. Transferring them rather than re-deriving them is deliberate: it keeps
the constants out of reach of the present evaluation. Three consequences are
reported rather than hidden:

1. The folds, seed and loss machinery are the **protocol's** (shared with the
   H1 statistical arms), not the prior bundle's, so the numbers produced here
   are a fresh run and not a reproduction of the prior bundle's figures.
2. The prior bundle promoted ``OH + SSDA`` to primary omic panel **on the
   observed loss**. That is a result-guided choice and is not admissible as a
   primary here, so the primary panel is ``OH-solanidine alone`` — the panel
   the prior bundle itself records as mechanistically prespecified, because
   4-OH-solanidine has an authentic synthesised standard and SSDA does not.
   The other panels run as declared sensitivities.
3. The Rüdesheim activity-score grid does not contain every activity score
   ARIANA reports. Interpolating onto it is a **decision**, and
   :func:`activity_score_grid_audit` measures how many patients it affects so
   the decision is reported instead of being made silently.
"""
from __future__ import annotations

from typing import Mapping, Sequence

#: Protocol §3 requires a literature source for every fixed constant. Each
#: entry is (module, the fixed coefficient, the publication). Nothing in this
#: table is estimated on ARIANA.
LITERATURE_CONSTANTS = (
    (
        "structural backbone",
        "one- or two-compartment oral steady state, tau = 12 h; ka, V/F and "
        "the distribution terms per backbone (see BACKBONES)",
        "Chai et al. 2026 Clin Transl Sci doi 10.1111/cts.70677; Eugene 2016 "
        "doi 10.19070/2167-910X-1600035; FDA 019962s033 (2006) Toprol-XL",
    ),
    (
        "CYP2D6 genotype",
        "fraction metabolised by CYP2D6 f_m = 0.80; relative catalytic "
        "activity 0 / 0.19 / 0.64 / 0.72 / 1.00 / 2.13 at activity score "
        "0 / 0.5 / 1.25 / 1.5 / 2 / 3",
        "Rüdesheim et al. 2020 Pharmaceutics 12:1200 "
        "doi 10.3390/pharmaceutics12121200, Table 3",
    ),
    (
        "phenoconversion",
        "activity score reduced by 0.54 per documented CYP2D6 "
        "substrate/inhibitor",
        "Müller/Sarömba et al. 2025 Br J Clin Pharmacol, external geriatric "
        "cohort n=88, PMC12272519",
    ),
    (
        "omic CYP2D6 phenotype",
        "relative activity = exp(gamma * (lnMR_i - lnMR_reference)), "
        "gamma = 1, capped at the published ultrarapid activity 2.13",
        "Medwid et al. 2024 Clin Pharmacol Ther doi 10.1002/cpt.3380; "
        "Müller et al. 2024 doi 10.1002/cpt.3106",
    ),
    (
        "allometry",
        "CL/F proportional to (weight/70)^0.75; V/F to (weight/70)^1",
        "Anderson & Holford, standard allometric scaling",
    ),
    (
        "sex",
        "CL/F male = CL/F female * exp(0.572)",
        "Eugene 2016 doi 10.19070/2167-910X-1600035",
    ),
)

# --- transferred constants -------------------------------------------------
#: Rüdesheim 2020 Table 3. The reference is activity score 2.
RUDESHEIM_AS_GRID = (0.0, 0.5, 1.25, 1.5, 2.0, 3.0)
RUDESHEIM_RELATIVE_ACTIVITY = (0.0, 0.19, 0.64, 0.72, 1.00, 2.13)
#: Fraction of metoprolol metabolism carried by CYP2D6 in a normal metaboliser.
FM_CYP2D6 = 0.80
#: Activity-score decrement per documented CYP2D6 substrate/inhibitor.
DELTA_AS_PER_INHIBITOR = 0.54
#: log-scale sex effect on apparent clearance (male versus female).
SEX_LOG_EFFECT = 0.572
ALLOMETRIC_EXPONENT_CL = 0.75
ALLOMETRIC_EXPONENT_V = 1.0
#: Structural exponent of the omic module: a metabolic ratio of a formation
#: clearance is to first order proportional to in vivo enzyme activity.
OMIC_GAMMA = 1.0
#: Genotype/omic weighting in the combined mechanistic activity. Fixed at 0.5
#: before the prior bundle's sensitivity sweep, which found the optimum near
#: 0.75: the fixed value is deliberately NOT the optimum.
BLEND_WEIGHT = 0.5
#: Activity-score window defining a genotypic normal metaboliser, used to set
#: the omic reference level inside the training fold.
NORMAL_METABOLISER_AS = (1.25, 2.25)
#: Minimum number of training normal metabolisers needed for a reference.
MIN_REFERENCE_PATIENTS = 3
#: Dosing interval and sampling time assumed throughout: two equal daily doses
#: and a trough drawn before the next. The source documents neither the clock
#: time nor the time since the last dose, so this is an assumption.
TAU_HOURS = 12.0
TIME_AFTER_DOSE_HOURS = 12.0
#: Salt equivalence: 47.5 mg succinate carries the metoprolol base of 50 mg
#: tartrate. A prespecified sensitivity, not applied by default.
SUCCINATE_TARTRATE_EQUIVALENCE = 1.0526

#: Published structural backbones. ER is the primary scenario, chosen on the
#: pharmaceutical argument that 47.5 and 95 mg unit doses are the German
#: extended-release succinate presentations, never on an observed loss.
BACKBONES: Mapping[str, Mapping[str, object]] = {
    "FDA2006_ER": {
        "name": "FDA 2006 Toprol-XL, two-compartment oral extended release",
        "structure": "2C",
        "ka": 0.0475,
        "tlag": 0.818,
        "v_ref": 180.0,
        "q": 741.0,
        "vp": 661.0,
        "source": "FDA 019962s033 (2006)",
        "role": "primary scenario",
    },
    "Chai2026_IR": {
        "name": "Chai 2026, one-compartment oral immediate release, elderly, troughs",
        "structure": "1C",
        "ka": 0.235,
        "tlag": 0.0,
        "v_ref": 562.0,
        "source": "doi 10.1111/cts.70677",
        "role": "structural sensitivity",
    },
    "Eugene2016_IR": {
        "name": "Eugene 2016, one-compartment oral immediate release, elderly",
        "structure": "1C",
        "ka": 0.235,
        "tlag": 0.469,
        "v_ref": 38.0,
        "source": "doi 10.19070/2167-910X-1600035",
        "role": "structural sensitivity",
    },
}
PRIMARY_BACKBONE = "FDA2006_ER"

#: Omic panels. The primary is OH alone: 4-OH-solanidine has an authentic
#: synthesised standard, SSDA is semi-quantitative without one. The prior
#: bundle's promotion of OH+SSDA to primary was made on the observed loss and
#: is therefore carried here as a sensitivity, not as the primary.
OMIC_PANELS: Mapping[str, Sequence[str]] = {
    "oh_only": ("lnMR OH-solanidine/solanidine",),
    "ssda_only": ("lnMR SSDA/solanidine",),
    "oh_and_ssda": (
        "lnMR OH-solanidine/solanidine",
        "lnMR SSDA/solanidine",
    ),
}
PRIMARY_OMIC_PANEL = "oh_only"

#: Free parameters fitted on ARIANA, per arm. Protocol §3 requires this count.
FREE_PARAMETERS = {
    "M4_mechanistic_pgx": 1,
    "M5_mechanistic_pgx_omic": 1,
    "M4c_backbone_only": 1,
    "M5c_mechanistic_omic_only": 1,
    "P1_poppk_published_coefficients": 1,
    "P2_poppk_estimated_coefficients": 3,
    "P3_poppk_estimated_plus_omic": 4,
    "HY_poppk_published_plus_ml_residual": "1 structural + the ridge penalty "
    "and coefficients of the residual stage, all fitted inside the training "
    "fold",
}

#: Which mechanistic activity each arm uses.
MECHANISTIC_ARMS = {
    "M4c_backbone_only": "none",
    "M4_mechanistic_pgx": "pgx",
    "M5c_mechanistic_omic_only": "omic",
    "M5_mechanistic_pgx_omic": "pgx_omic",
}
#: Protocol §4 key secondary contrast.
KEY_SECONDARY_CONTRAST = ("M5_mechanistic_pgx_omic", "M4_mechanistic_pgx")

#: Small offset keeping log(activity) finite for a poor metaboliser.
_ACTIVITY_EPS = 1e-6
#: Offset in the PopPK covariate design, same purpose.
POPPK_LOG_ACTIVITY_EPS = 0.05


class MechanisticError(RuntimeError):
    """Raised when a transferred constant or a prediction is out of domain."""


# --- the CYP2D6 genotype module -------------------------------------------
def activity_score_grid_audit(activity_scores) -> dict:
    """Report, rather than hide, the activity scores absent from the grid.

    The Rüdesheim grid predates the post-2021 CPIC change that scores ``*41``
    as 0.25, and it contains no entry at 1.0. Any ARIANA activity score off
    the grid receives a **linearly interpolated** relative activity, which is
    an alignment decision and not a published value. This function measures
    how much of the cohort that decision touches.
    """
    import numpy as np

    scores = np.asarray(activity_scores, dtype=float)
    on_grid = np.isin(scores, np.asarray(RUDESHEIM_AS_GRID, dtype=float))
    off = sorted({float(value) for value in scores[~on_grid]})
    return {
        "published_grid": list(RUDESHEIM_AS_GRID),
        "published_relative_activity": list(RUDESHEIM_RELATIVE_ACTIVITY),
        "observed_scores": sorted({float(value) for value in scores}),
        "scores_absent_from_the_published_grid": off,
        "n_patients_on_grid": int(on_grid.sum()),
        "n_patients_interpolated": int((~on_grid).sum()),
        "fraction_interpolated": float((~on_grid).mean()),
        "rule_applied": (
            "monotone linear interpolation between the published grid points, "
            "constant extrapolation beyond activity score 3"
        ),
        "this_is_a_decision_not_a_published_value": (
            "An interpolated relative activity is not a value Rüdesheim 2020 "
            "reports. It is declared here because protocol §8 forbids aligning "
            "activity-score conventions silently. Reading an arm's result as "
            "evidence about the published mapping is therefore not admissible "
            "for the interpolated share of the cohort."
        ),
    }


def relative_activity_from_activity_score(activity_scores):
    """Relative CYP2D6 catalytic activity from the published grid."""
    import numpy as np

    scores = np.asarray(activity_scores, dtype=float)
    if not np.all(np.isfinite(scores)):
        raise MechanisticError("non-finite activity score")
    if np.any(scores < 0):
        raise MechanisticError("negative activity score")
    capped = np.minimum(scores, max(RUDESHEIM_AS_GRID))
    return np.interp(capped, RUDESHEIM_AS_GRID, RUDESHEIM_RELATIVE_ACTIVITY)


def effective_activity_score(activity_scores, inhibitor,
                             delta_as=DELTA_AS_PER_INHIBITOR):
    """Phenoconversion: the activity score net of inhibitor intake."""
    import numpy as np

    scores = np.asarray(activity_scores, dtype=float)
    flag = np.asarray(inhibitor, dtype=float)
    return np.maximum(0.0, scores - delta_as * flag)


def clearance_multiplier(relative_activity, fm=FM_CYP2D6):
    """Apparent-clearance multiplier ``(1 - f_m) + f_m * a_rel``.

    At ``a_rel = 0`` this gives 0.20, i.e. a five-fold exposure increase in a
    poor metaboliser, which is the published order of magnitude. The poor
    metaboliser's higher bioavailability is therefore **already absorbed** in
    this apparent parameterisation and must not be added a second time.
    """
    import numpy as np

    activity = np.asarray(relative_activity, dtype=float)
    if not np.all(np.isfinite(activity)) or np.any(activity < 0):
        raise MechanisticError("invalid relative activity")
    if not 0.0 < fm < 1.0:
        raise MechanisticError("fraction metabolised out of domain")
    return (1.0 - fm) + fm * activity


# --- the omic module -------------------------------------------------------
def omic_reference(train_activity_scores, train_inhibitor, train_ratios):
    """Normal-metaboliser reference level per omic ratio, training fold only.

    The ratios are expressed in peak units and are not comparable between
    analytes, so each ratio keeps its own reference. The reference never uses
    the metoprolol concentration and never uses a test patient.
    """
    import numpy as np

    scores = np.asarray(train_activity_scores, dtype=float)
    flag = np.asarray(train_inhibitor, dtype=float)
    ratios = np.atleast_2d(np.asarray(train_ratios, dtype=float))
    if ratios.shape[0] != scores.size:
        ratios = ratios.T
    low, high = NORMAL_METABOLISER_AS
    genotypic_nm = (scores >= low) & (scores <= high)
    selected = genotypic_nm & (flag == 0)
    fallback_used = False
    if selected.sum() < MIN_REFERENCE_PATIENTS:
        selected = genotypic_nm
        fallback_used = True
    if selected.sum() < MIN_REFERENCE_PATIENTS:
        raise MechanisticError(
            "fewer than three training normal metabolisers: the omic "
            "reference is not definable on this fold"
        )
    reference = ratios[selected].mean(axis=0)
    if not np.all(np.isfinite(reference)):
        raise MechanisticError("non-finite omic reference")
    return {
        "reference": reference,
        "n_reference_patients": int(selected.sum()),
        "inhibitor_free_requirement_relaxed": bool(fallback_used),
    }


def relative_activity_from_omic(ratios, reference, gamma=OMIC_GAMMA,
                                cap=max(RUDESHEIM_RELATIVE_ACTIVITY)):
    """Relative CYP2D6 activity read from one or more metabolic ratios.

    With several ratios the activity is their **geometric mean**: these are
    repeated noisy readings of one quantity, in vivo CYP2D6 activity. No
    weight is fitted on the data.
    """
    import numpy as np

    matrix = np.atleast_2d(np.asarray(ratios, dtype=float))
    ref = np.asarray(reference, dtype=float).ravel()
    if matrix.shape[1] != ref.size:
        matrix = matrix.T
    if matrix.shape[1] != ref.size:
        raise MechanisticError("one reference per omic ratio is required")
    if not np.all(np.isfinite(matrix)):
        raise MechanisticError("non-finite omic ratio")
    log_activity = gamma * (matrix - ref[None, :])
    activity = np.exp(log_activity.mean(axis=1))
    return np.clip(activity, 0.0, cap)


def mechanistic_activity(frame_arrays, arm, reference, gamma=OMIC_GAMMA,
                         blend_weight=BLEND_WEIGHT,
                         delta_as=DELTA_AS_PER_INHIBITOR):
    """Relative CYP2D6 activity for one mechanistic arm.

    ``frame_arrays`` carries ``activity_score``, ``inhibitor`` and ``ratios``.
    """
    import numpy as np

    scores = np.asarray(frame_arrays["activity_score"], dtype=float)
    if arm == "none":
        return np.ones(scores.size)
    activity_pgx = relative_activity_from_activity_score(
        effective_activity_score(scores, frame_arrays["inhibitor"], delta_as)
    )
    if arm == "pgx":
        return activity_pgx
    activity_omic = relative_activity_from_omic(
        frame_arrays["ratios"], reference, gamma
    )
    if arm == "omic":
        return activity_omic
    if arm == "pgx_omic":
        return np.exp(
            (1.0 - blend_weight) * np.log(activity_pgx + _ACTIVITY_EPS)
            + blend_weight * np.log(activity_omic + _ACTIVITY_EPS)
        )
    raise MechanisticError(f"unknown mechanistic arm {arm!r}")


# --- the structural backbone ----------------------------------------------
def _steady_state_kernel(rate, ka, u, tau):
    """Steady-state first-order depot kernel, numerically stabilised."""
    import numpy as np

    def f(a):
        return np.exp(-a * u) / (-np.expm1(-a * tau))

    if abs(ka - rate) < 1e-5:
        h = 1e-5
        return -(f(rate + h) - f(rate - h)) / (2.0 * h)
    return (f(rate) - f(ka)) / (ka - rate)


def steady_state_trough(dose_per_administration_mg, clearance, volume,
                        backbone, time_after_dose_h=TIME_AFTER_DOSE_HOURS,
                        tau_h=TAU_HOURS):
    """Steady-state trough concentration, ng/mL, for a 1C or 2C oral model."""
    import numpy as np

    dose = np.atleast_1d(np.asarray(dose_per_administration_mg, dtype=float))
    cl = np.atleast_1d(np.asarray(clearance, dtype=float))
    v = np.atleast_1d(np.asarray(volume, dtype=float))
    size = max(dose.size, cl.size, v.size)
    dose = np.resize(dose, size)
    cl = np.resize(cl, size)
    v = np.resize(v, size)
    for name, values in (("dose", dose), ("clearance", cl), ("volume", v)):
        if not np.all(np.isfinite(values)):
            raise MechanisticError(f"non-finite {name} reached the backbone")
        if np.any(values <= 0.0):
            raise MechanisticError(f"non-positive {name} reached the backbone")
    ka = float(backbone["ka"])
    u = (time_after_dose_h - float(backbone["tlag"])) % tau_h

    values = np.empty(size, dtype=float)
    for index in range(size):
        if backbone["structure"] == "1C":
            values[index] = (
                dose[index] * 1000.0 * ka / v[index]
                * _steady_state_kernel(cl[index] / v[index], ka, u, tau_h)
            )
            continue
        k10 = cl[index] / v[index]
        k12 = float(backbone["q"]) / v[index]
        k21 = float(backbone["q"]) / float(backbone["vp"])
        total = k10 + k12 + k21
        product = k10 * k21
        discriminant = total * total - 4.0 * product
        if discriminant < 0.0:
            raise MechanisticError("two-compartment roots are not real")
        root = discriminant ** 0.5
        alpha = 0.5 * (total + root)
        beta = 0.5 * (total - root)
        if alpha <= beta:
            raise MechanisticError("degenerate two-compartment roots")
        a_coef = (alpha - k21) / (v[index] * (alpha - beta))
        b_coef = (k21 - beta) / (v[index] * (alpha - beta))
        values[index] = dose[index] * 1000.0 * ka * (
            a_coef * _steady_state_kernel(alpha, ka, u, tau_h)
            + b_coef * _steady_state_kernel(beta, ka, u, tau_h)
        )
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise MechanisticError(
            "invalid mechanistic concentration: check the transferred "
            "structural parameters"
        )
    return values


def structural_prediction(arrays, baseline_clearance, arm, backbone, reference,
                          gamma=OMIC_GAMMA, blend_weight=BLEND_WEIGHT,
                          fm=FM_CYP2D6, delta_as=DELTA_AS_PER_INHIBITOR,
                          dose_equivalence=1.0,
                          time_after_dose_h=TIME_AFTER_DOSE_HOURS,
                          tau_h=TAU_HOURS):
    """Predicted trough for one mechanistic arm, given a locked baseline."""
    import numpy as np

    activity = mechanistic_activity(arrays, arm, reference, gamma,
                                    blend_weight, delta_as)
    weight = np.asarray(arrays["weight"], dtype=float)
    clearance = (
        baseline_clearance
        * clearance_multiplier(activity, fm)
        * (weight / 70.0) ** ALLOMETRIC_EXPONENT_CL
        * np.exp(SEX_LOG_EFFECT * np.asarray(arrays["male"], dtype=float))
    )
    volume = float(backbone["v_ref"]) * (weight / 70.0) ** ALLOMETRIC_EXPONENT_V
    dose = np.asarray(arrays["daily_dose_mg"], dtype=float) / 2.0 * dose_equivalence
    return steady_state_trough(dose, clearance, volume, backbone,
                               time_after_dose_h, tau_h)


def lock_baseline_clearance(train_arrays, arm, backbone, reference, **kwargs):
    """Fit the single free parameter on the training patients, then lock it.

    The reference apparent clearance (activity score 2, 70 kg, female) is the
    **only** quantity estimated on the cohort: one degree of freedom, by least
    squares on the log residuals of the training patients of this outer fold.
    """
    import numpy as np
    from scipy.optimize import minimize_scalar

    observed_log = np.log(np.asarray(train_arrays["observed"], dtype=float))

    def objective(log_baseline):
        try:
            predicted = structural_prediction(
                train_arrays, float(np.exp(log_baseline)), arm, backbone,
                reference, **kwargs
            )
        except MechanisticError:
            return 1e10
        return float(np.sum((observed_log - np.log(predicted)) ** 2))

    fit = minimize_scalar(
        objective, bounds=(float(np.log(1.0)), float(np.log(2000.0))),
        method="bounded", options={"xatol": 1e-10},
    )
    baseline = float(np.exp(fit.x))
    if not np.isfinite(baseline) or baseline <= 0.0:
        raise MechanisticError("baseline clearance locking failed")
    return baseline


# --- the covariate PopPK arms ---------------------------------------------
def poppk_design(arrays, reference, delta_as=DELTA_AS_PER_INHIBITOR):
    """Covariate design shared by P1, P2 and P3.

    ``log CL = theta0 + 0.75*log(weight/70) + theta_sex*male
               + theta_as*log(a_rel + eps) [+ theta_omic*(lnMR - lnMR_ref)]``
    """
    import numpy as np

    scores = np.asarray(arrays["activity_score"], dtype=float)
    activity = relative_activity_from_activity_score(
        effective_activity_score(scores, arrays["inhibitor"], delta_as)
    )
    ratios = np.atleast_2d(np.asarray(arrays["ratios"], dtype=float))
    ref = np.asarray(reference, dtype=float).ravel()
    if ratios.shape[1] != ref.size:
        ratios = ratios.T
    return {
        "log_weight": np.log(np.asarray(arrays["weight"], dtype=float) / 70.0),
        "male": np.asarray(arrays["male"], dtype=float),
        "relative_activity": activity,
        "log_relative_activity": np.log(activity + POPPK_LOG_ACTIVITY_EPS),
        "omic_centred": (ratios - ref[None, :]).mean(axis=1),
    }


def poppk_clearance(design, theta, published=False, fm=FM_CYP2D6):
    """Apparent clearance from the covariate model."""
    import numpy as np

    base = np.exp(theta["log_cl0"]) * np.exp(
        ALLOMETRIC_EXPONENT_CL * design["log_weight"]
    )
    if published:
        return (
            base
            * np.exp(SEX_LOG_EFFECT * design["male"])
            * clearance_multiplier(design["relative_activity"], fm)
        )
    return (
        base
        * np.exp(theta["sex"] * design["male"])
        * np.exp(theta["activity_score"] * design["log_relative_activity"])
        * np.exp(theta.get("omic", 0.0) * design["omic_centred"])
    )


def poppk_predict(arrays, design, theta, backbone, published=False,
                  time_after_dose_h=TIME_AFTER_DOSE_HOURS, tau_h=TAU_HOURS,
                  fm=FM_CYP2D6):
    import numpy as np

    weight = np.asarray(arrays["weight"], dtype=float)
    clearance = poppk_clearance(design, theta, published, fm)
    volume = float(backbone["v_ref"]) * (weight / 70.0) ** ALLOMETRIC_EXPONENT_V
    dose = np.asarray(arrays["daily_dose_mg"], dtype=float) / 2.0
    return steady_state_trough(dose, clearance, volume, backbone,
                               time_after_dose_h, tau_h)


def poppk_fit(train_arrays, train_design, backbone, with_omic=False,
              time_after_dose_h=TIME_AFTER_DOSE_HOURS, tau_h=TAU_HOURS):
    """Estimate the covariate coefficients on the training patients only.

    This is the one-observation-per-subject equivalent of a population
    regression. No individual random effect is estimated, so there is no
    MAP-Bayesian step and no claim of an NLME fit.
    """
    import numpy as np
    from scipy.optimize import minimize

    names = ["log_cl0", "sex", "activity_score"]
    if with_omic:
        names.append("omic")
    start = {"log_cl0": float(np.log(110.0)), "sex": 0.5,
             "activity_score": 0.5, "omic": -0.3}
    observed_log = np.log(np.asarray(train_arrays["observed"], dtype=float))

    def objective(vector):
        theta = dict(zip(names, (float(value) for value in vector)))
        try:
            predicted = poppk_predict(train_arrays, train_design, theta,
                                      backbone, False, time_after_dose_h, tau_h)
        except MechanisticError:
            return 1e10
        return float(np.sum((observed_log - np.log(predicted)) ** 2))

    start_vector = [start[name] for name in names]
    fit = minimize(objective, start_vector, method="BFGS",
                   options={"maxiter": 400, "gtol": 1e-6})
    # A solver's success flag is a statement about its own stopping rule, not
    # about the estimate. The objective is differentiated by finite differences,
    # so a tight gradient tolerance is unreachable and the flag goes false on a
    # perfectly good optimum. Convergence is therefore established by refitting
    # with a derivative-free method and measuring the disagreement.
    refit = minimize(objective, fit.x, method="Nelder-Mead",
                     options={"maxiter": 20000, "xatol": 1e-10, "fatol": 1e-12})
    best = refit if refit.fun < fit.fun else fit
    theta = dict(zip(names, (float(value) for value in best.x)))
    discrepancy = float(np.max(np.abs(np.asarray(fit.x) - np.asarray(refit.x))))
    return {
        "theta": theta,
        "solver_success_flag": bool(fit.success),
        "cross_check_max_parameter_difference": discrepancy,
        "cross_check_objective_gain": float(fit.fun - best.fun),
        "converged": bool(discrepancy < 1e-4),
        "sum_of_squares": float(best.fun),
        "n_free_parameters": len(names),
    }
