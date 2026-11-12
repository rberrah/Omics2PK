"""The project's first clinical-pharmacokinetic contrast outside ARIANA.

Source: Lee J, et al. "Interindividual variability in CYP3A-mediated venetoclax
metabolism in vitro and in vivo." Clin Transl Sci 2024;17:e70106.
PMCID PMC11645475, CC BY-NC.

Twelve patients with chronic lymphocytic leukaemia, each with, per patient and
per visit:

  * venetoclax daily dose and venetoclax plasma concentration (main text);
  * plasma 4-beta-hydroxycholesterol, cholesterol and their ratio
    (supplement Table S3, reached through the BioStudies / EBI FTP channel).

That is a **clinical** drug concentration paired with a **plasma metabolite**
in the same individuals, which until now the project held only in ARIANA.

No new predictive architecture. The ridge, the inner penalty selection, the
loss triple, the stored-error bootstrap and the permutation rank are imported
from the frozen ``protein_clint_pilot``; the refitting out-of-bag interval is
imported from ``plasma_biomarker_validation``.

What this cohort is honest about
--------------------------------
**The baseline is the mean, not a clinical model.** The article publishes
patient characteristics -- age, weight, creatinine clearance, bilirubin,
transaminases -- **by sex group only**, so no per-patient covariate exists. The
comparator therefore cannot be "clinical plus pharmacogenetics", which is what
the funding question actually asks. A gain measured here would be a gain over
**no information at all**, which is a far weaker and far easier bar. That makes
this design **biased in favour of the biomarker**, and any non-detection
correspondingly harder to dismiss.

**The endpoint is a single timed concentration normalised by dose.** Not an
AUC, not a clearance, not a steady-state exposure metric. It is the reciprocal
of the quantity the article calls "estimated venetoclax CL", so it is a
concentration-scale reading of the same ratio and not an independent endpoint.

**The article already reports the marginal answer.** It states that the
dose-corrected trough was not significantly correlated with the plasma
4-beta-HC/cholesterol ratio. This module adds what the article does not
provide: an out-of-sample reading, a permutation null that refits, and an
interval that refits. A negative here is a confirmation with controls, not a
discovery.
"""
from __future__ import annotations

import numpy as np

from . import plasma_biomarker_validation as PBV
from . import protein_clint_pilot as PCP

FROZEN = {
    "seed": 20261110,
    "validation": "leave_one_patient_out",
    "ridge_alphas": PCP.FROZEN["ridge_alphas"],
    "permutation_replicates": 200,
    "n_bootstrap": 2000,
    "n_out_of_bag_bootstrap": 2000,
    "primary_analysis": "visit_1_only",
    "sensitivity_analysis": "both_visits_patient_held_together",
    "endpoint": "log10_dose_normalised_plasma_concentration",
    "primary_arm": "C1_plasma_metabolite",
    "primary_comparator": "C0_mean",
    "declared_family_size": 2,
    "status_ceiling": "exploratory_clinical_pk",
}

ARMS = {
    "C0_mean": (),
    "C1_plasma_metabolite": ("plasma_4b_hc_cholesterol_ratio",),
}

#: The three figures the article publishes about these patients' exposure, and
#: which any extraction of the dose and concentration columns must reproduce.
#: The article's "estimated venetoclax CL" is proportional to dose divided by
#: concentration -- its own method is in a Supporting Information file we do
#: not hold -- so the FOLD RANGE is reproducible exactly, and the median and
#: the bounds are reproducible once the single proportionality constant is
#: recovered from the median.
PUBLISHED_ANCHORS = {
    "fold_variation_in_dose_over_concentration": 28.8,
    "median_estimated_cl_l_per_h": 4.0,
    "minimum_estimated_cl_l_per_h": 1.02,
    "maximum_estimated_cl_l_per_h": 29.4,
    "highest_two_are_patients_3_and_7_at_visit_2": True,
    "tolerance_fold": 0.1,
    "tolerance_cl": 0.15,
}


class VenetoclaxContrastError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise VenetoclaxContrastError(message)


def dose_normalised_target(rows, keys):
    """log10(concentration / daily dose), the endpoint.

    Dividing by the dose is not a modelling choice dressed up as a transform:
    three patients had their dose reduced between visits, so an unnormalised
    concentration would carry the prescriber's response to toxicity as if it
    were pharmacokinetics.
    """
    values = []
    for key in keys:
        record = rows[key]
        _check(record["concentration_nm"] > 0.0,
               f"{key} has a non-positive concentration; a not-detected value "
               "must be handled by a declared rule, never by a silent drop")
        _check(record["dose_mg_per_day"] > 0.0, f"{key} has no dose")
        values.append(np.log10(record["concentration_nm"]
                               / record["dose_mg_per_day"]))
    return np.asarray(values, dtype=float)


def verify_published_anchors(rows):
    """Reproduce the article's own exposure figures, or raise."""
    ratios = {key: record["dose_mg_per_day"] / record["concentration_nm"]
              for key, record in rows.items()}
    values = sorted(ratios.values())
    fold = values[-1] / values[0]
    expected = PUBLISHED_ANCHORS["fold_variation_in_dose_over_concentration"]
    _check(abs(fold - expected) <= PUBLISHED_ANCHORS["tolerance_fold"],
           f"the dose-over-concentration fold range is {fold:.2f}, and the "
           f"article publishes {expected}. The extraction of the dose or the "
           "concentration column is wrong.")

    n = len(values)
    median = (values[n // 2] if n % 2
              else (values[n // 2 - 1] + values[n // 2]) / 2)
    scale = PUBLISHED_ANCHORS["median_estimated_cl_l_per_h"] / median
    low, high = values[0] * scale, values[-1] * scale
    for name, got in (("minimum_estimated_cl_l_per_h", low),
                      ("maximum_estimated_cl_l_per_h", high)):
        want = PUBLISHED_ANCHORS[name]
        _check(abs(got - want) <= PUBLISHED_ANCHORS["tolerance_cl"],
               f"{name} recomputes to {got:.2f} against a published {want}")

    top_two = sorted(ratios, key=lambda k: ratios[k], reverse=True)[:2]
    _check(set(top_two) == {("3", 2), ("7", 2)},
           f"the two highest clearances are {top_two}; the article names "
           "patients 3 and 7 at visit 2")
    return {
        "fold_variation_recomputed": round(fold, 2),
        "fold_variation_published": expected,
        "implied_proportionality_constant": round(scale, 3),
        "minimum_recomputed": round(low, 2),
        "maximum_recomputed": round(high, 2),
        "median_used_to_recover_the_constant": (
            PUBLISHED_ANCHORS["median_estimated_cl_l_per_h"]),
        "highest_two_recomputed": [list(k) for k in sorted(top_two)],
        "what_this_proves": (
            "that the dose and concentration columns were read correctly. It "
            "does NOT validate the article's clearance method, which is in a "
            "Supporting Information file this project does not hold -- the "
            "constant is recovered from the published median rather than "
            "derived."),
    }


def leave_one_patient_out(rows, keys, columns):
    """Out-of-sample prediction holding out ALL observations of one patient.

    With two visits per patient, holding out an observation rather than a
    patient would train on the same individual it scores. The resampling and
    permutation units are the patient for the same reason.
    """
    patients = sorted({key[0] for key in keys})
    _check(len(patients) >= 4,
           f"need at least four patients, got {len(patients)}")
    target = dose_normalised_target(rows, keys)
    block = (np.array([[rows[k][c] for c in columns] for k in keys],
                      dtype=float) if columns else np.zeros((len(keys), 0)))
    _check(np.isfinite(block).all(), f"the design for {columns} has a gap")
    predicted = np.zeros(len(keys))
    for patient in patients:
        held = np.array([i for i, k in enumerate(keys) if k[0] == patient])
        keep = np.array([i for i, k in enumerate(keys) if k[0] != patient])
        if columns:
            alpha = PCP.select_penalty(block[keep], target[keep],
                                       FROZEN["ridge_alphas"])
            predicted[held] = PCP._ridge_fit_predict(
                block[keep], target[keep], block[held], alpha)
        else:
            predicted[held] = float(target[keep].mean())
    return predicted, np.abs(predicted - target), target


def score_arm(rows, keys, arm):
    _check(arm in ARMS, f"{arm!r} is not a frozen arm")
    predicted, errors, target = leave_one_patient_out(rows, keys, ARMS[arm])
    native = np.abs(10.0 ** predicted - 10.0 ** target)
    out = PCP.losses(errors, native)
    out["arm"] = arm
    out["predictors"] = list(ARMS[arm])
    out["native_units"] = "nM of venetoclax per mg of daily dose"
    out["n_patients"] = len({key[0] for key in keys})
    out["n_observations"] = len(keys)
    return out, errors


def patient_level_bootstrap(rows, keys, *, n_boot=None, seed=None, alpha=0.05,
                            family_size=1):
    """Refitting out-of-bag interval, resampling PATIENTS and not observations."""
    n_boot = (FROZEN["n_out_of_bag_bootstrap"] if n_boot is None else n_boot)
    seed = FROZEN["seed"] if seed is None else seed
    patients = sorted({key[0] for key in keys})
    by_patient = {p: [i for i, k in enumerate(keys) if k[0] == p]
                  for p in patients}
    target = dose_normalised_target(rows, keys)
    block = np.array([[rows[k][c] for c in ARMS["C1_plasma_metabolite"]]
                      for k in keys], dtype=float)
    rng = np.random.default_rng(seed)
    draws, skipped = [], 0
    for _ in range(n_boot):
        drawn = rng.integers(0, len(patients), size=len(patients))
        in_bag = [patients[i] for i in drawn]
        out_of_bag = [p for p in patients if p not in set(in_bag)]
        train = [i for p in in_bag for i in by_patient[p]]
        test = [i for p in out_of_bag for i in by_patient[p]]
        if len(out_of_bag) < 2 or len(set(in_bag)) < 3:
            skipped += 1
            continue
        train, test = np.array(train), np.array(test)
        chosen = PCP.select_penalty(block[train], target[train],
                                    FROZEN["ridge_alphas"])
        arm = np.abs(PCP._ridge_fit_predict(block[train], target[train],
                                            block[test], chosen)
                     - target[test]).mean()
        base = np.abs(float(target[train].mean()) - target[test]).mean()
        draws.append(arm - base)
    _check(len(draws) >= 100,
           f"only {len(draws)} of {n_boot} draws left enough patients out")
    array = np.asarray(draws, dtype=float)
    half = 100.0 * alpha / (2.0 * family_size)
    low, high = (float(v) for v in np.percentile(array, [half, 100.0 - half]))
    return {
        "mean_difference_log10": float(array.mean()),
        "interval": [low, high],
        "interval_excludes_zero": bool(low * high > 0.0),
        "level": f"{100.0 * (1.0 - alpha / family_size):.2f} %",
        "family_size": int(family_size),
        "n_draws_used": int(array.size),
        "n_draws_skipped_too_small": int(skipped),
        "refits_the_model_inside_each_draw": True,
        "resampling_unit": "patient",
        "why_the_patient_and_not_the_observation": (
            "two visits of one patient are not two individuals. Resampling "
            "observations would treat them as such and narrow the interval "
            "without adding information."),
        "exploratory": True,
    }


def permutation_null(rows, keys, *, replicates=None, seed=None):
    """Permute the biomarker across PATIENTS, carrying each patient's visits.

    Permuting observations would destroy the within-patient structure as well
    as the between-patient association, which would make the null easier to
    beat than the design deserves. Carrying a patient's whole visit vector to
    another patient destroys only the association being tested.
    """
    replicates = (FROZEN["permutation_replicates"] if replicates is None
                  else replicates)
    seed = FROZEN["seed"] if seed is None else seed
    patients = sorted({key[0] for key in keys})
    rng = np.random.default_rng(seed + 1)
    column = ARMS["C1_plasma_metabolite"][0]
    values = []
    for _ in range(replicates):
        order = rng.permutation(len(patients))
        mapping = {patients[i]: patients[order[i]]
                   for i in range(len(patients))}
        shuffled = {}
        for key in keys:
            donor, visit = key
            source = (mapping[donor], visit)
            record = dict(rows[key])
            record[column] = rows.get(source, rows[key])[column]
            shuffled[key] = record
        _, arm_errors = score_arm(shuffled, keys, "C1_plasma_metabolite")
        _, base_errors = score_arm(shuffled, keys, "C0_mean")
        values.append(float((arm_errors - base_errors).mean()))
    return values


def contrast(rows, keys, *, family_size=1, seed=None):
    arm_score, arm_errors = score_arm(rows, keys, "C1_plasma_metabolite")
    base_score, base_errors = score_arm(rows, keys, "C0_mean")
    stored = PCP.paired_bootstrap(arm_errors, base_errors, seed=seed,
                                  family_size=family_size)
    stored["resampling_unit"] = (
        "observation, which OVERSTATES independence when a patient "
        "contributes two visits; the patient-level refitting interval is the "
        "one that governs")
    refitting = patient_level_bootstrap(rows, keys, seed=seed,
                                        family_size=family_size)
    null = permutation_null(rows, keys, seed=seed)
    control = PCP.permutation_rank(stored["mean_difference_log10"], null)
    control["null_deciles"] = PBV.deciles(null)
    control.pop("null_min", None)
    control.pop("null_max", None)
    supported = (refitting["interval_excludes_zero"]
                 and refitting["mean_difference_log10"] < 0.0
                 and control["rank_p_value"] <= 0.05)
    return {
        "arm_score": arm_score,
        "comparator_score": base_score,
        "stored_error_bootstrap": stored,
        "refitting_patient_level_bootstrap": refitting,
        "permutation_control": control,
        "verdict": "SUPPORTED" if supported else "NOT_SUPPORTED",
        "verdict_rule": (
            "SUPPORTED requires the refitting patient-level interval to "
            "exclude zero, its estimate to be an improvement, and the "
            "permutation rank p-value to be at most 0.05."),
    }


def what_this_cannot_claim():
    return {
        "the_baseline_is_not_clinical": (
            "no per-patient covariate is published -- age, weight and "
            "creatinine clearance appear by sex group only -- so the "
            "comparator is the mean. A gain here would be a gain over no "
            "information, not over clinical knowledge plus pharmacogenetics, "
            "which is what the funding question asks."),
        "not_an_auc_or_a_clearance": (
            "the endpoint is one timed concentration divided by the daily "
            "dose. It is the reciprocal of the article's 'estimated "
            "venetoclax CL' and is not an independent endpoint."),
        "twelve_patients": (
            "a cohort of twelve, two visits each, with two not-detected "
            "concentrations. Nothing here can establish an effect size."),
        "not_a_pharmacogenetic_comparison": (
            "no genotype is published for these patients, so the contrast "
            "that matters most to the funding case -- an omic layer against "
            "pharmacogenetics -- cannot be formed."),
        "one_drug_one_population": (
            "venetoclax in elderly patients with chronic lymphocytic "
            "leukaemia, median age 71.8, on a dose that three patients had "
            "reduced for toxicity."),
        "the_in_vitro_donors_are_different_people": (
            "the same article reports fifteen hepatocyte donors and twenty "
            "microsome lots. Those are not these patients, and no analysis "
            "may join them."),
    }
