"""ARIANA: the three-way statistical ablation, read against the mechanistic arms.

Executes the contrast family frozen in
``docs/OMICS2PK_ARIANA_STATISTICAL_VERSUS_MECHANISTIC_ABLATION_PRESPECIFICATION_20261110.md``.

This script fits **no new model and introduces no new arm**. Every arm is
imported from the two frozen runners -- the statistical arms from
``run_h1_ariana_incremental_omic_value.py`` and the mechanistic arms from
``run_h1_ariana_mechanistic_and_hybrid_arms.py`` -- so that all eleven
contrasts are computed on one fold list, which is the primary runner's own.
What this script adds is the declared family, its Bonferroni adjustment, and
the statistical-against-mechanistic contrasts no document has reported.

The endpoint is an observed pre-dose metoprolol concentration. It is not an
AUC, not a clearance and not a steady-state exposure metric. No patient
identifier, value or prediction is written to the report.

Usage:
    python -I scripts/run_ariana_statistical_versus_mechanistic_ablation.py \
        --source /path/to/Solanidine_metoprolol_ratios_PGx.xlsx \
        --out reports/data/ariana_statistical_versus_mechanistic.json
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omics2pk import ariana_mechanistic_arms as MECH  # noqa: E402

PRESPEC_PATH = (
    ROOT
    / "docs"
    / "OMICS2PK_ARIANA_STATISTICAL_VERSUS_MECHANISTIC_ABLATION_PRESPECIFICATION_20261110.md"
)
PROTOCOL_VERSION = "1.0.0"

#: The declared family, §4. Frozen here and checked against the document.
FAMILY_A = {
    "A1_genotype_vs_clinical": ("M1a_clinical_genotype", "M0_clinical"),
    "A2_comedication_vs_clinical": ("M1b_clinical_comedication", "M0_clinical"),
    "A3_metabolomics_vs_clinical": ("M2_clinical_omic", "M0_clinical"),
    "A4_both_vs_clinical": ("M3a_clinical_genotype_omic", "M0_clinical"),
    "A5_pgx_on_top_of_metabolomics": (
        "M3a_clinical_genotype_omic",
        "M2_clinical_omic",
    ),
    "A6_metabolomics_on_top_of_genotype": (
        "M3a_clinical_genotype_omic",
        "M1a_clinical_genotype",
    ),
}

#: (new arm, baseline arm, mechanistic panel for the baseline)
FAMILY_B = {
    "B1_statistical_vs_mechanistic_omic_matched": (
        "M2_clinical_omic",
        "M5c_mechanistic_omic_only",
        "oh_and_ssda",
    ),
    "B2_statistical_vs_mechanistic_with_genotype_matched": (
        "M3a_clinical_genotype_omic",
        "M5_mechanistic_pgx_omic",
        "oh_and_ssda",
    ),
    "B3_statistical_vs_mechanistic_genotype_only": (
        "M1a_clinical_genotype",
        "M4_mechanistic_pgx",
        "oh_and_ssda",
    ),
    "B4_the_two_baselines": (
        "M0_clinical",
        "M4c_backbone_only",
        "oh_and_ssda",
    ),
    "B5_the_unmatched_panel": (
        "M2_clinical_omic",
        "M5c_mechanistic_omic_only",
        "oh_only",
    ),
}

FAMILY_SIZE = len(FAMILY_A) + len(FAMILY_B)
NOTES_NAMED_CONTRAST = "A5_pgx_on_top_of_metabolomics"
MATCHED_PANEL = "oh_and_ssda"
UNMATCHED_PANEL = "oh_only"
STATISTICAL_MODEL = "ridge"


class RunError(RuntimeError):
    pass


def _load(name, filename):
    path = ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H1 = _load("run_h1_primary", "run_h1_ariana_incremental_omic_value.py")
MECHRUN = _load("run_h1_mech", "run_h1_ariana_mechanistic_and_hybrid_arms.py")


def check_imported_machinery():
    """Both frozen runners must still be the ones this family assumes."""
    expected = {
        "OUTER_FOLDS": 5,
        "OUTER_REPEATS": 2,
        "INNER_FOLDS": 3,
        "SEED": 20261107,
        "N_BOOTSTRAP": 2000,
        "RELATIVE_IMPROVEMENT_THRESHOLD": 0.10,
        "TRIM_FRACTION": 0.10,
    }
    mismatched = {
        name: (value, getattr(H1, name))
        for name, value in expected.items()
        if getattr(H1, name) != value
    }
    if mismatched:
        raise RunError(
            "the statistical runner's frozen constants have changed, so these "
            f"contrasts would not share its folds: {mismatched}"
        )
    needed_arms = {arm for pair in FAMILY_A.values() for arm in pair}
    needed_arms |= {new for new, _, _ in FAMILY_B.values()}
    missing = sorted(needed_arms - set(H1.ARMS))
    if missing:
        raise RunError(
            f"the statistical runner no longer defines the arms this family "
            f"needs: {missing}"
        )
    for panel in (MATCHED_PANEL, UNMATCHED_PANEL):
        if panel not in MECH.OMIC_PANELS:
            raise RunError(f"the mechanistic module no longer declares panel {panel!r}")
    if len(MECH.OMIC_PANELS[MATCHED_PANEL]) != len(H1.OMIC):
        raise RunError(
            f"panel {MATCHED_PANEL!r} has "
            f"{len(MECH.OMIC_PANELS[MATCHED_PANEL])} columns while the "
            f"statistical omic block has {len(H1.OMIC)}: the matched contrasts "
            "would not be matched"
        )
    if MECH.PRIMARY_OMIC_PANEL != UNMATCHED_PANEL:
        raise RunError(
            f"the mechanistic primary panel is {MECH.PRIMARY_OMIC_PANEL!r}, "
            f"not {UNMATCHED_PANEL!r}, so B5 no longer measures what it says"
        )
    MECHRUN.assert_no_prohibited_predictor()
    return {
        "statistical_constants_verified": sorted(expected),
        "mechanistic_primary_panel": MECH.PRIMARY_OMIC_PANEL,
        "matched_panel_columns": list(MECH.OMIC_PANELS[MATCHED_PANEL]),
        "statistical_omic_block": list(H1.OMIC),
        "prohibited_predictor_absent_everywhere": H1.COL_PROHIBITED,
    }


_FROZEN_IN_DOCUMENT = {
    "Family size": (lambda: FAMILY_SIZE, lambda s: int(s.strip())),
    "Fold seed": (lambda: H1.SEED, lambda s: int(s.strip("` "))),
    "Bootstrap draws": (lambda: H1.N_BOOTSTRAP,
                        lambda s: int(s.replace(",", "").strip())),
    "Relative MAE threshold": (lambda: H1.RELATIVE_IMPROVEMENT_THRESHOLD,
                               lambda s: float(s.strip())),
    "Matched mechanistic panel": (lambda: MATCHED_PANEL, lambda s: s.strip("` ")),
    "Unmatched mechanistic panel": (lambda: UNMATCHED_PANEL, lambda s: s.strip("` ")),
    "Note's named contrast": (lambda: NOTES_NAMED_CONTRAST, lambda s: s.strip("` ")),
    "Adjustment": (lambda: "Bonferroni", lambda s: s.strip()),
    "Status ceiling": (lambda: "exploratory", lambda s: s.strip("` ")),
}


def check_freeze_against_document(path=PRESPEC_PATH):
    if not path.is_file():
        raise RunError(f"prespecification absent: {path}")
    text = path.read_text(encoding="utf-8")
    checked = {}
    for label, (getter, cast) in _FROZEN_IN_DOCUMENT.items():
        pattern = re.compile(
            r"^\|\s*" + re.escape(label) + r"\s*\|\s*(.+?)\s*\|\s*$", re.MULTILINE
        )
        match = pattern.search(text)
        if not match:
            raise RunError(f"the prespecification no longer states {label!r}")
        stated = cast(match.group(1))
        expected = getter()
        if stated != expected:
            raise RunError(
                f"{label}: the document says {stated!r}, the runner says {expected!r}"
            )
        checked[label] = stated
    outer = re.search(
        r"^\|\s*Outer folds × repeats × inner\s*\|\s*(\d+) × (\d+) × (\d+)\s*\|\s*$",
        text,
        re.MULTILINE,
    )
    if not outer:
        raise RunError("the prespecification no longer states the fold structure")
    stated_folds = tuple(int(v) for v in outer.groups())
    if stated_folds != (H1.OUTER_FOLDS, H1.OUTER_REPEATS, H1.INNER_FOLDS):
        raise RunError(f"fold structure disagrees: document says {stated_folds}")
    checked["Outer folds × repeats × inner"] = list(stated_folds)

    for key in list(FAMILY_A) + list(FAMILY_B):
        if f"`{key}`" not in text and f"**`{key}`**" not in text:
            raise RunError(f"the prespecification does not name the contrast {key!r}")
    if "Family A below is not a blind prespecification" not in text:
        raise RunError(
            "the prespecification must state that Family A is not a blind test"
        )
    if "It does not isolate the representation" not in text:
        raise RunError(
            "the prespecification must state what Family B cannot measure"
        )
    return {
        "checked": checked,
        "prespecification_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }


# --------------------------------------------------------------------------
def statistical_errors(frame, observed, observed_log, folds):
    """Out-of-fold errors of the frozen statistical arms, ridge, same folds."""
    factory, grid = H1.model_factories()[STATISTICAL_MODEL]
    errors = {}
    losses = {}
    for arm, features in H1.ARMS.items():
        design = frame[features].to_numpy(dtype=float)
        predictions = H1.out_of_fold_predictions(
            design, observed_log, folds, factory, grid
        )
        errors[arm] = MECHRUN.arm_errors(observed, observed_log, predictions)
        losses[arm] = H1.loss_summary(observed, observed_log, predictions)
    return errors, losses


def mechanistic_errors(frame, observed, observed_log, folds, panel):
    """Out-of-fold errors of the frozen mechanistic arms at one declared panel.

    The call is the frozen runner's own, with the frozen transferred constants
    and the backbone that runner names primary. Only the omic panel varies, and
    both panels used here are already declared in the mechanistic module.
    """
    arrays = MECHRUN.build_arrays(frame, MECH.OMIC_PANELS[panel])
    backbone = MECH.BACKBONES[MECH.PRIMARY_BACKBONE]
    predictions, _ = MECHRUN.mechanistic_out_of_fold(
        arrays,
        folds,
        backbone,
        gamma=MECH.OMIC_GAMMA,
        blend_weight=MECH.BLEND_WEIGHT,
        fm=MECH.FM_CYP2D6,
        delta_as=MECH.DELTA_AS_PER_INHIBITOR,
        dose_equivalence=1.0,
    )
    errors = {
        arm: MECHRUN.arm_errors(observed, observed_log, matrix)
        for arm, matrix in predictions.items()
    }
    losses = {
        arm: H1.loss_summary(observed, observed_log, matrix)
        for arm, matrix in predictions.items()
    }
    return errors, losses


def one_contrast(new_errors, baseline_errors, observed, n_boot, family_size):
    contrast = MECHRUN.full_contrast(new_errors, baseline_errors, observed, n_boot)
    contrast["multiplicity_adjusted"] = MECHRUN.bonferroni_interval(
        new_errors["absolute"], baseline_errors["absolute"], n_boot, family_size
    )
    unadjusted = contrast["primary_loss"]
    adjusted = contrast["multiplicity_adjusted"]
    relative = unadjusted["relative_change_vs_baseline"]
    reaches = (
        relative is not None
        and -relative >= H1.RELATIVE_IMPROVEMENT_THRESHOLD
        and unadjusted["direction"] == "improvement"
    )
    contrast["decision"] = {
        "carries_a_claim": bool(reaches and adjusted["interval_excludes_zero"]),
        "reaches_the_relative_threshold": bool(reaches),
        "adjusted_interval_excludes_zero": bool(adjusted["interval_excludes_zero"]),
        "unadjusted_interval_excludes_zero": bool(
            unadjusted["interval_excludes_zero"]
        ),
        "inconclusive_within_the_family": bool(
            unadjusted["interval_excludes_zero"]
            and not adjusted["interval_excludes_zero"]
        ),
        "rule": (
            "A contrast carries a claim only if the ADJUSTED interval excludes "
            "zero and the relative improvement reaches "
            f"{H1.RELATIVE_IMPROVEMENT_THRESHOLD:.0%}. An unadjusted interval "
            "excluding zero while the adjusted one does not is inconclusive "
            "within the family, never a positive finding. An interval "
            "including zero means no effect was DETECTED, never that there is "
            "none."
        ),
    }
    return contrast


def run(source, n_boot=None):
    import numpy as np

    freeze = check_freeze_against_document()
    machinery = check_imported_machinery()
    n_boot = H1.N_BOOTSTRAP if n_boot is None else n_boot

    frame, checksum = H1.load_ariana(source)
    observed = frame[H1.COL_TROUGH].to_numpy(dtype=float)
    if np.any(observed <= 0):
        raise RunError("non-positive concentrations cannot be logged")
    observed_log = np.log(observed)
    folds = H1.make_folds(len(frame), H1.OUTER_FOLDS, H1.OUTER_REPEATS, H1.SEED)

    stat_errors, stat_losses = statistical_errors(
        frame, observed, observed_log, folds
    )
    mech = {}
    mech_losses = {}
    for panel in (MATCHED_PANEL, UNMATCHED_PANEL):
        mech[panel], mech_losses[panel] = mechanistic_errors(
            frame, observed, observed_log, folds, panel
        )

    report = {
        "objective": "H1 secondary: the three-way statistical ablation and the "
                     "statistical-against-mechanistic comparison",
        "protocol_version": PROTOCOL_VERSION,
        "cohort": "ariana_metoprolol_47",
        "n_patients": int(len(frame)),
        "source_sha256": checksum,
        "target": H1.COL_TROUGH,
        "target_is": (
            "observed pre-dose metoprolol concentration, ng/mL. NOT an AUC, "
            "NOT a clearance, NOT a steady-state exposure metric."
        ),
        "prohibited_predictor_absent": H1.COL_PROHIBITED,
        "map_be": "excluded by construction; one observation per patient",
        "status_ceiling": (
            "exploratory. This cohort has been analysed repeatedly by this "
            "project, so nothing here is an independent validation, a new "
            "cohort, or evidence transferable to another drug."
        ),
        "no_new_model_was_fitted": (
            "Every arm is imported from the two frozen runners. This protocol "
            "adds a declared family, its Bonferroni adjustment and the "
            "statistical-against-mechanistic contrasts, nothing else."
        ),
        "family_a_is_not_a_blind_prespecification": (
            "Family A's contrasts were already computed and visible in the "
            "committed reports/data/h1_ariana_trough.json before the protocol "
            "was written. Family A is an adjusted reading of an existing "
            "computation, never a newly prespecified test."
        ),
        "family_b_cannot_isolate_the_representation": (
            "A mechanistic arm fits one free parameter on these 47 patients "
            "and borrows the rest from the literature; a statistical arm fits "
            "a ridge penalty and six to eight coefficients on the same 47. "
            "Family B compares transferred structure plus one free parameter "
            "against estimation from 47 patients, not representation alone."
        ),
        "prespecification": freeze,
        "imported_machinery": machinery,
        "validation": {
            "outer_folds": H1.OUTER_FOLDS,
            "outer_repeats": H1.OUTER_REPEATS,
            "inner_folds": H1.INNER_FOLDS,
            "seed": H1.SEED,
            "fold_fingerprint_sha256": MECHRUN.fold_fingerprint(folds),
            "folds_identical_across_both_families": True,
            "how_identity_is_guaranteed": (
                "the fold list is the primary runner's own make_folds, "
                "imported rather than reimplemented"
            ),
            "bootstrap": {
                "n": int(n_boot),
                "unit": "patients; a fold repeat is never a second patient",
            },
            "adjustment": f"Bonferroni across a declared family of {FAMILY_SIZE}",
            "why_bonferroni_not_holm": (
                "Bonferroni is uniformly more conservative, so an interval "
                "excluding zero after Bonferroni would also exclude it after "
                "Holm. It also reuses tested code instead of adding an "
                "untested adjustment path."
            ),
        },
        "free_parameters_fitted_on_this_cohort": dict(MECH.FREE_PARAMETERS),
        "statistical_arms": {
            name: list(features) for name, features in H1.ARMS.items()
        },
        "losses": {
            "statistical_ridge": stat_losses,
            "mechanistic_matched_panel": mech_losses[MATCHED_PANEL],
            "mechanistic_unmatched_panel": mech_losses[UNMATCHED_PANEL],
        },
        "family_size": FAMILY_SIZE,
        "contrasts": {},
    }

    for key, (new, base) in FAMILY_A.items():
        report["contrasts"][key] = one_contrast(
            stat_errors[new], stat_errors[base], observed, n_boot, FAMILY_SIZE
        )
        report["contrasts"][key]["definition"] = (
            f"{new} versus {base}, both statistical, ridge, mean absolute error"
        )
        report["contrasts"][key]["family"] = "A"

    for key, (new, base, panel) in FAMILY_B.items():
        report["contrasts"][key] = one_contrast(
            stat_errors[new], mech[panel][base], observed, n_boot, FAMILY_SIZE
        )
        report["contrasts"][key]["definition"] = (
            f"{new} (statistical, ridge) versus {base} (mechanistic, panel "
            f"{panel}), mean absolute error, identical folds"
        )
        report["contrasts"][key]["family"] = "B"
        report["contrasts"][key]["mechanistic_panel"] = panel
        report["contrasts"][key]["panels_are_matched"] = bool(
            panel == MATCHED_PANEL
        )

    named = report["contrasts"][NOTES_NAMED_CONTRAST]
    difference = named["primary_loss"]["mean_difference"]
    report["the_question_the_note_named"] = {
        "contrast": NOTES_NAMED_CONTRAST,
        "question": (
            "Does pharmacogenetics add information the metabolomics does not "
            "already carry?"
        ),
        "point_estimate_mean_absolute_error_difference": difference,
        "direction": "improvement" if difference < 0 else "degradation",
        "answer": (
            "no detectable complementary information; the point estimate is a "
            "degradation"
            if difference > 0
            else "an improvement was detected"
            if named["decision"]["carries_a_claim"]
            else "an improvement in point estimate only, not carrying a claim"
        ),
        "adjusted_interval": named["multiplicity_adjusted"]["interval"],
        "adjusted_interval_excludes_zero": named["multiplicity_adjusted"][
            "interval_excludes_zero"
        ],
    }
    report["claims_carried"] = sorted(
        key
        for key, value in report["contrasts"].items()
        if value["decision"]["carries_a_claim"]
    )
    # A contrast whose adjusted interval excludes zero in the WRONG direction
    # carries no claim under §6, but it is a detected difference and must not
    # be filed alongside the contrasts where nothing was detected.
    report["differences_detected_in_favour_of_the_baseline"] = sorted(
        key
        for key, value in report["contrasts"].items()
        if value["multiplicity_adjusted"]["interval_excludes_zero"]
        and value["primary_loss"]["direction"] == "degradation"
    )
    report["nothing_detected"] = sorted(
        key
        for key, value in report["contrasts"].items()
        if not value["multiplicity_adjusted"]["interval_excludes_zero"]
    )
    report["how_to_read_the_three_lists"] = (
        "claims_carried: adjusted interval excludes zero AND the improvement "
        "reaches the prespecified size. differences_detected_in_favour_of_the_"
        "baseline: the adjusted interval excludes zero but the new arm is "
        "WORSE, which is a finding about the baseline and not an absence of "
        "effect. nothing_detected: the adjusted interval includes zero, which "
        "means no difference was detected and never that there is none."
    )
    report["inconclusive_within_the_family"] = sorted(
        key
        for key, value in report["contrasts"].items()
        if value["decision"]["inconclusive_within_the_family"]
    )
    report["no_clinical_pk_proof"] = (
        "No result here is a clinical-pharmacokinetic proof. It is a "
        "comparison of representations on one observed trough concentration in "
        "47 patients from one study."
    )
    report["no_person_identifiers_values_or_predictions_exported"] = True
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, help="ARIANA workbook")
    parser.add_argument("--out", default=None)
    parser.add_argument("--bootstrap", type=int, default=None)
    args = parser.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        report = run(args.source, n_boot=args.bootstrap)
    except (RunError, H1.AnalysisError) as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    text = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
        print(f"written: {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
