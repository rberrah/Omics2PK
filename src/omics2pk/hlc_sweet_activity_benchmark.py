"""Site-disjoint HLC RNA -> measured CYP activity benchmark for inductive SWEET.

This module is designed for the Human Liver Cohort (HLC; Synapse syn4499)
after an authenticated same-donor expression/phenotype join has already been
verified. It never performs that identity join itself.

Primary validation holds out one tissue-collection/study site at a time because
site was an important covariate in the source HLC study. No site indicator is
fed to the predictor: the held-out site is intentionally unseen.

Research evidence level E2 only: measured microsomal CYP activity, NOT clinical
plasma PK, dose selection, PopPK or a medical device.
"""
from __future__ import annotations

from math import expm1, isfinite, log1p
from statistics import mean
from typing import Any

import numpy as np

from omics2pk.epa_liver_cyp_validation import _ridge
from omics2pk.inductive_sweet import CONSTRUCTION_V1
from omics2pk.transcriptome_preprocessing import (
    EDGE_SELECTION_TRAIN_VARIANCE,
    MODE_PASSTHROUGH,
    build_preprocessed_inductive_outer_split,
    preprocess_outer_split,
)


class HLCSweetBenchmarkError(ValueError):
    pass


def _vector(values: Any, *, name: str, n: int, numeric: bool = True):
    if len(values) != n:
        raise HLCSweetBenchmarkError(f"{name} must align one-to-one with subjects")
    if numeric:
        try:
            out = np.asarray(values, dtype=np.float64)
        except (TypeError, ValueError) as exc:
            raise HLCSweetBenchmarkError(f"{name} must be numeric") from exc
        if out.shape != (n,) or not np.isfinite(out).all():
            raise HLCSweetBenchmarkError(f"{name} must be a finite vector")
        return out
    out = tuple(str(x).strip() for x in values)
    if any(not x for x in out):
        raise HLCSweetBenchmarkError(f"{name} contains blank labels")
    return out


def _pca_train_test(train_x: np.ndarray, test_x: np.ndarray, n_components: int = 4):
    if train_x.ndim != 2 or test_x.ndim != 2 or train_x.shape[1] != test_x.shape[1]:
        raise HLCSweetBenchmarkError("Invalid PCA matrices")
    n_components = min(n_components, train_x.shape[1], train_x.shape[0] - 1)
    if n_components < 1:
        raise HLCSweetBenchmarkError("Insufficient training data for PCA comparator")
    center = train_x.mean(axis=0)
    centered = train_x - center
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    components = vt[:n_components]
    return centered @ components.T, (test_x - center) @ components.T


def _predict_arm(train_x, train_y, test_x):
    pred = _ridge(
        [list(map(float, row)) for row in np.asarray(train_x)],
        [float(v) for v in train_y],
        list(map(float, np.asarray(test_x))),
    )
    value = max(0.0, expm1(min(30.0, float(pred))))
    if not isfinite(value):
        raise HLCSweetBenchmarkError("Nonfinite held-site prediction")
    return value


def evaluate_hlc_sweet_activity(
    expression: Any,
    *,
    genes,
    subject_ids,
    site_ids,
    ages,
    sexes,
    activities,
    direct_cyp_gene,
    activity_label: str,
    activity_unit: str,
    edge_top_k: int = 32,
    internal_crossfit_folds: int = 5,
) -> dict[str, Any]:
    """Compare simple RNA representations with inductive SWEET v1 by held-out site.

    Preconditions:
    - expression is source-normalized liver RNA, genes x biological donors;
    - every subject has authenticated matched microsomal CYP activity;
    - site labels are original source collection/study-site labels;
    - direct_cyp_gene is a source-annotated gene/probe identity fixed before fit.

    Gene selection for the SWEET edge panel is outcome-free and performed on
    each outer training split only. The complete expression universe is used for
    sample-weight similarity after the source-normalized passthrough transform.
    """
    try:
        x = np.asarray(expression, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise HLCSweetBenchmarkError("Expression must be numeric") from exc
    genes = tuple(str(g).strip() for g in genes)
    subjects = tuple(str(s).strip() for s in subject_ids)
    if (
        x.ndim != 2
        or x.shape != (len(genes), len(subjects))
        or not np.isfinite(x).all()
        or len(genes) < max(8, edge_top_k)
        or len(subjects) < 30
        or len(set(genes)) != len(genes)
        or len(set(subjects)) != len(subjects)
        or any(not g for g in genes)
        or any(not s for s in subjects)
    ):
        raise HLCSweetBenchmarkError("Invalid source-normalized expression matrix or identities")
    direct_genes = (
        (direct_cyp_gene,)
        if isinstance(direct_cyp_gene, str)
        else tuple(str(g).strip() for g in direct_cyp_gene)
    )
    if (
        not direct_genes
        or len(set(direct_genes)) != len(direct_genes)
        or any(not g or g not in genes for g in direct_genes)
    ):
        raise HLCSweetBenchmarkError(
            "Direct CYP gene/probe(s) must be unique and source-annotated in expression"
        )
    if not activity_label.strip() or not activity_unit.strip():
        raise HLCSweetBenchmarkError("Measured CYP activity label and physical unit are required")
    if isinstance(edge_top_k, bool) or not isinstance(edge_top_k, int) or not 8 <= edge_top_k <= 128:
        raise HLCSweetBenchmarkError("edge_top_k must be an integer in 8..128")
    if isinstance(internal_crossfit_folds, bool) or not isinstance(internal_crossfit_folds, int):
        raise HLCSweetBenchmarkError("internal_crossfit_folds must be an integer")

    n = len(subjects)
    sites = _vector(site_ids, name="site_ids", n=n, numeric=False)
    age = _vector(ages, name="ages", n=n)
    sex = _vector(sexes, name="sexes", n=n)
    activity = _vector(activities, name="activities", n=n)
    if np.any(age < 0) or np.any(age > 120):
        raise HLCSweetBenchmarkError("Age outside prespecified human range")
    if not np.all(np.isin(sex, (0.0, 1.0))):
        raise HLCSweetBenchmarkError("Sex must be encoded as 0/1 from the original source")
    if np.any(activity < 0):
        raise HLCSweetBenchmarkError("Measured CYP activity cannot be negative")

    unique_sites = tuple(sorted(set(sites)))
    if len(unique_sites) < 3:
        raise HLCSweetBenchmarkError("Primary HLC benchmark requires at least three original sites")
    site_counts = {site: sites.count(site) for site in unique_sites}
    if any(v < 5 for v in site_counts.values()):
        raise HLCSweetBenchmarkError("Each held-out site needs at least five complete donors")

    gene_pos = {g: i for i, g in enumerate(genes)}
    direct_idx = [gene_pos[g] for g in direct_genes]
    arms = (
        "age_sex",
        "age_sex_plus_direct_cyp_rna",
        "age_sex_plus_train_selected_rna",
        "age_sex_plus_train_selected_pca4",
        "age_sex_plus_inductive_sweet_v1",
    )
    all_errors = {arm: [] for arm in arms}
    site_results = {}

    for site in unique_sites:
        test_idx = [i for i, value in enumerate(sites) if value == site]
        train_idx = [i for i, value in enumerate(sites) if value != site]
        if len(train_idx) < 20:
            raise HLCSweetBenchmarkError("Too few outer-training donors after site holdout")
        train_ids = tuple(subjects[i] for i in train_idx)
        test_ids = tuple(subjects[i] for i in test_idx)

        prep = preprocess_outer_split(
            x,
            genes=genes,
            subject_ids=subjects,
            train_ids=train_ids,
            test_ids=test_ids,
            mode=MODE_PASSTHROUGH,
        )
        network = build_preprocessed_inductive_outer_split(
            x,
            genes=genes,
            subject_ids=subjects,
            train_ids=train_ids,
            test_ids=test_ids,
            tissue="liver",
            mode=MODE_PASSTHROUGH,
            folds=internal_crossfit_folds,
            construction=CONSTRUCTION_V1,
            edge_selection_mode=EDGE_SELECTION_TRAIN_VARIANCE,
            edge_top_k=edge_top_k,
        )
        selected = tuple(network["edge_genes"])
        selected_idx = [gene_pos[g] for g in selected]

        train_expression = np.asarray(prep["training"], dtype=float).T
        test_expression = np.asarray(prep["test"], dtype=float).T
        train_selected = train_expression[:, selected_idx]
        test_selected = test_expression[:, selected_idx]
        train_pca, test_pca = _pca_train_test(train_selected, test_selected, 4)

        train_network_records = {r["subject_id"]: r for r in network["training"]}
        test_network_records = {r["subject_id"]: r for r in network["test"]}
        train_sweet = np.asarray(
            [[train_network_records[s]["features"][g] for g in selected] for s in train_ids],
            dtype=float,
        )
        test_sweet = np.asarray(
            [[test_network_records[s]["features"][g] for g in selected] for s in test_ids],
            dtype=float,
        )

        train_age = age[train_idx]
        test_age = age[test_idx]
        train_sex = sex[train_idx]
        test_sex = sex[test_idx]
        train_activity = activity[train_idx]
        test_activity = activity[test_idx]
        y_train = np.log1p(train_activity)

        errors = {arm: [] for arm in arms}
        for j in range(len(test_idx)):
            base_train = np.column_stack((train_age, train_sex))
            base_test = np.asarray([test_age[j], test_sex[j]], dtype=float)
            designs = {
                "age_sex": (base_train, base_test),
                "age_sex_plus_direct_cyp_rna": (
                    np.column_stack((base_train, train_expression[:, direct_idx])),
                    np.r_[base_test, test_expression[j, direct_idx]],
                ),
                "age_sex_plus_train_selected_rna": (
                    np.column_stack((base_train, train_selected)),
                    np.r_[base_test, test_selected[j]],
                ),
                "age_sex_plus_train_selected_pca4": (
                    np.column_stack((base_train, train_pca)),
                    np.r_[base_test, test_pca[j]],
                ),
                "age_sex_plus_inductive_sweet_v1": (
                    np.column_stack((base_train, train_sweet)),
                    np.r_[base_test, test_sweet[j]],
                ),
            }
            truth = float(test_activity[j])
            for arm, (train_design, test_design) in designs.items():
                pred = _predict_arm(train_design, y_train, test_design)
                err = abs(pred - truth)
                errors[arm].append(err)
                all_errors[arm].append(err)

        test_qc = [test_network_records[s] for s in test_ids]
        site_results[site] = {
            "n_train": len(train_idx),
            "n_test": len(test_idx),
            "mae_by_model": {arm: round(mean(errors[arm]), 6) for arm in arms},
            "edge_selection_fingerprint": network["edge_selection_fingerprint"],
            "outer_reference_fingerprint": network["outer_reference_fingerprint"],
            "n_edge_genes": len(selected),
            "heldout_weight_clamp_fraction": round(
                sum(bool(r["weight_clamped_outside_reference"]) for r in test_qc) / len(test_qc),
                6,
            ),
            "heldout_abs_similarity_z_gt3_fraction": round(
                sum(abs(float(r["sample_similarity_z_vs_reference"])) > 3 for r in test_qc)
                / len(test_qc),
                6,
            ),
        }

    overall = {arm: round(mean(all_errors[arm]), 6) for arm in arms}
    return {
        "status": "HLC_SITE_DISJOINT_E2_BENCHMARK_EXECUTED",
        "evidence_level": "E2_measured_human_liver_microsomal_CYP_activity",
        "activity_label": activity_label,
        "activity_unit": activity_unit,
        "n_complete_same_donor_subjects": n,
        "n_original_sites": len(unique_sites),
        "validation": "leave_one_original_site_out",
        "internal_training_network_features": (
            f"{internal_crossfit_folds}-fold subject-disjoint crossfit"
        ),
        "preprocessing": MODE_PASSTHROUGH,
        "edge_selection": EDGE_SELECTION_TRAIN_VARIANCE,
        "edge_top_k": edge_top_k,
        "construction": CONSTRUCTION_V1,
        "n_direct_cyp_expression_features": len(direct_genes),
        "direct_cyp_expression_features": direct_genes,
        "site_used_as_predictor": False,
        "site_used_only_for_outer_holdout": True,
        "pk_pd_target_used_to_build_network": False,
        "heldout_expression_used_to_fit_preprocessing_or_gene_selection": False,
        "mae_by_model": overall,
        "delta_MAE_sweet_minus_direct_cyp_rna": round(
            overall["age_sex_plus_inductive_sweet_v1"]
            - overall["age_sex_plus_direct_cyp_rna"],
            6,
        ),
        "delta_MAE_sweet_minus_selected_rna": round(
            overall["age_sex_plus_inductive_sweet_v1"]
            - overall["age_sex_plus_train_selected_rna"],
            6,
        ),
        "delta_MAE_sweet_minus_pca4": round(
            overall["age_sex_plus_inductive_sweet_v1"]
            - overall["age_sex_plus_train_selected_pca4"],
            6,
        ),
        "site_results": site_results,
        "individual_predictions_exported": False,
        "not_clinical_PK_or_PopPK": True,
        "not_external_independent_cohort_validation": True,
        "equivalent_to_official_sweet": False,
    }
