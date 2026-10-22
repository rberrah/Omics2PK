"""Bridge the qualified private HLC ZIP parser to inductive SWEET v1.

No original donor identifiers or individual molecular/activity measurements are
returned by the public-facing benchmark result. The donor IDs are used only
in-memory for exact source joins and deterministic centre assignment.

The centre codebook is source-derived from Schadt et al. 2008 Supplementary
Table S1:
  111.... -> Vanderbilt
  222.... -> Pittsburgh
  333.... -> Merck
Any other private donor prefix fails closed; no centre is guessed.
"""
from __future__ import annotations

import hashlib
from pathlib import Path
import re
from typing import Any

import numpy as np

from omics2pk.hlc_real_cyp_validation import (
    ACTIVITIES,
    OriginalHLCQualificationError,
    load_original_hlc,
)
from omics2pk.hlc_sweet_activity_benchmark import evaluate_hlc_sweet_activity

EXPRESSION_ARCHIVE_SHA256 = (
    "35a80d532b967fdfbc73eaff571ddade17b90b30121ebe48465870245fe5e376"
)
PHENOTYPE_ARCHIVE_SHA256 = (
    "141187cd0e226ba43389bd10ab6973661ad935fa8bad51f9984f1833a2cf15a7"
)
CENTER_PREFIX = {
    "111": "Vanderbilt",
    "222": "Pittsburgh",
    "333": "Merck",
}
ACTIVITY_UNIT = "nmol/min/mg protein"


def _sha256(path) -> str:
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def source_center_from_hlc_donor_id(donor_id: str) -> str:
    """Use the public HLC source codebook; never infer an unknown prefix."""
    donor = str(donor_id).strip()
    if len(donor) != 7 or not donor.isdecimal():
        raise OriginalHLCQualificationError(
            "HLC source donor ID must be an original seven-digit identifier"
        )
    center = CENTER_PREFIX.get(donor[:3])
    if center is None:
        raise OriginalHLCQualificationError(
            "HLC donor centre prefix is not present in the public source codebook"
        )
    return center


def prepare_authorized_hlc_sweet_input(
    parsed: dict[str, Any],
    *,
    activity_key: str,
) -> dict[str, Any]:
    """Create an in-memory complete-case HLC source view for SWEET.

    This consumes the already-qualified load_original_hlc representation.
    It does not parse or re-join the original ZIPs independently.
    """
    if activity_key not in ACTIVITIES:
        raise OriginalHLCQualificationError("Unknown HLC CYP activity endpoint")
    required = {"sample_ids", "clinical", "rna", "feature_ids", "gene_rows", "outcomes"}
    if set(parsed) != required:
        raise OriginalHLCQualificationError(
            "Unexpected HLC parsed-source contract; requalify before SWEET"
        )

    sample_ids = tuple(str(x) for x in parsed["sample_ids"])
    clinical = np.asarray(parsed["clinical"], dtype=float)
    rna = np.asarray(parsed["rna"], dtype=float)
    feature_ids = tuple(str(x) for x in parsed["feature_ids"])
    y = np.asarray(parsed["outcomes"][activity_key], dtype=float)

    if (
        clinical.shape != (len(sample_ids), 2)
        or rna.shape != (len(sample_ids), len(feature_ids))
        or y.shape != (len(sample_ids),)
        or len(set(feature_ids)) != len(feature_ids)
    ):
        raise OriginalHLCQualificationError("HLC SWEET source arrays are misaligned")

    mask = np.isfinite(y) & (y >= 0) & np.isfinite(clinical).all(axis=1)
    if int(mask.sum()) < 300:
        raise OriginalHLCQualificationError(
            "Insufficient complete HLC donors for site-disjoint CYP benchmark"
        )
    eligible_rna = rna[mask]
    if not np.isfinite(eligible_rna).all():
        raise OriginalHLCQualificationError(
            "HLC expression contains missing/nonfinite values; qualify train-only "
            "source-scale imputation before running SWEET"
        )

    eligible_ids = tuple(sid for sid, keep in zip(sample_ids, mask) if keep)
    sites = tuple(source_center_from_hlc_donor_id(sid) for sid in eligible_ids)
    if set(sites) != set(CENTER_PREFIX.values()):
        raise OriginalHLCQualificationError(
            "Complete-case CYP cohort does not retain all three source centres"
        )

    target_gene = "CYP3A4" if activity_key.startswith("CYP3A4") else activity_key.split("_")[0]
    positions = parsed["gene_rows"].get(target_gene)
    if positions is None:
        raise OriginalHLCQualificationError("Target CYP probe annotation missing")
    direct_features = tuple(feature_ids[int(i)] for i in np.asarray(positions).tolist())
    if not direct_features:
        raise OriginalHLCQualificationError("No source-annotated target CYP probes")

    eligible_clinical = clinical[mask]
    return {
        "expression": eligible_rna.T,
        "genes": feature_ids,
        "subject_ids": eligible_ids,
        "site_ids": sites,
        "ages": eligible_clinical[:, 0],
        "sexes": eligible_clinical[:, 1],
        "activities": y[mask],
        "direct_cyp_features": direct_features,
        "activity_label": activity_key,
        "activity_unit": ACTIVITY_UNIT,
        "n_complete": int(mask.sum()),
        "center_counts": {
            center: sites.count(center) for center in sorted(set(sites))
        },
    }


def _load_qualified_authorized_hlc(expression_zip, phenotype_zip):
    """Hash and parse the two exact private source archives once."""
    expression_sha = _sha256(expression_zip)
    phenotype_sha = _sha256(phenotype_zip)
    if expression_sha != EXPRESSION_ARCHIVE_SHA256:
        raise OriginalHLCQualificationError(
            "Authorized HLC expression archive SHA-256 differs from qualified source"
        )
    if phenotype_sha != PHENOTYPE_ARCHIVE_SHA256:
        raise OriginalHLCQualificationError(
            "Authorized HLC phenotype archive SHA-256 differs from qualified source"
        )
    return (
        load_original_hlc(expression_zip, phenotype_zip),
        expression_sha,
        phenotype_sha,
    )


def _annotate_authorized_result(
    result: dict[str, Any],
    prepared: dict[str, Any],
    *,
    expression_sha: str,
    phenotype_sha: str,
) -> dict[str, Any]:
    result.update({
        "status": "AUTHORIZED_ORIGINAL_HLC_SITE_DISJOINT_SWEET_E2",
        "source_expression_archive_sha256": expression_sha,
        "source_phenotype_archive_sha256": phenotype_sha,
        "source_expression_probes": len(prepared["genes"]),
        "source_complete_case_donors": prepared["n_complete"],
        "source_center_complete_case_counts": prepared["center_counts"],
        "center_codebook": "Schadt_2008_Table_S1_prefix_111_222_333",
        "source_processed_expression_precedes_outer_CV": True,
        "source_article_reports_center_age_sex_control_in_analysis": True,
        "archive_exact_covariate_residualization_verified": False,
        "end_to_end_raw_expression_train_frozen": False,
        "individual_source_rows_or_ids_exported": False,
        "evidence_boundary": (
            "E2 human liver microsomal CYP activity; not clinical PK/PopPK"
        ),
    })
    return result


def _evaluate_prepared_authorized_hlc_sweet(
    prepared: dict[str, Any],
    *,
    edge_top_k: int,
    internal_crossfit_folds: int,
) -> dict[str, Any]:
    return evaluate_hlc_sweet_activity(
        prepared["expression"],
        genes=prepared["genes"],
        subject_ids=prepared["subject_ids"],
        site_ids=prepared["site_ids"],
        ages=prepared["ages"],
        sexes=prepared["sexes"],
        activities=prepared["activities"],
        direct_cyp_gene=prepared["direct_cyp_features"],
        activity_label=prepared["activity_label"],
        activity_unit=prepared["activity_unit"],
        edge_top_k=edge_top_k,
        internal_crossfit_folds=internal_crossfit_folds,
    )


def benchmark_authorized_hlc_sweet(
    expression_zip,
    phenotype_zip,
    *,
    activity_key: str = "CYP3A4M_activity",
    edge_top_k: int = 32,
    internal_crossfit_folds: int = 5,
) -> dict[str, Any]:
    """Run one prespecified site-disjoint SWEET E2 endpoint."""
    parsed, expression_sha, phenotype_sha = _load_qualified_authorized_hlc(
        expression_zip, phenotype_zip
    )
    prepared = prepare_authorized_hlc_sweet_input(parsed, activity_key=activity_key)
    result = _evaluate_prepared_authorized_hlc_sweet(
        prepared,
        edge_top_k=edge_top_k,
        internal_crossfit_folds=internal_crossfit_folds,
    )
    return _annotate_authorized_result(
        result,
        prepared,
        expression_sha=expression_sha,
        phenotype_sha=phenotype_sha,
    )


def benchmark_authorized_hlc_sweet_all(
    expression_zip,
    phenotype_zip,
    *,
    edge_top_k: int = 32,
    internal_crossfit_folds: int = 5,
) -> dict[str, Any]:
    """Run the complete prespecified 10-endpoint HLC SWEET family.

    Endpoints remain separate because their activity scales differ. The summary
    contains only descriptive win counts, never pooled/averaged MAE.
    """
    parsed, expression_sha, phenotype_sha = _load_qualified_authorized_hlc(
        expression_zip, phenotype_zip
    )
    endpoints = {}
    for activity_key in ACTIVITIES:
        prepared = prepare_authorized_hlc_sweet_input(
            parsed, activity_key=activity_key
        )
        result = _evaluate_prepared_authorized_hlc_sweet(
            prepared,
            edge_top_k=edge_top_k,
            internal_crossfit_folds=internal_crossfit_folds,
        )
        endpoints[activity_key] = _annotate_authorized_result(
            result,
            prepared,
            expression_sha=expression_sha,
            phenotype_sha=phenotype_sha,
        )

    def wins(delta_key: str) -> int:
        return sum(float(row[delta_key]) < 0 for row in endpoints.values())

    return {
        "status": "AUTHORIZED_ORIGINAL_HLC_ALL_10_SITE_DISJOINT_SWEET_E2",
        "prespecified_endpoints": list(ACTIVITIES),
        "n_prespecified_endpoints": len(ACTIVITIES),
        "endpoint_results": endpoints,
        "descriptive_endpoint_win_counts": {
            "SWEET_better_than_direct_CYP_RNA": wins(
                "delta_MAE_sweet_minus_direct_cyp_rna"
            ),
            "SWEET_better_than_train_selected_RNA": wins(
                "delta_MAE_sweet_minus_selected_rna"
            ),
            "SWEET_better_than_PCA4": wins(
                "delta_MAE_sweet_minus_pca4"
            ),
            "SWEET_better_than_topology_permuted_control": wins(
                "delta_MAE_sweet_minus_topology_permuted_control"
            ),
        },
        "MAE_not_pooled_across_endpoints": True,
        "endpoint_scales_may_differ": True,
        "endpoint_win_counts_are_descriptive_not_independent_tests": True,
        "familywise_confirmatory_inference_performed": False,
        "source_expression_archive_sha256": expression_sha,
        "source_phenotype_archive_sha256": phenotype_sha,
        "individual_source_rows_or_ids_exported": False,
        "evidence_boundary": (
            "E2 human liver microsomal CYP activity; not clinical PK/PopPK"
        ),
    }
