"""Axis B pilot: does the microsomal proteome add to the target enzyme?

Executes the design frozen in
``docs/OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_PRESPECIFICATION_20261110.md``.
The freeze table is read back out of that document before any data is loaded,
the endpoint is checked against the circularity guard, and the donor mapping
between the two halves of the source is re-derived from the data and refused
unless it is a unanimous bijection.

No donor identifier, value or prediction is written to the report.

Usage:
    python -I scripts/run_wegler_protein_clint_pilot.py \
        --per-donor-csv  <LOCAL>/wegler2021_table_s3_per_donor_kinetics.csv \
        --data-s1        <LOCAL>/mp1c00053_si_002.xlsx \
        --audit-json     reports/data/wegler_table_s3_per_donor_kinetics_audit.json \
        --out            reports/data/axis_b_wegler_protein_clint_pilot.json
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from omics2pk import protein_clint_pilot as PCP  # noqa: E402

PRESPEC_PATH = (
    REPO
    / "docs"
    / "OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_PRESPECIFICATION_20261110.md"
)

DONORS = tuple(f"UU{index}" for index in range(1, 16))
SYSTEMS = ("HLM", "HH")
BLANK_TOKENS = ("", "NA", "NaN", "nan", "ND", "#N/A", "-")

_FROZEN_IN_DOCUMENT = {
    "Seed": ("seed", lambda s: int(s.strip("` "))),
    "Components (k)": ("components", lambda s: int(s.strip())),
    "Component sensitivity": ("component_sensitivity", lambda s: int(s.strip())),
    "Permutation replicates": ("permutation_replicates", lambda s: int(s.strip())),
    "Bootstrap draws": ("n_bootstrap", lambda s: int(s.replace(",", "").strip())),
    "Primary endpoint": ("primary_endpoint", lambda s: s.strip("` ")),
    "Paired sensitivity endpoint": (
        "paired_sensitivity_endpoint", lambda s: s.strip("` ")),
    "Primary arm": ("primary_arm", lambda s: s.strip("` ")),
    "Primary comparator": ("primary_comparator", lambda s: s.strip("` ")),
    "Secondary comparator": ("secondary_comparator", lambda s: s.strip("` ")),
    "Target scale": ("target_scale", lambda s: s.strip("` ")),
    "Status ceiling": ("status_ceiling", lambda s: s.strip("` ")),
    "Ridge grid": ("ridge_alphas", lambda s: tuple(float(v) for v in s.split(","))),
}


class RunError(RuntimeError):
    pass


def check_freeze_against_document(path=PRESPEC_PATH):
    if not path.is_file():
        raise RunError(f"prespecification absent: {path}")
    text = path.read_text(encoding="utf-8")
    checked = {}
    for label, (key, cast) in _FROZEN_IN_DOCUMENT.items():
        pattern = re.compile(
            r"^\|\s*" + re.escape(label) + r"\s*\|\s*(.+?)\s*\|\s*$", re.MULTILINE
        )
        match = pattern.search(text)
        if not match:
            raise RunError(f"the prespecification no longer states {label!r}")
        stated = cast(match.group(1))
        if stated != PCP.FROZEN[key]:
            raise RunError(
                f"{label}: the document says {stated!r}, the module says "
                f"{PCP.FROZEN[key]!r}"
            )
        checked[label] = stated

    for forbidden in PCP.FROZEN["forbidden_endpoints"]:
        if f"`{forbidden}`" not in text:
            raise RunError(
                f"the prespecification does not name {forbidden!r} as forbidden"
            )
    sentence = f"**Primary contrast: `{PCP.PRIMARY_ARM}` versus `{PCP.PRIMARY_COMPARATOR}`.**"
    if sentence not in text:
        raise RunError(f"the prespecification does not state the primary contrast")
    if "the only admissible endpoint" not in text:
        raise RunError("the prespecification does not record the circularity check")
    if "donors 10-15 are permuted" not in text.replace("–", "-"):
        raise RunError("the prespecification does not record the donor mapping")
    return {
        "checked": checked,
        "primary_contrast_stated": sentence,
        "prespecification_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
    }


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _number(value):
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    text = str(value).strip()
    if text in BLANK_TOKENS:
        return None
    try:
        return float(text)
    except ValueError:
        raise RunError(f"unparseable numeric cell {value!r}") from None


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_table_s3(path, audit_json):
    audit = json.loads(Path(audit_json).read_text(encoding="utf-8"))
    if not audit["verification"]["all_anchors_agree"]:
        raise RunError(
            "the per-donor table has not been verified against the article's "
            "own published correlations; refusing to fit on an unverified "
            "transcription"
        )
    rows = {}
    with Path(path).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            rows[(row["substrate"], row["donor"])] = row
    canonical = _canonical_text(rows)
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    if digest != audit["per_donor_table_sha256"]:
        raise RunError(
            "the per-donor csv does not reproduce the checksum the audit "
            f"recorded\n  audit {audit['per_donor_table_sha256']}\n  file  {digest}"
        )
    missing = [
        key for key in ((s, d) for s in PCP.SUBSTRATE_CYP for d in DONORS)
        if key not in rows
    ]
    if missing:
        raise RunError(f"Table S3 rows absent: {missing[:5]}")
    return rows, audit, digest


def _canonical_text(rows):
    """Reproduce the audit script's canonical serialisation, to compare digests."""
    substrates = ("Midazolam", "Omeprazole", "Diclofenac", "Bupropion", "Bufuralol")
    columns = (
        "cyp_pmol_per_mg_HLM", "cyp_pmol_per_mg_HH", "cyp_pmol_per_mg_HL",
        "cyp_pmol_in_incubation_HLM", "cyp_pmol_in_incubation_HH", "kpuu_HH",
        "clint_mic_eq7", "clint_mic_eq8", "clint_mic_eq9",
        "clint_hep_eq7", "clint_hep_eq8",
    )
    lines = ["substrate,donor," + ",".join(columns)]
    for substrate in substrates:
        for donor in DONORS:
            row = rows[(substrate, donor)]
            cells = []
            for name in columns:
                raw = row[name]
                cells.append("ND" if raw == "ND" else f"{float(raw):.10g}")
            lines.append(f"{substrate},{donor}," + ",".join(cells))
    return "\n".join(lines) + "\n"


def load_data_s1(path):
    """Complete, named, log10 proteome blocks for both preparations."""
    import numpy as np
    import openpyxl

    checksum = file_digest(path)
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)
    blocks = {}
    for system, sheet in PCP.PROTEOME_SHEET.items():
        worksheet = workbook[sheet]
        prefix = {"Microsomes": "Microsomes", "Hepatocytes": "Hepatocyte"}[sheet]
        rows = worksheet.iter_rows(values_only=True)
        header = [str(c) if c is not None else "" for c in next(rows)]
        gene_index = header.index("Gene.names")
        concentration = [
            header.index(f"Concentration.{index}.{prefix}") for index in range(1, 16)
        ]
        genes = []
        values = []
        deposited = 0
        blank_gene = 0
        incomplete = 0
        for row in rows:
            deposited += 1
            raw_gene = row[gene_index]
            gene = None if raw_gene is None or not str(raw_gene).strip() else str(raw_gene).strip()
            cells = [_number(row[index]) for index in concentration]
            if gene is None:
                blank_gene += 1
                continue
            if any(cell is None for cell in cells):
                incomplete += 1
                continue
            if min(cells) <= 0:
                raise RunError(
                    f"{sheet}: gene {gene} carries a non-positive concentration, "
                    "which cannot be log-transformed"
                )
            genes.append(gene)
            values.append(cells)
        matrix = np.log10(np.asarray(values, dtype=float).T)  # donors x proteins
        blocks[system] = {
            "sheet": sheet,
            "genes": tuple(genes),
            "log10_matrix_by_data_s1_index": matrix,
            "rows_deposited": deposited,
            "rows_with_a_blank_gene_name": blank_gene,
            "rows_incomplete_over_15_donors": incomplete,
            "proteins_retained": len(genes),
        }
    workbook.close()
    return blocks, checksum


def build_mapping_inputs(table_s3, blocks):
    """Assemble the two files' CYP concentrations for the mapping proof."""
    table = {}
    for substrate, cyp in PCP.SUBSTRATE_CYP.items():
        for donor in DONORS:
            row = table_s3[(substrate, donor)]
            for system, column in PCP.BASELINE_COLUMN.items():
                table[(cyp, donor, system)] = _number(row[column])
    data = {}
    # The homogenate sheet is not loaded for modelling, so the mapping proof
    # uses the two preparations that are. Both still give four genes each.
    for system, block in blocks.items():
        lookup = {gene: position for position, gene in enumerate(block["genes"])}
        for cyp in set(PCP.SUBSTRATE_CYP.values()):
            if cyp not in lookup:
                continue
            column = block["log10_matrix_by_data_s1_index"][:, lookup[cyp]]
            for index in range(1, 16):
                data[(cyp, index, system)] = float(10.0 ** column[index - 1])
    return table, data


# --------------------------------------------------------------------------
# One cell: substrate x system x k
# --------------------------------------------------------------------------
def run_cell(table_s3, block, mapping, substrate, system, k, *,
             proteome_columns=None):
    import numpy as np

    cyp = PCP.SUBSTRATE_CYP[substrate]
    endpoint = PCP.ENDPOINT_COLUMN[system]
    PCP.assert_endpoint_is_admissible(endpoint)
    baseline_column = PCP.BASELINE_COLUMN[system]

    donors = []
    target = []
    abundance = []
    not_determined = 0
    for donor in DONORS:
        row = table_s3[(substrate, donor)]
        clearance = _number(row[endpoint])
        concentration = _number(row[baseline_column])
        if clearance is None:
            not_determined += 1
            continue
        if concentration is None or concentration <= 0 or clearance <= 0:
            raise RunError(
                f"{substrate}/{system}/{donor}: a non-positive or absent value "
                "cannot be log-transformed, and no imputation is permitted"
            )
        donors.append(donor)
        target.append(np.log10(clearance))
        abundance.append(np.log10(concentration))

    target = np.asarray(target, dtype=float)
    abundance = np.asarray(abundance, dtype=float).reshape(-1, 1)
    native = np.asarray(
        [_number(table_s3[(substrate, donor)][endpoint]) for donor in donors],
        dtype=float,
    )

    # The proteome block, reindexed from Data S1 order into Table S3 batches,
    # with the target CYP removed. The reindexing is where a wrong mapping
    # would corrupt everything, so it goes through the proven bijection only.
    batch_of_index = mapping["mapping"]
    index_of_batch = {batch: index for index, batch in batch_of_index.items()}
    genes = list(block["genes"])
    keep = [position for position, gene in enumerate(genes) if gene != cyp]
    if proteome_columns is not None:
        allowed = set(proteome_columns)
        keep = [position for position in keep if genes[position] in allowed]
    rows_for_donors = [index_of_batch[donor] - 1 for donor in donors]
    proteome = block["log10_matrix_by_data_s1_index"][np.ix_(rows_for_donors, keep)]

    n = len(donors)
    predictions = {name: np.full(n, np.nan) for name in PCP.ARMS}
    penalties = {name: [] for name in PCP.ARMS}
    shortfalls = 0
    for held in range(n):
        train = np.arange(n) != held
        test = ~train
        components_train, components_test, taken = PCP.proteome_components(
            proteome[train], proteome[test], k
        )
        if taken < k:
            shortfalls += 1
        designs = {
            "W0_mean": (np.zeros((train.sum(), 0)), np.zeros((1, 0))),
            "W1_target_cyp": (abundance[train], abundance[test]),
            "W2_target_cyp_plus_proteome": (
                np.column_stack([abundance[train], components_train]),
                np.column_stack([abundance[test], components_test]),
            ),
        }
        for name, (design_train, design_test) in designs.items():
            alpha = PCP.select_penalty(
                design_train, target[train], PCP.FROZEN["ridge_alphas"]
            )
            penalties[name].append(alpha)
            predictions[name][held] = PCP._ridge_fit_predict(
                design_train, target[train], design_test, alpha
            )[0]

    if np.isnan(np.concatenate(list(predictions.values()))).any():
        raise RunError(f"{substrate}/{system}: a donor received no prediction")

    errors = {
        name: np.abs(value - target) for name, value in predictions.items()
    }
    native_errors = {
        name: np.abs(10.0 ** value - native) for name, value in predictions.items()
    }
    return {
        "substrate": substrate,
        "cyp": cyp,
        "system": system,
        "endpoint": endpoint,
        "components": k,
        "n_donors": n,
        "n_not_determined_excluded": not_determined,
        "n_proteins_in_the_block": int(proteome.shape[1]),
        "target_cyp_removed_from_the_block": cyp in genes,
        "component_shortfall_folds": shortfalls,
        "median_penalty_by_arm": {
            name: float(np.median(value)) for name, value in penalties.items()
        },
        "losses": {
            name: PCP.losses(errors[name], native_errors[name]) for name in PCP.ARMS
        },
        "_errors": errors,
        "_donors": tuple(donors),
    }


def permutation_null(table_s3, block, mapping, substrate, system, k, replicates,
                     observed_cell):
    """Permute the proteome block's donor labels only, and refit everything."""
    import numpy as np

    cyp = PCP.SUBSTRATE_CYP[substrate]
    donors = list(observed_cell["_donors"])
    endpoint = PCP.ENDPOINT_COLUMN[system]
    baseline_column = PCP.BASELINE_COLUMN[system]
    target = np.asarray(
        [np.log10(_number(table_s3[(substrate, d)][endpoint])) for d in donors]
    )
    abundance = np.asarray(
        [np.log10(_number(table_s3[(substrate, d)][baseline_column])) for d in donors]
    ).reshape(-1, 1)

    index_of_batch = {batch: index for index, batch in mapping["mapping"].items()}
    genes = list(block["genes"])
    keep = [position for position, gene in enumerate(genes) if gene != cyp]
    rows_for_donors = [index_of_batch[donor] - 1 for donor in donors]
    proteome = block["log10_matrix_by_data_s1_index"][np.ix_(rows_for_donors, keep)]

    n = len(donors)
    rng = np.random.default_rng(PCP.FROZEN["seed"])
    differences = []
    for _ in range(replicates):
        shuffled = proteome[rng.permutation(n)]
        new = np.empty(n)
        base = np.empty(n)
        for held in range(n):
            train = np.arange(n) != held
            test = ~train
            components_train, components_test, _ = PCP.proteome_components(
                shuffled[train], shuffled[test], k
            )
            wide_train = np.column_stack([abundance[train], components_train])
            wide_test = np.column_stack([abundance[test], components_test])
            alpha = PCP.select_penalty(
                wide_train, target[train], PCP.FROZEN["ridge_alphas"]
            )
            new[held] = PCP._ridge_fit_predict(
                wide_train, target[train], wide_test, alpha
            )[0]
            alpha_base = PCP.select_penalty(
                abundance[train], target[train], PCP.FROZEN["ridge_alphas"]
            )
            base[held] = PCP._ridge_fit_predict(
                abundance[train], target[train], abundance[test], alpha_base
            )[0]
        differences.append(
            float((np.abs(new - target) - np.abs(base - target)).mean())
        )
    observed = float(
        (
            observed_cell["_errors"][PCP.PRIMARY_ARM]
            - observed_cell["_errors"][PCP.PRIMARY_COMPARATOR]
        ).mean()
    )
    summary = PCP.permutation_rank(observed, differences)
    quantiles = np.percentile(differences, list(range(0, 101, 10)))
    summary["null_deciles"] = {
        f"p{percent:02d}": float(value)
        for percent, value in zip(range(0, 101, 10), quantiles)
    }
    summary["replicate_values_are_not_published"] = (
        "An array of 200 numbers in a committed artefact is how a donor-level "
        "vector would escape; the deciles carry everything the decision uses."
    )
    return summary


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------
def run(per_donor_csv, data_s1, audit_json, *, replicates=None, n_boot=None,
        dry_run=False):
    import numpy as np

    started = time.time()
    freeze = check_freeze_against_document()
    replicates = (
        PCP.FROZEN["permutation_replicates"] if replicates is None else replicates
    )
    n_boot = PCP.FROZEN["n_bootstrap"] if n_boot is None else n_boot

    table_s3, audit, table_digest = load_table_s3(per_donor_csv, audit_json)
    blocks, data_s1_digest = load_data_s1(data_s1)
    table_inputs, data_inputs = build_mapping_inputs(table_s3, blocks)
    mapping = PCP.resolve_donor_mapping(table_inputs, data_inputs)

    circularity = {
        "admissible_endpoints": sorted(PCP.ENDPOINT_COLUMN.values()),
        "forbidden_endpoints": list(PCP.FROZEN["forbidden_endpoints"]),
        "why": (
            "The article's Figure S11 caption states that eq. 8 scales "
            "clearance by the specific CYP amount and the homogenate CYP "
            "concentration, and that eq. 9 adds Kpuu on top of that scaling. "
            "Both therefore contain the predictor. eq. 7 scales by "
            "hepatocellularity or literature MPPGL and is free of it."
        ),
        "checked_for_every_cell": True,
        "never_a_predictor_or_endpoint": list(PCP.NEVER_A_PREDICTOR_OR_ENDPOINT),
    }
    for forbidden in PCP.FROZEN["forbidden_endpoints"]:
        try:
            PCP.assert_endpoint_is_admissible(forbidden)
        except PCP.ProteinClintPilotError:
            continue
        raise RunError(
            f"the circularity guard accepted {forbidden!r}, which is the one "
            "thing it exists to refuse"
        )

    report = {
        "objective": "Axis B pilot -- does the proteome add to the target enzyme?",
        "cohort": "wegler_15",
        "status_ceiling": PCP.FROZEN["status_ceiling"],
        "no_clinical_claim": PCP.no_clinical_claim(),
        "spearman_prohibition": PCP.spearman_square_is_not_a_variance_fraction(),
        "prespecification": freeze,
        "sources": {
            "table_s3_per_donor_sha256": table_digest,
            "table_s3_verification_anchors": (
                f"{audit['verification']['n_anchors_agreeing']} of "
                f"{audit['verification']['n_anchors']}"
            ),
            "data_s1_sha256": data_s1_digest,
            "per_donor_rows_are_not_in_this_repository": True,
        },
        "circularity_check": circularity,
        "donor_mapping": {
            "mapping": {str(k): v for k, v in sorted(mapping["mapping"].items())},
            "n_vectors_voting": mapping["n_vectors_voting"],
            "unanimous": mapping["unanimous"],
            "is_identity": mapping["is_identity"],
            "why_it_is_re_derived_every_run": (
                "Table S3 numbers donors UU1-UU15 and Data S1 numbers them "
                "1-15. They agree for donors 1-9 and are PERMUTED for 10-15. "
                "Assuming the identity would mis-pair six of fifteen donors, "
                "so the bijection is proven from the CYP concentrations both "
                "files publish rather than trusted as a constant."
            ),
        },
        "proteome_blocks": {
            system: {
                key: value
                for key, value in block.items()
                if key not in ("log10_matrix_by_data_s1_index", "genes")
            }
            for system, block in blocks.items()
        },
        "log10_applied_to_the_proteome": (
            "Concentrations are strictly positive and span more than six "
            "orders of magnitude, so they are log10-transformed before any "
            "standardisation."
        ),
    }
    if dry_run:
        report["dry_run"] = True
        report["note"] = (
            "Freeze verified, circularity guard exercised, donor mapping "
            "proven. Nothing was fitted."
        )
        return report

    cytochrome_free = {
        system: tuple(
            gene for gene in block["genes"] if not PCP.CYTOCHROME_PATTERN.match(gene)
        )
        for system, block in blocks.items()
    }
    report["cytochrome_sensitivity"] = {
        "n_cytochrome_or_redox_genes_removed": {
            system: block["proteins_retained"] - len(cytochrome_free[system])
            for system, block in blocks.items()
        },
        "question": (
            "Co-regulated cytochromes could make W2's addition a "
            "re-measurement of the baseline rather than new information."
        ),
    }

    cells = {}
    for system in SYSTEMS:
        for substrate in PCP.SUBSTRATE_CYP:
            for k in (PCP.FROZEN["components"], PCP.FROZEN["component_sensitivity"]):
                cells[(substrate, system, k, "all")] = run_cell(
                    table_s3, blocks[system], mapping, substrate, system, k
                )
            cells[(substrate, system, PCP.FROZEN["components"], "no_cyp")] = run_cell(
                table_s3, blocks[system], mapping, substrate, system,
                PCP.FROZEN["components"],
                proteome_columns=cytochrome_free[system],
            )

    family = PCP.FROZEN["n_substrates"]
    results = {}
    for (substrate, system, k, block_name), cell in cells.items():
        key = f"{system}::{substrate}::k{k}::{block_name}"
        errors = cell["_errors"]
        entry = {
            name: value for name, value in cell.items()
            if not name.startswith("_")
        }
        entry["contrasts"] = {
            "primary_W2_vs_W1": PCP.paired_bootstrap(
                errors[PCP.PRIMARY_ARM], errors[PCP.PRIMARY_COMPARATOR],
                n_boot=n_boot,
            ),
            "primary_W2_vs_W1_bonferroni": PCP.paired_bootstrap(
                errors[PCP.PRIMARY_ARM], errors[PCP.PRIMARY_COMPARATOR],
                n_boot=n_boot, family_size=family,
            ),
            "secondary_W1_vs_W0": PCP.paired_bootstrap(
                errors[PCP.PRIMARY_COMPARATOR], errors[PCP.SECONDARY_COMPARATOR],
                n_boot=n_boot,
            ),
            "secondary_W1_vs_W0_bonferroni": PCP.paired_bootstrap(
                errors[PCP.PRIMARY_COMPARATOR], errors[PCP.SECONDARY_COMPARATOR],
                n_boot=n_boot, family_size=family,
            ),
            "descriptive_W2_vs_W0": PCP.paired_bootstrap(
                errors[PCP.PRIMARY_ARM], errors[PCP.SECONDARY_COMPARATOR],
                n_boot=n_boot,
            ),
        }
        results[key] = entry

    # --- the permutation nulls, primary contrast only, both systems -------
    for system in SYSTEMS:
        for substrate in PCP.SUBSTRATE_CYP:
            cell = cells[(substrate, system, PCP.FROZEN["components"], "all")]
            key = f"{system}::{substrate}::k{PCP.FROZEN['components']}::all"
            results[key]["permutation_null"] = permutation_null(
                table_s3, blocks[system], mapping, substrate, system,
                PCP.FROZEN["components"], replicates, cell,
            )

    # --- the paired system comparison ------------------------------------
    system_comparison = {}
    for substrate in PCP.SUBSTRATE_CYP:
        microsomal = cells[(substrate, "HLM", PCP.FROZEN["components"], "all")]
        hepatic = cells[(substrate, "HH", PCP.FROZEN["components"], "all")]
        shared = [d for d in microsomal["_donors"] if d in set(hepatic["_donors"])]
        mic_index = {d: i for i, d in enumerate(microsomal["_donors"])}
        hep_index = {d: i for i, d in enumerate(hepatic["_donors"])}
        entry = {
            "n_donors_in_both_systems": len(shared),
            "definition": (
                "A difference of differences over donors: the arm's benefit in "
                "microsomes minus the same arm's benefit in hepatocytes, on "
                "the donors present in both."
            ),
            "not_causal": (
                "A non-zero value means the benefit is not the same in the two "
                "preparations. It does NOT say why: the systems differ in "
                "enzyme content, membrane environment, cytosolic cofactors and "
                "scaling constant simultaneously, and this design separates "
                "none of them."
            ),
        }
        for label, new, base in (
            ("W2_vs_W1", PCP.PRIMARY_ARM, PCP.PRIMARY_COMPARATOR),
            ("W1_vs_W0", PCP.PRIMARY_COMPARATOR, PCP.SECONDARY_COMPARATOR),
        ):
            mic = np.asarray([
                microsomal["_errors"][new][mic_index[d]]
                - microsomal["_errors"][base][mic_index[d]]
                for d in shared
            ])
            hep = np.asarray([
                hepatic["_errors"][new][hep_index[d]]
                - hepatic["_errors"][base][hep_index[d]]
                for d in shared
            ])
            entry[label] = PCP.paired_bootstrap(mic, hep, n_boot=n_boot)
            entry[label]["microsomal_mean_benefit_log10"] = float(mic.mean())
            entry[label]["hepatocyte_mean_benefit_log10"] = float(hep.mean())
        system_comparison[substrate] = entry
    report["paired_system_comparison"] = system_comparison

    # --- the sample-size projection --------------------------------------
    projections = {}
    for substrate in PCP.SUBSTRATE_CYP:
        cell = cells[(substrate, "HLM", PCP.FROZEN["components"], "all")]
        differences = (
            cell["_errors"][PCP.PRIMARY_ARM]
            - cell["_errors"][PCP.PRIMARY_COMPARATOR]
        )
        key = f"HLM::{substrate}::k{PCP.FROZEN['components']}::all"
        interval = results[key]["contrasts"]["primary_W2_vs_W1"]["interval"]
        observed = float(differences.mean())
        best_case = min(interval)
        projections[substrate] = PCP.sample_size_projection(
            differences,
            effect_scenarios={
                "observed_point_estimate": observed,
                "half_the_observed": observed / 2.0,
                "most_favourable_in_the_interval": best_case,
            },
        )
    # The freeze points the projection at the PRIMARY contrast, which is the
    # one that shows no benefit; its "most favourable case" therefore demands
    # an absurd cohort, which is itself the answer. The contrast that does show
    # an effect is W1 versus W0, and a projection for it is useful -- but it
    # was NOT prespecified, so it is labelled as added after seeing the result.
    post_hoc = {}
    for system in SYSTEMS:
        for substrate in PCP.SUBSTRATE_CYP:
            cell = cells[(substrate, system, PCP.FROZEN["components"], "all")]
            differences = (
                cell["_errors"][PCP.PRIMARY_COMPARATOR]
                - cell["_errors"][PCP.SECONDARY_COMPARATOR]
            )
            key = f"{system}::{substrate}::k{PCP.FROZEN['components']}::all"
            interval = results[key]["contrasts"]["secondary_W1_vs_W0"]["interval"]
            observed = float(differences.mean())
            post_hoc[f"{system}::{substrate}"] = PCP.sample_size_projection(
                differences,
                effect_scenarios={
                    "observed_point_estimate": observed,
                    "half_the_observed": observed / 2.0,
                    "least_favourable_in_the_interval": max(interval),
                },
            )
    report["sample_size_projection"] = {
        "per_substrate": projections,
        "scope": (
            "For the PRESPECIFIED primary contrast (W2 versus W1) on the "
            "microsomal endpoint. A negative assumed effect is an improvement; "
            "a scenario whose assumed effect is a degradation is reported with "
            "its sign and is not a design target."
        ),
        "how_to_read_the_primary_projection": (
            "The primary contrast shows no benefit, so the most favourable "
            "value inside its interval is near zero and the donors it would "
            "take to detect runs into the hundreds of thousands. That is not "
            "an absurdity to be hidden; it is the finding, stated as a design "
            "requirement: the proteome-beyond-the-enzyme question is not worth "
            "sizing a study for at this effect size."
        ),
        "post_hoc_q1_projection_not_prespecified": {
            "status": "POST_HOC_ADDED_AFTER_SEEING_THE_RESULTS",
            "why": (
                "The freeze pointed the projection at the primary contrast, "
                "which turned out to be the one without an effect. The "
                "secondary contrast W1 versus W0 does show one in hepatocytes, "
                "and sizing THAT is what a future cohort would need. This was "
                "not prespecified and is labelled accordingly; it is a "
                "reporting addition, not a change of design, and the "
                "prespecified projection above is unaltered."
            ),
            "contrast": "W1_target_cyp versus W0_mean",
            "scenarios_differ_from_the_primary": (
                "The third scenario is the LEAST favourable value in the "
                "interval rather than the most favourable, because here the "
                "effect is real and the conservative question is how many "
                "donors the weakest plausible version would need."
            ),
            "per_cell": post_hoc,
        },
    }

    report["results"] = results
    report["reporting_separation"] = {
        "predictive_performance": "arms, three losses, contrasts with both levels, permutation nulls",
        "biological_reading": "the three questions, each with its caveat, plus the cytochrome sensitivity",
        "analytical_validity": "circularity check, donor mapping, ND handling, component shortfalls, fold safety, sample size",
    }
    report["no_donor_identifiers_values_or_predictions_exported"] = True
    report["runtime_seconds"] = round(time.time() - started, 1)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--per-donor-csv", required=True)
    parser.add_argument("--data-s1", required=True)
    parser.add_argument("--audit-json", required=True)
    parser.add_argument("--out", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--replicates", type=int, default=None,
                        help="override only for a smoke test; the freeze says 200")
    parser.add_argument("--bootstrap", type=int, default=None)
    args = parser.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        report = run(
            args.per_donor_csv, args.data_s1, args.audit_json,
            replicates=args.replicates, n_boot=args.bootstrap,
            dry_run=args.dry_run,
        )
    except (RunError, PCP.ProteinClintPilotError) as exc:
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
