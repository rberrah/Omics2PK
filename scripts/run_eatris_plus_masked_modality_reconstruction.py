"""Axis A: reconstruct one masked molecular layer in EATRIS-Plus v4.

Executes the design frozen in
``docs/OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_PRESPECIFICATION_20261110.md``.
The freeze table is read back out of that document before any data is loaded,
and the run refuses to start if the document and the module disagree.

This cohort has no drug, no dose, no concentration, no enzyme activity and no
clinical endpoint. Nothing this script prints is evidence about drug
metabolism or pharmacokinetics. No individual identifier, no individual value
and no individual prediction is written to the report.

Usage:
    python -I scripts/run_eatris_plus_masked_modality_reconstruction.py \
        --source /path/to/EATRIS_PLUS_mudata.h5mu --out reports/data/axis_a.json

    python -I scripts/run_eatris_plus_masked_modality_reconstruction.py \
        --source ... --dry-run
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from omics2pk import multiomic_layer_complementarity as MLC  # noqa: E402

PRESPEC_PATH = (
    REPO
    / "docs"
    / "OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_PRESPECIFICATION_20261110.md"
)
EXPECTED_SOURCE_SHA256 = (
    "a764a6cd7246edbd5c3a4c4ce166cd2517dc93eaa11834f2d1e8530245581c7c"
)

#: Freeze-table labels in §13 of the document, and how to parse each.
_FROZEN_IN_DOCUMENT = {
    "Fold seed": ("seed", lambda s: int(s.strip("` "))),
    "Bootstrap draws": ("n_bootstrap", lambda s: int(s.replace(",", "").strip())),
    "Components per layer (k)": ("components_per_layer", lambda s: int(s.strip())),
    "Permutation replicates": ("permutation_replicates", lambda s: int(s.strip())),
    "Parsimony tolerance (median ΔR²)": (
        "parsimony_tolerance",
        lambda s: float(s.strip()),
    ),
    "Primary target": ("primary_target", lambda s: s.strip("` ")),
    "Secondary target": ("secondary_target", lambda s: s.strip("` ")),
    "Primary arm": ("primary_arm", lambda s: s.strip("` ")),
    "Primary comparator": ("primary_comparator", lambda s: s.strip("` ")),
    "Complementarity comparator": (
        "complementarity_comparator",
        lambda s: s.strip("` "),
    ),
    "Status ceiling": ("status_ceiling", lambda s: s.strip("` ")),
    "Standardised-value clip": (
        "standardised_clip",
        lambda s: float(s.replace("±", "").replace("training-fold SD", "").strip()),
    ),
    "Ridge grid": (
        "ridge_alphas",
        lambda s: tuple(float(v) for v in s.split(",")),
    ),
    "Component sensitivities": (
        "component_sensitivities",
        lambda s: tuple(int(v) for v in s.split(",")),
    ),
}


class RunError(RuntimeError):
    pass


def check_freeze_against_document(path=PRESPEC_PATH):
    """A freeze in code only moves silently; a freeze in prose only is unchecked."""
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
        if stated != MLC.FROZEN[key]:
            raise RunError(
                f"{label}: the document says {stated!r}, the module says "
                f"{MLC.FROZEN[key]!r}"
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
    module_folds = (
        MLC.FROZEN["outer_folds"],
        MLC.FROZEN["outer_repeats"],
        MLC.FROZEN["inner_folds"],
    )
    if stated_folds != module_folds:
        raise RunError(
            f"fold structure: the document says {stated_folds}, the module "
            f"says {module_folds}"
        )
    checked["Outer folds × repeats × inner"] = list(stated_folds)

    sentence = (
        f"`{MLC.PRIMARY_ARM}` versus `{MLC.PRIMARY_COMPARATOR}`"
    )
    if sentence not in text:
        raise RunError(f"the prespecification does not state the primary contrast ({sentence})")
    for name, spec in MLC.CONFIGURATIONS.items():
        if spec["individual_set_sha256_prefix"] not in text:
            raise RunError(
                f"the prespecification does not record the individual set of {name}"
            )
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


def set_checksum(identifiers):
    payload = "\n".join(sorted(str(x) for x in identifiers)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def _decode(values):
    return [v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values]


def _index(group):
    if "_index" in group:
        return _decode(group["_index"][...])
    attr = group.attrs.get("_index")
    key = attr.decode("utf-8") if isinstance(attr, bytes) else str(attr)
    return _decode(group[key][...])


def _obs_column(obs_group, name):
    node = obs_group[name]
    if hasattr(node, "keys") and "codes" in node:
        codes = node["codes"][...]
        categories = _decode(node["categories"][...])
        return [categories[c] if c >= 0 else None for c in codes]
    if node.dtype.kind in "SOU":
        return _decode(node[...])
    return list(node[...])


def load_container(source):
    import h5py
    import numpy as np

    checksum = file_digest(source)
    if checksum != EXPECTED_SOURCE_SHA256:
        raise RunError(
            f"source checksum mismatch\n  expected {EXPECTED_SOURCE_SHA256}\n"
            f"  actual   {checksum}"
        )
    with h5py.File(source, "r") as handle:
        layer_names = sorted(handle["mod"].keys())
        MLC.assert_exact_partition(layer_names)
        values = {}
        identifiers = {}
        logged = {}
        for name in layer_names:
            group = handle["mod"][name]
            identifiers[name] = _index(group["obs"])
            block = np.asarray(group["X"][...], dtype=float)
            # §8.2.1 rule 1 is a property of the deposited layer, so it is
            # evaluated here, once, on every row -- never per fold.
            if MLC.log_rule_applies(block):
                logged[name] = True
                block = np.log10(block)
            else:
                logged[name] = False
            values[name] = block
        person_ids = _index(handle["obs"])
        covariates = {}
        for name in MLC.COVARIATES_K0 + MLC.COVARIATES_H:
            if name not in handle["obs"]:
                raise RunError(f"declared covariate {name!r} is absent from the person table")
            covariates[name] = _obs_column(handle["obs"], name)
    return {
        "source_sha256": checksum,
        "layer_names": layer_names,
        "values": values,
        "identifiers": identifiers,
        "log10_applied": logged,
        "person_ids": person_ids,
        "covariates": covariates,
    }


def _numeric_covariate(name, column, order):
    import numpy as np

    out = np.full(len(order), np.nan)
    for position, index in enumerate(order):
        raw = column[index]
        if raw is None:
            continue
        if name == "Sex":
            text = str(raw).strip().lower()
            if text in ("f", "female", "0"):
                out[position] = 0.0
            elif text in ("m", "male", "1"):
                out[position] = 1.0
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value):
            out[position] = value
    return out


def build_configuration(container, name):
    import numpy as np

    spec = MLC.CONFIGURATIONS[name]
    kept = [x for x in container["layer_names"] if x not in spec["drop_layers"]]
    sets = [set(container["identifiers"][layer]) for layer in kept]
    complete = set.intersection(*sets)
    individuals = sorted(complete)
    if len(individuals) != spec["n"]:
        raise RunError(
            f"{name}: complete case is {len(individuals)}, the freeze says {spec['n']}"
        )
    digest = set_checksum(individuals)
    if not digest.startswith(spec["individual_set_sha256_prefix"]):
        raise RunError(
            f"{name}: individual-set checksum {digest[:24]} does not match the "
            f"frozen {spec['individual_set_sha256_prefix']}"
        )

    layers = {}
    screened = {}
    for layer in kept:
        row_of = {x: i for i, x in enumerate(container["identifiers"][layer])}
        rows = [row_of[x] for x in individuals]
        block = container["values"][layer][rows, :]
        mask = MLC.complete_columns(block)
        screened[layer] = {
            "columns": int(block.shape[1]),
            "dropped_for_missingness": int((~mask).sum()),
            "retained": int(mask.sum()),
            "log10_applied": bool(container["log10_applied"][layer]),
        }
        layers[layer] = np.ascontiguousarray(block[:, mask])

    person_of = {x: i for i, x in enumerate(container["person_ids"])}
    order = [person_of[x] for x in individuals]
    k0 = np.column_stack([
        _numeric_covariate(cov, container["covariates"][cov], order)
        for cov in MLC.COVARIATES_K0
    ])
    if not np.isfinite(k0).all():
        missing = {
            cov: int(np.isnan(_numeric_covariate(
                cov, container["covariates"][cov], order)).sum())
            for cov in MLC.COVARIATES_K0
        }
        raise RunError(
            f"{name}: the baseline block K0 must be complete, missing {missing}"
        )
    haematology = np.column_stack([
        _numeric_covariate(cov, container["covariates"][cov], order)
        for cov in MLC.COVARIATES_H
    ])
    return {
        "name": name,
        "role": spec["role"],
        "layers_kept": kept,
        "n": len(individuals),
        "individual_set_sha256": digest,
        "layer_blocks": layers,
        "screened": screened,
        "K0": k0,
        "H": haematology,
        "haematology_missing_cells": int(np.isnan(haematology).sum()),
        "families_present": {
            family: list(layers_)
            for family, layers_ in MLC.families_present(kept).items()
        },
    }


# --------------------------------------------------------------------------
# Fitting
# --------------------------------------------------------------------------
def make_folds(n):
    from sklearn.model_selection import RepeatedKFold

    splitter = RepeatedKFold(
        n_splits=MLC.FROZEN["outer_folds"],
        n_repeats=MLC.FROZEN["outer_repeats"],
        random_state=MLC.FROZEN["seed"],
    )
    return [(train, test) for train, test in splitter.split(range(n))]


def component_cache(config, folds, k_max):
    """Training-fold components for every (layer, fold), extracted once.

    The first j columns of a k-component extraction are exactly the j-component
    extraction, because the directions are ordered, so the sensitivities at
    k = 2 and k = 10 are slices of one computation rather than three runs.
    """
    cache = {}
    shortfalls = {}
    for fold_index, (train, test) in enumerate(folds):
        for layer, block in config["layer_blocks"].items():
            scores_train, scores_test, taken = MLC.layer_components(
                block[train], block[test], k_max
            )
            cache[(layer, fold_index)] = (scores_train, scores_test)
            if taken < k_max:
                entry = shortfalls.setdefault(
                    layer,
                    {"folds_affected": 0, "min_available": taken,
                     "max_available": taken, "requested": int(k_max)},
                )
                entry["folds_affected"] += 1
                entry["min_available"] = min(entry["min_available"], taken)
                entry["max_available"] = max(entry["max_available"], taken)
    return cache, shortfalls


def _impute_haematology(train_block, test_block):
    import numpy as np

    median = np.nanmedian(train_block, axis=0)
    median = np.where(np.isfinite(median), median, 0.0)
    filled_train = np.where(np.isnan(train_block), median, train_block)
    filled_test = np.where(np.isnan(test_block), median, test_block)
    return filled_train, filled_test


def assemble(config, cache, fold_index, train, test, layers, k, with_haematology):
    import numpy as np

    blocks_train = [config["K0"][train]]
    blocks_test = [config["K0"][test]]
    if with_haematology:
        filled_train, filled_test = _impute_haematology(
            config["H"][train], config["H"][test]
        )
        blocks_train.append(filled_train)
        blocks_test.append(filled_test)
    for layer in layers:
        scores_train, scores_test = cache[(layer, fold_index)]
        take = min(k, scores_train.shape[1])
        blocks_train.append(scores_train[:, :take])
        blocks_test.append(scores_test[:, :take])
    return np.column_stack(blocks_train), np.column_stack(blocks_test)


def total_sums_of_squares(targets, folds, n_repeats):
    """Per-individual squared deviation from that individual's training mean."""
    import numpy as np

    n_features = targets.shape[1]
    n = targets.shape[0]
    total = np.zeros((n_features, n))
    for train, test in folds:
        train_mean = targets[train].mean(axis=0)
        total[:, test] += ((targets[test] - train_mean) ** 2).T
    return total / n_repeats


def arm_errors(config, cache, folds, targets, layers, k, with_haematology):
    import numpy as np

    n_features = targets.shape[1]
    n = targets.shape[0]
    error = np.zeros((n_features, n))
    penalties = []
    width = None
    for fold_index, (train, test) in enumerate(folds):
        repeat = fold_index // MLC.FROZEN["outer_folds"]
        design_train, design_test = assemble(
            config, cache, fold_index, train, test, layers, k, with_haematology
        )
        width = design_train.shape[1]
        predictions, chosen, _ = MLC.ridge_per_target(
            design_train,
            targets[train],
            design_test,
            seed=MLC.FROZEN["seed"] + repeat,
        )
        error[:, test] += ((predictions - targets[test]) ** 2).T
        penalties.append(float(np.median(chosen)))
    error /= MLC.FROZEN["outer_repeats"]
    # Written as a mapping, not a list: the committed-artefact guard refuses
    # any numeric array of ten or more entries, because such an array in a
    # published report is how an individual-level vector would escape.
    return error, {
        "design_columns": int(width),
        "median_penalty_by_fold": {
            f"fold_{index + 1:02d}": value for index, value in enumerate(penalties)
        },
    }


# --------------------------------------------------------------------------
# One configuration, one target
# --------------------------------------------------------------------------
def run_target(config, cache, folds, target, *, k, run_permutations,
               run_parsimony, extra_arms=None, n_boot=None):
    import numpy as np

    n_boot = MLC.FROZEN["n_bootstrap"] if n_boot is None else n_boot
    targets = config["layer_blocks"][target]
    total = total_sums_of_squares(targets, folds, MLC.FROZEN["outer_repeats"])

    arms = MLC.arms_for(target, config["layers_kept"])
    subsets = MLC.family_subsets(target, config["layers_kept"])
    if extra_arms:
        arms.update(extra_arms)

    designs = {}
    for name, layers in arms.items():
        designs[name] = (layers, name == "A4_covariates_haematology")
    for name, (_, layers) in subsets.items():
        designs.setdefault(name, (layers, False))

    errors = {}
    shapes = {}
    for name, (layers, with_h) in designs.items():
        errors[name], shapes[name] = arm_errors(
            config, cache, folds, targets, layers, k, with_h
        )

    r2 = {
        name: MLC.r2_per_feature(block, total) for name, block in errors.items()
    }
    arm_table = {}
    for name in sorted(errors):
        values = r2[name]
        layers = designs[name][0]
        arm_table[name] = {
            "layers": list(layers),
            "n_layers": len(layers),
            "n_source_features": int(sum(
                config["layer_blocks"][layer].shape[1] for layer in layers
            )),
            "design_columns": shapes[name]["design_columns"],
            "median_r2": float(np.nanmedian(values)),
            "mean_r2": float(np.nanmean(values)),
            "fraction_r2_above_zero": float(np.nanmean(values > 0.0)),
            "r2_quartiles": [float(v) for v in np.nanpercentile(values, [25, 50, 75])],
            "median_penalty_by_fold": shapes[name]["median_penalty_by_fold"],
        }

    # Every reported contrast is named here as a (new, baseline) pair, so that
    # the permutation control below can be computed for ALL of them uniformly
    # rather than for whichever ones looked interesting afterwards.
    pairs = {
        "primary_A1_vs_A0": (MLC.PRIMARY_ARM, MLC.PRIMARY_COMPARATOR),
        "complementarity_A1_vs_A2": (MLC.PRIMARY_ARM, MLC.COMPLEMENTARITY_COMPARATOR),
        "sibling_A2_vs_A0": (MLC.COMPLEMENTARITY_COMPARATOR, MLC.PRIMARY_COMPARATOR),
        "all_but_target_A3_vs_A1": ("A3_all_but_target", MLC.PRIMARY_ARM),
        "haematology_A4_vs_A0": (
            "A4_covariates_haematology",
            MLC.PRIMARY_COMPARATOR,
        ),
    }
    target_family = MLC.family_of(target)
    for family in sorted(config["families_present"]):
        if f"B_{family}" in designs:
            pairs[f"single_family_B_{family}_vs_A0"] = (
                f"B_{family}",
                MLC.PRIMARY_COMPARATOR,
            )
            pairs[f"leave_out_L_{family}_vs_A1"] = (f"L_{family}", MLC.PRIMARY_ARM)
    if extra_arms:
        for name in extra_arms:
            pairs[f"{MLC.PRIMARY_ARM}_vs_{name}"] = (MLC.PRIMARY_ARM, name)
    # The nine-layer subset inside the primary set: the eleven-layer arm
    # against the same arm without the two RNA families, on the SAME people and
    # the SAME folds, which is the only clean way to price the transcriptome.
    rna_free = tuple(
        sorted(
            family
            for family in config["families_present"]
            if family not in (target_family, "transcriptome", "small_rna")
        )
    )
    rna_free_arm = "S_" + ("+".join(rna_free) if rna_free else "none")
    # In a configuration that holds no RNA layer at all, this subset IS the
    # primary arm, and the contrast would be an arm against itself: a
    # tautological zero dressed as a comparison. Compare the layer content,
    # not the arm names, so the degenerate case is actually caught.
    degenerate = rna_free_arm in designs and set(designs[rna_free_arm][0]) == set(
        designs[MLC.PRIMARY_ARM][0]
    )
    if rna_free_arm in designs and not degenerate:
        pairs["two_rna_families_A1_vs_" + rna_free_arm] = (
            MLC.PRIMARY_ARM,
            rna_free_arm,
        )

    def contrast(new, base):
        if new not in errors or base not in errors:
            return {"status": "NOT_EVALUABLE",
                    "reason": f"missing arm(s): {[x for x in (new, base) if x not in errors]}"}
        out = MLC.paired_bootstrap_median_delta(
            errors[new], errors[base], total, n_boot=n_boot
        )
        out["new_arm"] = new
        out["baseline_arm"] = base
        out["new_median_r2"] = arm_table[new]["median_r2"]
        out["baseline_median_r2"] = arm_table[base]["median_r2"]
        return out

    report = {
        "target_layer": target,
        "target_features_deposited": int(
            config["screened"][target]["columns"]
        ),
        "target_features_scored": int(targets.shape[1]),
        "target_same_family_sibling": list(
            x for x in MLC.siblings_of(target) if x in config["layers_kept"]
        ),
        "biochemically_adjacent_families": list(
            MLC.BIOCHEMICALLY_ADJACENT.get(MLC.family_of(target), ())
        ),
        "components_per_layer": k,
        "arms": arm_table,
        "contrasts": {key: contrast(*value) for key, value in pairs.items()},
    }

    if run_parsimony:
        table = {}
        for name, (families, layers) in subsets.items():
            table[name] = {
                "families": tuple(families),
                "n_features": int(sum(
                    config["layer_blocks"][layer].shape[1] for layer in layers
                )),
                "median_delta_r2": float(
                    arm_table[name]["median_r2"]
                    - arm_table[MLC.PRIMARY_COMPARATOR]["median_r2"]
                ),
            }
        report["parsimony"] = MLC.parsimonious_subset(table)
        ranked = sorted(table.items(), key=lambda kv: -kv[1]["median_delta_r2"])
        report["parsimony"]["ranked_best_to_worst"] = [name for name, _ in ranked]
        report["parsimony"]["why_a_ranked_list_is_given_separately"] = (
            "The committed report is serialised with sorted keys, so the "
            "enumeration mapping is alphabetical. The ranking is therefore "
            "recorded explicitly rather than implied by key order."
        )
        report["parsimony"]["enumeration"] = {
            name: {
                "families": list(row["families"]),
                "n_features": row["n_features"],
                "median_delta_r2": row["median_delta_r2"],
            }
            for name, row in ranked
        }

    if run_permutations:
        # The null is computed for EVERY reported contrast, not only the
        # primary one. Doing it uniformly is what keeps it a control: a null
        # added afterwards to whichever contrast looked significant would be
        # selection, while a null added to all of them can only temper claims.
        needed = sorted({name for pair in pairs.values() for name in pair})
        rng = np.random.default_rng(MLC.FROZEN["seed"])
        null_deltas = {key: [] for key in pairs}
        for _ in range(MLC.FROZEN["permutation_replicates"]):
            order = rng.permutation(targets.shape[0])
            shuffled = targets[order]
            null_total = total_sums_of_squares(
                shuffled, folds, MLC.FROZEN["outer_repeats"]
            )
            null_r2 = {}
            for name in needed:
                layers, with_h = designs[name]
                block, _ = arm_errors(
                    config, cache, folds, shuffled, layers, k, with_h
                )
                null_r2[name] = MLC.r2_per_feature(block, null_total)
            for key, (new, base) in pairs.items():
                null_deltas[key].append(
                    float(np.nanmedian(null_r2[new] - null_r2[base]))
                )

        definition = (
            "The individual-to-target pairing is permuted before the whole "
            "cross-validation. Predictors, covariates and folds are untouched, "
            "so the component bases are identical across replicates -- an exact "
            "identity, not an approximation. The null is computed for every "
            "reported contrast, uniformly."
        )
        for key, values in null_deltas.items():
            entry = report["contrasts"][key]
            if entry.get("status") == "NOT_EVALUABLE":
                continue
            entry["permutation_null"] = MLC.permutation_rank(
                entry["median_delta_r2"], values
            )
            entry["permutation_null"]["replicates"] = {
                f"replicate_{index + 1:02d}": value
                for index, value in enumerate(values)
            }
        report["permutation_control"] = dict(
            report["contrasts"]["primary_A1_vs_A0"]["permutation_null"]
        )
        report["permutation_control"]["definition"] = definition
        report["permutation_control"]["arms_permuted"] = needed
        report["permutation_control"]["why_the_null_is_not_centred_on_zero"] = (
            "A wider arm shrinks harder under an inner-CV penalty, so it lands "
            "closer to the training mean than a narrow arm does and scores a "
            "slightly higher R-squared with no information at all. The null "
            "measures the size of that artefact in this exact design, which is "
            "why a small positive delta cannot be read as a finding without it."
        )

    primary = report["contrasts"]["primary_A1_vs_A0"]
    complementarity = report["contrasts"]["complementarity_A1_vs_A2"]
    permutation = report.get("permutation_control")
    criteria = {
        "median_delta_positive": bool(primary.get("median_delta_r2", 0.0) > 0.0),
        "interval_excludes_zero": bool(primary.get("interval_excludes_zero", False)),
        "exceeds_every_permutation_replicate": (
            bool(permutation["exceeds_every_replicate"]) if permutation else None
        ),
        "beats_the_same_family_sibling": bool(
            complementarity.get("median_delta_r2", 0.0) > 0.0
            and complementarity.get("interval_excludes_zero", False)
        ),
    }
    supported = (
        criteria["median_delta_positive"]
        and criteria["interval_excludes_zero"]
        and criteria["exceeds_every_permutation_replicate"] is True
    )
    report["prespecified_decision"] = {
        "criteria": criteria,
        "primary_claim_supported": bool(supported),
        "complementarity_claim_supported": bool(
            supported and criteria["beats_the_same_family_sibling"]
        ),
        "verdict": (
            "SUPPORTED_WITH_COMPLEMENTARITY"
            if supported and criteria["beats_the_same_family_sibling"]
            else "SUPPORTED_BUT_ATTRIBUTABLE_TO_THE_SAME_FAMILY_SIBLING"
            if supported
            else "NOT_SUPPORTED"
        ),
        "sentence_fixed_before_the_result": (
            "If criteria 1 to 3 hold but the sibling comparison fails, the "
            "conclusion is that the ensemble predicts the masked layer and the "
            "prediction is attributable to a second assay of the same "
            "molecules, not to inter-family information."
        ),
    }
    return report


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------
def run(source, *, dry_run=False, n_boot=None):
    import numpy as np

    started = time.time()
    freeze = check_freeze_against_document()
    container = load_container(source)
    configs = {
        name: build_configuration(container, name) for name in MLC.CONFIGURATIONS
    }

    report = {
        "objective": "Axis A -- multi-omic method, EATRIS-Plus v4",
        "cohort": "eatris_plus_v4",
        "source_sha256": container["source_sha256"],
        "prespecification": freeze,
        "status_ceiling": MLC.FROZEN["status_ceiling"],
        "no_pharmacological_claim": MLC.no_pharmacological_claim(),
        "compartment": (
            "Predictors are leukocyte transcriptome and methylome and plasma "
            "small molecules; the target is a plasma protein panel. Blood is "
            "not liver."
        ),
        "single_site": (
            "All individuals carry the CZC prefix: one biobank, so no centre "
            "confound and no internal replication either."
        ),
        "no_identifier_value_or_prediction_exported": True,
        "configurations": {
            name: {
                "role": config["role"],
                "layers_kept": config["layers_kept"],
                "n_layers": len(config["layers_kept"]),
                "n_individuals": config["n"],
                "individual_set_sha256": config["individual_set_sha256"],
                "families_present": config["families_present"],
                "column_screening": config["screened"],
                "haematology_missing_cells": config["haematology_missing_cells"],
            }
            for name, config in configs.items()
        },
        "log10_applied_to": sorted(
            name for name, flag in container["log10_applied"].items() if flag
        ),
    }
    if dry_run:
        report["dry_run"] = True
        report["note"] = "Freeze verified and configurations built. Nothing was fitted."
        return report

    primary_name = MLC.FROZEN["primary_configuration"]
    k = MLC.FROZEN["components_per_layer"]
    k_max = max((k,) + tuple(MLC.FROZEN["component_sensitivities"]))

    caches = {}
    shortfalls = {}
    folds = {}
    for name, config in configs.items():
        folds[name] = make_folds(config["n"])
        caches[name], shortfalls[name] = component_cache(config, folds[name], k_max)
    report["component_shortfalls"] = {
        name: value for name, value in shortfalls.items() if value
    }

    results = {}

    # --- primary configuration, primary target ---------------------------
    config = configs[primary_name]
    results[f"{primary_name}::proteins"] = run_target(
        config,
        caches[primary_name],
        folds[primary_name],
        MLC.FROZEN["primary_target"],
        k=k,
        run_permutations=True,
        run_parsimony=True,
        n_boot=n_boot,
    )
    nine_layer_key = next(
        (
            key
            for key in results[f"{primary_name}::proteins"]["contrasts"]
            if key.startswith("two_rna_families_")
        ),
        None,
    )
    results[f"{primary_name}::proteins"]["nine_layer_comparison_inside_the_set"] = {
        "definition": (
            "The eleven-layer arm against the nine-layer arm on the SAME 79 "
            "individuals and the SAME folds, which isolates mRNA and "
            "miRNA-mature at fixed n."
        ),
        "contrast_key": nine_layer_key,
    }

    # --- secondary target, primary configuration -------------------------
    results[f"{primary_name}::acylcarnitines"] = run_target(
        config,
        caches[primary_name],
        folds[primary_name],
        MLC.FROZEN["secondary_target"],
        k=k,
        run_permutations=True,
        run_parsimony=False,
        n_boot=n_boot,
    )

    # --- component sensitivities -----------------------------------------
    sensitivity = {}
    for alternative in MLC.FROZEN["component_sensitivities"]:
        targets = config["layer_blocks"][MLC.FROZEN["primary_target"]]
        total = total_sums_of_squares(targets, folds[primary_name],
                                      MLC.FROZEN["outer_repeats"])
        arms = MLC.arms_for(MLC.FROZEN["primary_target"], config["layers_kept"])
        block = {}
        for name in (MLC.PRIMARY_COMPARATOR, MLC.PRIMARY_ARM,
                     MLC.COMPLEMENTARITY_COMPARATOR):
            block[name], _ = arm_errors(
                config, caches[primary_name], folds[primary_name], targets,
                arms[name], alternative, False,
            )
        entry = {
            "median_r2": {
                name: float(np.nanmedian(MLC.r2_per_feature(value, total)))
                for name, value in block.items()
            },
            "primary_contrast": MLC.paired_bootstrap_median_delta(
                block[MLC.PRIMARY_ARM], block[MLC.PRIMARY_COMPARATOR], total,
                n_boot=n_boot,
            ),
        }
        sensitivity[f"k={alternative}"] = entry
    results[f"{primary_name}::proteins"]["component_sensitivity"] = sensitivity

    # --- C12: what the qRT-PCR layer adds, at fixed n --------------------
    twelve = configs["C12_all_layers"]
    all_other = MLC.arms_for(MLC.FROZEN["primary_target"], twelve["layers_kept"])[
        MLC.PRIMARY_ARM
    ]
    minus_qrtpcr = tuple(x for x in all_other if x != "miRNA-qRT-PCR")
    results["C12_all_layers::proteins"] = run_target(
        twelve,
        caches["C12_all_layers"],
        folds["C12_all_layers"],
        MLC.FROZEN["primary_target"],
        k=k,
        run_permutations=True,
        run_parsimony=True,
        extra_arms={"A1_minus_miRNA_qRT_PCR": minus_qrtpcr},
        n_boot=n_boot,
    )
    results["C12_all_layers::proteins"]["what_the_qrtpcr_layer_adds"] = {
        "definition": (
            "All non-target families including miRNA-qRT-PCR against the same "
            "arm with that single layer removed, on the SAME 55 individuals "
            "and the SAME folds. The difference is attributable to the layer."
        ),
        "contrast_key": f"{MLC.PRIMARY_ARM}_vs_A1_minus_miRNA_qRT_PCR",
    }

    # --- C9: wider coverage, fewer modalities ----------------------------
    results["C9_no_rna::proteins"] = run_target(
        configs["C9_no_rna"],
        caches["C9_no_rna"],
        folds["C9_no_rna"],
        MLC.FROZEN["primary_target"],
        k=k,
        run_permutations=True,
        run_parsimony=True,
        n_boot=n_boot,
    )
    results["C9_no_rna::proteins"]["forbidden_comparison"] = (
        "The nine-layer arm on 114 individuals and the nine-layer arm on 79 "
        "individuals are DIFFERENT PEOPLE in DIFFERENT FOLDS. No paired "
        "interval is computed between them and no difference between them may "
        "be attributed to a modality. Performance differences between distinct "
        "sample sizes are not modality effects."
    )

    report["results"] = results
    report["reporting_separation"] = {
        "predictive_performance": "arms, contrasts, permutation null",
        "biological_relevance": (
            "families and target features, each statement carrying the "
            "compartment caveat, and none carrying a pharmacological reading"
        ),
        "analytical_validity": (
            "fold safety, component shortfalls, individual-set checksums, "
            "covariate completeness, the optimism of the parsimony selection, "
            "and the single-site limitation"
        ),
    }
    report["runtime_seconds"] = round(time.time() - started, 1)
    report["how_to_read_runtime"] = (
        "Wall clock on one machine. A suspended process inflates it without "
        "computing anything, so it is never read as a cost model."
    )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, help="EATRIS_PLUS_mudata.h5mu")
    parser.add_argument("--out", default=None)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--bootstrap", type=int, default=None,
                        help="override only for a smoke test; the freeze says 2000")
    args = parser.parse_args(argv)
    # The freeze-table labels carry Greek and superscript characters, and a
    # Windows console defaults to cp1252. Reconfigure rather than drop them,
    # so the printed report is the same text as the written one.
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        report = run(args.source, dry_run=args.dry_run, n_boot=args.bootstrap)
    except (RunError, MLC.LayerComplementarityError) as exc:
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
