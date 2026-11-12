"""Enumerate the EATRIS-Plus v4 layer configurations, before any model is fitted.

This script fits nothing and predicts nothing. It exists so that the three
individual sets used by the Axis A prespecification -- eleven layers on the
larger complete case, twelve layers on the smaller one, and nine layers
without either RNA layer -- are *derived from the container* and recorded with
a checksum, rather than asserted from the ISA workbook (which over-counts: it
claims 127 individuals for every assay sheet and the deposited data does not
agree).

No individual identifier is printed. Each individual set is summarised by its
size and by the SHA-256 of its sorted identifiers, which is enough to prove
later that the runner used the same people without publishing who they are.

Usage:
    python -I scripts/audit_eatris_plus_layer_configurations.py \
        --source /path/to/EATRIS_PLUS_mudata.h5mu
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

EXPECTED_SHA256_PREFIX = "a764a6cd7246edbd"

#: Twelve deposited layers grouped into six biological families. The grouping
#: is declared HERE, before any fit, because the commissioning note requires
#: family-level ablation: reconstructing one plasma-proteome measurement from
#: another plasma-proteome measurement is not inter-omic complementarity.
FAMILIES = {
    "transcriptome": ("mRNA",),
    "small_rna": ("miRNA-mature", "miRNA-qRT-PCR"),
    "methylome": ("EM-seq",),
    "proteome": ("proteins", "proteomics-ipp"),
    "lipidome": (
        "lipidomics-neg",
        "lipidomics-pos",
        "lipidomics-spline",
        "fattyacids",
    ),
    "metabolome_targeted": ("acylcarnitines", "aminoacids"),
}

#: Configurations, named by (layers, intended complete case). The layer lists
#: are the decision the commissioning note froze; the individual counts are
#: measured here and must agree with the note.
CONFIGURATIONS = {
    "C11_primary": {
        "drop_layers": ("miRNA-qRT-PCR",),
        "expected_n": 79,
        "role": "primary",
    },
    "C12_all_layers": {
        "drop_layers": (),
        "expected_n": 55,
        "role": "secondary_what_the_qrtpcr_layer_adds",
    },
    "C9_no_rna": {
        "drop_layers": ("miRNA-qRT-PCR", "miRNA-mature", "mRNA"),
        "expected_n": 114,
        "role": "secondary_wider_coverage_without_the_transcriptome",
    },
}

#: Person-level covariates admitted as the baseline block. Age.group and
#: BMI.group are deliberately excluded: they are coarsenings of Age and BMI
#: and would enter the same information twice.
COVARIATES_NUMERIC = (
    "Age",
    "BMI",
    "Height",
    "Weight",
    "Leukocytes",
    "Neutrophils",
    "Lymphocytes",
    "Hematocrit",
    "Erythrocytes",
    "Hemoglobin",
    "MCV",
    "MCH",
    "MCHC",
    "MXD",
    "Platelets",
)
COVARIATES_BINARY = ("Sex",)


class AuditError(RuntimeError):
    pass


def sha256_of_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def set_checksum(identifiers):
    """SHA-256 of the sorted identifiers, so a set can be proved identical
    across runs without any identifier leaving this machine."""
    payload = "\n".join(sorted(str(x) for x in identifiers)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _decode(values):
    out = []
    for value in values:
        out.append(value.decode("utf-8") if isinstance(value, bytes) else str(value))
    return out


def read_index(group):
    import h5py  # noqa: F401  (imported for the type contract only)

    if "_index" in group:
        return _decode(group["_index"][...])
    attr = group.attrs.get("_index")
    key = attr.decode("utf-8") if isinstance(attr, bytes) else str(attr)
    return _decode(group[key][...])


def read_obs_column(obs_group, name):
    """Return a column of the global obs table, decoding h5ad categoricals."""
    node = obs_group[name]
    if hasattr(node, "keys") and "codes" in node:
        codes = node["codes"][...]
        categories = _decode(node["categories"][...])
        return [categories[c] if c >= 0 else None for c in codes]
    return _decode(node[...]) if node.dtype.kind in "SOU" else list(node[...])


def audit(source):
    import h5py
    import numpy as np

    checksum = sha256_of_file(source)
    if not checksum.startswith(EXPECTED_SHA256_PREFIX):
        raise AuditError(
            f"source checksum {checksum} does not start with the registered "
            f"prefix {EXPECTED_SHA256_PREFIX}"
        )

    report = {
        "source": Path(source).name,
        "source_sha256": checksum,
        "fits_nothing": True,
        "identifiers_are_not_printed": True,
        "layers": {},
        "families": {k: list(v) for k, v in FAMILIES.items()},
        "configurations": {},
    }

    with h5py.File(source, "r") as handle:
        modalities = sorted(handle["mod"].keys())
        per_layer_ids = {}
        for name in modalities:
            group = handle["mod"][name]
            ids = read_index(group["obs"])
            features = read_index(group["var"])
            shape = group["X"].shape
            if (len(ids), len(features)) != tuple(shape):
                raise AuditError(
                    f"layer {name}: obs/var lengths {(len(ids), len(features))} "
                    f"disagree with X shape {tuple(shape)}"
                )
            if len(set(ids)) != len(ids):
                raise AuditError(f"layer {name}: duplicated individual identifiers")
            per_layer_ids[name] = set(ids)
            report["layers"][name] = {
                "individuals": len(ids),
                "features": len(features),
            }

        declared = sorted(x for layers in FAMILIES.values() for x in layers)
        if declared != modalities:
            raise AuditError(
                "the declared family grouping is not an exact partition of the "
                f"deposited layers: declared {declared}, deposited {modalities}"
            )

        global_ids = read_index(handle["obs"])
        union = set()
        for ids in per_layer_ids.values():
            union |= ids
        report["union_individuals"] = len(union)
        report["global_obs_rows"] = len(global_ids)
        report["every_layer_individual_has_a_person_row"] = union <= set(global_ids)
        if not report["every_layer_individual_has_a_person_row"]:
            raise AuditError("a layer carries individuals absent from the person table")

        # How many layers each individual carries -- the coverage profile that
        # makes the three configurations the only sensible ones.
        coverage = {}
        for identifier in union:
            count = sum(1 for ids in per_layer_ids.values() if identifier in ids)
            coverage[count] = coverage.get(count, 0) + 1
        report["individuals_by_number_of_layers_held"] = {
            str(k): coverage[k] for k in sorted(coverage, reverse=True)
        }

        covariate_rows = {}
        for name in COVARIATES_NUMERIC + COVARIATES_BINARY:
            if name not in handle["obs"]:
                raise AuditError(f"declared covariate {name!r} is not in the person table")
            covariate_rows[name] = read_obs_column(handle["obs"], name)
        index_of = {identifier: i for i, identifier in enumerate(global_ids)}

        for config_name, spec in CONFIGURATIONS.items():
            kept = [m for m in modalities if m not in spec["drop_layers"]]
            complete = set(union)
            for layer in kept:
                complete &= per_layer_ids[layer]
            families_present = {
                family: [layer for layer in layers if layer in kept]
                for family, layers in FAMILIES.items()
            }
            families_present = {k: v for k, v in families_present.items() if v}

            missing_covariates = {}
            for name, column in covariate_rows.items():
                absent = 0
                for identifier in complete:
                    value = column[index_of[identifier]]
                    if value is None or (
                        isinstance(value, str) and value.strip() in ("", "NA", "nan")
                    ):
                        absent += 1
                    elif isinstance(value, float) and not np.isfinite(value):
                        absent += 1
                if absent:
                    missing_covariates[name] = absent

            entry = {
                "role": spec["role"],
                "layers_kept": kept,
                "n_layers": len(kept),
                "layers_dropped": list(spec["drop_layers"]),
                "complete_case_n": len(complete),
                "expected_n": spec["expected_n"],
                "agrees_with_commissioned_count": len(complete) == spec["expected_n"],
                "individual_set_sha256": set_checksum(complete),
                "families_present": families_present,
                "n_families": len(families_present),
                "features_total": sum(
                    report["layers"][layer]["features"] for layer in kept
                ),
                "covariate_missing_counts": missing_covariates,
            }
            report["configurations"][config_name] = entry

        # Nesting matters for the factorial reading the note demands: the
        # miRNA-qRT-PCR question is only answerable at fixed n if the twelve
        # layer set is a subset of the eleven layer set's people.
        s12 = set(union)
        for layer in CONFIGURATIONS["C12_all_layers"]["drop_layers"]:
            pass
        for layer in modalities:
            s12 &= per_layer_ids[layer]
        s11 = set(union)
        for layer in [m for m in modalities if m != "miRNA-qRT-PCR"]:
            s11 &= per_layer_ids[layer]
        s9 = set(union)
        for layer in [
            m for m in modalities
            if m not in ("miRNA-qRT-PCR", "miRNA-mature", "mRNA")
        ]:
            s9 &= per_layer_ids[layer]
        report["nesting"] = {
            "C12_subset_of_C11": s12 <= s11,
            "C11_subset_of_C9": s11 <= s9,
            "C12_subset_of_C9": s12 <= s9,
            "why_it_matters": (
                "Nested sets let the layer question and the sample-size "
                "question be separated: the qRT-PCR layer is judged on the "
                "SAME people, and n is judged with the SAME layers. A "
                "performance difference between sets of different size is "
                "never read as a modality effect."
            ),
        }

        report["candidate_target_layers"] = {
            "proteins": {
                "features": report["layers"]["proteins"]["features"],
                "individuals": report["layers"]["proteins"]["individuals"],
                "same_family_sibling": "proteomics-ipp",
                "biochemically_adjacent_families": [],
            },
            "acylcarnitines": {
                "features": report["layers"]["acylcarnitines"]["features"],
                "individuals": report["layers"]["acylcarnitines"]["individuals"],
                "same_family_sibling": "aminoacids",
                "biochemically_adjacent_families": ["lipidome"],
            },
        }

    report["no_pharmacological_endpoint"] = (
        "This cohort holds no drug, dose, concentration, time, enzyme activity "
        "or clinical endpoint. Nothing computed on it may be reported as "
        "evidence about drug metabolism or pharmacokinetics."
    )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, help="EATRIS_PLUS_mudata.h5mu")
    parser.add_argument("--out", default=None, help="write JSON here as well as stdout")
    args = parser.parse_args(argv)
    try:
        report = audit(args.source)
    except AuditError as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    text = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
