#!/usr/bin/env python
"""F-series: do signed regulatory activities add anything to the target transcript?

Frozen by `docs/OMICS2PK_F_SERIES_FUNCTIONAL_OMICS_PRESPECIFICATION_20261109.md`
version 1.0.0. The runner refuses to load any data until it has checked that
every value in that document's section 12 still matches
`omics2pk.functional_omics.FROZEN`.

The estimator is not reimplemented. The ridge path, the inner-fold penalty
choice, the trimmed mean and the paired donor bootstrap are **imported by path**
from `scripts/run_h2a_hlc_selection_comparison.py`, so the F-series is H2a's
estimator with different columns and nothing else. Importing rather than copying
means the two cannot drift apart without a test noticing.

Three blocks come out, and they are never merged:

* `predictive_performance` -- the only block that answers the primary question;
* `biological_relevance` -- whether fitted coefficients point the way CollecTRI
  says they should, which is a question about biology, not about prediction;
* `analytical_validity` -- what each guard removed, what is `NOT_EVALUABLE` and
  why, and the algebraic identity check.

No donor identifier, no per-donor value and no expression row is written out.

Run from the repository root:

    python scripts/run_f_series_hlc_functional_omics.py \
        --hlc "<directory with curatedExpressionLiver.zip and the phenotype zip>" \
        --network "<directory>/collectri.tsv" \
        --identity "<directory>/probe_identity_full_table.json" \
        --out reports/data/f_series_hlc_functional_omics.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import time
from datetime import date, timezone, datetime
from pathlib import Path
from zipfile import ZipFile

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk import functional_omics as FO  # noqa: E402


def _import_h2a():
    """Load H2a's estimator from its own file, so the two cannot drift apart."""
    import importlib.util

    path = REPO_ROOT / "scripts" / "run_h2a_hlc_selection_comparison.py"
    if not path.is_file():
        raise FO.FunctionalOmicsError(f"H2a runner absent: {path}")
    spec = importlib.util.spec_from_file_location("_h2a_estimator", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H2A = _import_h2a()

PHENOTYPE_MEMBER = "phenotype.txt"
EXPRESSION_MEMBER = "expression.txt"
FEATURES_MEMBER = "features.txt"
POOL_NAME = "A_qualified_human"
POOL_PATH = REPO_ROOT / "registry" / "QUALIFIED_PROBE_POOLS.json"
PRESPEC_PATH = REPO_ROOT / "docs" / (
    "OMICS2PK_F_SERIES_FUNCTIONAL_OMICS_PRESPECIFICATION_20261109.md")

ENDPOINTS = tuple(FO.ENDPOINT_ENZYME)
AGE_ROW = "AGE_(YRS)"
SEX_ROW = "GENDER"
SEX_CODES = {"Male": 1.0, "Female": 0.0}


# --------------------------------------------------------------------------
# Governance: the frozen table has to still be the frozen table
# --------------------------------------------------------------------------
def check_estimator_machinery_unchanged():
    """The F-series must use H2a's estimator, not a second implementation of it."""
    expected = {
        "OUTER_FOLDS": FO.FROZEN["outer_folds"],
        "OUTER_REPEATS": FO.FROZEN["outer_repeats"],
        "INNER_FOLDS": FO.FROZEN["inner_folds"],
        "N_BOOTSTRAP": FO.FROZEN["n_bootstrap"],
        "TRIM_FRACTION": FO.FROZEN["trim_fraction"],
        "MAX_PROBE_MISSING_FRACTION": FO.FROZEN["max_probe_missing_fraction"],
    }
    for name, value in expected.items():
        found = getattr(H2A, name, None)
        if found != value:
            raise FO.FunctionalOmicsError(
                f"H2a's {name} is {found!r}, the F-series freeze says {value!r}")
    if tuple(H2A.RIDGE_ALPHAS) != tuple(FO.FROZEN["ridge_alphas"]):
        raise FO.FunctionalOmicsError("the ridge grid differs from the freeze")
    if H2A.SEED == FO.FROZEN["seed"]:
        raise FO.FunctionalOmicsError(
            "the F-series fold seed equals H2a's; the freeze states a distinct seed")
    return True


_FROZEN_IN_DOCUMENT = {
    "Fold seed": ("seed", lambda text: int(text.strip("`"))),
    "Bootstrap draws": ("n_bootstrap", lambda text: int(text.replace(",", ""))),
    "Relative MAE threshold": ("relative_mae_threshold", float),
    "Trim fraction": ("trim_fraction", float),
    "Regulon floor (admissible members)": ("regulon_floor", int),
    "Regulator floor per endpoint": ("regulator_floor", int),
    "Donor floor per endpoint": ("donor_floor", int),
    "Control replicates": ("control_replicates",
                           lambda text: int(text.split()[0])),
    "Status ceiling": ("status_ceiling", lambda text: text.strip("`")),
}


def check_freeze_against_document(path=PRESPEC_PATH):
    """Read section 12 back out of the prespecification and compare it.

    A freeze that lives only in code is a freeze one edit away from moving. A
    freeze that lives only in prose is one nobody checks. This reads the
    document's own table and refuses to run if the two disagree.
    """
    if not path.is_file():
        raise FO.FunctionalOmicsError(f"prespecification absent: {path}")
    text = path.read_text(encoding="utf-8")
    checked = {}
    for label, (key, cast) in _FROZEN_IN_DOCUMENT.items():
        pattern = re.compile(
            r"^\|\s*" + re.escape(label) + r"\s*\|\s*(.+?)\s*\|\s*$", re.MULTILINE)
        match = pattern.search(text)
        if not match:
            raise FO.FunctionalOmicsError(
                f"the prespecification no longer states {label!r}")
        stated = cast(match.group(1))
        if stated != FO.FROZEN[key]:
            raise FO.FunctionalOmicsError(
                f"{label}: the document says {stated!r}, the module says "
                f"{FO.FROZEN[key]!r}")
        checked[label] = stated
    # The primary contrast is named in prose, not in the table, so it is checked
    # as prose: the document must state which arm is compared against which on
    # which endpoint, or the freeze does not identify the experiment.
    primary_sentence = (f"`{FO.PRIMARY_ARM}` versus `{FO.PRIMARY_COMPARATOR}` on "
                        f"`{FO.PRIMARY_ENDPOINT}`")
    if primary_sentence not in text:
        raise FO.FunctionalOmicsError(
            f"the prespecification does not state the primary contrast "
            f"({primary_sentence})")
    return {"checked": checked, "primary_contrast_stated": primary_sentence,
            "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest()}


def file_digest(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_hlc(hlc_dir: Path):
    import pandas as pd

    expression_zip = hlc_dir / "curatedExpressionLiver.zip"
    phenotype_zip = hlc_dir / "curatedPhenotypesynapse.zip"
    for path in (expression_zip, phenotype_zip):
        if not path.is_file():
            raise FO.FunctionalOmicsError(f"expected archive absent: {path}")
    with ZipFile(expression_zip) as archive:
        with archive.open(FEATURES_MEMBER) as handle:
            features = pd.read_csv(handle, sep="\t", dtype=str)
        with archive.open(EXPRESSION_MEMBER) as handle:
            expression = pd.read_csv(handle, sep="\t", index_col=0)
    with ZipFile(phenotype_zip) as archive:
        with archive.open(PHENOTYPE_MEMBER) as handle:
            phenotype = pd.read_csv(handle, sep="\t", index_col=0, dtype=str)

    if expression.shape != (H2A.EXPECTED_PROBES, H2A.EXPECTED_EXPRESSION_DONORS):
        raise FO.FunctionalOmicsError(
            f"expression dimensions changed: {expression.shape}")
    if not np.array_equal(expression.index.to_numpy(dtype=str),
                          features["feature_id"].to_numpy(dtype=str)):
        raise FO.FunctionalOmicsError("probe identifiers no longer align")

    donors = expression.columns.astype(str).to_numpy()
    missing = sorted(set(donors) - set(phenotype.columns.astype(str)))
    if missing:
        raise FO.FunctionalOmicsError(
            f"{len(missing)} expression donors absent from the phenotype table")
    phenotype = phenotype.loc[:, donors]
    age = pd.to_numeric(phenotype.loc[AGE_ROW], errors="coerce").to_numpy(float)
    sex = phenotype.loc[SEX_ROW].map(SEX_CODES).to_numpy(float)
    return {
        "values": expression.to_numpy(dtype=np.float64),
        "symbols": features["genesymbol"].fillna("").to_numpy(dtype=str),
        "feature_ids": features["feature_id"].astype(str).to_numpy(dtype=str),
        "donors": donors,
        "age": age,
        "sex": sex,
        "phenotype": phenotype,
    }


def load_network(path: Path):
    import pandas as pd

    frame = pd.read_csv(path, sep="\t", dtype=str)
    frame.columns = [str(name).strip().lower() for name in frame.columns]
    if "sources" in frame.columns:
        with_collectri = int(frame["sources"].astype(str)
                             .str.contains("CollecTRI", case=False).sum())
    else:
        with_collectri = -1
    return frame, {
        "file": path.name,
        "sha256": file_digest(path),
        "rows": int(len(frame)),
        "rows_naming_collectri": with_collectri,
    }


# --------------------------------------------------------------------------
# The two universes: scoring genes (no cytochrome) and the target transcript
# --------------------------------------------------------------------------
def build_universes(data, pool_path=POOL_PATH, pool_name=POOL_NAME):
    """Admissible probes, grouped twice.

    `scoring` holds every admissible probe whose gene is **not** a cytochrome,
    by symbol: this is the universe the activity scores are computed over, and
    the three guards of prespecification section 4.3 all end here.

    `cytochrome` holds the admissible probes of the cytochromes, by symbol.
    They are available as a *predictor column* for the arms that carry the
    target transcript, and they never enter a score. Keeping the two sets in
    separate dictionaries is what makes that separation checkable instead of
    merely intended.
    """
    pool = json.loads(Path(pool_path).read_text(encoding="utf-8"))
    admissible = pool.get("admissible_features", {}).get(pool_name)
    if not admissible:
        raise FO.FunctionalOmicsError(f"no admissible set named {pool_name!r}")
    admissible = set(map(str, admissible))

    scoring: dict[str, list[int]] = {}
    cytochrome: dict[str, list[int]] = {}
    refused_alias = 0
    for position, identifier in enumerate(data["feature_ids"]):
        if str(identifier) not in admissible:
            continue
        symbol = str(data["symbols"][position])
        if not symbol:
            continue
        resolved = str(data["resolved"].get(str(identifier), "") or "")
        if FO.is_cytochrome(symbol):
            cytochrome.setdefault(symbol, []).append(position)
            continue
        if FO.is_cytochrome(resolved):
            refused_alias += 1
            continue
        scoring.setdefault(symbol, []).append(position)
    if not scoring:
        raise FO.FunctionalOmicsError("the scoring universe is empty")

    FO.assert_no_cytochrome_reaches_scoring(scoring.keys())
    return {
        "scoring": scoring,
        "cytochrome": cytochrome,
        "admissible": admissible,
        "pool": {"file": Path(pool_path).name, "pool_name": pool_name,
                 "rule_version": pool.get("rule_version"),
                 "admissible_probes": len(admissible),
                 "scoring_symbols": len(scoring),
                 "cytochrome_symbols": len(cytochrome),
                 "refused_by_resolved_alias": refused_alias},
    }


def choose_probes(values, groups, train):
    """One probe per symbol, picked on the training fold by H2a's rule."""
    chosen: dict[str, int] = {}
    floor = 1.0 - FO.FROZEN["max_probe_missing_fraction"]
    for symbol, positions in groups.items():
        best, best_mean = None, -np.inf
        for position in positions:
            observed = values[position, train]
            present = np.isfinite(observed)
            if present.mean() < floor:
                continue
            mean = float(observed[present].mean())
            if mean > best_mean:
                best, best_mean = position, mean
        if best is not None:
            chosen[symbol] = best
    return chosen


def assemble(values, chosen, donor_index, train):
    """Rows for the chosen probes over `donor_index`, training-mean imputed."""
    genes = sorted(chosen)
    block = np.empty((len(genes), len(donor_index)), dtype=np.float64)
    for row, symbol in enumerate(genes):
        series = values[chosen[symbol], donor_index].copy()
        absent = ~np.isfinite(series)
        if absent.any():
            observed = series[train][np.isfinite(series[train])]
            series[absent] = float(observed.mean()) if observed.size else 0.0
        block[row] = series
    return genes, block


# --------------------------------------------------------------------------
# Per-endpoint specification, fixed before any fitting
# --------------------------------------------------------------------------
def specify_endpoints(network_frame, universes, data, log):
    """Candidate regulators, evaluable regulators and evaluability, per endpoint."""
    signed = FO.resolve_signed_network(network_frame)
    regulons = signed.regulons()
    unsigned = FO.resolve_signed_network(network_frame, drop_cytochrome_targets=False)
    scoring_symbols = set(universes["scoring"])

    import pandas as pd

    expression_donors = set(map(str, data["donors"]))
    specification = {}
    for endpoint in ENDPOINTS:
        enzyme = FO.ENDPOINT_ENZYME[endpoint]
        candidates = unsigned.regulators_of(enzyme)
        evaluable, refused = FO.evaluable_regulators(
            regulons, candidates, scoring_symbols)
        series = pd.to_numeric(data["phenotype"].loc[endpoint], errors="coerce")
        with_value = {str(name) for name in series.dropna().index}
        donors_with_both = len(with_value & expression_donors)
        target_probes = universes["cytochrome"].get(enzyme, [])
        verdict = FO.endpoint_evaluability(
            endpoint=endpoint,
            target_admissible_probes=len(target_probes),
            evaluable_regulator_count=len(evaluable),
            donors_with_both=donors_with_both)
        # F2's columns: the regulators' own transcripts, complexes expanded.
        transcripts = sorted({member
                              for regulator in evaluable
                              for member in FO.complex_members(regulator)
                              if member in scoring_symbols})
        recorded = {regulator: unsigned.edges.get((regulator, enzyme))
                    for regulator in evaluable}
        specification[endpoint] = {
            "enzyme": enzyme,
            "candidate_regulators": candidates,
            "evaluable_regulators": evaluable,
            "refused_regulators": refused,
            "regulator_transcripts": transcripts,
            "recorded_edge_sign": recorded,
            "target_admissible_probes": len(target_probes),
            "donors_with_endpoint_and_expression": donors_with_both,
            "evaluability": verdict,
            "regulon_sizes": {
                regulator: {
                    "members": len(regulons.get(regulator, {})),
                    "admissible_members": sum(
                        1 for gene in regulons.get(regulator, {})
                        if gene in scoring_symbols)}
                for regulator in evaluable},
        }
        log(f"  {endpoint}: enzyme {enzyme}, candidates {len(candidates)}, "
            f"evaluable {len(evaluable)}, donors {donors_with_both}, "
            f"target probes {len(target_probes)}, "
            f"F5 {verdict['verdict']['F5']}")
    return specification, signed, regulons, unsigned


# --------------------------------------------------------------------------
# Control settings, enumerated before any of them is run
# --------------------------------------------------------------------------
def control_settings(replicates):
    """Every activity variant to compute, named and seeded in advance."""
    settings = [{"name": "real", "kind": "real", "fraction": None, "replicate": 0,
                 "seed": None}]
    for fraction in FO.FROZEN["perturbation_fractions"]:
        for replicate in range(replicates):
            settings.append({
                "name": f"sign_perturbation_{fraction:.2f}_{replicate:02d}",
                "kind": "sign_perturbation", "fraction": float(fraction),
                "replicate": replicate, "seed": FO.FROZEN["seed"] + replicate})
    for replicate in range(replicates):
        settings.append({
            "name": f"rewired_null_{replicate:02d}", "kind": "rewired_null",
            "fraction": None, "replicate": replicate,
            "seed": FO.FROZEN["seed"] + 100 + replicate})
    for replicate in range(replicates):
        settings.append({
            "name": f"donor_permuted_{replicate:02d}", "kind": "donor_permuted",
            "fraction": None, "replicate": replicate,
            "seed": FO.FROZEN["seed"] + 200 + replicate})
    settings.append({
        "name": "algebraic_identity_check", "kind": "sign_perturbation",
        "fraction": FO.FROZEN["identity_check_fraction"], "replicate": 0,
        "seed": FO.FROZEN["seed"]})
    return settings


def activities_for_setting(setting, *, regulons, regulators, genes, standardised,
                           universe_genes):
    """The donors-by-regulators activity block for one control setting."""
    kind = setting["kind"]
    if kind == "real":
        weights = FO.weight_matrix(regulons, regulators, genes)
        return FO.ulm_activities(standardised, weights)
    rng = np.random.default_rng(setting["seed"])
    if kind == "sign_perturbation":
        altered = FO.perturb_regulon_signs(regulons, setting["fraction"], rng)
        weights = FO.weight_matrix(altered, regulators, genes)
        return FO.ulm_activities(standardised, weights)
    if kind == "rewired_null":
        altered = FO.rewire_regulons(regulons, universe_genes, rng)
        weights = FO.weight_matrix(altered, regulators, genes)
        return FO.ulm_activities(standardised, weights)
    if kind == "donor_permuted":
        weights = FO.weight_matrix(regulons, regulators, genes)
        scores = FO.ulm_activities(standardised, weights)
        return scores[rng.permutation(scores.shape[0]), :]
    raise FO.FunctionalOmicsError(f"unknown control kind {kind!r}")


# --------------------------------------------------------------------------
# One run
# --------------------------------------------------------------------------
def run(args, log):
    check_estimator_machinery_unchanged()
    freeze = check_freeze_against_document()
    log(f"freeze verified against the prespecification "
        f"({len(freeze['checked'])} values)")

    identity = json.loads(Path(args.identity).read_text(encoding="utf-8"))
    resolved = {str(key): str((value or {}).get("human_symbol") or "")
                for key, value in identity.items()}

    data = load_hlc(Path(args.hlc))
    data["resolved"] = resolved
    network_frame, network_audit = load_network(Path(args.network))
    log(f"network {network_audit['rows']} rows, "
        f"{network_audit['rows_naming_collectri']} naming CollecTRI")

    universes = build_universes(data)
    guards = FO.cytochrome_guard_audit(
        array_symbols=data["symbols"], feature_ids=data["feature_ids"],
        resolved_symbols=resolved, admissible_features=universes["admissible"])
    log(f"scoring universe {len(universes['scoring'])} genes; cytochrome guards "
        f"refused {guards['refused_by_array_symbol']} by symbol and "
        f"{guards['refused_by_resolved_alias']} by resolved alias")

    log("per-endpoint specification:")
    specification, signed, regulons, unsigned = specify_endpoints(
        network_frame, universes, data, log)

    replicates = int(args.control_replicates)
    settings = control_settings(replicates)
    log(f"{len(settings)} activity settings per fold "
        f"({replicates} replicates per control)")

    # One shared activity computation per fold serves every endpoint: the union
    # of candidate regulators is scored once and each endpoint takes its own
    # columns, so a control replicate costs one matrix product rather than ten.
    union = sorted({regulator
                    for spec in specification.values()
                    for regulator in spec["evaluable_regulators"]})
    log(f"union of evaluable regulators across endpoints: {len(union)}")

    if getattr(args, "dry_run", False):
        # Everything above is a measurement of the sources, already stated in the
        # prespecification. Stopping here exercises the loading, the three
        # cytochrome guards, the sign rules and the evaluability verdicts without
        # fitting a single model, so the machinery can be validated before the
        # freeze without any contrast existing.
        log("dry run: stopping before any model is fitted")
        return {"dry_run": True, "specification": specification,
                "universes": universes, "guards": guards,
                "network_audit": network_audit, "freeze": freeze,
                "signed_network_audit": signed.audit, "union": union,
                "settings": settings}

    donors = data["donors"]
    folds = H2A.repeated_folds(len(donors), n_splits=FO.FROZEN["outer_folds"],
                               n_repeats=FO.FROZEN["outer_repeats"],
                               seed=FO.FROZEN["seed"])
    repeats, repeat_of = H2A.repeat_assignment(
        len(folds), n_splits=FO.FROZEN["outer_folds"], repeated=True)

    import pandas as pd

    endpoint_rows = {}
    for endpoint in ENDPOINTS:
        series = pd.to_numeric(data["phenotype"].loc[endpoint], errors="coerce")
        values = series.to_numpy(float)
        keep = np.nonzero(np.isfinite(values))[0]
        endpoint_rows[endpoint] = {"keep": keep, "y": values[keep]}

    # errors[endpoint][setting][arm] -> (repeats, n_donors) absolute errors
    errors: dict = {endpoint: {} for endpoint in ENDPOINTS}
    coefficients: dict = {endpoint: {} for endpoint in ENDPOINTS}
    identity_predictions: dict = {}
    chosen_probe_trace: dict = {}

    started = time.time()
    fold_seconds: list[float] = []
    for fold_number, (train, test) in enumerate(folds):
        repeat = repeat_of[fold_number]
        chosen_scoring = choose_probes(data["values"], universes["scoring"], train)
        genes, block = assemble(data["values"], chosen_scoring,
                                np.arange(len(donors)), train)
        FO.assert_no_cytochrome_reaches_scoring(genes, resolved_symbols=None)
        standardised = FO.standardise_on_training(block, train)
        chosen_cyp = choose_probes(data["values"], universes["cytochrome"], train)
        if fold_number == 0:
            chosen_probe_trace = {
                "fold": 0, "scoring_genes": len(genes),
                "cytochrome_genes": len(chosen_cyp),
                "scoring_probe_examples": {
                    gene: str(data["feature_ids"][chosen_scoring[gene]])
                    for gene in sorted(chosen_scoring)[:3]}}

        activity_cache = {}
        for setting in settings:
            activity_cache[setting["name"]] = activities_for_setting(
                setting, regulons=regulons, regulators=union, genes=genes,
                standardised=standardised, universe_genes=genes)
        union_position = {regulator: index for index, regulator in enumerate(union)}

        for endpoint in ENDPOINTS:
            spec = specification[endpoint]
            verdict = spec["evaluability"]["verdict"]
            rows = endpoint_rows[endpoint]
            keep, y = rows["keep"], rows["y"]
            local_train = np.searchsorted(keep, np.intersect1d(keep, train))
            local_test = np.searchsorted(keep, np.intersect1d(keep, test))
            if local_train.size < 3 or local_test.size == 0:
                continue

            clinical = H2A.baseline_matrix(data, keep, local_train)
            target = None
            if spec["enzyme"] in chosen_cyp:
                series = data["values"][chosen_cyp[spec["enzyme"]], keep].copy()
                absent = ~np.isfinite(series)
                if absent.any():
                    observed = series[local_train][np.isfinite(series[local_train])]
                    series[absent] = float(observed.mean()) if observed.size else 0.0
                target = series[:, None]
            transcripts = [gene for gene in spec["regulator_transcripts"]
                           if gene in chosen_scoring]
            transcript_block = None
            if transcripts:
                position = {gene: index for index, gene in enumerate(genes)}
                rows_wanted = [position[gene] for gene in transcripts
                               if gene in position]
                if rows_wanted:
                    transcript_block = block[np.asarray(rows_wanted)][:, keep].T

            def fit(design, label, setting_name):
                predicted, _ = H2A.ridge_prediction(
                    design[local_train], y[local_train], design[local_test],
                    inner_folds=FO.FROZEN["inner_folds"],
                    seed=FO.FROZEN["seed"] + fold_number)
                store = errors[endpoint].setdefault(setting_name, {})
                if label not in store:
                    store[label] = np.full((repeats, y.size), np.nan)
                store[label][repeat, local_test] = np.abs(predicted - y[local_test])
                return predicted

            # Arms without activities are computed once, under "real".
            for label in ("F0", "F1", "F2"):
                if verdict[label] == "NOT_EVALUABLE":
                    continue
                if label == "F0":
                    design = clinical
                elif label == "F1":
                    if target is None:
                        continue
                    design = np.column_stack([clinical, target])
                else:
                    if target is None:
                        continue
                    pieces = [clinical, target]
                    if transcript_block is not None:
                        pieces.append(transcript_block)
                    design = np.column_stack(pieces)
                fit(design, label, "real")

            selected = [r for r in spec["evaluable_regulators"]
                        if r in union_position]
            if not selected:
                continue
            columns = [union_position[r] for r in selected]
            for setting in settings:
                scores = activity_cache[setting["name"]][keep][:, columns]
                names = list(selected)
                if not np.all(np.isfinite(scores)):
                    finite = np.isfinite(scores).all(axis=0)
                    scores = scores[:, finite]
                    names = [name for name, ok in zip(names, finite) if ok]
                    if scores.shape[1] == 0:
                        continue
                if verdict["F4p"] != "NOT_EVALUABLE":
                    fit(np.column_stack([clinical, scores]), "F4p", setting["name"])
                if verdict["F5"] != "NOT_EVALUABLE" and target is not None:
                    design = np.column_stack([clinical, target, scores])
                    predicted = fit(design, "F5", setting["name"])
                    if setting["name"] in ("real", "algebraic_identity_check"):
                        identity_predictions.setdefault(
                            endpoint, {}).setdefault(setting["name"], {})[
                                fold_number] = predicted.copy()
                    if (setting["name"] == "real"
                            and endpoint == FO.PRIMARY_ENDPOINT):
                        coefficients[endpoint][fold_number] = _standardised_weights(
                            design[local_train], y[local_train],
                            ["age", "sex", "target_transcript"] + names)
        fold_seconds.append(round(time.time() - started, 1))
        log(f"  fold {fold_number + 1}/{len(folds)} done "
            f"({fold_seconds[-1]:.0f}s elapsed)")

    return {
        "errors": errors, "endpoint_rows": endpoint_rows,
        "specification": specification, "universes": universes,
        "guards": guards, "network_audit": network_audit, "freeze": freeze,
        "settings": settings, "repeats": repeats, "union": union,
        "coefficients": coefficients,
        "identity_predictions": identity_predictions,
        "chosen_probe_trace": chosen_probe_trace,
        "folds": folds,
        "signed_network_audit": signed.audit,
        "fold_seconds": fold_seconds,
    }


def _standardised_weights(x_train, y_train, names):
    """Ridge weights at the inner-chosen penalty, for the biology block only."""
    predicted, alpha = H2A.ridge_prediction(
        x_train, y_train, x_train, inner_folds=FO.FROZEN["inner_folds"],
        seed=FO.FROZEN["seed"])
    mean = x_train.mean(axis=0)
    spread = x_train.std(axis=0)
    spread = np.where(spread > 0.0, spread, 1.0)
    standardised = (x_train - mean) / spread
    centred = y_train - y_train.mean()
    gram = standardised.T @ standardised
    weights = np.linalg.solve(
        gram + alpha * np.eye(gram.shape[0]), standardised.T @ centred)
    return {"alpha": float(alpha),
            "weights": {str(name): float(value)
                        for name, value in zip(names, weights)}}


# --------------------------------------------------------------------------
# Reporting: three blocks, never merged
# --------------------------------------------------------------------------
def per_donor(block):
    """One absolute error per donor, averaged over repeats. A repeat is not a donor.

    A donor absent from every repeat of an arm produces no error, which is a
    legitimate state (an endpoint's donor set differs from the expression
    donors) and not a numerical accident, so the all-missing warning is
    suppressed rather than worked around.
    """
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        averaged = np.nanmean(block, axis=0)
    return averaged[np.isfinite(averaged)], np.isfinite(averaged)


def bonferroni_interval(difference, reference, *, seed, family_size):
    """A paired bootstrap interval widened for the family it belongs to."""
    difference = np.asarray(difference, dtype=np.float64)
    reference = np.asarray(reference, dtype=np.float64)
    rng = np.random.default_rng(seed)
    n = difference.size
    draws = np.empty(FO.FROZEN["n_bootstrap"], dtype=np.float64)
    relative = np.empty(FO.FROZEN["n_bootstrap"], dtype=np.float64)
    for index in range(FO.FROZEN["n_bootstrap"]):
        pick = rng.integers(0, n, size=n)
        draws[index] = difference[pick].mean()
        denominator = reference[pick].mean()
        relative[index] = draws[index] / denominator if denominator > 0 else np.nan
    alpha = 0.05 / max(1, int(family_size))
    low, high = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    rel_low, rel_high = np.nanpercentile(
        relative, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return {"alpha": float(alpha), "family_size": int(family_size),
            "ci_low": float(low), "ci_high": float(high),
            "relative_ci_low": float(rel_low), "relative_ci_high": float(rel_high),
            "interval_excludes_zero": bool(low > 0 or high < 0)}


def contrast(errors, endpoint, arm, comparator, y, *, family_size, setting="real"):
    """Paired `arm - comparator`, on the donors both arms predicted."""
    # The arm may come from a control setting; the comparator always comes from
    # the real network, because a control perturbs the arm under test and not
    # the baseline it is measured against. Two separate lookups, therefore --
    # reading the comparator out of `store` would ask the settings dictionary
    # for an arm name, which silently finds nothing.
    store = errors[endpoint].get(setting, {})
    real = errors[endpoint].get("real", {})
    if arm not in store or comparator not in real:
        missing = [name for name, present in
                   ((arm, arm in store), (comparator, comparator in real))
                   if not present]
        return {"label": f"{arm}_vs_{comparator}", "endpoint": endpoint,
                "setting": setting, "status": "NOT_EVALUABLE",
                "reason": f"no prediction for {', '.join(missing)}"}
    a_block, b_block = store[arm], real[comparator]
    a_values, a_mask = per_donor(a_block)
    b_values, b_mask = per_donor(b_block)
    shared = a_mask & b_mask
    a_full = np.nanmean(a_block, axis=0)[shared]
    b_full = np.nanmean(b_block, axis=0)[shared]
    predicted_y = np.asarray(y)[shared]
    difference = a_full - b_full
    base = float(b_full.mean())
    interval = bonferroni_interval(difference, b_full,
                                   seed=FO.FROZEN["seed"], family_size=family_size)
    arm_loss = {"mae": float(a_full.mean()),
                "trimmed_mae": H2A.trimmed_mean(a_full, FO.FROZEN["trim_fraction"])}
    comparator_loss = {
        "mae": base,
        "trimmed_mae": H2A.trimmed_mean(b_full, FO.FROZEN["trim_fraction"])}
    relative = (float(difference.mean()) / base) if base > 0 else None
    out = {
        "label": f"{arm}_vs_{comparator}", "endpoint": endpoint, "setting": setting,
        "n_donors": int(shared.sum()),
        "arm": arm_loss, "comparator": comparator_loss,
        "delta_mae": float(difference.mean()), "relative": relative,
        "meets_relative_threshold": bool(
            relative is not None
            and relative <= -FO.FROZEN["relative_mae_threshold"]),
        "status_ceiling": FO.FROZEN["status_ceiling"],
    }
    out.update(interval)
    reversal_trimmed = (arm_loss["trimmed_mae"] > comparator_loss["trimmed_mae"])
    out["reversal_on_trimmed_mae"] = bool(reversal_trimmed)
    out["how_to_read"] = (
        "delta_mae negative means the first arm has the smaller error. "
        "interval_excludes_zero false means no difference was DETECTED, never "
        "that there is none. Success needs all three of: relative reduction at "
        "least the frozen threshold, interval excluding zero, and no reversal on "
        "both robustness co-criteria.")
    out["predicted_endpoint_mean"] = float(predicted_y.mean())
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", required=True)
    parser.add_argument("--network", required=True)
    parser.add_argument("--identity", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--control-replicates", type=int,
                        default=FO.FROZEN["control_replicates"])
    parser.add_argument("--dry-run", action="store_true",
                        help=("validate loading, the cytochrome guards, the sign "
                              "rules and the evaluability verdicts, then stop "
                              "before fitting anything"))
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    def log(message):
        if not args.quiet:
            print(message, flush=True)

    started = time.time()
    state = run(args, log)
    if state.get("dry_run"):
        log("\nno contrast was computed, nothing was written")
        return 0
    errors, specification = state["errors"], state["specification"]

    # --- block 1: predictive performance -------------------------------
    arm_losses = {}
    for endpoint in ENDPOINTS:
        y = state["endpoint_rows"][endpoint]["y"]
        arm_losses[endpoint] = {}
        for arm, block in errors[endpoint].get("real", {}).items():
            values, mask = per_donor(block)
            predictions_y = y[mask]
            arm_losses[endpoint][arm] = {
                "mae": float(values.mean()),
                "trimmed_mae": H2A.trimmed_mean(values, FO.FROZEN["trim_fraction"]),
                "n_donors": int(values.size),
                "endpoint_mean": float(predictions_y.mean()),
            }

    primary = contrast(errors, FO.PRIMARY_ENDPOINT, FO.PRIMARY_ARM,
                       FO.PRIMARY_COMPARATOR,
                       state["endpoint_rows"][FO.PRIMARY_ENDPOINT]["y"],
                       family_size=1)
    key_secondary = [
        contrast(errors, FO.PRIMARY_ENDPOINT, arm, comparator,
                 state["endpoint_rows"][FO.PRIMARY_ENDPOINT]["y"],
                 family_size=len(FO.KEY_SECONDARY_CONTRASTS))
        for arm, comparator in FO.KEY_SECONDARY_CONTRASTS]

    other_endpoints = [e for e in ENDPOINTS
                       if e not in (FO.PRIMARY_ENDPOINT, FO.SAME_ENZYME_REPLICATE)]
    other_secondary = []
    for arm, comparator in (("F5", "F1"),) + FO.KEY_SECONDARY_CONTRASTS:
        family = [e for e in other_endpoints
                  if specification[e]["evaluability"]["verdict"][arm]
                  != "NOT_EVALUABLE"]
        for endpoint in family:
            other_secondary.append(contrast(
                errors, endpoint, arm, comparator,
                state["endpoint_rows"][endpoint]["y"], family_size=len(family)))
    replicate_contrast = contrast(
        errors, FO.SAME_ENZYME_REPLICATE, FO.PRIMARY_ARM, FO.PRIMARY_COMPARATOR,
        state["endpoint_rows"][FO.SAME_ENZYME_REPLICATE]["y"], family_size=1)
    replicate_contrast["independence"] = (
        "NOT an independent validation: the same microsomes from the same donors, "
        "a second reaction of the same enzyme")

    # --- the controls, summarised by kind ------------------------------
    control_summary = {}
    for endpoint in ENDPOINTS:
        per_kind: dict[str, dict] = {}
        for setting in state["settings"]:
            name = setting["name"]
            if name == "real":
                continue
            store = errors[endpoint].get(name, {})
            for arm in ("F4p", "F5"):
                if arm not in store:
                    continue
                values, _ = per_donor(store[arm])
                key = (f"{setting['kind']}"
                       + (f"_{setting['fraction']:.2f}"
                          if setting["fraction"] is not None else ""))
                if name == "algebraic_identity_check":
                    key = "algebraic_identity_check"
                bucket = per_kind.setdefault(key, {}).setdefault(arm, [])
                bucket.append(float(values.mean()))
        summary = {}
        for key, arms in per_kind.items():
            summary[key] = {
                arm: {"replicates": len(values),
                      "mae_mean": float(np.mean(values)),
                      "mae_min": float(np.min(values)),
                      "mae_max": float(np.max(values))}
                for arm, values in arms.items()}
        if summary:
            control_summary[endpoint] = summary

    identity_check = {}
    for endpoint, settings_seen in state["identity_predictions"].items():
        real = settings_seen.get("real", {})
        flipped = settings_seen.get("algebraic_identity_check", {})
        shared = sorted(set(real) & set(flipped))
        if not shared:
            continue
        worst = max(float(np.max(np.abs(real[f] - flipped[f]))) for f in shared)
        identity_check[endpoint] = {
            "folds_compared": len(shared),
            "max_absolute_prediction_difference": worst,
            "tolerance": FO.FROZEN["identity_check_tolerance"],
            "identity_holds": bool(worst <= FO.FROZEN["identity_check_tolerance"]),
            "what_it_means": (
                "A full sign inversion negates every activity and ridge answers "
                "by negating one coefficient, so the predictions must be "
                "identical. This is an algebraic identity, not a control, and it "
                "is reported here to show the implementation is the model the "
                "prespecification describes."),
        }

    # --- block 2: biological relevance ---------------------------------
    sign_agreement = {}
    fold_weights = state["coefficients"].get(FO.PRIMARY_ENDPOINT, {})
    if fold_weights:
        recorded = specification[FO.PRIMARY_ENDPOINT]["recorded_edge_sign"]
        agree = disagree = absent = 0
        per_regulator: dict[str, list[float]] = {}
        for fold in sorted(fold_weights):
            for name, weight in fold_weights[fold]["weights"].items():
                if name in ("age", "sex", "target_transcript"):
                    continue
                per_regulator.setdefault(name, []).append(weight)
        for regulator, weights in per_regulator.items():
            mean_weight = float(np.mean(weights))
            stated = recorded.get(regulator)
            if stated is None:
                absent += 1
            elif np.sign(mean_weight) == np.sign(stated):
                agree += 1
            else:
                disagree += 1
        total = agree + disagree
        sign_agreement = {
            "endpoint": FO.PRIMARY_ENDPOINT,
            "regulators_with_a_recorded_sign": total,
            "coefficient_sign_agrees": agree,
            "coefficient_sign_disagrees": disagree,
            "no_recorded_sign": absent,
            "agreement_fraction": (agree / total) if total else None,
            "binomial_p_two_sided": _binomial_p(agree, total) if total else None,
            "what_it_is_not": (
                "Sign agreement is a statement about direction, not about "
                "prediction. It never upgrades the predictive block, and a model "
                "can agree with every recorded sign while predicting no better "
                "than the transcript alone."),
        }

    # --- block 3: analytical validity ----------------------------------
    not_evaluable = {
        endpoint: {"verdict": specification[endpoint]["evaluability"]["verdict"],
                   "reasons": specification[endpoint]["evaluability"]["reasons"]}
        for endpoint in ENDPOINTS
        if any(v == "NOT_EVALUABLE"
               for v in specification[endpoint]["evaluability"]["verdict"].values())
    }
    fold_fingerprint = hashlib.sha256(
        json.dumps([[int(i) for i in test] for _, test in state["folds"]],
                   sort_keys=True).encode("utf-8")).hexdigest()

    report = {
        "report": "F-series: signed regulatory activities versus the target transcript",
        "prespecification": FO.FROZEN["prespecification"],
        "prespecification_version": FO.FROZEN["prespecification_version"],
        "prespecification_sha256": state["freeze"]["sha256"],
        "built_on": date.today().isoformat(),
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "runtime_seconds": round(time.time() - started, 1),
        # Written as a mapping, not a list. The repository's committed-artefact
        # guard refuses any numeric array of ten or more entries, as a crude
        # proxy for participant-level data; ten fold timings would trip it. The
        # guard is deliberately crude and is worth more than the convenience of
        # a list, so the shape changes here rather than the guard there.
        "fold_cumulative_seconds": {
            f"fold_{index + 1:02d}": value
            for index, value in enumerate(state.get("fold_seconds", []))},
        "how_to_read_runtime": (
            "runtime_seconds is wall clock and therefore includes any time the "
            "process spent suspended, which on a laptop can exceed the "
            "computation by two orders of magnitude. fold_cumulative_seconds is "
            "the cumulative time at the end of each outer fold and is the figure "
            "to read for cost: a fold that took far longer than its neighbours "
            "did not do more work. Neither number affects any result -- the run "
            "is fully determined by the frozen seeds."),
        "status_ceiling": FO.FROZEN["status_ceiling"],
        "status_ceiling_reason": (
            "The Human Liver Cohort has been inspected repeatedly by this project "
            "(H2a, H2d, the probe-pool construction, the destroyed-signal "
            "diagnostic). Every F-series result is exploratory whatever its "
            "direction."),
        "endpoint_class": "functional_in_vitro",
        "clinical_pk_transfer": "NOT_EVALUATED",
        "clinical_pk_transfer_reason": (
            "This cohort has no dose, no concentration and no time. Nothing here "
            "was scaled to a clearance or an exposure, and no simulated "
            "concentration appears in this artefact. Closing the chain needs "
            "individual transcriptome and individual concentration-time "
            "pharmacokinetics in the same publicly downloadable subjects, which "
            "the project's registry does not contain."),
        "frozen": {key: (list(value) if isinstance(value, tuple) else value)
                   for key, value in FO.FROZEN.items()},
        "arms": {
            "F0": "age and sex",
            "F1": "F0 + the target enzyme's own transcript (primary comparator)",
            "F2": "F1 + the evaluable regulators' own transcripts",
            "F4p": "F0 + regulator activities from cytochrome-free regulons",
            "F5": "F1 + regulator activities (primary arm)",
        },
        "reserved_unimplemented_arms": FO.RESERVED_UNIMPLEMENTED_ARMS,
        "predictive_performance": {
            "primary_contrast": primary,
            "key_secondary_contrasts": key_secondary,
            "same_enzyme_replicate": replicate_contrast,
            "other_secondary_contrasts": other_secondary,
            "arm_losses": arm_losses,
            "controls": control_summary,
            "control_reading": (
                "A gain the rewired null reproduces is refuted by its own "
                "control. Flat performance across sign_perturbation fractions "
                "means any gain came from the gene list, not from the signs."),
        },
        "biological_relevance": {
            "coefficient_sign_agreement": sign_agreement,
            "regulon_sizes": {
                endpoint: specification[endpoint]["regulon_sizes"]
                for endpoint in (FO.PRIMARY_ENDPOINT,)},
            "candidate_and_evaluable_counts": {
                endpoint: {
                    "enzyme": specification[endpoint]["enzyme"],
                    "candidate_regulators": len(
                        specification[endpoint]["candidate_regulators"]),
                    "evaluable_regulators": len(
                        specification[endpoint]["evaluable_regulators"]),
                    "regulator_transcripts_for_F2": len(
                        specification[endpoint]["regulator_transcripts"]),
                }
                for endpoint in ENDPOINTS},
            "refused_regulators": {
                endpoint: specification[endpoint]["refused_regulators"]
                for endpoint in ENDPOINTS
                if specification[endpoint]["refused_regulators"]},
        },
        "analytical_validity": {
            "network": state["network_audit"],
            "signed_network_audit": state["signed_network_audit"],
            "probe_pool": state["universes"]["pool"],
            "cytochrome_guards": state["guards"],
            "not_evaluable": not_evaluable,
            "donors_with_endpoint_and_expression": {
                endpoint: specification[endpoint][
                    "donors_with_endpoint_and_expression"]
                for endpoint in ENDPOINTS},
            "fold_fingerprint_sha256": fold_fingerprint,
            "repeats": state["repeats"],
            "union_of_evaluable_regulators": len(state["union"]),
            "activity_settings_per_fold": len(state["settings"]),
            "probe_choice_trace_fold_0": state["chosen_probe_trace"],
            "algebraic_identity_check": identity_check,
            "no_individual_data_written": True,
        },
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    log(f"\nwritten {out}")
    log(f"primary {primary.get('label')}: delta {primary.get('delta_mae')!r} "
        f"relative {primary.get('relative')!r} "
        f"excludes zero {primary.get('interval_excludes_zero')!r}")
    return 0


def _binomial_p(successes, trials):
    """Two-sided exact binomial p against a half, without SciPy."""
    from math import comb

    if trials <= 0:
        return None
    probabilities = [comb(trials, k) * 0.5 ** trials for k in range(trials + 1)]
    observed = probabilities[successes]
    return float(sum(p for p in probabilities if p <= observed + 1e-15))


if __name__ == "__main__":
    raise SystemExit(main())
