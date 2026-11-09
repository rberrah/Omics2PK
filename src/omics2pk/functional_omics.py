"""Signed regulatory activities as pharmacokinetic-adjacent features (F-series).

Frozen by `docs/OMICS2PK_F_SERIES_FUNCTIONAL_OMICS_PRESPECIFICATION_20261109.md`,
version 1.0.0. Every constant in `FROZEN` is read back from that document by
`tests/test_functional_omics.py`; a value that moves in one place and not the
other fails the suite.

The module does four things and refuses to do a fifth:

1. turns a CollecTRI download into signed regulons, with the sign rule, the
   duplicate rule and the cytochrome removal all explicit and audited;
2. computes the published univariate-linear-model activity score, implemented
   from its closed form rather than imported, and standardised inside the
   training fold only;
3. builds the three controls the prespecification fixes -- partial sign
   perturbation, size-and-sign-preserving rewiring, donor permutation -- and
   labels the full sign inversion as the algebraic identity it is;
4. reports `NOT_EVALUABLE` with a reason instead of returning a number it
   cannot support.

The fifth thing, which it does not do: nothing here converts an activity into a
clearance, a concentration or an exposure. The cohort this feeds has no
clinical pharmacokinetics, so there is no function to call.

**Circularity.** Each of the nine study enzymes is a member of its own
candidate regulators' regulons -- CYP3A4 of 51, CYP2D6 of 6. Scoring a
transcription factor from a regulon that still contains the predicted enzyme
and then predicting that enzyme is predicting a variable from itself. The
enzymes are therefore removed at the edge level, at the probe level and at the
resolved-alias level, and `cytochrome_guard_audit` reports what each removal
took.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

# --------------------------------------------------------------------------
# Frozen parameters -- prespecification section 12
# --------------------------------------------------------------------------
FROZEN = {
    "prespecification": (
        "docs/OMICS2PK_F_SERIES_FUNCTIONAL_OMICS_PRESPECIFICATION_20261109.md"),
    "prespecification_version": "1.0.0",
    "seed": 20261109,
    "outer_folds": 5,
    "outer_repeats": 2,
    "inner_folds": 3,
    "n_bootstrap": 2000,
    "ridge_alphas": (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0),
    "relative_mae_threshold": 0.10,
    "trim_fraction": 0.10,
    "regulon_floor": 5,
    "regulator_floor": 5,
    "donor_floor": 50,
    "max_probe_missing_fraction": 0.20,
    "control_replicates": 20,
    "perturbation_fractions": (0.25, 0.50, 0.75),
    "identity_check_fraction": 1.00,
    "identity_check_tolerance": 1e-8,
    "status_ceiling": "exploratory",
}

# The arms. `F4p` is the prespecification's `F4'`; a prime does not survive
# every tool that reads these artefacts, so the code spells it `p`.
ARMS = ("F0", "F1", "F2", "F4p", "F5")
PRIMARY_ARM = "F5"
PRIMARY_COMPARATOR = "F1"
PRIMARY_ENDPOINT = "CYP3A4M_activity"
SAME_ENZYME_REPLICATE = "CYP3A4T_activity"
KEY_SECONDARY_CONTRASTS = (("F4p", "F1"), ("F5", "F2"))

# F3 is reserved and deliberately absent: PROGENy's fourteen footprints contain
# no xenobiotic-metabolism pathway, so the arm would carry no column with a
# stated mechanistic relation to a cytochrome.
RESERVED_UNIMPLEMENTED_ARMS = {
    "F3": ("pathway footprints (PROGENy): none of its 14 footprints is a "
           "xenobiotic-metabolism pathway, so no column would have a stated "
           "mechanistic relation to a cytochrome"),
}

ENDPOINT_ENZYME = {
    "CYP1A2_activity": "CYP1A2",
    "CYP2A6_activity": "CYP2A6",
    "CYP2B6_activity": "CYP2B6",
    "CYP2C8_activity": "CYP2C8",
    "CYP2C9_activity": "CYP2C9",
    "CYP2C19_activity": "CYP2C19",
    "CYP2D6_activity": "CYP2D6",
    "CYP2E1_activity": "CYP2E1",
    "CYP3A4M_activity": "CYP3A4",
    "CYP3A4T_activity": "CYP3A4",
}

# Matches a cytochrome P450 gene symbol, including the pseudogene suffix, so
# that a paralogue or a pseudogene of the predicted enzyme cannot re-enter a
# regulon under a name the narrow list would miss.
CYP_PATTERN = re.compile(r"^CYP[0-9]+[A-Z]*[0-9]*(P[0-9]+)?$", re.IGNORECASE)

# CollecTRI writes a protein complex as its members joined by an underscore.
# A human gene symbol does not contain one, so the split is unambiguous.
COMPLEX_SEPARATOR = "_"

STIMULATION = 1
INHIBITION = -1


class FunctionalOmicsError(RuntimeError):
    """Raised instead of returning a number the data does not support."""


def is_cytochrome(symbol) -> bool:
    return bool(symbol) and bool(CYP_PATTERN.match(str(symbol)))


def complex_members(regulator: str) -> list[str]:
    """The named members of a complex regulator, or the regulator itself."""
    text = str(regulator)
    if COMPLEX_SEPARATOR not in text:
        return [text]
    return [part for part in text.split(COMPLEX_SEPARATOR) if part]


# --------------------------------------------------------------------------
# Network: sign resolution, duplicates, cytochrome removal
# --------------------------------------------------------------------------
@dataclass
class SignedNetwork:
    """Resolved signed edges, plus the audit of every row that did not survive."""

    edges: dict[tuple[str, str], int]
    audit: dict = field(default_factory=dict)

    def regulons(self) -> dict[str, dict[str, int]]:
        out: dict[str, dict[str, int]] = {}
        for (source, target), weight in self.edges.items():
            out.setdefault(source, {})[target] = weight
        return out

    def regulators_of(self, target: str) -> list[str]:
        return sorted({source for (source, candidate) in self.edges
                       if candidate == target})


def _boolean_column(frame, name):
    if name not in frame.columns:
        raise FunctionalOmicsError(f"the network has no column named {name!r}")
    return frame[name].astype(str).str.strip().str.lower().eq("true").to_numpy()


def resolve_signed_network(frame, *, drop_cytochrome_targets=True) -> SignedNetwork:
    """Turn a downloaded CollecTRI table into one signed weight per edge.

    Three rules, all frozen in prespecification section 4 and all audited here
    rather than applied silently:

    * **sign** -- the consensus columns decide. `+1` for stimulation alone,
      `-1` for inhibition alone, and the edge *leaves the regulon* when the
      consensus says both or neither. It is not kept with weight zero, because
      a zero-weight member still counts towards the regulon size and towards
      the evaluability floor.
    * **duplicates** -- collapsed to one edge; dropped outright when two rows
      for the same pair resolve to opposite signs.
    * **cytochromes** -- every edge whose target is any CYP gene is removed
      before a regulon exists, so that no score can be computed from the
      variable it is meant to predict.
    """
    required = ("source_genesymbol", "target_genesymbol", "consensus_stimulation",
                "consensus_inhibition")
    for column in required:
        if column not in frame.columns:
            raise FunctionalOmicsError(f"the network has no column named {column!r}")

    source = frame["source_genesymbol"].astype(str).to_numpy()
    target = frame["target_genesymbol"].astype(str).to_numpy()
    stimulation = _boolean_column(frame, "consensus_stimulation")
    inhibition = _boolean_column(frame, "consensus_inhibition")

    rows_total = int(len(frame))
    dropped_sign_indeterminate = 0
    dropped_cytochrome_target = 0
    dropped_conflicting_duplicate = 0
    collapsed_duplicate_rows = 0

    resolved: dict[tuple[str, str], int] = {}
    conflicting: set[tuple[str, str]] = set()
    for index in range(rows_total):
        up, down = bool(stimulation[index]), bool(inhibition[index])
        if up == down:
            dropped_sign_indeterminate += 1
            continue
        if drop_cytochrome_targets and is_cytochrome(target[index]):
            dropped_cytochrome_target += 1
            continue
        key = (source[index], target[index])
        weight = STIMULATION if up else INHIBITION
        if key in resolved:
            collapsed_duplicate_rows += 1
            if resolved[key] != weight:
                conflicting.add(key)
        else:
            resolved[key] = weight

    for key in conflicting:
        resolved.pop(key, None)
        dropped_conflicting_duplicate += 1

    if not resolved:
        raise FunctionalOmicsError("sign resolution left no edge at all")

    return SignedNetwork(
        edges=resolved,
        audit={
            "rows_read": rows_total,
            "edges_retained": len(resolved),
            "dropped_sign_indeterminate": dropped_sign_indeterminate,
            "dropped_cytochrome_target": dropped_cytochrome_target,
            "duplicate_rows_collapsed": collapsed_duplicate_rows,
            "dropped_conflicting_duplicate": dropped_conflicting_duplicate,
            "cytochrome_targets_removed": bool(drop_cytochrome_targets),
            "how_to_read": (
                "dropped_cytochrome_target is the circularity guard at the edge "
                "level; it must be positive on the real network, because every "
                "study enzyme has recorded regulators."),
        },
    )


def candidate_regulators(frame, enzyme: str) -> list[str]:
    """The transcription factors whose CollecTRI regulon contains `enzyme`.

    Read from the network **before** the cytochrome removal, which is the whole
    point: the enzyme's identity selects the regulators, and the enzyme is then
    removed from every regulon used for scoring. Choosing the regulators is a
    mechanistic prior; keeping the enzyme inside them would be circular. The
    two steps are separate and only the first uses the enzyme's name.
    """
    if not enzyme:
        raise FunctionalOmicsError("an enzyme symbol is required")
    unsigned = resolve_signed_network(frame, drop_cytochrome_targets=False)
    return unsigned.regulators_of(enzyme)


def recorded_edge_sign(frame, regulator: str, enzyme: str):
    """The sign CollecTRI records from `regulator` to `enzyme`, or None.

    Kept for the biological-relevance block alone. This edge is removed from
    every regulon used for scoring; it is read back only to ask whether a
    fitted coefficient points the way the resource says it should, which is a
    question about biology and not about prediction.
    """
    unsigned = resolve_signed_network(frame, drop_cytochrome_targets=False)
    return unsigned.edges.get((str(regulator), str(enzyme)))


# --------------------------------------------------------------------------
# Evaluability
# --------------------------------------------------------------------------
def evaluable_regulators(regulons, regulators, admissible_symbols, *,
                         floor=None):
    """Which candidate regulators have enough admissible, cytochrome-free members.

    Returns the evaluable regulators in a fixed order and, for every regulator
    refused, the reason and the counts behind it. A regulator is never dropped
    silently: `NOT_EVALUABLE` with a reason is the reportable outcome.
    """
    floor = FROZEN["regulon_floor"] if floor is None else int(floor)
    admissible = set(map(str, admissible_symbols))
    kept, refused = [], {}
    for regulator in sorted(map(str, regulators)):
        regulon = regulons.get(regulator)
        if not regulon:
            refused[regulator] = {
                "reason": "no regulon survived the sign and cytochrome rules",
                "members": 0, "admissible_members": 0}
            continue
        members = [gene for gene in regulon if not is_cytochrome(gene)]
        usable = [gene for gene in members if gene in admissible]
        if len(usable) < floor:
            refused[regulator] = {
                "reason": f"fewer than {floor} admissible members",
                "members": len(members), "admissible_members": len(usable)}
            continue
        kept.append(regulator)
    return kept, refused


def endpoint_evaluability(*, endpoint, target_admissible_probes,
                          evaluable_regulator_count, donors_with_both):
    """Per-arm evaluability for one endpoint, with the reason attached.

    The three floors are independent, so an endpoint can be evaluable for the
    functional arm and not for the transcript arms. `CYP1A2_activity` is
    exactly that case on this cohort: none of its four probes is admissible, so
    `F1`, `F2` and `F5` have no target column, while `F4p` needs none.
    """
    reasons: dict[str, str] = {}
    verdict = {arm: "EVALUABLE" for arm in ARMS}

    if donors_with_both < FROZEN["donor_floor"]:
        for arm in ARMS:
            verdict[arm] = "NOT_EVALUABLE"
        reasons["all_arms"] = (
            f"{donors_with_both} donors carry both the endpoint and expression, "
            f"below the frozen floor of {FROZEN['donor_floor']}")
        return {"endpoint": endpoint, "verdict": verdict, "reasons": reasons}

    if int(target_admissible_probes) <= 0:
        for arm in ("F1", "F2", "F5"):
            verdict[arm] = "NOT_EVALUABLE"
        reasons["F1/F2/F5"] = (
            "the target enzyme has no admissible probe, so no arm carrying its "
            "transcript can be built; this is not a negative result")
    if int(evaluable_regulator_count) < FROZEN["regulator_floor"]:
        for arm in ("F4p", "F5"):
            verdict[arm] = "NOT_EVALUABLE"
        reasons["F4p/F5"] = (
            f"{evaluable_regulator_count} evaluable regulators, below the frozen "
            f"floor of {FROZEN['regulator_floor']}")
    return {"endpoint": endpoint, "verdict": verdict, "reasons": reasons}


# --------------------------------------------------------------------------
# The univariate linear model, written from its closed form
# --------------------------------------------------------------------------
def standardise_on_training(values, train):
    """Gene-wise z-score using the training donors' mean and spread only.

    The transformation is unsupervised -- it never sees the endpoint -- but it
    is still fitted inside the fold. A gene with no spread in the training
    donors becomes a constant zero rather than a division by zero.
    """
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 2:
        raise FunctionalOmicsError("expression must be genes by donors")
    train = np.asarray(train, dtype=np.int64)
    if train.size < 2:
        raise FunctionalOmicsError("standardisation needs at least two training donors")
    centre = values[:, train].mean(axis=1, keepdims=True)
    spread = values[:, train].std(axis=1, keepdims=True)
    safe = np.where(spread > 0.0, spread, 1.0)
    out = (values - centre) / safe
    out[np.repeat(spread <= 0.0, values.shape[1], axis=1)] = 0.0
    return out


def weight_matrix(regulons, regulators, genes):
    """The signed weight of every gene in every regulon, as a genes-by-TF matrix."""
    position = {str(gene): index for index, gene in enumerate(genes)}
    matrix = np.zeros((len(position), len(regulators)), dtype=np.float64)
    for column, regulator in enumerate(regulators):
        for gene, weight in regulons.get(str(regulator), {}).items():
            index = position.get(str(gene))
            if index is not None:
                matrix[index, column] = float(weight)
    return matrix


def ulm_activities(standardised, weights):
    """decoupleR's `run_ulm` statistic, vectorised over donors and regulators.

    For donor `j` and regulator `e` the score is the t-statistic of the slope in
    the ordinary least squares fit of that donor's standardised profile on the
    regulon's weight vector, over the whole scoring universe with zeros outside
    the regulon:

        beta = (w - mean w) . x_j / sum (w - mean w)^2
        RSS  = sum (x_j - mean x_j)^2 - beta^2 sum (w - mean w)^2
        t    = beta / sqrt( (RSS / (G - 2)) / sum (w - mean w)^2 )

    Returned as donors by regulators. A regulator whose weight vector has no
    spread over the universe gets `nan`, never a silent zero: no score is not
    the same statement as a score of zero.
    """
    standardised = np.asarray(standardised, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if standardised.ndim != 2 or weights.ndim != 2:
        raise FunctionalOmicsError("both inputs must be two-dimensional")
    genes = standardised.shape[0]
    if weights.shape[0] != genes:
        raise FunctionalOmicsError(
            f"weights cover {weights.shape[0]} genes, expression covers {genes}")
    if genes < 3:
        raise FunctionalOmicsError(
            "the t-statistic needs at least three genes to have a residual degree "
            "of freedom")

    centred_weights = weights - weights.mean(axis=0, keepdims=True)
    sum_squares = (centred_weights ** 2).sum(axis=0)
    degenerate = sum_squares <= 0.0
    safe = np.where(degenerate, 1.0, sum_squares)

    centred_profile = standardised - standardised.mean(axis=0, keepdims=True)
    slope = (centred_weights.T @ centred_profile) / safe[:, None]
    total = (centred_profile ** 2).sum(axis=0)
    residual = total[None, :] - (slope ** 2) * safe[:, None]
    residual = np.clip(residual, 1e-300, None)
    standard_error = np.sqrt((residual / (genes - 2)) / safe[:, None])
    scores = slope / standard_error
    scores[degenerate, :] = np.nan
    return scores.T


# --------------------------------------------------------------------------
# Controls
# --------------------------------------------------------------------------
def perturb_regulon_signs(regulons, fraction, rng):
    """Flip the sign of a random `fraction` of each regulon's edges.

    A **partial** flip is the control that bites. A full flip is not a control
    at all: under this statistic, replacing `w` by `-w` leaves the residuals
    untouched and negates the slope, so every activity is exactly negated and a
    ridge on standardised columns answers by negating one coefficient. The
    predictions are then identical, which is an algebraic identity and not a
    biological finding. `fraction=1.0` is accepted only so the identity can be
    *checked*; the prespecification labels that run
    `algebraic_identity_check`, never `control`.

    A partial flip changes the relative weighting inside the regulon, so the
    score genuinely changes and is no longer the published statistic. That is
    what makes it informative.
    """
    fraction = float(fraction)
    if not 0.0 <= fraction <= 1.0:
        raise FunctionalOmicsError("the perturbation fraction must lie in [0, 1]")
    out: dict[str, dict[str, int]] = {}
    for regulator, regulon in regulons.items():
        genes = sorted(regulon)
        count = int(round(fraction * len(genes)))
        flipped = set()
        if count > 0:
            picks = rng.choice(len(genes), size=count, replace=False)
            flipped = {genes[int(p)] for p in picks}
        out[regulator] = {gene: (-weight if gene in flipped else weight)
                          for gene, weight in regulon.items()}
    return out


def rewire_regulons(regulons, universe, rng):
    """A null network with the same shape and the same signed composition.

    Each regulon keeps its exact size and its exact count of activating and
    repressing edges; only the identity of the members changes, resampled
    uniformly from the scoring universe. That isolates the one claim worth
    testing -- that the network knows *which* genes report on a regulator --
    from the much weaker fact that a weighted average of two hundred arbitrary
    transcripts can track liver enzyme activity.
    """
    pool = [str(gene) for gene in universe if not is_cytochrome(gene)]
    if not pool:
        raise FunctionalOmicsError("the rewiring universe is empty")
    out: dict[str, dict[str, int]] = {}
    for regulator, regulon in regulons.items():
        size = len(regulon)
        if size > len(pool):
            raise FunctionalOmicsError(
                f"regulon {regulator!r} has {size} members but the universe has "
                f"{len(pool)} genes")
        positives = sum(1 for weight in regulon.values() if weight > 0)
        picks = rng.choice(len(pool), size=size, replace=False)
        genes = [pool[int(p)] for p in picks]
        out[regulator] = {gene: (STIMULATION if index < positives else INHIBITION)
                          for index, gene in enumerate(genes)}
    return out


def permute_donors(values, rng):
    """Destroy the expression-to-endpoint link, keeping network and endpoint intact."""
    values = np.asarray(values, dtype=np.float64)
    order = rng.permutation(values.shape[1])
    return values[:, order], order


# --------------------------------------------------------------------------
# The cytochrome guards, audited rather than asserted
# --------------------------------------------------------------------------
def cytochrome_guard_audit(*, array_symbols, feature_ids, resolved_symbols,
                           admissible_features):
    """What each of the three cytochrome guards removed, counted on real data.

    Guard 1 lives in `resolve_signed_network` (edge level). Guards 2 and 3 are
    computed here: a probe is refused when its array symbol is a cytochrome,
    and refused again when the gene it *resolves* to is one even though its
    array symbol is not. The third guard is the one a narrow symbol list
    misses; on this array it catches `FLJ20359`, `FLJ39501` and `FLJ16008`.

    The audit reports the number found. It does not assert that the number is
    zero, because a later annotation release could make it non-zero.
    """
    array_symbols = [str(symbol or "") for symbol in array_symbols]
    feature_ids = [str(identifier) for identifier in feature_ids]
    if len(array_symbols) != len(feature_ids):
        raise FunctionalOmicsError("symbols and feature identifiers differ in length")
    admissible = set(map(str, admissible_features))

    by_symbol, by_alias = [], []
    for symbol, identifier in zip(array_symbols, feature_ids):
        if identifier not in admissible:
            continue
        resolved = str(resolved_symbols.get(identifier, "") or "")
        if is_cytochrome(symbol):
            by_symbol.append(identifier)
        elif is_cytochrome(resolved):
            by_alias.append({"feature_id": identifier, "array_symbol": symbol,
                             "resolved_symbol": resolved})
    return {
        "refused_by_array_symbol": len(by_symbol),
        "refused_by_resolved_alias": len(by_alias),
        "resolved_alias_examples": by_alias[:10],
        "refused_feature_ids": sorted(set(by_symbol)
                                      | {row["feature_id"] for row in by_alias}),
        "how_to_read": (
            "refused_by_resolved_alias counts probes whose array symbol is not a "
            "cytochrome but whose resolved gene is. Excluding by array symbol "
            "alone would leave those probes in the regulons."),
    }


def assert_no_cytochrome_reaches_scoring(symbols, resolved_symbols=None):
    """Stop the run rather than score a cytochrome into its own predictor."""
    offenders = []
    for symbol in symbols:
        text = str(symbol or "")
        if is_cytochrome(text):
            offenders.append(text)
        elif resolved_symbols is not None:
            resolved = str(resolved_symbols.get(text, "") or "")
            if is_cytochrome(resolved):
                offenders.append(f"{text}->{resolved}")
    if offenders:
        raise FunctionalOmicsError(
            f"{len(offenders)} cytochrome measurements reached the scoring "
            f"universe, starting with {sorted(set(offenders))[:5]}; the F-series "
            "scores would then be computed from the variable they predict")
    return True


__all__ = [
    "ARMS", "COMPLEX_SEPARATOR", "CYP_PATTERN", "ENDPOINT_ENZYME", "FROZEN",
    "FunctionalOmicsError", "INHIBITION", "KEY_SECONDARY_CONTRASTS",
    "PRIMARY_ARM", "PRIMARY_COMPARATOR", "PRIMARY_ENDPOINT",
    "RESERVED_UNIMPLEMENTED_ARMS", "SAME_ENZYME_REPLICATE", "STIMULATION",
    "SignedNetwork", "assert_no_cytochrome_reaches_scoring",
    "candidate_regulators", "complex_members", "cytochrome_guard_audit",
    "endpoint_evaluability", "evaluable_regulators", "is_cytochrome",
    "permute_donors", "perturb_regulon_signs", "recorded_edge_sign",
    "resolve_signed_network", "rewire_regulons", "standardise_on_training",
    "ulm_activities", "weight_matrix",
]
