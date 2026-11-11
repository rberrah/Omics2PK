"""Does the presumed mechanism of a variant tell you which measurement to trust?

The question, and why it is not the Fusing question
---------------------------------------------------
`nested_fusion_evaluation` asked whether combining several **predictions** pays.
It does not, on HLC. This module asks a different question, and a negative
answer there does not transfer: given what is known about a variant — where it
sits, what it changes, what the protein's annotated machinery is — can we
choose **which biological measurement carries the information**?

The case that makes the question real: a variant can keep protein abundance
and still destroy activity. Measured here on public deep-mutational-scanning
data, the abundance-activity rank correlation is **+0.76 for CYP2C9 but only
+0.34 for NUDT15**, so the measurement closest to the gene is not reliably the
one that reports its function, and *how* unreliable depends on the protein.

What is predicted, and what must never be
-----------------------------------------
The target is an **activity score from a functional screen**. It is not a
clearance, not an AUC, not a concentration, and not an in-vivo phenotype. A
functional screen is also **not** a multi-omic measurement in individuals:
these are engineered single variants in cell or yeast systems. The relevance
of a *measurement type* and the relevance of a *specific biomarker* are two
questions, and this module addresses only the first.

The two assays also differ in system. For CYP2C9 the abundance assay is
VAMP-seq in HEK293T while the activity assay is Click-seq in humanised yeast,
so a discordance there may be an assay-system artefact rather than biology.
NUDT15's two assays are both in human cells. **This confound is reported with
every CYP2C9 result and is the reason NUDT15 is the primary gene.**

Why the folds are grouped by position
-------------------------------------
Several substitutions at one residue are not independent observations: the
NUDT15 discordant set contains four different substitutions at Arg34. A random
split would train on `p.Arg34Ala` and test on `p.Arg34Gly` and call it
out-of-sample. Every split here is therefore **grouped by protein position**,
and a second, harder split groups **contiguous sequence regions**, because
neighbouring residues share structural context too.

The mechanism rules propose, they do not decide
-----------------------------------------------
`mechanism_hypotheses` returns a ranked list of candidate mechanisms with
prior weights, never a single deterministic label. Nothing downstream is
allowed to treat a hypothesis as a fact, and the rules' usefulness is measured,
not assumed: `M3` (mechanism features alone) is a comparator that can lose.
"""
from __future__ import annotations

import math
import re

import numpy as np

from .fusion_headroom_diagnostics import paired_gap, ridge_predict

#: Frozen before execution.
FROZEN = {
    "seed": 20261112,
    "position_folds": 5,
    "region_blocks": 8,
    "n_bootstrap": 2000,
    "alpha": 0.05,
    "minimum_interesting_relative_gain": 0.05,
    "ridge_penalty": 1.0,
    "loss": "absolute_error",
    "target": "activity_score",
    "grouping": "protein_position, and contiguous region as a harder split",
}

#: Kyte-Doolittle hydropathy.
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
#: Residue volume, cubic angstroms.
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}
CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.5}
POLAR = set("NQSTYCGH")

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}

_HGVS = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter|=|\*)$")


class VariantSelectionError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise VariantSelectionError(message)


def parse_protein_variant(hgvs):
    """`p.Arg34Ala` -> (`R`, 34, `A`). Returns None for anything else.

    Nonsense is returned with `None` as the replacement amino acid and is kept,
    because a premature stop is one of the mechanism classes under test.
    Synonymous (`p.Xxx123=`) is returned with `"="`.
    """
    if not isinstance(hgvs, str):
        return None
    match = _HGVS.match(hgvs.strip())
    if match is None:
        return None
    reference = THREE_TO_ONE.get(match.group(1))
    if reference is None:
        return None
    position = int(match.group(2))
    raw = match.group(3)
    if raw in ("Ter", "*"):
        return reference, position, None
    if raw == "=":
        return reference, position, "="
    return reference, position, THREE_TO_ONE.get(raw)


# --------------------------------------------------------------------------
# the mechanism rules: they rank hypotheses, they do not label
# --------------------------------------------------------------------------
#: Each rule is (name, predicate, prior weight, the measurement it argues for).
#: The weights are **declared priors, not probabilities of predictive
#: usefulness** — the distinction the commissioning note insists on for
#: interaction scores applies here too.
MECHANISM_CLASSES = (
    "nonsense_truncation",
    "functional_site_proximal",
    "motif_internal",
    "buried_core_candidate",
    "surface_tolerant_candidate",
    "terminal_anchor",
)


def mechanism_hypotheses(reference, position, replacement, annotation):
    """A ranked list of candidate mechanisms. Never a single verdict.

    `annotation` carries what public annotation actually provides: the sequence
    length, the annotated functional-site positions and any annotated motif
    span. Where annotation is sparse — CYP2C9 has **one** annotated binding
    residue in UniProt — the rules say so by returning low weights, rather than
    inventing structure.
    """
    length = annotation["length"]
    sites = annotation.get("functional_sites") or ()
    motifs = annotation.get("motifs") or ()
    out = []

    if replacement is None:
        out.append(("nonsense_truncation", 0.90,
                    "transcript decay or a truncated protein: both abundance "
                    "and activity are expected to fall, so neither "
                    "measurement is more informative than the other"))

    nearest = min((abs(position - s) for s in sites), default=None)
    if nearest is not None and nearest <= 8:
        out.append(("functional_site_proximal",
                    0.70 if nearest <= 3 else 0.40,
                    "protein may be present but inactive: ACTIVITY is the "
                    "informative measurement and abundance may mislead"))
    for start, end in motifs:
        if start <= position <= end:
            out.append(("motif_internal", 0.55,
                        "inside an annotated functional motif: activity "
                        "preferred over abundance"))
            break

    if replacement not in (None, "="):
        delta_hydropathy = HYDROPATHY[replacement] - HYDROPATHY[reference]
        delta_volume = VOLUME[replacement] - VOLUME[reference]
        if HYDROPATHY[reference] > 1.5 and (delta_hydropathy < -3.0
                                            or abs(delta_volume) > 50.0):
            out.append(("buried_core_candidate", 0.45,
                        "a hydrophobic residue replaced by something polar or "
                        "much larger may destabilise: ABUNDANCE is the "
                        "informative measurement"))
        if (HYDROPATHY[reference] < -2.0 and abs(delta_volume) < 40.0
                and nearest is not None and nearest > 12):
            out.append(("surface_tolerant_candidate", 0.30,
                        "a conservative change at an exposed, unannotated "
                        "position: little effect expected on either"))

    if position <= 20 or position >= length - 10:
        out.append(("terminal_anchor", 0.35,
                    "a terminus: for a microsomal CYP the N-terminal segment "
                    "is the membrane anchor, so abundance may change without "
                    "catalytic change"))

    if not out:
        out.append(("surface_tolerant_candidate", 0.10,
                    "no rule fired: the annotation available does not support "
                    "any mechanistic hypothesis for this variant"))
    out.sort(key=lambda row: -row[1])
    return [{"mechanism": name, "prior_weight": weight, "argument": why}
            for name, weight, why in out]


# --------------------------------------------------------------------------
# feature blocks, each of which is a comparator's information set
# --------------------------------------------------------------------------
def substitution_features(reference, position, replacement):
    """Physicochemistry of the change. Available for every variant."""
    if replacement in (None, "="):
        base = [0.0] * 7
        return base + [1.0 if replacement is None else 0.0]
    delta_hydropathy = HYDROPATHY[replacement] - HYDROPATHY[reference]
    delta_volume = VOLUME[replacement] - VOLUME[reference]
    charge_from = CHARGE.get(reference, 0.0)
    charge_to = CHARGE.get(replacement, 0.0)
    return [
        delta_hydropathy,
        abs(delta_hydropathy),
        delta_volume / 100.0,
        abs(delta_volume) / 100.0,
        charge_to - charge_from,
        1.0 if replacement == "P" else 0.0,
        1.0 if reference == "G" else 0.0,
        0.0,
    ]


SUBSTITUTION_NAMES = ("d_hydropathy", "abs_d_hydropathy", "d_volume",
                      "abs_d_volume", "d_charge", "proline_introduced",
                      "glycine_lost", "is_nonsense")


def mechanism_features(reference, position, replacement, annotation):
    """Where the variant sits relative to what annotation actually knows."""
    length = annotation["length"]
    sites = annotation.get("functional_sites") or ()
    motifs = annotation.get("motifs") or ()
    nearest = min((abs(position - s) for s in sites), default=float(length))
    in_motif = any(start <= position <= end for start, end in motifs)
    ranked = mechanism_hypotheses(reference, position, replacement, annotation)
    top = {row["mechanism"]: row["prior_weight"] for row in ranked}
    return ([
        position / float(length),
        math.log1p(nearest),
        1.0 if nearest <= 3 else 0.0,
        1.0 if in_motif else 0.0,
        min(position, length - position) / float(length),
    ] + [top.get(name, 0.0) for name in MECHANISM_CLASSES])


MECHANISM_NAMES = (("relative_position", "log1p_distance_to_site",
                    "within_3_of_site", "inside_motif",
                    "distance_to_nearest_terminus")
                   + MECHANISM_CLASSES)


def build_design(variants, annotation):
    """Assemble the per-variant blocks. `variants` is a list of parsed tuples."""
    substitution, mechanism = [], []
    for reference, position, replacement in variants:
        substitution.append(
            substitution_features(reference, position, replacement))
        mechanism.append(
            mechanism_features(reference, position, replacement, annotation))
    return (np.asarray(substitution, dtype=float),
            np.asarray(mechanism, dtype=float))


# --------------------------------------------------------------------------
# grouped splits, because substitutions at one residue are not independent
# --------------------------------------------------------------------------
def position_groups(positions, folds, seed):
    """Assign whole protein POSITIONS to folds, never individual variants."""
    unique = np.unique(np.asarray(positions, dtype=int))
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(unique)
    assignment = {int(p): index % folds
                  for index, p in enumerate(shuffled)}
    return np.asarray([assignment[int(p)] for p in positions], dtype=int)


def region_groups(positions, blocks, length):
    """Assign CONTIGUOUS sequence regions to folds. The harder split.

    Neighbouring residues share structural context, so a position-grouped split
    can still train on the same structural neighbourhood it tests on. Blocking
    the sequence removes that, at the cost of a harsher extrapolation.
    """
    positions = np.asarray(positions, dtype=int)
    edges = np.linspace(1, length + 1, blocks + 1)
    return np.clip(np.digitize(positions, edges) - 1, 0, blocks - 1)


def contiguous_blocks(positions, width):
    """Label each row by the contiguous sequence block of `width` it falls in.

    This is the **resampling** unit, not a fold assignment, and the two want
    different granularities. Folds want few, large, well-separated blocks;
    a percentile bootstrap wants enough units to have a distribution at all —
    resampling the eight region folds themselves would leave eight. Blocks of
    ten residues keep the local structural dependence inside a cluster while
    leaving tens of clusters to resample.
    """
    positions = np.asarray(positions, dtype=int)
    _check(width >= 1, "a block must be at least one residue wide")
    return positions // int(width)


# --------------------------------------------------------------------------
# the comparators
# --------------------------------------------------------------------------
#: `M5` is the only new idea. Everything else exists to make it refutable.
COMPARATORS = ("M0_intercept", "M1_abundance_only", "M2_substitution_only",
               "M3_mechanism_only", "M4_fixed_union",
               "M5_conditional_interaction", "M6_trust_gate")


def evaluate(activity, abundance, substitution, mechanism, groups, positions,
             *, seed=None, penalty=None, permute_abundance=False,
             permute_mechanism=False):
    """Every comparator on identical grouped folds, nothing fitted on a test row.

    Two conditional forms are run, because they can disagree and the
    disagreement is informative:

    * `M5_conditional_interaction` adds abundance x mechanism products to the
      union. It is the minimal sound way to use abundance conditionally, and
      it reduces **exactly** to `M4` when those coefficients vanish.
    * `M6_trust_gate` is the interpretable form: in each training fold it
      fits, from the non-abundance blocks only, a model of how wrong the
      abundance-based prediction is, and uses it to blend an abundance
      prediction with a non-abundance one. It reports a per-variant trust
      weight, at the cost of blending two predictors that were each already
      pooled across regimes — which is a real handicap, not a bug.
    """
    seed = FROZEN["seed"] if seed is None else seed
    penalty = FROZEN["ridge_penalty"] if penalty is None else penalty
    activity = np.asarray(activity, dtype=float)
    abundance = np.asarray(abundance, dtype=float).reshape(-1, 1)
    substitution = np.asarray(substitution, dtype=float)
    mechanism = np.asarray(mechanism, dtype=float)
    groups = np.asarray(groups, dtype=int)
    positions = np.asarray(positions, dtype=int)
    n = activity.size
    _check(positions.size == n, "the position vector disagrees on n")
    for name, block in (("abundance", abundance), ("substitution", substitution),
                        ("mechanism", mechanism)):
        _check(block.shape[0] == n, f"the {name} block disagrees on n")
        _check(np.isfinite(block).all(), f"the {name} block has a missing value")
    _check(np.isfinite(activity).all(), "an activity score is missing")
    _check(np.unique(groups).size >= 3,
           f"only {np.unique(groups).size} groups; a grouped split needs 3")

    rng = np.random.default_rng(seed)
    if permute_abundance:
        abundance = abundance[rng.permutation(n)]
    if permute_mechanism:
        mechanism = mechanism[rng.permutation(n)]

    scored = {name: np.full(n, np.nan) for name in COMPARATORS}
    blend = np.full(n, np.nan)
    for fold in np.unique(groups):
        held = np.flatnonzero(groups == fold)
        keep = np.flatnonzero(groups != fold)
        if keep.size < 10 or held.size == 0:
            continue
        target = activity[keep]
        centre = float(target.mean())

        def fit(block):
            return ridge_predict(block[keep], target, block[held],
                                 penalty=penalty)

        scored["M0_intercept"][held] = centre
        scored["M1_abundance_only"][held] = fit(abundance)
        scored["M2_substitution_only"][held] = fit(substitution)
        scored["M3_mechanism_only"][held] = fit(mechanism)
        union = np.column_stack([abundance, substitution, mechanism])
        scored["M4_fixed_union"][held] = fit(union)

        # --- M5: abundance used CONDITIONALLY on the mechanism block.
        # Adding abundance x mechanism products is the minimal honest way to
        # say "rely on abundance here but not there": it stays a linear solve,
        # and it collapses back to M4 exactly when those coefficients are
        # zero, so M5 vs M4 isolates the conditional part and nothing else.
        # A blend of two POOLED predictors cannot express this, which is why
        # the first version of M5 failed its own positive control: each
        # component had already been attenuated by averaging over both
        # regimes before the blend could choose between them.
        interaction = abundance * mechanism
        scored["M5_conditional_interaction"][held] = fit(
            np.column_stack([union, interaction]))

        # --- M6: learn, WITHOUT abundance, how much abundance can be trusted
        without = np.column_stack([substitution, mechanism])
        without_keep, abundance_keep = without[keep], abundance[keep]
        # Group the inner split by POSITION WITHIN the training fold too, so
        # the trust model is not fitted and evaluated on the same residue.
        inner = position_groups(positions[keep], 3, seed + 1)
        residual = np.full(keep.size, np.nan)
        for inner_fold in np.unique(inner):
            inner_held = np.flatnonzero(inner == inner_fold)
            inner_keep = np.flatnonzero(inner != inner_fold)
            if inner_keep.size < 5:
                continue
            inner_target = target[inner_keep]
            from_abundance = ridge_predict(
                abundance_keep[inner_keep], inner_target,
                abundance_keep[inner_held], penalty=penalty)
            from_without = ridge_predict(
                without_keep[inner_keep], inner_target,
                without_keep[inner_held], penalty=penalty)
            residual[inner_held] = np.abs(
                from_abundance - target[inner_held]) - np.abs(
                    from_without - target[inner_held])
        usable = np.isfinite(residual)
        if usable.sum() < 10:
            scored["M6_trust_gate"][held] = scored["M4_fixed_union"][held]
            continue
        # Positive residual means abundance did worse than the alternative.
        advantage = ridge_predict(without_keep[usable], residual[usable],
                                  without[held], penalty=penalty)
        spread = float(np.std(residual[usable])) or 1.0
        weight = 1.0 / (1.0 + np.exp(advantage / spread))
        blend[held] = weight
        scored["M6_trust_gate"][held] = (
            weight * fit(abundance) + (1.0 - weight) * fit(without))

    complete = np.all([np.isfinite(scored[name]) for name in COMPARATORS],
                      axis=0)
    _check(complete.sum() >= 20,
           f"only {int(complete.sum())} variants were scored by every arm")
    losses = {name: np.abs(scored[name][complete] - activity[complete])
              for name in COMPARATORS}
    return {
        "n_scored": int(complete.sum()),
        "n_groups": int(np.unique(groups).size),
        # The mask is returned so a caller can (a) attach group labels to the
        # rows that were actually scored, for a cluster bootstrap, and (b)
        # assert that two runs it wants to pair scored the SAME rows.
        "scored_rows": complete,
        "permuted_abundance": bool(permute_abundance),
        "permuted_mechanism": bool(permute_mechanism),
        "losses": losses,
        "mean_loss": {name: float(values.mean())
                      for name, values in losses.items()},
        "mean_abundance_blend_weight": float(
            np.nanmean(blend[complete])) if np.isfinite(blend[complete]).any()
            else float("nan"),
        "blend_weight_spread": float(np.nanstd(blend[complete]))
            if np.isfinite(blend[complete]).any() else float("nan"),
        "what_the_blend_weight_means": (
            "1.0 means the conditional rule relied entirely on abundance for "
            "that variant and 0.0 means it ignored abundance. A spread near "
            "zero means the rule chose the same thing everywhere, so it is "
            "M4 under another name."),
    }


def _ladder(interval, threshold):
    """The project's six rungs, applied to the interval's LOWER bound.

    Written once and shared, because the last time a second contrast grew its
    own copy of this ladder the copy announced a gain of +0.62 % as present
    against a 5 % margin. A rung is a claim about where the interval sits
    relative to the smallest gain worth having, and `threshold` is that gain.
    """
    low, high = float(interval[0]), float(interval[1])
    if low >= threshold:
        return "GAIN_REACHES_MARGIN"
    if low > 0.0 and high < threshold:
        return "GAIN_PRESENT_BUT_BELOW_MARGIN"
    if low > 0.0:
        return "GAIN_PRESENT_MARGIN_UNDECIDED"
    if high < 0.0:
        return "DETECTABLY_WORSE"
    if high < threshold:
        return "REFUTED_GAIN_BELOW_MARGIN"
    return "NO_DETECTABLE_GAIN_UNDERPOWERED"


def ablation(richer, plainer, arm, *, label, family_size=1, margin=None,
             seed=None, groups=None):
    """Paired contrast of ONE arm between two knowledge layers.

    `M4` beating `M1` says the union beats the best single measurement. It
    does **not** say what the newest block contributed, because the union
    already contained everything that came before. The direct question is
    whether *these* features, added to *that* arm, on *those* rows, reduce the
    loss — and it has to be asked as a paired interval, not as a difference of
    two mean losses reported side by side.

    Both runs must have scored exactly the same rows, which is asserted rather
    than assumed: a fold where one layer produced a non-finite prediction
    would otherwise silently compare different cohorts.
    """
    margin = (FROZEN["minimum_interesting_relative_gain"] if margin is None
              else margin)
    for name, block in (("richer", richer), ("plainer", plainer)):
        _check(arm in block["losses"], f"{name} has no arm {arm!r}")
        _check(block.get("scored_rows") is not None,
               f"{name} predates scored_rows; re-run evaluate")
    mask_r = np.asarray(richer["scored_rows"], dtype=bool)
    mask_p = np.asarray(plainer["scored_rows"], dtype=bool)
    _check(mask_r.size == mask_p.size and bool((mask_r == mask_p).all()),
           "the two layers scored different rows, so the arms are not paired")
    cluster = None
    if groups is not None:
        labels = np.asarray(groups)
        _check(labels.size == mask_r.size,
               "the group labels must cover every row evaluate saw")
        cluster = labels[mask_r]
    worse, better = plainer["losses"][arm], richer["losses"][arm]
    gap = paired_gap(worse, better, family_size=family_size,
                     n_boot=FROZEN["n_bootstrap"], alpha=FROZEN["alpha"],
                     seed=FROZEN["seed"] if seed is None else seed,
                     groups=cluster)
    reference = float(np.mean(worse))
    _check(reference > 0.0, "the reference loss must be positive")
    threshold = margin * reference
    return {"contrast": f"{arm}_with_{label}_minus_{arm}_without",
            "arm": arm, "added_layer": label,
            "verdict": _ladder(gap["interval"], threshold),
            "relative_gain": float(gap["gap"] / reference),
            "minimum_interesting_gain": float(threshold),
            "what_this_isolates": (
                f"the loss of {arm} with the {label} block minus the same arm "
                "without it, on identical rows and identical folds. Nothing "
                "else differs between the two fits."),
            **gap}


def contrast(result, better, worse, *, family_size=1, margin=None, seed=None,
             groups=None):
    """Paired interval and the project's six-rung verdict ladder.

    `groups` is a **full-length** per-row label array — the protein position,
    or a contiguous-region block. It is subset to the rows this run actually
    scored and handed to `paired_gap` as the resampling unit. Supply it
    whenever rows are not independent, which for a deep-mutational-scanning
    cohort means always: the folds being grouped does not make the interval
    grouped, and only the interval decides a verdict.
    """
    margin = (FROZEN["minimum_interesting_relative_gain"] if margin is None
              else margin)
    losses = result["losses"]
    _check(better in losses and worse in losses,
           f"unknown comparator in {better!r} vs {worse!r}")
    cluster = None
    if groups is not None:
        labels = np.asarray(groups)
        mask = result.get("scored_rows")
        _check(mask is not None,
               "this result predates scored_rows; re-run evaluate")
        _check(labels.size == np.asarray(mask).size,
               "the group labels must cover every row evaluate saw")
        cluster = labels[np.asarray(mask, dtype=bool)]
    gap = paired_gap(losses[worse], losses[better], family_size=family_size,
                     n_boot=FROZEN["n_bootstrap"], alpha=FROZEN["alpha"],
                     seed=FROZEN["seed"] if seed is None else seed,
                     groups=cluster)
    reference = float(np.mean(losses[worse]))
    _check(reference > 0.0, "the reference loss must be positive")
    threshold = margin * reference
    name = _ladder(gap["interval"], threshold)
    return {"contrast": f"{better}_minus_{worse}", "verdict": name,
            "relative_gain": float(gap["gap"] / reference),
            "minimum_interesting_gain": float(threshold),
            "an_interval_containing_zero_is_not_a_refutation": True, **gap}


# --------------------------------------------------------------------------
# positive control: the harness must see a conditional rule that exists
# --------------------------------------------------------------------------
def conditionally_gateable_cohort(n_positions=180, seed=0):
    """A protein where abundance is informative on one half and useless on the
    other, and a mechanism feature says which.

    `M1` is mediocre everywhere, `M4` cannot fix it because the fault is
    conditional, and a rule that reads the mechanism feature recovers most of
    it. Without this control a null `M5` vs `M4` would be unreadable.
    """
    rng = np.random.default_rng(seed)
    reference_positions = np.repeat(np.arange(1, n_positions + 1), 4)
    n = reference_positions.size
    near_site = (reference_positions <= n_positions // 2)
    truth = rng.normal(0.0, 1.0, size=n)
    # Abundance tracks activity away from the site and is pure noise near it.
    abundance = np.where(near_site, rng.normal(0.0, 1.0, size=n),
                         truth + rng.normal(0.0, 0.15, size=n))
    # One substitution feature carries the activity signal near the site.
    carrier = np.where(near_site, truth + rng.normal(0.0, 0.15, size=n),
                       rng.normal(0.0, 1.0, size=n))
    substitution = np.column_stack([carrier, rng.normal(size=(n, 3))])
    mechanism = np.column_stack([near_site.astype(float),
                                 rng.normal(size=(n, 2))])
    return (truth, abundance, substitution, mechanism, reference_positions)


def self_test():
    """Assert the harness recovers a conditional rule, and loses it when the
    mechanism feature is permuted."""
    truth, abundance, substitution, mechanism, positions = \
        conditionally_gateable_cohort()
    groups = position_groups(positions, FROZEN["position_folds"],
                             FROZEN["seed"])
    live = evaluate(truth, abundance, substitution, mechanism, groups,
                    positions)
    dead = evaluate(truth, abundance, substitution, mechanism, groups,
                    positions, permute_mechanism=True)
    gain = contrast(live, "M5_conditional_interaction", "M4_fixed_union")
    control = contrast(dead, "M5_conditional_interaction", "M4_fixed_union")
    gate = contrast(live, "M6_trust_gate", "M4_fixed_union")
    # The cluster bootstrap must still SEE an effect that is there. Widening
    # every interval until nothing is ever detected would not be rigour, it
    # would be trading one error for the opposite one, so the positive control
    # is run under the resampling unit the real analysis uses.
    clustered = contrast(live, "M5_conditional_interaction", "M4_fixed_union",
                         groups=positions)
    _check(gain["verdict"] == "GAIN_REACHES_MARGIN",
           "the harness must detect a conditional rule that exists; got "
           f"{gain['verdict']} ({gain['relative_gain']:+.3f})")
    _check(clustered["verdict"] == "GAIN_REACHES_MARGIN",
           "the POSITION-CLUSTERED bootstrap must also detect a conditional "
           f"rule that exists; got {clustered['verdict']}")
    # No assertion on the DIRECTION of the width change. A cluster bootstrap
    # widens an interval when the statistic is positively correlated within a
    # cluster and can narrow it otherwise, and in this synthetic cohort only
    # the mechanism indicator is shared within a position while the activity,
    # the abundance and the carrier are all drawn per row — so the gap has
    # little intra-cluster correlation and the clustered interval comes out
    # slightly NARROWER. This assertion was written the other way round and
    # fired on its own control. The ratio is measured and reported; it is not
    # a law.
    _check(live["blend_weight_spread"] > 0.01,
           "the trust gate must actually vary its weight; spread "
           f"{live['blend_weight_spread']:.4f}")
    _check(not control["verdict"].startswith("GAIN_REACHES"),
           "permuting the mechanism feature must remove the advantage; got "
           f"{control['verdict']}")
    return {
        "with_real_mechanism": {"mean_loss": live["mean_loss"],
                                "M5_vs_M4": gain["verdict"],
                                "relative_gain": gain["relative_gain"],
                                "M6_vs_M4": gate["verdict"],
                                "M6_relative_gain": gate["relative_gain"],
                                "blend_spread": live["blend_weight_spread"]},
        "with_permuted_mechanism": {"M5_vs_M4": control["verdict"],
                                    "relative_gain": control["relative_gain"]},
        "position_clustered_bootstrap": {
            "M5_vs_M4": clustered["verdict"],
            "interval": clustered["interval"],
            "n_clusters": clustered["n_resampling_units"],
            "width_ratio_against_row_wise": float(
                (clustered["interval"][1] - clustered["interval"][0])
                / (gain["interval"][1] - gain["interval"][0])),
            "why": ("a resampling unit that can no longer detect a real "
                    "effect is not a stricter analysis, it is a broken one."),
        },
        "what_this_proves": (
            "a mechanism-conditional choice of which measurement to trust is "
            "recoverable by this harness when it is really there, and is lost "
            "when the mechanism feature is permuted. A null on real data is "
            "therefore about the data."),
    }


def what_this_module_does_not_do():
    return {
        "is_not_a_pharmacokinetic_result": (
            "the target is a functional-screen activity score. Not a "
            "clearance, not an AUC, not a concentration."),
        "is_not_a_multi_omic_measurement_in_individuals": (
            "these are engineered single variants in cell or yeast systems. "
            "No individual is measured on several layers."),
        "does_not_confirm_a_mechanism": (
            "mechanism_hypotheses ranks candidates and nothing treats a "
            "candidate as established."),
        "does_not_treat_an_annotation_as_causal": (
            "a UniProt binding-site annotation is an annotation. Proximity to "
            "one is a hypothesis about function, not a causal claim, and its "
            "predictive usefulness is measured by M3 rather than assumed."),
        "measurement_type_versus_biomarker": (
            "this asks which measurement TYPE is informative. Which specific "
            "biomarker to measure in a cohort is a separate question."),
    }
