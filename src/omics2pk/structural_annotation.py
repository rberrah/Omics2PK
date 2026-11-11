"""Three-dimensional distance to a protein's catalytic machinery, from a PDB file.

Why this module exists
----------------------
`variant_mechanism_omic_selection` measured a real mechanistic effect — a
variant within 3 residues of an annotated functional site is 6.74x more likely
to keep protein abundance and lose activity — and then found that the
annotation reaches fewer than half the discordant variants. The residues it
missed are 32 to 134 positions away from any annotation **in sequence**.

Sequence distance is the wrong metric for a folded protein. CYP2C9's most
frequently discordant residue is 301, which is 134 positions from the only
UniProt annotation (Cys435) and sits in the I-helix that packs directly against
the heme. A one-dimensional ruler cannot see that. This module replaces it with
a three-dimensional one.

**That makes the earlier diagnosis refutable.** If annotation coverage was
really the binding constraint, a structural ruler should reach the discordant
residues a sequence ruler misses. If it does not, the diagnosis was wrong, and
this module is the instrument that says so.

What a distance here does and does not mean
-------------------------------------------
Proximity to a cofactor or a bound substrate is a **hypothesis** that a variant
affects catalysis rather than folding. It is not a causal claim and not a
prediction of predictive usefulness — the same discipline this project applies
to a STRING confidence score. Its value is measured by whether the arms that
use it beat the arms that do not.

The numbering guard
-------------------
A crystal structure's residue numbers need not agree with the reference
sequence: constructs are truncated, and a silent offset would corrupt every
distance while still producing plausible-looking numbers. `verify_catalytic_
geometry` therefore checks a fact that is only true when the numbering is
right: **UniProt calls Cys435 CYP2C9's axial heme ligand, so its sulphur must
sit about 2.3 A from the iron.** An offset of even one residue breaks that, and
the check fails loudly rather than returning a plausible wrong answer.
"""
from __future__ import annotations

import math

import numpy as np

#: Frozen before execution.
FROZEN = {
    # Heavy atoms only: hydrogens are absent from most crystal structures, so
    # including them would make distances depend on whether a depositor
    # modelled them.
    "exclude_hydrogen": True,
    "exclude_water": True,
    # A residue is "at" the machinery when any of its heavy atoms comes within
    # this distance of any ligand heavy atom. 5 A is contact range; 8 A is the
    # second shell. Both are reported so the choice is visible.
    "contact_angstrom": 5.0,
    "second_shell_angstrom": 8.0,
    # The axial ligand check. A coordinate bond is ~2.3 A; 3.0 allows for
    # resolution without allowing a one-residue offset to pass.
    "axial_bond_max_angstrom": 3.0,
}

#: Catalytic machinery per protein, by PDB ligand code. Declared, not guessed.
CATALYTIC_CENTRES = {
    "CYP2C9": {
        "pdb": "1R9O",
        "resolution_angstrom": 2.0,
        "title": "CYP2C9 with flurbiprofen bound",
        "cofactor_ligands": ("HEM",),
        "substrate_ligands": ("FLP",),
        "metal_atom": "FE",
        # UniProt P11712: "axial binding residue".
        "axial_residue": 435,
        "axial_atom": "SG",
        "uniprot": "P11712",
    },
    "NUDT15": {
        "pdb": "5LPG",
        "resolution_angstrom": 1.7,
        "title": "NUDT15 in complex with 6-thio-GMP",
        "cofactor_ligands": ("MG",),
        "substrate_ligands": ("71V",),
        "metal_atom": "MG",
        # The Nudix box coordinates the metal; no single axial residue is
        # annotated, so the guard below uses the metal-to-motif distance.
        "axial_residue": None,
        "axial_atom": None,
        "motif_span": (48, 69),
        "uniprot": "Q9NV35",
    },
}

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


class StructuralAnnotationError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise StructuralAnnotationError(message)


def parse_pdb(text):
    """Minimal fixed-column PDB reader. No new dependency, no silent coercion.

    Returns `(residues, ligands)`. `residues` maps a residue number to
    `{"name": three-letter, "atoms": {atom_name: (x, y, z)}}` for the first
    polymer chain encountered, because a variant position refers to one chain.
    `ligands` maps a ligand code to a list of `(atom_name, (x, y, z))`.
    """
    residues, ligands, chain = {}, {}, None
    for line in text.splitlines():
        record = line[:6]
        if record not in ("ATOM  ", "HETATM"):
            continue
        element = line[76:78].strip().upper()
        if FROZEN["exclude_hydrogen"] and element in ("H", "D"):
            continue
        name = line[17:20].strip().upper()
        if FROZEN["exclude_water"] and name in ("HOH", "DOD", "WAT"):
            continue
        atom = line[12:16].strip().upper()
        try:
            point = (float(line[30:38]), float(line[38:46]),
                     float(line[46:54]))
        except ValueError:
            continue
        if record == "HETATM" and name not in _THREE_TO_ONE:
            ligands.setdefault(name, []).append((atom, point))
            continue
        this_chain = line[21]
        if chain is None:
            chain = this_chain
        if this_chain != chain:
            continue
        # Skip insertion codes: they are not reference-sequence positions.
        if line[26].strip():
            continue
        try:
            number = int(line[22:26])
        except ValueError:
            continue
        entry = residues.setdefault(number, {"name": name, "atoms": {}})
        entry["atoms"][atom] = point
    _check(residues, "no polymer residue was parsed from this file")
    return residues, ligands


def _distance(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                     + (a[2] - b[2]) ** 2)


def _min_distance(residue, points):
    if not points:
        return float("nan")
    return min(_distance(atom, point)
               for atom in residue["atoms"].values() for point in points)


def residue_distances(residues, ligands, centre):
    """Per-residue minimum heavy-atom distance to cofactor, metal and substrate."""
    cofactor = [point for code in centre["cofactor_ligands"]
                for _, point in ligands.get(code, ())]
    substrate = [point for code in centre["substrate_ligands"]
                 for _, point in ligands.get(code, ())]
    metal = [point for code in centre["cofactor_ligands"]
             for atom, point in ligands.get(code, ())
             if atom == centre["metal_atom"]]
    _check(cofactor, f"no cofactor atoms found for "
                     f"{centre['cofactor_ligands']} in {centre['pdb']}")
    _check(substrate, f"no substrate atoms found for "
                      f"{centre['substrate_ligands']} in {centre['pdb']}")
    out = {}
    for number, residue in residues.items():
        out[number] = {
            "to_cofactor": _min_distance(residue, cofactor),
            "to_metal": _min_distance(residue, metal) if metal
                        else float("nan"),
            "to_substrate": _min_distance(residue, substrate),
            "residue": _THREE_TO_ONE.get(residue["name"], "X"),
        }
    return out


def verify_catalytic_geometry(residues, ligands, centre):
    """Fail loudly if the structure's numbering does not match the reference.

    The check is a fact that holds only when the numbering is right. For
    CYP2C9, UniProt annotates Cys435 as the **axial** heme ligand, so its
    sulphur must sit within a coordinate-bond distance of the iron. A
    one-residue offset puts a different residue there and breaks it.

    For NUDT15 no single axial residue is annotated, so the guard instead
    requires the catalytic metal to be inside the Nudix box's contact shell,
    which is equally false under an offset of any size that would matter.
    """
    report = {"pdb": centre["pdb"], "uniprot": centre["uniprot"],
              "n_residues_parsed": len(residues),
              "residue_range": [min(residues), max(residues)]}
    metal = [point for code in centre["cofactor_ligands"]
             for atom, point in ligands.get(code, ())
             if atom == centre["metal_atom"]]
    _check(metal, f"no {centre['metal_atom']} atom in {centre['pdb']}")

    axial = centre.get("axial_residue")
    if axial is not None:
        _check(axial in residues,
               f"{centre['pdb']} has no residue {axial}; the construct may be "
               "truncated or renumbered")
        found = residues[axial]["name"]
        _check(found == "CYS",
               f"{centre['pdb']} residue {axial} is {found}, not CYS. The "
               "numbering does not match the reference sequence, so every "
               "distance computed from it would be wrong")
        atom = residues[axial]["atoms"].get(centre["axial_atom"])
        _check(atom is not None,
               f"residue {axial} has no {centre['axial_atom']} atom")
        bond = min(_distance(atom, point) for point in metal)
        _check(bond <= FROZEN["axial_bond_max_angstrom"],
               f"Cys{axial} SG is {bond:.2f} A from the "
               f"{centre['metal_atom']}, beyond the "
               f"{FROZEN['axial_bond_max_angstrom']} A a coordinate bond "
               "allows. The numbering is almost certainly offset")
        report["axial_bond_angstrom"] = round(bond, 3)
        report["axial_residue_checked"] = f"Cys{axial}"
    else:
        start, end = centre["motif_span"]
        inside = [n for n in residues if start <= n <= end]
        _check(inside, f"{centre['pdb']} has no residue in the motif span "
                       f"{start}-{end}")
        nearest = min(_min_distance(residues[n], metal) for n in inside)
        _check(nearest <= FROZEN["contact_angstrom"],
               f"the nearest motif residue is {nearest:.2f} A from the "
               f"{centre['metal_atom']}, beyond contact range. The numbering "
               "does not match the reference sequence")
        report["motif_to_metal_angstrom"] = round(nearest, 3)
        report["motif_span_checked"] = list(centre["motif_span"])
    report["verified"] = True
    return report


def structural_features(position, distances):
    """The block a comparator uses. Unresolved positions are marked, not faked.

    A residue absent from the crystal structure — a disordered loop, or the
    truncated membrane anchor of a microsomal CYP — gets the median distance
    of the resolved residues and an explicit `unresolved` flag. Imputing
    silently would let "not in the structure" masquerade as "far from the
    machinery", and for a CYP's N-terminal anchor those are opposite claims.
    """
    resolved = distances.get(int(position))
    if resolved is None:
        return [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
    cofactor = resolved["to_cofactor"]
    substrate = resolved["to_substrate"]
    return [
        math.log1p(cofactor),
        math.log1p(substrate),
        1.0 if cofactor <= FROZEN["contact_angstrom"] else 0.0,
        1.0 if substrate <= FROZEN["contact_angstrom"] else 0.0,
        1.0 if min(cofactor, substrate) <= FROZEN["second_shell_angstrom"]
        else 0.0,
        0.0,
    ]


STRUCTURAL_NAMES = ("log1p_angstrom_to_cofactor", "log1p_angstrom_to_substrate",
                    "in_cofactor_contact", "in_substrate_contact",
                    "in_second_shell", "unresolved_in_structure")


def build_structural_block(positions, distances):
    rows = [structural_features(int(p), distances) for p in positions]
    block = np.asarray(rows, dtype=float)
    median = np.median([d["to_cofactor"] for d in distances.values()
                        if np.isfinite(d["to_cofactor"])] or [0.0])
    unresolved = block[:, STRUCTURAL_NAMES.index("unresolved_in_structure")]
    # Fill the marked rows with the median, keeping the flag so a model can
    # learn that "unresolved" is its own category.
    mask = unresolved > 0.0
    if mask.any():
        block[mask, 0] = math.log1p(float(median))
        block[mask, 1] = math.log1p(float(median))
    return block


def coverage(positions, discordant, distances, *, cutoff=None):
    """How many discordant variants a structural ruler reaches.

    This is the number that decides whether the earlier diagnosis was right.
    A sequence ruler reached 48.8 % of NUDT15's discordant variants and 12.5 %
    of CYP2C9's; if 3D proximity does not do better, coverage was not the
    binding constraint.
    """
    cutoff = FROZEN["second_shell_angstrom"] if cutoff is None else cutoff
    positions = np.asarray(positions, dtype=int)
    discordant = np.asarray(discordant, dtype=bool)
    reach = np.array([
        min(distances[int(p)]["to_cofactor"], distances[int(p)]["to_substrate"])
        <= cutoff if int(p) in distances else False
        for p in positions])
    near, far = reach, ~reach
    out = {
        "cutoff_angstrom": float(cutoff),
        "n_within": int(near.sum()),
        "n_beyond": int(far.sum()),
        "n_unresolved": int(sum(1 for p in positions
                                if int(p) not in distances)),
        "discordance_rate_within": float(discordant[near].mean())
            if near.any() else float("nan"),
        "discordance_rate_beyond": float(discordant[far].mean())
            if far.any() else float("nan"),
        "share_of_discordant_variants_reached": float(near[discordant].mean())
            if discordant.any() else float("nan"),
    }
    if far.any() and discordant[far].mean() > 0 and near.any():
        out["enrichment_ratio"] = float(discordant[near].mean()
                                        / discordant[far].mean())
        out["enrichment_unavailable_because"] = None
    else:
        # A ratio with a zero denominator is unbounded, not large. Saying so
        # beats returning a bare NaN that a reader has to guess about.
        out["enrichment_ratio"] = float("nan")
        if not near.any():
            reason = "no variant lies within the cutoff"
        elif not far.any():
            reason = "every variant lies within the cutoff"
        else:
            reason = ("no discordant variant lies beyond the cutoff, so the "
                      "ratio is unbounded rather than large")
        out["enrichment_unavailable_because"] = reason
    return out


def enrichment_interval(positions, discordant, distances, *, cutoff=None,
                        n_boot=2000, seed=20261112):
    cutoff = FROZEN["second_shell_angstrom"] if cutoff is None else cutoff
    positions = np.asarray(positions, dtype=int)
    discordant = np.asarray(discordant, dtype=bool)
    reach = np.array([
        min(distances[int(p)]["to_cofactor"], distances[int(p)]["to_substrate"])
        <= cutoff if int(p) in distances else False
        for p in positions])
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n_boot):
        picked = rng.integers(0, positions.size, size=positions.size)
        near, far = reach[picked], ~reach[picked]
        if near.sum() < 5 or far.sum() < 5:
            continue
        denominator = discordant[picked][far].mean()
        if denominator <= 0:
            continue
        draws.append(discordant[picked][near].mean() / denominator)
    if len(draws) < 100:
        return [float("nan"), float("nan")]
    return [float(np.percentile(draws, 2.5)),
            float(np.percentile(draws, 97.5))]


def what_this_module_does_not_do():
    return {
        "a_distance_is_not_a_cause": (
            "proximity to a cofactor or a bound substrate is a hypothesis "
            "that a variant affects catalysis rather than folding. Its "
            "predictive usefulness is measured by the arms that use it, never "
            "assumed."),
        "one_structure_is_one_conformation": (
            "a crystal structure is one conformation of one construct, often "
            "with a flexible region unmodelled and, for a microsomal CYP, the "
            "membrane anchor truncated. Unresolved positions are flagged, not "
            "quietly imputed."),
        "a_bound_ligand_is_not_every_substrate": (
            "flurbiprofen in 1R9O and 6-thio-GMP in 5LPG define one pocket "
            "each. A different substrate can occupy a different volume, so "
            "`to_substrate` is specific to the complex that was crystallised."),
        "numbering_is_verified_not_trusted": (
            "verify_catalytic_geometry fails unless a fact that depends on "
            "correct numbering holds, because an offset would produce "
            "plausible wrong distances rather than an error."),
    }
