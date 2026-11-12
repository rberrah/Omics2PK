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
    # Pre-registered in OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md
    # §3.3, selected by the same rule as BRAF: the best-resolution RCSB entry
    # for UniProt P10398 carrying an ATP analogue in the ARAF chain. The rule
    # picks 9AXM uniquely -- it is the ONLY ARAF entry with an ATP analogue at
    # all -- so dropping the "legacy PDB format" clause could not steer the
    # choice. See §4.2 for that deviation.
    "ARAF": {
        "pdb": "9AXM",
        "resolution_angstrom": 2.42,
        "title": ("ARAF:MEK1 with AMP-PNP and Mg, inactive aC-out "
                  "conformation"),
        # MEK1 is chain A and carries its own AMP-PNP and Mg, exactly as in
        # 6PP9, so the chain restriction is load-bearing here too.
        "chain": "B",
        "cofactor_ligands": ("ANP", "MG"),
        "substrate_ligands": (),
        "substrate_chain": "A",
        "substrate_is_a_protein": "MEK1",
        "metal_atom": "MG",
        "axial_residue": None,
        # Verified against the coordinates BEFORE being declared, which is
        # what BRAF taught: 336 is UniProt's catalytic lysine, 429 its
        # proton-acceptor aspartate, 447-449 the DFG motif found by scanning
        # the activation segment rather than assumed, 354 the aC glutamate
        # and 382 the gatekeeper by alignment with BRAF's E501 and T529.
        "residue_identities": {336: "LYS", 354: "GLU", 382: "THR",
                               429: "ASP", 447: "ASP", 448: "PHE",
                               449: "GLY"},
        "nucleotide_contact": {"residue": 336, "atom": "NZ", "ligand": "ANP",
                               "max_angstrom": 4.0},
        "uniprot": "P10398",
    },
    # Pre-registered in OMICS2PK_BRAF_PREREGISTRATION_20261111.md, selected
    # by rule rather than by hand: the best-resolution RCSB entry for
    # UniProt P15056 carrying an ATP-analogue ligand in the BRAF chain.
    "BRAF": {
        "pdb": "6PP9",
        "resolution_angstrom": 2.59,
        "title": ("BRAF:MEK1 with AMP-PNP and Mg, inactive aC-out "
                  "conformation"),
        # MEK1 is chain B and carries its own AMP-PNP and Mg, so the chain
        # restriction is load-bearing, not cosmetic.
        "chain": "A",
        "cofactor_ligands": ("ANP", "MG"),
        "substrate_ligands": (),
        "substrate_chain": "B",
        "substrate_is_a_protein": "MEK1",
        "metal_atom": "MG",
        "axial_residue": None,
        "residue_identities": {483: "LYS", 501: "GLU", 529: "THR",
                               594: "ASP", 595: "PHE", 596: "GLY",
                               600: "VAL"},
        "nucleotide_contact": {"residue": 483, "atom": "NZ", "ligand": "ANP",
                               "max_angstrom": 4.0},
        "uniprot": "P15056",
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


def parse_pdb(text, chain=None):
    """Minimal fixed-column PDB reader. No new dependency, no silent coercion.

    Returns `(residues, ligands)`. `residues` maps a residue number to
    `{"name": three-letter, "atoms": {atom_name: (x, y, z)}}` for one polymer
    chain, because a variant position refers to one chain.
    `ligands` maps a ligand code to a list of `(atom_name, (x, y, z))`.

    **`chain` restricts the ligands too, and for a complex it must.** In
    6PP9, BRAF is chain A and its substrate MEK1 is chain B — and MEK1 carries
    its own AMP-PNP and its own magnesium. Taking ligands from the whole file
    would silently measure every BRAF residue's distance to whichever
    nucleotide happened to be closer, including the other kinase's. That is
    the same failure mode as a numbering offset: plausible wrong distances
    rather than an error. With `chain=None` the first polymer chain
    encountered is used and every ligand is kept, which is what a
    single-chain file needs and what the two earlier proteins were read with.
    """
    residues, ligands = {}, {}
    chain = chain.strip().upper() if isinstance(chain, str) else None
    wanted = chain
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
        this_chain = line[21].strip().upper()
        if record == "HETATM" and name not in _THREE_TO_ONE:
            if wanted is None or this_chain == wanted:
                ligands.setdefault(name, []).append((atom, point))
            continue
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
    _check(residues, "no polymer residue was parsed from this file"
                     + (f" for chain {wanted}" if wanted else ""))
    return residues, ligands


def parse_mmcif(text, chain=None):
    """Same contract as `parse_pdb`, for entries RCSB serves only as mmCIF.

    Recent depositions have no legacy PDB file at all -- 9AXM, the only ARAF
    entry carrying an ATP analogue, returns 404 for `.pdb` and 200 for
    `.cif`. Refusing it would declare a protocol unevaluable because of this
    module's reader rather than because of the data, which is the wrong kind
    of strictness when the brief requires using public data that is actually
    accessible.

    **`auth_asym_id` and `auth_seq_id` are used, not the label_ fields.** The
    author numbering is the one a UniProt position refers to and the one a PDB
    file would have carried; `label_seq_id` restarts at 1 per entity and would
    silently shift every residue number. Only the first model is read.
    """
    wanted = chain.strip().upper() if isinstance(chain, str) else None
    residues, ligands = {}, {}
    lines = text.splitlines()
    index = 0
    columns = None
    while index < len(lines):
        line = lines[index].strip()
        if line == "loop_":
            headers, cursor = [], index + 1
            while cursor < len(lines) and lines[cursor].strip().startswith("_"):
                headers.append(lines[cursor].strip())
                cursor += 1
            if headers and headers[0].startswith("_atom_site."):
                columns = {name.split(".", 1)[1]: position
                           for position, name in enumerate(headers)}
                index = cursor
                break
            index = cursor
            continue
        index += 1
    _check(columns is not None, "no _atom_site loop found in this mmCIF file")
    for required in ("group_PDB", "label_atom_id", "label_comp_id",
                     "Cartn_x", "Cartn_y", "Cartn_z", "auth_seq_id",
                     "auth_asym_id"):
        _check(required in columns,
               f"the mmCIF _atom_site loop lacks {required}")
    while index < len(lines):
        raw = lines[index]
        index += 1
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            break
        if stripped.startswith("_") or stripped == "loop_":
            break
        fields = stripped.split()
        if len(fields) < len(columns):
            continue
        if "pdbx_PDB_model_num" in columns and \
                fields[columns["pdbx_PDB_model_num"]] not in ("1", "."):
            continue
        record = fields[columns["group_PDB"]]
        name = fields[columns["label_comp_id"]].strip().upper()
        if FROZEN["exclude_water"] and name in ("HOH", "DOD", "WAT"):
            continue
        element = (fields[columns["type_symbol"]].strip().upper()
                   if "type_symbol" in columns else "")
        if FROZEN["exclude_hydrogen"] and element in ("H", "D"):
            continue
        atom = fields[columns["label_atom_id"]].strip().strip('"').upper()
        try:
            point = (float(fields[columns["Cartn_x"]]),
                     float(fields[columns["Cartn_y"]]),
                     float(fields[columns["Cartn_z"]]))
        except ValueError:
            continue
        this_chain = fields[columns["auth_asym_id"]].strip().upper()
        if record == "HETATM" and name not in _THREE_TO_ONE:
            if wanted is None or this_chain == wanted:
                ligands.setdefault(name, []).append((atom, point))
            continue
        if wanted is not None and this_chain != wanted:
            continue
        if "pdbx_PDB_ins_code" in columns and \
                fields[columns["pdbx_PDB_ins_code"]] not in (".", "?"):
            continue
        try:
            number = int(fields[columns["auth_seq_id"]])
        except ValueError:
            continue
        entry = residues.setdefault(number, {"name": name, "atoms": {}})
        entry["atoms"][atom] = point
    _check(residues, "no polymer residue was parsed from this mmCIF file"
                     + (f" for chain {wanted}" if wanted else ""))
    return residues, ligands


def parse_structure(text, chain=None):
    """Read either format. mmCIF is detected by its own loop header."""
    if "_atom_site." in text:
        return parse_mmcif(text, chain=chain)
    return parse_pdb(text, chain=chain)


def mmcif_chain_points(text, chain):
    """`polymer_chain_points` for mmCIF: a protein substrate's heavy atoms."""
    wanted = chain.strip().upper()
    residues, _ = parse_mmcif(text, chain=wanted)
    points = [point for entry in residues.values()
              for point in entry["atoms"].values()]
    _check(points, f"no polymer atom found for chain {wanted}")
    return points


def chain_points(text, chain):
    """Either format."""
    if "_atom_site." in text:
        return mmcif_chain_points(text, chain)
    return polymer_chain_points(text, chain)


def polymer_chain_points(text, chain):
    """Heavy-atom coordinates of one polymer chain, for a PROTEIN substrate.

    A kinase's substrate is another protein, so `to_substrate` cannot be the
    distance to a bound small molecule. It becomes the distance to the nearest
    heavy atom of the substrate chain — an **interface** distance, not a
    pocket distance. The two are different quantities and are reported as
    such; calling them both "to substrate" without saying so would be the
    kind of quiet redefinition this project exists to avoid.
    """
    wanted = chain.strip().upper()
    points = []
    for line in text.splitlines():
        if line[:6] != "ATOM  " or line[21].strip().upper() != wanted:
            continue
        if FROZEN["exclude_hydrogen"] and line[76:78].strip().upper() in ("H", "D"):
            continue
        try:
            points.append((float(line[30:38]), float(line[38:46]),
                           float(line[46:54])))
        except ValueError:
            continue
        _ = None
    _check(points, f"no polymer atom found for chain {wanted}")
    return points


def _distance(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                     + (a[2] - b[2]) ** 2)


def _min_distance(residue, points):
    if not points:
        return float("nan")
    return min(_distance(atom, point)
               for atom in residue["atoms"].values() for point in points)


def residue_distances(residues, ligands, centre, *, substrate_points=None):
    """Per-residue minimum heavy-atom distance to cofactor, metal and substrate.

    `substrate_points` supplies the substrate when it is a **protein chain**
    rather than a ligand, as it is for a kinase. A centre that declares no
    substrate ligand must then be given one, or the call fails: a missing
    substrate silently becoming `inf` would turn `min(cofactor, substrate)`
    into a cofactor-only ruler without saying so.
    """
    cofactor = [point for code in centre["cofactor_ligands"]
                for _, point in ligands.get(code, ())]
    substrate = [point for code in centre["substrate_ligands"]
                 for _, point in ligands.get(code, ())]
    if substrate_points is not None:
        _check(not substrate,
               f"{centre['pdb']} declares substrate ligands AND a substrate "
               "chain; only one can define the substrate")
        substrate = list(substrate_points)
    metal = [point for code in centre["cofactor_ligands"]
             for atom, point in ligands.get(code, ())
             if atom == centre["metal_atom"]]
    _check(cofactor, f"no cofactor atoms found for "
                     f"{centre['cofactor_ligands']} in {centre['pdb']}")
    _check(substrate, f"no substrate atoms found for "
                      f"{centre['substrate_ligands']} in {centre['pdb']}; a "
                      "protein substrate needs substrate_points")
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

    identities = centre.get("residue_identities")
    if identities:
        # A kinase has no axial ligand and no short catalytic motif, so the
        # guard is a set of residue identities at specific numbers plus one
        # distance that depends on them. Five identities at five numbers do
        # not survive an offset of any size.
        checked = {}
        for number, expected in sorted(identities.items()):
            _check(number in residues,
                   f"{centre['pdb']} has no residue {number}; the construct "
                   "may be truncated or renumbered")
            found = residues[number]["name"]
            _check(found == expected,
                   f"{centre['pdb']} residue {number} is {found}, not "
                   f"{expected}. The numbering does not match the reference "
                   "sequence, so every distance from it would be wrong")
            checked[str(number)] = expected
        report["residue_identities_checked"] = checked
        contact = centre.get("nucleotide_contact")
        _check(contact is not None,
               "an identity guard needs one distance that depends on it")
        number, atom_name = contact["residue"], contact["atom"]
        atom = residues[number]["atoms"].get(atom_name)
        _check(atom is not None,
               f"residue {number} has no {atom_name} atom")
        points = [point for a, point in ligands.get(contact["ligand"], ())]
        _check(points,
               f"no {contact['ligand']} atom in {centre['pdb']} chain "
               f"{centre.get('chain')}")
        bond = min(_distance(atom, point) for point in points)
        _check(bond <= contact["max_angstrom"],
               f"residue {number} {atom_name} is {bond:.2f} A from "
               f"{contact['ligand']}, beyond the {contact['max_angstrom']} A "
               "this contact allows. The numbering is almost certainly offset")
        report["nucleotide_contact_angstrom"] = round(bond, 3)
        report["nucleotide_contact_checked"] = (
            f"{identities[number]}{number} {atom_name} to "
            f"{contact['ligand']}")
        report["why_not_the_salt_bridge"] = (
            "the obvious guard for a kinase, the catalytic lysine to aC "
            "glutamate salt bridge, is BROKEN in an inactive aC-out "
            "conformation such as 6PP9 (12.69 A measured). Declaring it would "
            "have failed on correctly numbered coordinates. A numbering guard "
            "must catch a plausible wrong answer without crying wolf on a "
            "right one.")
        report["verified"] = True
        return report

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
                        n_boot=2000, seed=20261112, groups=None):
    """Interval for the enrichment ratio, resampling rows or clusters.

    **The resampling unit matters more here than anywhere else in the
    project.** `reach` is a function of the *position*, so every substitution
    at one residue carries an identical value of the explanatory variable. A
    row bootstrap therefore treats nineteen substitutions at one near-the-heme
    residue as nineteen independent pieces of evidence that near-the-heme
    residues are discordant, when they are one. Pass `groups` — the position
    label, or a coarser region label — to resample clusters instead.
    """
    cutoff = FROZEN["second_shell_angstrom"] if cutoff is None else cutoff
    positions = np.asarray(positions, dtype=int)
    discordant = np.asarray(discordant, dtype=bool)
    reach = np.array([
        min(distances[int(p)]["to_cofactor"], distances[int(p)]["to_substrate"])
        <= cutoff if int(p) in distances else False
        for p in positions])
    return discordance_ratio_interval(reach, discordant, groups=groups,
                                      n_boot=n_boot, seed=seed)


def discordance_ratio_interval(reach, discordant, *, groups=None, n_boot=2000,
                               seed=20261112):
    """Interval for P(discordant | reached) / P(discordant | not reached).

    One implementation serves the sequence ruler and the three-dimensional
    one, so the two enrichment figures cannot be resampled by different rules
    and then printed in the same table.

    With `groups`, clusters are resampled with replacement. The ratio over a
    resampled cluster set is a ratio of summed cluster totals, which is exact,
    so nothing is materialised per draw.
    """
    reach = np.asarray(reach, dtype=bool)
    discordant = np.asarray(discordant, dtype=bool)
    if reach.shape != discordant.shape:
        raise ValueError("reach and discordant must align")
    rng = np.random.default_rng(seed)
    if groups is None:
        codes = np.arange(reach.size)
        n_units = reach.size
    else:
        labels = np.asarray(groups)
        if labels.shape != reach.shape:
            raise ValueError("one group label per row is required")
        keys, codes = np.unique(labels, return_inverse=True)
        n_units = keys.size
        if n_units < 3:
            return [float("nan"), float("nan")]
    n_near = np.bincount(codes, weights=reach.astype(float), minlength=n_units)
    d_near = np.bincount(codes, weights=(reach & discordant).astype(float),
                         minlength=n_units)
    n_far = np.bincount(codes, weights=(~reach).astype(float),
                        minlength=n_units)
    d_far = np.bincount(codes, weights=(~reach & discordant).astype(float),
                        minlength=n_units)
    picked = rng.integers(0, n_units, size=(n_boot, n_units))
    near_total, far_total = n_near[picked].sum(axis=1), n_far[picked].sum(axis=1)
    d_near_total = d_near[picked].sum(axis=1)
    d_far_total = d_far[picked].sum(axis=1)
    usable = (near_total >= 5) & (far_total >= 5) & (d_far_total > 0)
    if usable.sum() < 100:
        return [float("nan"), float("nan")]
    draws = ((d_near_total[usable] / near_total[usable])
             / (d_far_total[usable] / far_total[usable]))
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
        "a_protein_substrate_is_an_interface_not_a_pocket": (
            "for BRAF the substrate is MEK1, a protein chain, so "
            "`to_substrate` is the distance to the nearest substrate-chain "
            "atom. That is an interface distance and not the same quantity as "
            "the distance to flurbiprofen in a CYP pocket. The two are not "
            "comparable across proteins and are not claimed to be."),
        "a_bound_ligand_is_not_every_substrate": (
            "flurbiprofen in 1R9O and 6-thio-GMP in 5LPG define one pocket "
            "each. A different substrate can occupy a different volume, so "
            "`to_substrate` is specific to the complex that was crystallised."),
        "numbering_is_verified_not_trusted": (
            "verify_catalytic_geometry fails unless a fact that depends on "
            "correct numbering holds, because an offset would produce "
            "plausible wrong distances rather than an error."),
    }
