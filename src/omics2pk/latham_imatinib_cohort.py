"""Reader for the Latham PHH table, written after an extraction bug cost a day.

What went wrong, kept here because the fix is the lesson
--------------------------------------------------------
The published table marks a missing value with a **dash**, and the dash it uses
is `U+2212 MINUS SIGN`, not the ASCII hyphen. A completeness check that tested
membership of a hand-written set of tokens reported three absent values as
present, and the conclusion drawn from that was that the article contradicted
itself: its text says protein quantitation came from 13 donors, its table
appeared to carry 14. **The article was right and the reader was wrong.**
Donor GKP has no CYP3A4, CYP3A5 or total CYP3A concentration.

So missingness here is decided by **Unicode category**, not by a list someone
has to remember to extend. Every character class that renders as a dash has
`unicodedata.category` of `Pd`, which covers the hyphen, the en and em dashes,
the figure dash and the non-breaking hyphen; `U+2212` is `Sm`, a maths symbol,
and is named explicitly beside them because it is the one that caused the bug.

And a missing value is **never** a zero. A concentration of zero is a
measurement; a dash is the absence of one. `parse_number` returns `None`, and
an unrecognised token **raises** rather than being quietly treated as missing,
because silently widening the missing set is how the first error happened.
"""
from __future__ import annotations

import csv
import unicodedata

import numpy as np

#: Tokens that mean "not measured". Dashes are detected by category as well,
#: so this list does not have to be exhaustive to be correct.
MISSING_WORDS = frozenset({"", "na", "n/a", "nd", "n.d.", "nr", "none",
                           "not measured", "not determined", "."})

#: The character that caused the bug, named so a reader sees why it is here.
MINUS_SIGN = "−"

#: Columns of Table 4, by role. Indices are checked against the header text.
COLUMNS = {
    "donor": "Donor ID",
    "age": "Age (yrs)",
    "sex": "Sex",
    "race": "Reported Race",
    "bmi": "BMI",
    "activity_cyp2c8": "CYP2C8 Activity",
    "activity_cyp3a_midazolam": "hydroxy-midazolam",
    "activity_cyp3a_testosterone": "hydroxy-testoster",
    "genotype_cyp2c8": "CYP2C8 Genotype",
    "genotype_cyp3a4": "CYP3A4 Genotype",
    "genotype_cyp3a5": "CYP3A5 Genotype",
    "metabolizer_cyp3a5": "CYP3A5 Metabolizer Status",
    "protein_cyp2c8": "CYP2C8 Conc",
    "protein_cyp3a4": "CYP3A4 Conc",
    "protein_cyp3a5": "CYP3A5 Conc",
    "protein_cyp3a_total": "CYP3A Conc",
    "clearance": "Imatinib Clearance",
    "metabolite": "desmethyl imatinib formation",
}

#: CYP3A5 phenotype, ordered a priori and never learned: poor < intermediate
#: < normal expressor. Fixing it in advance is what keeps it leak-free.
CYP3A5_ORDINAL = {"PM": 0.0, "IM": 1.0, "NM": 2.0}

#: Reference diplotype per gene. A non-reference indicator is 1 otherwise.
REFERENCE_DIPLOTYPE = {"genotype_cyp2c8": "*1/*1", "genotype_cyp3a4": "*1/*1"}


class LathamReadError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise LathamReadError(message)


def is_missing(value):
    """True when a cell means "not measured". Dashes count, zero does not."""
    if value is None:
        return True
    text = str(value).strip()
    if text.lower() in MISSING_WORDS:
        return True
    if not text:
        return True
    # A cell made only of dash-like characters is a missing marker. U+2212 is
    # categorised Sm rather than Pd, which is exactly why a hand-written token
    # list missed it.
    return all(unicodedata.category(ch) == "Pd" or ch == MINUS_SIGN
               for ch in text)


def parse_number(value, *, field=""):
    """Float, or None when missing. Never zero for a missing value.

    An unrecognised token raises. Treating anything unparseable as missing is
    how three absent concentrations were read as present.
    """
    if is_missing(value):
        return None
    text = str(value).strip().replace(",", "")
    # A published value may carry a footnote marker or a thousands separator;
    # anything else unparseable is an error, not a silent missing.
    stripped = text.rstrip("abcdefgh*†‡ ")
    try:
        return float(stripped)
    except ValueError:
        raise LathamReadError(
            f"unparseable value {value!r} in {field or 'a numeric column'}; "
            "declare it in MISSING_WORDS if it means 'not measured', but do "
            "not let it become a missing value by accident")


def load_table4(path):
    """Per-donor rows plus a completeness count for every block.

    Returns parsed columns and `completeness`, which is the number the protocol
    needs: n is per block, not per table.
    """
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle))
    _check(len(rows) >= 2, "the table has no data rows")
    headers, data = rows[0], [r for r in rows[1:] if any(c.strip() for c in r)]

    index = {}
    for role, fragment in COLUMNS.items():
        hits = [i for i, h in enumerate(headers)
                if fragment.lower() in h.lower()]
        _check(hits, f"no column matching {fragment!r} for role {role}")
        # CYP3A Conc also matches CYP3A4/CYP3A5; take the shortest header.
        index[role] = min(hits, key=lambda i: len(headers[i]))
    _check(len(set(index.values())) == len(index),
           f"two roles resolved to one column: {index}")

    numeric_roles = ("age", "bmi", "activity_cyp2c8",
                     "activity_cyp3a_midazolam",
                     "activity_cyp3a_testosterone", "protein_cyp2c8",
                     "protein_cyp3a4", "protein_cyp3a5",
                     "protein_cyp3a_total", "clearance", "metabolite")
    parsed = {"donor": [r[index["donor"]].strip() for r in data]}
    for role in numeric_roles:
        parsed[role] = [parse_number(r[index[role]], field=role) for r in data]
    for role in ("sex", "race", "genotype_cyp2c8", "genotype_cyp3a4",
                 "genotype_cyp3a5", "metabolizer_cyp3a5"):
        parsed[role] = [r[index[role]].strip() for r in data]

    _check(len(set(parsed["donor"])) == len(parsed["donor"]),
           "donor identifiers are not unique")

    blocks = {
        "clinical": ("age", "bmi"),
        "activities": ("activity_cyp2c8", "activity_cyp3a_midazolam",
                       "activity_cyp3a_testosterone"),
        "proteomics": ("protein_cyp2c8", "protein_cyp3a4", "protein_cyp3a5"),
        "clearance": ("clearance",),
        "metabolite": ("metabolite",),
    }
    completeness = {}
    for name, roles in blocks.items():
        complete = [parsed["donor"][i] for i in range(len(data))
                    if all(parsed[role][i] is not None for role in roles)]
        completeness[name] = {
            "n_complete": len(complete),
            "donors_missing": [d for d in parsed["donor"]
                               if d not in set(complete)],
        }
    return {
        "n_rows": len(data),
        "columns": index,
        "parsed": parsed,
        "completeness": completeness,
        "why_n_is_per_block": (
            "the published text says activity came from 14 donors and protein "
            "quantitation from 13. That is a statement about blocks, not about "
            "the table, and it is correct: one donor has no CYP3A4, CYP3A5 or "
            "total CYP3A concentration."),
    }


def primary_cohort(loaded):
    """Donors complete on the proteomic block AND the clearance endpoint."""
    parsed = loaded["parsed"]
    keep = [i for i in range(loaded["n_rows"])
            if all(parsed[role][i] is not None
                   for role in ("protein_cyp2c8", "protein_cyp3a4",
                                "protein_cyp3a5", "clearance", "age", "bmi"))]
    return keep, [parsed["donor"][i] for i in keep]


def encode_blocks(loaded, rows, *, train_rows):
    """Design blocks, with every learned quantity learned on `train_rows` only.

    What is fixed a priori and therefore leak-free: the CYP3A5 phenotype
    ordering, and which diplotype is the reference. What is learned in-fold:
    the centring and scaling, and nothing else.

    A category absent from the training rows produces a **constant** column
    rather than an error. At thirteen donors with a singleton `*1/*15` carrier
    and a singleton `NM`, the fold that holds that donor out has no variation
    left in its column, and ridge gives it a zero coefficient. That is the
    correct behaviour and it is why the encoding is not one-hot.
    """
    parsed = loaded["parsed"]
    rows = list(rows)
    train_rows = list(train_rows)
    _check(set(train_rows) <= set(rows), "training rows must be a subset")

    def numeric(role):
        return np.array([parsed[role][i] for i in rows], dtype=float)

    def standardise(values):
        train = np.array([parsed_value for index, parsed_value
                          in zip(rows, values) if index in set(train_rows)],
                         dtype=float)
        centre = float(train.mean())
        spread = float(train.std()) or 1.0
        return (values - centre) / spread

    clinical = np.column_stack([
        standardise(numeric("age")),
        standardise(numeric("bmi")),
        np.array([1.0 if parsed["sex"][i].upper().startswith("M") else 0.0
                  for i in rows]),
    ])
    pgx = np.column_stack([
        np.array([0.0 if parsed["genotype_cyp2c8"][i]
                  == REFERENCE_DIPLOTYPE["genotype_cyp2c8"] else 1.0
                  for i in rows]),
        np.array([0.0 if parsed["genotype_cyp3a4"][i]
                  == REFERENCE_DIPLOTYPE["genotype_cyp3a4"] else 1.0
                  for i in rows]),
        np.array([CYP3A5_ORDINAL.get(
            parsed["metabolizer_cyp3a5"][i].upper(), np.nan) for i in rows]),
    ])
    _check(np.isfinite(pgx).all(),
           "an unmapped CYP3A5 metaboliser status; declare it in "
           "CYP3A5_ORDINAL rather than letting it become a missing value")
    proteomics = np.column_stack([
        standardise(np.log(numeric(role)))
        for role in ("protein_cyp2c8", "protein_cyp3a4", "protein_cyp3a5")
    ])
    return {"clinical": clinical, "pgx": pgx, "proteomics": proteomics,
            "learned_in_fold": ("centring and scaling only",),
            "fixed_a_priori": ("CYP3A5 PM<IM<NM ordering",
                               "the reference diplotype per gene")}


def what_this_module_does_not_do():
    return {
        "a_dash_is_not_a_zero": (
            "U+2212 MINUS SIGN marks three absent concentrations in this "
            "table. A completeness check built on a hand-written token list "
            "read them as present and produced the claim that the article "
            "contradicted itself. Missingness is decided by Unicode category "
            "here, and an unrecognised token raises."),
        "n_is_per_block_not_per_table": (
            "activity comes from 14 donors and protein quantitation from 13. "
            "Reporting one n for the dataset would hide which contrast rests "
            "on which donors."),
        "an_in_vitro_clearance_is_not_a_clearance": (
            "the endpoint is hepatocyte intrinsic clearance, measured with no "
            "absorption, no plasma binding, no hepatic blood flow and no "
            "renal excretion."),
        "a_singleton_category_is_not_an_effect": (
            "one *1/*15 carrier and one normal CYP3A5 expressor. In the fold "
            "that holds either out, the column is constant in training. The "
            "encoding tolerates that by design; it does not make the category "
            "estimable."),
    }
