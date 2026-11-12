"""The catalogue: what can be measured, what it would answer, and what it is not.

Why this module exists
----------------------
The programme's rule is to start from a pharmacokinetic question, identify the
mechanisms, decide what information is needed, and only then choose a
statistical method. Nothing enforced that order. This module is the enforcing
part: a structured vocabulary of inputs, mechanisms, sources and hypotheses,
with the three confusions that would quietly invalidate a matrix built from
them written in as refusals rather than as warnings.

The three refusals
------------------
1. **A measured activity is not a pharmacokinetic parameter.** An enzyme
   activity in microsomes is a proxy for *metabolic* clearance, never for
   total clearance, and never for bioavailability, volume or absorption.
   `PK_PARAMETERS` and `FUNCTIONAL_ENDPOINTS` are disjoint, and a hypothesis
   that requires one refuses the other.
2. **Pairing is asserted, never assumed.** Two layers measured in different
   publications are not paired, whatever they are about. A dataset must
   declare `same_individuals` or `same_variants` before any hypothesis
   needing joint data is called evaluable.
3. **The legacy objectives H1-H5 are not this programme's H1-H5.** The
   registry's labels mean different things -- legacy H4 is "when omics become
   harmful", the new H4 is longitudinal intra-individual variability -- and
   reusing the tags would claim capabilities that do not exist. The crosswalk
   is explicit and mixing the two vocabularies raises.

What a source is, and is not
----------------------------
A knowledge base, an annotation service and an integration method are three
different things, and calling DIABLO a data source would put a statistical
choice where a biological fact belongs. `SOURCES` carries a `kind` for exactly
that reason.
"""
from __future__ import annotations

#: Pharmacokinetic parameters. A hypothesis about PK needs one of these.
PK_PARAMETERS = {
    "CL": "total clearance",
    "CL_metabolic": "metabolic clearance, a COMPONENT of CL and not CL",
    "F": "oral bioavailability",
    "V": "volume of distribution",
    "ka": "absorption rate",
    "AUC": "exposure",
    "Cmax": "peak concentration",
    "Cmin": "trough concentration",
    "concentration": "a measured concentration at a stated time",
    "half_life": "terminal half-life",
}

#: Functional and mechanistic endpoints. Informative about PK, never PK.
FUNCTIONAL_ENDPOINTS = {
    "enzyme_activity": "probe turnover by an enzyme, in vitro",
    "clint_in_vitro": ("intrinsic clearance in microsomes or hepatocytes. A "
                       "component of metabolic clearance measured outside a "
                       "body, with no absorption, binding or blood flow"),
    "transporter_activity": "uptake or efflux rate",
    "protein_abundance": "quantified protein concentration",
    "transcript_abundance": "quantified transcript level",
    "metabolite_ratio": "a metabolic ratio used as a phenotype proxy",
    "variant_functional_score": "a deep-mutational-scanning score",
    "binding": "plasma protein or tissue binding, fraction unbound",
}

assert not (set(PK_PARAMETERS) & set(FUNCTIONAL_ENDPOINTS)), \
    "a parameter cannot be both a PK quantity and a functional endpoint"

#: A. Hospital and clinical inputs, by what they act on mechanistically.
CLINICAL_INPUTS = {
    "age": ("covariate", "maturation and decline of enzyme and renal function"),
    "sex": ("covariate", "enzyme expression differences, body composition"),
    "weight": ("covariate", "V and allometric scaling of CL"),
    "renal_function": ("organ function", "renal CL of parent and metabolites"),
    "hepatic_function": ("organ function", "metabolic CL and first pass"),
    "inflammation": ("state", "CYP down-regulation, acute-phase binding"),
    "albumin": ("binding", "fraction unbound, hence CL and V of bound drugs"),
    "hematocrit": ("matrix", "blood-to-plasma ratio, assay interpretation"),
    "comedication": ("perturbation", "induction, inhibition, transport"),
    "ddi": ("perturbation", "a named interaction with a mechanism"),
    "diet": ("perturbation", "CYP1A2 and CYP3A modulation, absorption"),
    "dose": ("design", "linearity, saturation"),
    "formulation": ("design", "ka and F"),
    "timing": ("design", "what a concentration means"),
    "adherence": ("design", "whether the dose history is real"),
    "concentration": ("endpoint", "the PK observation itself"),
}

#: B. Biological inputs, with the PK parameter each can plausibly inform.
BIOLOGICAL_INPUTS = {
    "snp": ("genomic", ("CL_metabolic", "F")),
    "diplotype": ("genomic", ("CL_metabolic", "F")),
    "rare_variant": ("genomic", ("CL_metabolic",)),
    "transcript": ("transcriptomic", ("CL_metabolic", "F")),
    "splicing": ("transcriptomic", ("CL_metabolic",)),
    "protein": ("proteomic", ("CL_metabolic", "F", "V")),
    "ptm": ("proteomic", ("CL_metabolic",)),
    "metabolite": ("metabolomic", ("CL_metabolic", "AUC")),
    "lipid": ("lipidomic", ("V", "F")),
    "enzyme_activity": ("functional", ("CL_metabolic",)),
    "transporter_activity": ("functional", ("F", "CL", "V")),
    "binding_protein": ("proteomic", ("V", "CL")),
}

#: C. The context without which a measurement cannot be interpreted.
MEASUREMENT_CONTEXT = (
    "organ", "tissue", "compartment", "cell_type", "sampling_time",
    "experimental_condition", "drug_or_substrate", "analytical_method",
    "quantification", "temporal_stability",
)

#: D. Mechanistic knowledge classes. Deliberately not CYP-only.
MECHANISTIC_KNOWLEDGE = {
    "phase_i_enzymes": ("CYP, FMO, ADH, ALDH, esterases", ("CL_metabolic",)),
    "phase_ii_enzymes": ("UGT, SULT, NAT, GST, TPMT", ("CL_metabolic",)),
    "efflux_transporters": ("ABCB1, ABCG2, ABCC family", ("F", "CL", "V")),
    "uptake_transporters": ("SLCO, SLC22, SLC47 family", ("F", "CL", "V")),
    "binding_proteins": ("albumin, AAG, lipoproteins", ("V", "CL")),
    "redox_partners": ("POR, cytochrome b5", ("CL_metabolic",)),
    "nuclear_receptors": ("PXR, CAR, AhR, induction cascades", ("CL_metabolic",)),
    "absorption_machinery": ("gut transporters, gut wall CYP3A", ("F", "ka")),
    "targets": ("the pharmacological target, for PD not PK", ()),
    "complexes_and_pathways": ("protein complexes, signalling cascades", ()),
}

#: E. Sources, separated by KIND. A method is not a source of biology.
SOURCES = {
    "STRING": ("database", "protein interaction confidence, NOT causality"),
    "Reactome": ("database", "curated pathways"),
    "KEGG": ("database", "curated pathways and drug maps"),
    "UniProt": ("database", "protein annotation, functional sites"),
    "PharmVar": ("database", "star-allele definitions"),
    "ClinPGx": ("database", "pharmacogenomic annotation, successor to PharmGKB"),
    "CPIC": ("database", "genotype-to-phenotype guidelines"),
    "RCSB": ("database", "experimental structures and ligands"),
    "MaveDB": ("database", "deep-mutational-scanning scores"),
    "GEO": ("repository", "transcriptomic series, indexed by ASSAY not phenotype"),
    "PRIDE": ("repository", "proteomic deposits, often beside DMPK papers"),
    "MetaboLights": ("repository", "metabolomic studies"),
    "MetabolomicsWorkbench": ("repository", "metabolomic studies"),
    "EnsemblVEP": ("annotation_tool", "variant consequence prediction"),
    "DIABLO": ("method", "multi-block integration. NOT a biological source"),
    "mixOmics": ("method", "the package DIABLO belongs to"),
    "ElasticNet": ("method", "penalised selection"),
    "XGBoost_SHAP": ("method", "boosted trees plus an attribution"),
    "Spearman": ("method", "rank correlation, the baseline selector"),
}

#: This programme's hypotheses, with what each one REQUIRES to be evaluable.
HYPOTHESES = {
    "H1": {
        "question": ("do omics add information beyond PGx, clinical and "
                     "routine biology for PK variability?"),
        "priority": "primary",
        "endpoint_kind": "pk",
        "minimum_comparator": "clinical + biology + PGx, same individuals",
        "needs_pairing": "same_individuals",
        "needs_abstention_arm": True,
    },
    "H2": {
        "question": ("can the informative omics be selected automatically, "
                     "ideally BEFORE acquisition?"),
        "priority": "primary",
        "endpoint_kind": "pk_or_functional",
        "minimum_comparator": "Spearman or elastic net on the same pool",
        "needs_pairing": "same_individuals_or_variants",
        "needs_abstention_arm": True,
        # Choosing two candidates out of two is not selection: no selector can
        # distinguish itself on a pool that size, and the near-tie it produces
        # says nothing. Established on ARIANA, which has exactly two.
        "needs_candidate_pool_of_at_least": 5,
    },
    "H3": {
        "question": ("do omics predict the biological functions involved in "
                     "PK better than PGx alone? NOT limited to enzymes: "
                     "transporters, binding proteins, absorption"),
        "priority": "intermediate",
        "endpoint_kind": "functional",
        "minimum_comparator": "PGx alone, same individuals or variants",
        "needs_pairing": "same_individuals_or_variants",
        "needs_abstention_arm": True,
    },
    "H4": {
        "question": ("do omic changes explain INTRA-individual PK "
                     "variability between occasions?"),
        "priority": "longitudinal",
        "endpoint_kind": "pk",
        "minimum_comparator": "within-individual change in clinical state",
        "needs_pairing": "same_individuals_repeated",
        "needs_abstention_arm": True,
    },
    "H5": {
        "question": "does omic integration improve dose personalisation?",
        "priority": "final_application",
        "endpoint_kind": "pk",
        "minimum_comparator": "the dosing rule already in use",
        "needs_pairing": "same_individuals",
        "needs_abstention_arm": True,
    },
}

#: The registry's OWN H1-H5, which mean different things. Never interchange.
LEGACY_OBJECTIVES = {
    "H1": "incremental value of omics (ARIANA)",
    "H2": "reliability of mechanism-guided selection",
    "H3": "parsimony and minimal information",
    "H4": "when omics become HARMFUL",
    "H5": "biological, analytical and temporal robustness",
}

#: Where the two vocabularies are near enough to be confused, and where not.
LEGACY_CROSSWALK = {
    "H1": ("H1", "compatible: the legacy objective is one drug's instance of "
                 "the new question"),
    "H2": ("H2", "compatible in spirit; the new H2 adds the prospective, "
                 "before-acquisition task the legacy one never had"),
    "H3": (None, "legacy H3 is parsimony, a PROPERTY of a selection. The new "
                 "H3 is a question about functional endpoints. No mapping"),
    "H4": (None, "legacy H4 is when omics become harmful. The new H4 is "
                 "longitudinal. No mapping, and reusing the tag would claim "
                 "longitudinal data that does not exist"),
    "H5": (None, "legacy H5 is analytical robustness. The new H5 is dose "
                 "personalisation. No mapping"),
}

#: The four mechanistic axes the programme is short of data on, named so that
#: coverage can be counted rather than asserted. Each says which PK parameters
#: it informs, and metabolic clearance is deliberately absent from three of
#: them: that is the point of listing them.
AXES = {
    "metabolic_clearance": {
        "mechanisms": ("phase_i_enzymes", "phase_ii_enzymes",
                       "redox_partners"),
        "informs": ("CL_metabolic",),
        "note": ("the only axis with data. Four experiments and counting on "
                 "CYP2C8, CYP3A4, CYP3A5, NUDT15 and two RAF kinases"),
    },
    "transport": {
        "mechanisms": ("efflux_transporters", "uptake_transporters"),
        "informs": ("F", "CL", "V"),
        "note": ("a transporter abundance does not inform F or ka by itself: "
                 "the interpretation depends on the mechanism and on the PK "
                 "model that uses it"),
    },
    "binding_distribution": {
        "mechanisms": ("binding_proteins",),
        "informs": ("V", "CL"),
        "note": "fraction unbound, albumin, AAG, tissue distribution",
    },
    "absorption": {
        "mechanisms": ("absorption_machinery",),
        "informs": ("F", "ka"),
        "note": "intestinal expression, gut-wall first pass, gut transporters",
    },
    "temporal": {
        "mechanisms": ("nuclear_receptors",),
        "informs": ("CL", "F", "V"),
        "note": ("inflammation, infection, transplantation, comedication. "
                 "Requires repeated occasions, not a between-group contrast"),
    },
}

#: The four statuses a dataset can have against one hypothesis.
STATUSES = ("READY", "PARTIAL", "ACCESS_REQUIRED", "NOT_EVALUABLE")

#: Admissible pairing claims. Anything else makes joint analysis inadmissible.
PAIRINGS = ("same_individuals", "same_individuals_repeated", "same_variants",
            "different_publications", "unknown")


class CatalogueError(ValueError):
    pass


def _check(condition, message):
    if not condition:
        raise CatalogueError(message)


def endpoint_kind(endpoint):
    """`pk`, `functional`, or a refusal. The first refusal, enforced."""
    if endpoint in PK_PARAMETERS:
        return "pk"
    if endpoint in FUNCTIONAL_ENDPOINTS:
        return "functional"
    raise CatalogueError(
        f"unknown endpoint {endpoint!r}; declare it in PK_PARAMETERS or in "
        "FUNCTIONAL_ENDPOINTS, which are disjoint on purpose")


def assert_not_a_pk_parameter(endpoint):
    """Raise if a functional measurement is being treated as a PK parameter.

    `clint_in_vitro` is the one that tempts: it is called a clearance, it has
    clearance units, and it is measured without absorption, binding or blood
    flow. It is a component of metabolic clearance in a tube, and metabolic
    clearance is itself only a component of CL.
    """
    kind = endpoint_kind(endpoint)
    _check(kind == "functional",
           f"{endpoint!r} is a PK parameter; this guard is for functional "
           "measurements being passed off as one")
    return {"endpoint": endpoint, "is_pk_parameter": False,
            "what_it_is": FUNCTIONAL_ENDPOINTS[endpoint],
            "may_inform": sorted({
                pk for _, pks in BIOLOGICAL_INPUTS.values() for pk in pks
                if endpoint in BIOLOGICAL_INPUTS}) or ["CL_metabolic"]}


def source_kind(name):
    _check(name in SOURCES, f"unknown source {name!r}")
    return SOURCES[name][0]


def assert_is_a_database(name):
    """Refuse to treat a method or a tool as a source of biology."""
    kind = source_kind(name)
    _check(kind in ("database", "repository"),
           f"{name} is a {kind}, not a source of biological information. "
           "Calling it one puts a statistical choice where a biological fact "
           "belongs")
    return True


def assert_hypothesis_vocabulary(label, *, programme):
    """Keep the two H1-H5 vocabularies apart.

    `programme` is `"current"` or `"legacy"`. Asking for a label without
    saying which vocabulary it belongs to is the mistake this prevents.
    """
    _check(programme in ("current", "legacy"),
           "say which vocabulary: 'current' or 'legacy'")
    table = HYPOTHESES if programme == "current" else LEGACY_OBJECTIVES
    _check(label in table, f"{label} is not in the {programme} vocabulary")
    if programme == "legacy":
        mapped, why = LEGACY_CROSSWALK[label]
        return {"label": label, "programme": "legacy",
                "means": LEGACY_OBJECTIVES[label],
                "maps_to_current": mapped, "why": why}
    return {"label": label, "programme": "current",
            "question": HYPOTHESES[label]["question"]}


def feasibility(hypothesis, *, endpoint, pairing, access, layers_available,
                comparator_available, n_units, longitudinal=False,
                n_candidate_omics=None):
    """Status of one dataset against one hypothesis, computed not asserted.

    Returns `(status, reasons)`. The reasons are the point: a `NOT_EVALUABLE`
    that does not say which requirement failed cannot be acted on.

    **`comparator_available` is per hypothesis, not per dataset.** One flag
    cannot answer five different questions: ARIANA has the clinical and PGx
    baseline H1 needs and no dosing rule for H5 to beat, and a single boolean
    reported it as ready for both.

    **`n_candidate_omics` guards H2.** Choosing two candidates out of two is
    not selection.
    """
    _check(hypothesis in HYPOTHESES, f"unknown hypothesis {hypothesis!r}")
    _check(pairing in PAIRINGS, f"unknown pairing {pairing!r}")
    spec = HYPOTHESES[hypothesis]
    reasons = []

    kind = endpoint_kind(endpoint)
    wanted = spec["endpoint_kind"]
    if wanted == "pk" and kind != "pk":
        reasons.append(
            f"{hypothesis} needs a PK parameter; {endpoint} is a functional "
            "measurement and must not be read as one")
    elif wanted == "functional" and kind != "functional":
        reasons.append(
            f"{hypothesis} is about functional endpoints; {endpoint} is a PK "
            "parameter, which belongs to H1, H4 or H5")

    need = spec["needs_pairing"]
    if need == "same_individuals" and pairing != "same_individuals":
        reasons.append(f"needs same_individuals; declared {pairing}")
    elif need == "same_individuals_repeated" and \
            pairing != "same_individuals_repeated":
        reasons.append(f"needs repeated measures; declared {pairing}")
    elif need == "same_individuals_or_variants" and pairing not in (
            "same_individuals", "same_individuals_repeated", "same_variants"):
        reasons.append(f"needs per-individual or per-variant pairing; "
                       f"declared {pairing}")
    if pairing in ("different_publications", "unknown"):
        reasons.append("joint analysis is inadmissible: two layers from "
                       "different publications are not paired")

    if hypothesis == "H4" and not longitudinal:
        reasons.append("H4 needs repeated occasions; this dataset has none")
    if not layers_available:
        reasons.append("no omic layer is available")
    if not comparator_available:
        reasons.append(f"the minimum comparator is absent: "
                       f"{spec['minimum_comparator']}")
    floor = spec.get("needs_candidate_pool_of_at_least")
    if floor is not None and n_candidate_omics is not None \
            and n_candidate_omics < floor:
        reasons.append(
            f"only {n_candidate_omics} candidate omic variables; selection "
            f"needs at least {floor}. Choosing {n_candidate_omics} out of "
            f"{n_candidate_omics} is not selection")
    if n_units is not None and n_units < 20:
        reasons.append(f"only {n_units} units, below any useful power")

    if access == "unavailable":
        return "NOT_EVALUABLE", reasons + ["the data are not obtainable"]
    if reasons:
        return "NOT_EVALUABLE", reasons
    if access in ("author_request", "governed_application"):
        return "ACCESS_REQUIRED", ["every scientific requirement is met; "
                                   f"access is {access}"]
    if access == "public_registered":
        return "PARTIAL", ["public but behind a registration, and the brief "
                           "admits only channels needing no authorisation"]
    return "READY", ["every requirement met and the channel is open"]


def axis_of(mechanism):
    """Which mechanistic axis a knowledge class belongs to."""
    _check(mechanism in MECHANISTIC_KNOWLEDGE,
           f"unknown mechanism {mechanism!r}")
    for axis, spec in AXES.items():
        if mechanism in spec["mechanisms"]:
            return axis
    return None


def axis_coverage(rows):
    """How many datasets reach each axis, and with what kind of endpoint.

    Counting rather than asserting, because "we have no transporter data" is a
    claim that should be re-derivable from the matrix at any time.
    """
    out = {}
    for axis, spec in AXES.items():
        here = [r for r in rows if axis in (r.get("axes") or ())]
        # Rows are hypothesis x dataset, so a dataset appears once per
        # hypothesis. Counting rows would multiply every dataset by five and
        # report a coverage that does not exist.
        datasets = {r["dataset"] for r in here}
        with_pk = {r["dataset"] for r in here
                   if r.get("endpoint_kind") == "pk"}
        ready = {r["dataset"] for r in here if r.get("status") == "READY"}
        out[axis] = {
            "informs": list(spec["informs"]),
            "n_datasets": len(datasets),
            "n_with_a_pk_endpoint": len(with_pk),
            "n_ready_for_any_hypothesis": len(ready),
            "datasets": sorted(datasets),
            "note": spec["note"],
        }
    return out


def what_this_module_does_not_do():
    return {
        "a_catalogue_is_not_a_result": (
            "listing what could be measured says nothing about whether it "
            "carries information. Every entry here is a candidate until a "
            "pre-registered contrast says otherwise."),
        "an_activity_is_not_a_pk_parameter": (
            "clint_in_vitro has clearance units and is measured with no "
            "absorption, no binding and no blood flow. It informs metabolic "
            "clearance, which is itself a component of CL."),
        "pairing_is_never_inferred": (
            "two layers reported in different publications are not paired, "
            "however similar their cohorts. The pairing field must be "
            "asserted from the source, and `unknown` blocks joint analysis."),
        "the_legacy_labels_are_a_different_vocabulary": (
            "the registry's H4 is 'when omics become harmful' and its H5 is "
            "analytical robustness. Reusing those tags for the longitudinal "
            "and dosing hypotheses would claim data that does not exist."),
        "a_method_is_not_a_source": (
            "DIABLO integrates blocks; it supplies no biology. The kind field "
            "exists so a matrix cannot put a statistical choice where a "
            "biological fact belongs."),
    }
