#!/usr/bin/env python
"""Hypothesis x Dataset x Endpoint x Comparator x Feasibility, computed.

    python -I scripts/build_hypothesis_dataset_matrix.py

No experiment is run. This builds the Phase 1 and Phase 2 deliverables: the
catalogue of inputs, mechanisms and sources, and the matrix saying which
hypothesis each dataset can actually carry.

Two things it refuses to do
---------------------------
**It does not infer pairing.** The registry carries no pairing field, and
whether two layers were measured in the same individuals is exactly the thing
that must never be assumed. Every dataset below declares its pairing with a
`provenance` string saying where that claim comes from; anything unsubstantiated
is `unknown`, which makes joint analysis inadmissible and says so. A matrix
full of `unknown` is not a failure -- it is a list of what to verify.

**It does not reuse the registry's H1-H5 tags.** Those are the legacy
objectives and they mean different things: legacy H4 is "when omics become
harmful", this programme's H4 is longitudinal intra-individual variability.
Thirteen datasets carry a legacy `H5` tag; mapping them onto dose
personalisation would claim data that does not exist. Every classification
here is computed from the dataset's contents.
"""
from __future__ import annotations

import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from omics2pk import input_catalogue as CAT  # noqa: E402

REGISTRY = ROOT / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"

#: What each dataset actually offers, declared with where the claim comes from.
#: `endpoint` must be a key of PK_PARAMETERS or FUNCTIONAL_ENDPOINTS, so a
#: functional measurement cannot be smuggled in as a PK parameter.
DECLARED = {
    "ariana_metoprolol_47": {
        "drug": "metoprolol",
        "endpoint": "concentration",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("metabolomics_targeted",),
        "comparator_available": {"H1": True, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": 2,
        "provenance": ("H1_ARIANA_INCREMENTAL_OMIC_VALUE_RESULTS_20261107.md: "
                       "47 patients with a per-patient concentration and a "
                       "per-patient targeted metabolite panel"),
        "note": ("the frozen H1 contrast of this project. Two candidate omic "
                 "variables only, which is why it cannot carry selection"),
    },
    "hlc_syn4499": {
        "drug": "probe substrates, ten CYPs",
        "endpoint": "enzyme_activity",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("transcriptomics", "genotype"),
        "comparator_available": {"H1": False, "H2": True, "H3": True, "H4": False, "H5": False},
        "n_candidate_omics": 31,
        "provenance": ("HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS: "
                       "~390 donors with transcript, activity and two "
                       "source-verified SNPs on the same donor ids"),
        "note": ("genotype covers 171-172 donors only, so PGx enters as a "
                 "stratum. Access is Synapse-registered"),
    },
    "wegler_15": {
        "drug": "not applicable",
        "endpoint": "protein_abundance",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("proteomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "provenance": "registry entry wegler_15, n=15 liver donors",
        "note": "proteome without a paired functional endpoint in the deposit",
    },
    "lee_jackson_cyp3a_15": {
        "drug": "midazolam and testosterone probes",
        "endpoint": "enzyme_activity",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("proteomics", "targeted_metabolomics"),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "provenance": "registry entry lee_jackson_cyp3a_15, n=15",
        "note": "no genotype in the source, so the PGx comparator is absent",
    },
    "lee_venetoclax_patients_12": {
        "drug": "venetoclax",
        "endpoint": "concentration",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("targeted_metabolomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "provenance": "registry entry lee_venetoclax_patients_12, n=12",
        "note": "twelve patients; below any useful power for a gain contrast",
    },
    "jarmusch_midazolam": {
        "drug": "midazolam",
        "endpoint": "concentration",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("metabolomics_untargeted",),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "provenance": ("OMICS2PK_FIVE_OBJECTIVES_VALIDATION_PROTOCOL.md "
                       "v1.3.0: the GNPS table joined per donor gives 13 "
                       "matched subjects at plasma / day 1 / time 0"),
        "note": ("13 matched subjects. A GNPS-named feature is not a "
                 "confirmed structure"),
    },
    "epa_liver_100": {
        "drug": "eight CYP probes",
        "endpoint": "enzyme_activity",
        "pairing": "unknown",
        "longitudinal": False,
        "layers": (),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "provenance": ("no archive is present on this machine; the readers "
                       "exist but a reader without data is not a cohort"),
        "note": "NOT_EVALUABLE for want of local data, recorded 2026-11-11",
    },
    "mavedb_nudt15_paired_stability_activity": {
        "drug": "thiopurine",
        "endpoint": "variant_functional_score",
        "pairing": "same_variants",
        "longitudinal": False,
        "layers": ("proteomic", "functional"),
        "comparator_available": {"H1": False, "H2": True, "H3": True, "H4": False, "H5": False},
        "n_candidate_omics": 21,
        "provenance": "registry entry, same variant library in both assays",
        "note": "a screen score is not a PK quantity and not an individual",
    },
    "mavedb_cyp2c9_paired_abundance_activity": {
        "drug": "not applicable",
        "endpoint": "variant_functional_score",
        "pairing": "same_variants",
        "longitudinal": False,
        "layers": ("proteomic", "functional"),
        "comparator_available": {"H1": False, "H2": True, "H3": True, "H4": False, "H5": False},
        "n_candidate_omics": 21,
        "provenance": "registry entry, same variant library in both assays",
        "note": "abundance in HEK293T, activity in humanised yeast",
    },
    "mavedb_braf_paired_abundance_activity": {
        "drug": "not applicable",
        "endpoint": "variant_functional_score",
        "pairing": "same_variants",
        "longitudinal": False,
        "layers": ("proteomic", "functional"),
        "comparator_available": {"H1": False, "H2": True, "H3": True, "H4": False, "H5": False},
        "n_candidate_omics": 21,
        "provenance": "registry entry, one C-terminal construct",
        "note": "gain-of-function library; the structural rule is refuted here",
    },
    "mavedb_araf_abundance_conditions_and_activity": {
        "drug": "not applicable",
        "endpoint": "variant_functional_score",
        "pairing": "same_variants",
        "longitudinal": False,
        "layers": ("proteomic", "functional"),
        "comparator_available": {"H1": False, "H2": True, "H3": True, "H4": False, "H5": False},
        "n_candidate_omics": 21,
        "provenance": "registry entry, one C-terminal construct",
        "note": ("two abundance CONDITIONS, not two omic layers. The "
                 "structural rule is refuted here too"),
    },
}

#: Candidates found today in public repositories. None is downloaded.
CANDIDATES = {
    "pride_pxd066482_imatinib_clearance": {
        "label": ("PRIDE PXD066482 - interindividual variability in imatinib "
                  "metabolism in human liver microsomes and primary human "
                  "hepatocytes, CYP2C8 and CYP3A phenotypes"),
        "repository": "PRIDE",
        "accession": "PXD066482",
        "licence": "CC0",
        "drug": "imatinib",
        "endpoint": "clint_in_vitro",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("proteomics", "functional", "genotype"),
        "comparator_available": {"H1": False, "H2": True, "H3": True, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "access": "open_public",
        "n_units": None,
        "provenance": ("PRIDE project description, retrieved 2026-11-12: "
                       "individual-donor HLM and PHH incubated with imatinib "
                       "for clearance and N-desmethyl-imatinib formation, "
                       "with CYP2C8 and CYP3A activity, protein concentration "
                       "AND genotype on the same donors"),
        "why_it_matters": ("the only inventoried source carrying genotype, "
                           "protein, activity and a clearance on the same "
                           "units. That is the H3 design, and the PGx "
                           "comparator H3 needs is present"),
        "what_is_missing": ("the PRIDE deposit holds the proteomics; the "
                            "clearance, activity and genotype tables are in "
                            "the publication's supplement and are NOT yet "
                            "verified as public. n per donor is not stated in "
                            "the deposit metadata"),
    },
    "pride_pxd057282_tucatinib_clearance": {
        "label": ("PRIDE PXD057282 - impact of variation in CYP3A and CYP2C8 "
                  "on tucatinib metabolic clearance in human liver "
                  "microsomes"),
        "repository": "PRIDE",
        "accession": "PXD057282",
        "licence": "CC0",
        "drug": "tucatinib",
        "endpoint": "clint_in_vitro",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("proteomics", "functional"),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "access": "open_public",
        "n_units": 21,
        "provenance": ("PRIDE project description, retrieved 2026-11-12: 21 "
                       "donors, 11 male and 10 female, tucatinib depletion "
                       "plus CYP2C8, CYP3A4 and CYP3A5 concentrations by QTAP "
                       "and measured CYP activities"),
        "why_it_matters": ("same design as PXD066482 with a different drug, so "
                           "it is a replication candidate rather than a "
                           "second discovery set"),
        "what_is_missing": ("no genotype stated, so the PGx comparator H3 "
                            "needs is absent. Commercial microsomes, so no "
                            "clinical covariates beyond sex"),
    },
    "pride_pxd060653_alcohol_dme_ensemble": {
        "label": ("PRIDE PXD060653 - effects of chronic alcohol intake on the "
                  "ensemble of drug-metabolising enzymes"),
        "repository": "PRIDE", "accession": "PXD060653", "licence": None,
        "drug": "not applicable",
        "endpoint": "protein_abundance",
        "pairing": "unknown",
        "longitudinal": False,
        "layers": ("proteomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "access": "open_public",
        "n_units": None,
        "provenance": "PRIDE search, retrieved 2026-11-12",
        "why_it_matters": ("a PERTURBATION of the enzyme proteome, which is "
                           "the shape H4 needs even though this design is "
                           "between groups rather than within individuals"),
        "what_is_missing": ("no paired PK, no repeated occasions, pairing not "
                            "asserted. Cannot carry H4 as it stands"),
    },
    "pride_pxd014131_liver_compartments": {
        "label": ("PRIDE PXD014131 - human liver homogenate, microsomes, "
                  "low-speed pellet and hepatocytes"),
        "repository": "PRIDE", "accession": "PXD014131", "licence": None,
        "drug": "not applicable",
        "endpoint": "protein_abundance",
        "pairing": "unknown",
        "longitudinal": False,
        "layers": ("proteomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False, "H4": False, "H5": False},
        "n_candidate_omics": None,
        "access": "open_public",
        "n_units": None,
        "provenance": "PRIDE search, retrieved 2026-11-12",
        "why_it_matters": ("four COMPARTMENTS of the same organ, which is the "
                           "'which compartment to measure' axis of the brief "
                           "and the only inventoried source addressing it"),
        "what_is_missing": "no functional endpoint and no PK; pairing unstated",
    },
}

#: Found in the autonomous search of 2026-11-12, across the four named gaps.
#: Nothing is downloaded. Every `axes` tag is from the deposit's own
#: description, and every endpoint is the strongest one the deposit carries.
CANDIDATES_2026_11_12 = {
    "pride_pxd020987_human_small_intestine": {
        "label": ("PRIDE PXD020987 - quantitative proteomics of human small "
                  "intestine: 29 enzymes and 32 transporters in ileum and "
                  "jejunum mucosa from 16 surgical patients"),
        "repository": "PRIDE", "accession": "PXD020987", "licence": "CC0",
        "drug": "not applicable",
        "endpoint": "protein_abundance",
        "pairing": "same_individuals",
        "longitudinal": False,
        "layers": ("proteomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False,
                                 "H4": False, "H5": False},
        "n_candidate_omics": 61,
        "axes": ("absorption", "transport"),
        "access": "open_public",
        "n_units": 16,
        "provenance": ("PRIDE project description, retrieved 2026-11-12: "
                       "ileum and jejunum biopsies from 16 patients, mucosal "
                       "fractions, 29 enzymes, 32 transporters, 6 tight "
                       "junction proteins, alkaline phosphatase and markers "
                       "quantified, 60 for the first time"),
        "why_it_matters": ("the first inventoried source on the absorption "
                           "and transport axes at all, and the only one "
                           "carrying TWO intestinal regions per patient, "
                           "which is the compartment question the brief asks"),
        "what_is_missing": ("no PK, no clearance, no PGx and no independent "
                            "functional endpoint. An abundance predicting "
                            "another abundance is axis A, not H3. This is a "
                            "REFERENCE layer for what intestinal abundance "
                            "looks like, not a hypothesis-testing cohort"),
    },
    "pride_pxd011498_jejunum_before_after_gastrectomy": {
        "label": ("PRIDE PXD011498 - peptidomics of human jejunum biopsies "
                  "before and after gastrectomy"),
        "repository": "PRIDE", "accession": "PXD011498", "licence": None,
        "drug": "not applicable",
        "endpoint": "protein_abundance",
        "pairing": "same_individuals_repeated",
        "longitudinal": True,
        "layers": ("proteomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False,
                                 "H4": False, "H5": False},
        "n_candidate_omics": None,
        "axes": ("absorption", "temporal"),
        "access": "open_public",
        "n_units": None,
        "provenance": "PRIDE search, retrieved 2026-11-12",
        "why_it_matters": ("the only inventoried source with a genuine "
                           "WITHIN-SUBJECT repeated design on an axis the "
                           "programme needs. H4 has had no dataset with that "
                           "shape at all"),
        "what_is_missing": ("peptidomics rather than transporter "
                            "quantification, and no PK on either occasion, "
                            "so the design is right and the measurements are "
                            "the wrong ones"),
    },
    "epmc_pmc13683814_placental_transporters": {
        "label": ("Europe PMC PMC13683814 - proteomic profiling of xenobiotic "
                  "and nutrient transporters in human placenta"),
        "repository": "EuropePMC", "accession": "PMC13683814",
        "licence": "open access",
        "drug": "not applicable",
        "endpoint": "protein_abundance",
        "pairing": "unknown",
        "longitudinal": False,
        "layers": ("proteomics",),
        "comparator_available": {"H1": False, "H2": False, "H3": False,
                                 "H4": False, "H5": False},
        "n_candidate_omics": None,
        "axes": ("transport", "binding_distribution"),
        "access": "open_public",
        "n_units": None,
        "provenance": "Europe PMC search, retrieved 2026-11-12",
        "why_it_matters": ("placental transport is a distribution barrier, so "
                           "it touches V in a compartment nothing else in the "
                           "inventory reaches"),
        "what_is_missing": ("pairing not asserted, no PK, and a placental "
                            "barrier informs fetal exposure rather than the "
                            "maternal parameters the programme predicts"),
    },
}

#: Searches run today, with what they returned. Negative results included.
SEARCHES = [
    {"repository": "GEO", "query": "drug clearance AND expression profiling",
     "hits": 2, "on_target": 0,
     "verdict": ("both hits are off target: a rat Bcrp/Pgp study and a "
                 "vitamin B2 / AhR study")},
    {"repository": "GEO",
     "query": "liver AND cytochrome P450 activity AND Homo sapiens",
     "hits": 2, "on_target": 0,
     "verdict": ("iPSC-derived hepatocyte bioreactor yields and a "
                 "toxicogenomics classifier; neither has per-donor PK")},
    {"repository": "GEO", "query": "pharmacokinetics AND transcriptome",
     "hits": 3867, "on_target": None,
     "verdict": ("the term matches loosely and the surface is not usable: "
                 "GEO indexes the ASSAY, not the phenotype, so a per-donor "
                 "clearance almost never appears in the metadata")},
    {"repository": "PRIDE", "query": "pharmacokinetics", "hits": 25,
     "on_target": 1,
     "verdict": "PXD057282, tucatinib metabolic clearance"},
    {"repository": "PRIDE",
     "query": "human liver microsomes / intrinsic clearance / hepatocyte CYP",
     "hits": 13, "on_target": 4,
     "verdict": ("PXD066482, PXD057282, PXD060653 and PXD014131. PRIDE is "
                 "the productive seam because proteomics deposits accompany "
                 "DMPK papers, which carry the PK in their supplements")},
    {"repository": "MetabolomicsWorkbench",
     "query": "study title contains 'pharmacokinetic'", "hits": 1,
     "on_target": 0,
     "verdict": "ST003436, fasnall in NSG mice. Not human"},
]

#: The autonomous search of 2026-11-12, aimed at the four named gaps.
#: The negatives are the substance here: two channels are structurally unable
#: to reach what the programme needs, and that is worth knowing once.
SEARCHES_2026_11_12 = [
    {"repository": "PRIDE", "query": "transporter", "hits": 60,
     "on_target": 2,
     "verdict": ("both are in-vitro models -- an MDR1 efflux study in cancer "
                 "cells and a MASH cell model -- and neither has a PK or a "
                 "transport-function endpoint")},
    {"repository": "PRIDE", "query": "jejunum", "hits": 37, "on_target": 5,
     "verdict": ("PXD020987 is the find: 16 patients, ileum AND jejunum, 29 "
                 "enzymes and 32 transporters. PXD011498 adds a within-"
                 "subject before/after design")},
    {"repository": "PRIDE", "query": "intestinal", "hits": 42,
     "on_target": 5, "verdict": "overlaps the jejunum set"},
    {"repository": "PRIDE", "query": "albumin", "hits": 60, "on_target": 4,
     "verdict": ("albumin appears as an abundant contaminant or a disease "
                 "marker, never as a binding measurement paired with a free "
                 "fraction")},
    {"repository": "PRIDE", "query": "unbound", "hits": 60, "on_target": 1,
     "verdict": "no deposit measures a fraction unbound"},
    {"repository": "PRIDE", "query": "tacrolimus", "hits": 7, "on_target": 0,
     "verdict": ("no tacrolimus omics in PRIDE at all, which closes the most "
                 "obvious route to a longitudinal TDM cohort")},
    {"repository": "PRIDE", "query": "longitudinal", "hits": 60,
     "on_target": 10,
     "verdict": ("plasma proteomic time courses exist -- critical illness, "
                 "radiotherapy, sarcopenia -- and none carries a drug "
                 "concentration, so the H4 design exists without the H4 "
                 "endpoint")},
    {"repository": "PRIDE", "query": "ABCB1", "hits": 11, "on_target": 1,
     "verdict": "no human ABCB1 quantification paired with an endpoint"},
    {"repository": "EuropePMC",
     "query": ("\"fraction unbound\" AND pharmacokinetic* AND OPEN_ACCESS:Y "
               "AND HAS_SUPPL:Y"),
     "hits": 2628, "on_target": 0,
     "verdict": ("STRUCTURAL NEGATIVE: the hits are dominated by PBPK and "
                 "IVIVE modelling papers, which CONSUME a fraction unbound as "
                 "a model input rather than measure it beside omics in the "
                 "same subjects. The keyword surface finds the consumers, not "
                 "the producers")},
    {"repository": "EuropePMC",
     "query": ("transporter AND (abundance OR proteomic*) AND "
               "pharmacokinetic* AND OPEN_ACCESS:Y AND HAS_SUPPL:Y"),
     "hits": 3150, "on_target": 2,
     "verdict": ("same pattern. The two on-axis hits are placental and "
                 "skeletal-muscle transporter atlases, neither paired with "
                 "a PK parameter")},
    {"repository": "EuropePMC",
     "query": "tacrolimus AND (transcriptom* OR proteom* OR metabolom*) AND trough",
     "hits": 508, "on_target": 0,
     "verdict": ("reviews, association studies and a BXD mouse systems-"
                 "genetics paper. No deposited human cohort with repeated "
                 "omics beside repeated troughs")},
    {"repository": "GEO", "query": "jejunum AND Homo sapiens", "hits": 2146,
     "on_target": None,
     "verdict": ("the tissue surface is large and unusable for this purpose, "
                 "which confirms the earlier finding: GEO indexes the ASSAY, "
                 "so a per-subject drug concentration never appears in its "
                 "metadata whatever the query")},
]


def main(argv=None):
    registry = json.loads(REGISTRY.read_text(encoding="utf-8"))
    by_key = {d["key"]: d for d in registry["datasets"]}

    rows = []
    for key, declared in DECLARED.items():
        entry = by_key.get(key)
        access = (entry or {}).get("access_tier", "unknown")
        n_units = (entry or {}).get("n_analysed")
        for hypothesis in CAT.HYPOTHESES:
            status, reasons = CAT.feasibility(
                hypothesis, endpoint=declared["endpoint"],
                pairing=declared["pairing"], access=access,
                layers_available=declared["layers"],
                comparator_available=declared["comparator_available"][
                    hypothesis],
                n_units=n_units, longitudinal=declared["longitudinal"],
                n_candidate_omics=declared.get("n_candidate_omics"))
            rows.append({
                "source": "registry", "dataset": key,
                "label": (entry or {}).get("label", key),
                "hypothesis": hypothesis,
                "question": CAT.HYPOTHESES[hypothesis]["question"],
                "drug": declared["drug"], "endpoint": declared["endpoint"],
                "endpoint_kind": CAT.endpoint_kind(declared["endpoint"]),
                "comparator": CAT.HYPOTHESES[hypothesis]["minimum_comparator"],
                "pairing": declared["pairing"],
                "pairing_provenance": declared["provenance"],
                "access": access, "n_units": n_units,
                "status": status, "reasons": reasons,
                "note": declared["note"],
                "axes": list(declared.get("axes")
                             or ("metabolic_clearance",)),
            })

    for key, candidate in {**CANDIDATES, **CANDIDATES_2026_11_12}.items():
        for hypothesis in CAT.HYPOTHESES:
            status, reasons = CAT.feasibility(
                hypothesis, endpoint=candidate["endpoint"],
                pairing=candidate["pairing"], access=candidate["access"],
                layers_available=candidate["layers"],
                comparator_available=(
                    candidate["comparator_available"][hypothesis]
                    if isinstance(candidate["comparator_available"], dict)
                    else candidate["comparator_available"]),
                n_units=candidate["n_units"],
                longitudinal=candidate["longitudinal"],
                n_candidate_omics=candidate.get("n_candidate_omics"))
            rows.append({
                "source": "new_candidate", "dataset": key,
                "label": candidate["label"], "hypothesis": hypothesis,
                "question": CAT.HYPOTHESES[hypothesis]["question"],
                "drug": candidate["drug"], "endpoint": candidate["endpoint"],
                "endpoint_kind": CAT.endpoint_kind(candidate["endpoint"]),
                "comparator": CAT.HYPOTHESES[hypothesis]["minimum_comparator"],
                "pairing": candidate["pairing"],
                "pairing_provenance": candidate["provenance"],
                "access": candidate["access"], "n_units": candidate["n_units"],
                "status": status, "reasons": reasons,
                "note": candidate["what_is_missing"],
                "why_it_matters": candidate["why_it_matters"],
                "repository": candidate["repository"],
                "accession": candidate["accession"],
                "licence": candidate["licence"],
                "axes": list(candidate.get("axes") or ()),
            })

    summary = {}
    for hypothesis in CAT.HYPOTHESES:
        here = [r for r in rows if r["hypothesis"] == hypothesis]
        summary[hypothesis] = {
            status: [r["dataset"] for r in here if r["status"] == status]
            for status in CAT.STATUSES}

    payload = {
        "run": "hypothesis_dataset_matrix",
        "generated": "2026-11-12",
        "no_experiment_was_run": True,
        "catalogue": {
            "A_clinical_inputs": {k: {"role": v[0], "mechanism": v[1]}
                                  for k, v in CAT.CLINICAL_INPUTS.items()},
            "B_biological_inputs": {
                k: {"layer": v[0], "may_inform": list(v[1])}
                for k, v in CAT.BIOLOGICAL_INPUTS.items()},
            "C_measurement_context": list(CAT.MEASUREMENT_CONTEXT),
            "D_mechanistic_knowledge": {
                k: {"members": v[0], "may_inform": list(v[1])}
                for k, v in CAT.MECHANISTIC_KNOWLEDGE.items()},
            "E_sources": {k: {"kind": v[0], "what_it_is": v[1]}
                          for k, v in CAT.SOURCES.items()},
            "pk_parameters": CAT.PK_PARAMETERS,
            "functional_endpoints": CAT.FUNCTIONAL_ENDPOINTS,
        },
        "hypotheses": CAT.HYPOTHESES,
        "legacy_objectives": CAT.LEGACY_OBJECTIVES,
        "legacy_crosswalk": {k: {"maps_to_current": v[0], "why": v[1]}
                             for k, v in CAT.LEGACY_CROSSWALK.items()},
        "searches_run_today": SEARCHES,
        "searches_2026_11_12": SEARCHES_2026_11_12,
        "axes": CAT.AXES,
        "axis_coverage": CAT.axis_coverage(rows),
        "new_candidates": {**CANDIDATES, **CANDIDATES_2026_11_12},
        "matrix": rows,
        "summary_by_hypothesis": summary,
        "what_this_is_not": CAT.what_this_module_does_not_do(),
    }
    out = ROOT / "reports" / "data" / "hypothesis_dataset_matrix.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False,
                              default=str) + "\n", encoding="utf-8")

    print(f"artefact: {out}\n")
    print("axis coverage -- datasets reaching each mechanism, and how many "
          "carry an endpoint at all:")
    for axis, block in payload["axis_coverage"].items():
        print("  %-22s informs %-14s datasets %2d  with a PK endpoint %2d  "
              "READY anywhere %2d" % (axis, ",".join(block["informs"]),
                                      block["n_datasets"],
                                      block["n_with_a_pk_endpoint"],
                                      block["n_ready_for_any_hypothesis"]))
    print()
    for hypothesis, block in summary.items():
        counts = {s: len(v) for s, v in block.items() if v}
        print(f"{hypothesis}  {CAT.HYPOTHESES[hypothesis]['priority']:<18} "
              f"{counts}")
        for dataset in block["READY"]:
            print(f"     READY           {dataset}")
        for dataset in block["ACCESS_REQUIRED"]:
            print(f"     ACCESS_REQUIRED {dataset}")
        for dataset in block["PARTIAL"]:
            print(f"     PARTIAL         {dataset}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
