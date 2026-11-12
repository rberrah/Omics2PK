#!/usr/bin/env python
"""Execute the same-donor readiness gate against every registered source.

`src/omics2pk/hepatic_multiview_readiness.py` has existed with unit tests since
an earlier phase, but it had never been run against a registered source, so the
project had no recorded verdict on which cohorts can carry a same-donor
RNA-plus-protein-plus-outcome comparison. This runner produces that verdict.

It does not fit a model and it computes no performance number. For each source
it supplies the gate with **real donor identifiers** wherever the individual
rows are on disk, and with `source_rows_acquired=False` wherever they are not —
which is the honest input, because the gate is explicitly documented to require
authenticated individual rows rather than "article statements, scatterplot
pixels, group summaries or row-order guesses".

Two cases deserve naming in advance:

* **The cross-study fusion the module exists to prevent** is run deliberately,
  as `INVALID_FUSION_hlc_rna_plus_wegler_protein`. It must be refused. If it is
  ever allowed, the gate is broken.
* **EATRIS-Plus is the only source whose RNA and protein are same-donor.** It
  fails on the outcome, not on the donors — and it is blood, not liver, which
  the gate does not test. That compartment limit is reported beside the verdict
  so a reader cannot mistake a donor-structure pass for biological relevance.

Run from the repository root. Every path is optional; a source whose rows are
absent is reported as not acquired rather than skipped:

    python scripts/run_hepatic_multiview_readiness_on_registered_sources.py \
        --hlc "<HLC directory>" \
        --wegler-design "<directory with the four experimentaldesigntemplate.txt>" \
        --eatris "<EATRIS_PLUS_mudata.h5mu>" \
        --out reports/data/hepatic_multiview_readiness_verdict.json
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from datetime import date
from pathlib import Path
from zipfile import ZipFile

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from omics2pk.hepatic_multiview_readiness import (  # noqa: E402
    assess_same_donor_hepatic_multiview,
)

REGISTRY = REPO_ROOT / "registry" / "PUBLIC_DATA_EVIDENCE_REGISTRY.json"

#: The gate accepts exactly these three outcome kinds and refuses anything else
#: with `ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME`. The first version of this runner
#: invented its own vocabulary ("measured_microsomal_cyp_activity" and friends)
#: and every source came back with that code and zero donor counts -- a report
#: that looked plausible and said nothing. The values are therefore named here
#: and `check_outcome_vocabulary` probes the gate at start-up rather than
#: trusting this comment.
CLINT = "drug_specific_clint"
ACTIVITY = "drug_specific_metabolic_activity"
CLINICAL_PK = "clinical_pk"
ACCEPTED_OUTCOME_KINDS = (CLINT, ACTIVITY, CLINICAL_PK)

#: EATRIS-Plus has no drug-specific outcome of any kind, so no accepted value
#: describes it. `ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME` is its true verdict, not an
#: artefact of a wrong label.
NO_OUTCOME = "none_measured"


def check_outcome_vocabulary():
    """Fail loudly if the gate no longer accepts the kinds this runner passes."""
    probe = dict(source_id="vocabulary_probe", rna_donors=["a"],
                 protein_donors=["a"], outcome_donors=["a"],
                 source_rows_acquired=True)
    for kind in ACCEPTED_OUTCOME_KINDS:
        verdict = assess_same_donor_hepatic_multiview(outcome_kind=kind, **probe)
        if verdict.status == "ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME":
            raise SystemExit(
                f"the gate no longer accepts outcome_kind {kind!r}; this runner "
                "would report a uniform refusal that means nothing")
    rejected = assess_same_donor_hepatic_multiview(
        outcome_kind=NO_OUTCOME, **probe)
    if rejected.status != "ABSTAIN_NO_DRUG_SPECIFIC_OUTCOME":
        raise SystemExit(
            "the gate now accepts an unnamed outcome kind; the EATRIS verdict "
            "would no longer mean what this runner reports")
    return True

HLC_ACTIVITY_ROWS = ("CYP1A2_activity", "CYP2A6_activity", "CYP2B6_activity",
                     "CYP2C8_activity", "CYP2C9_activity", "CYP2C19_activity",
                     "CYP2D6_activity", "CYP2E1_activity", "CYP3A4M_activity",
                     "CYP3A4T_activity")


def hlc_donors(hlc_dir: Path):
    """Expression donors and activity donors, read from the archives."""
    import pandas as pd

    expression_zip = hlc_dir / "curatedExpressionLiver.zip"
    phenotype_zip = hlc_dir / "curatedPhenotypesynapse.zip"
    if not (expression_zip.is_file() and phenotype_zip.is_file()):
        return None, None
    with ZipFile(expression_zip) as archive:
        with archive.open("expression.txt") as handle:
            rna = [str(c) for c in
                   pd.read_csv(handle, sep="\t", index_col=0, nrows=1).columns]
    with ZipFile(phenotype_zip) as archive:
        with archive.open("phenotype.txt") as handle:
            phenotype = pd.read_csv(handle, sep="\t", index_col=0, dtype=str)
    outcome = set()
    for row in HLC_ACTIVITY_ROWS:
        if row in phenotype.index:
            series = pd.to_numeric(phenotype.loc[row], errors="coerce")
            outcome |= {str(i) for i in series.dropna().index}
    return rna, sorted(outcome)


def wegler_donors(design_dir: Path):
    """Protein donors per preparation, read from the deposited design templates."""
    import csv

    out = {}
    for prep in ("Hom", "Mic", "Hep", "Pellet"):
        path = design_dir / f"{prep}_experimentaldesigntemplate.txt"
        if not path.is_file():
            continue
        with open(path, encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        samples = sorted({str(r["Sample_ID"]) for r in rows
                          if r.get("Sample_ID")})
        out[prep] = samples
    return out


def eatris_donors(h5mu: Path):
    """Per-modality individual identifiers, read from the deposited container."""
    try:
        import h5py
    except ImportError:
        return None
    if not h5mu.is_file():
        return None
    out = {}
    with h5py.File(h5mu, "r") as handle:
        for name in handle["mod"]:
            group = handle["mod"][name]
            obs = group.get("obs")
            if obs is None:
                continue
            key = obs.attrs.get("_index", "_index")
            key = key.decode() if isinstance(key, bytes) else str(key)
            if key not in obs:
                key = "_index"
            if key not in obs:
                continue
            raw = obs[key][:]
            out[name] = sorted(
                x.decode() if isinstance(x, bytes) else str(x) for x in raw)
    return out


def case(**kwargs):
    """One gate call, with its inputs recorded beside its verdict."""
    extra = {k: kwargs.pop(k) for k in list(kwargs)
             if k in ("note", "compartment", "registry_key")}
    verdict = assess_same_donor_hepatic_multiview(**kwargs)
    row = asdict(verdict)
    row["missing_requirements"] = list(verdict.missing_requirements)
    row["inputs"] = {
        "rna_donors_supplied": kwargs["rna_donors"] is not None,
        "protein_donors_supplied": kwargs["protein_donors"] is not None,
        "outcome_donors_supplied": kwargs["outcome_donors"] is not None,
        "source_rows_acquired": kwargs["source_rows_acquired"],
    }
    row.update(extra)
    return row


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hlc", default=None)
    parser.add_argument("--wegler-design", default=None)
    parser.add_argument("--eatris", default=None)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    check_outcome_vocabulary()

    hlc_rna, hlc_outcome = (hlc_donors(Path(args.hlc))
                            if args.hlc else (None, None))
    wegler = wegler_donors(Path(args.wegler_design)) if args.wegler_design else {}
    eatris = eatris_donors(Path(args.eatris)) if args.eatris else None

    cases = []

    # --- sources whose individual rows are on disk ---------------------
    cases.append(case(
        registry_key="hlc_syn4499",
        source_id="hlc_syn4499",
        rna_donors=hlc_rna,
        protein_donors=None,
        outcome_donors=hlc_outcome,
        outcome_kind=ACTIVITY,
        source_rows_acquired=hlc_rna is not None,
        compartment="liver tissue RNA; activity in microsomes from the same livers",
        note=("The project's only transcriptome. No protein layer exists in "
              "these donors, so a same-donor RNA-versus-protein comparison is "
              "impossible here however large the cohort."),
    ))

    wegler_protein = wegler.get("Mic")
    cases.append(case(
        registry_key="wegler_15",
        source_id="wegler_15_microsomes",
        rna_donors=None,
        protein_donors=wegler_protein,
        outcome_donors=None,
        outcome_kind=CLINT,
        source_rows_acquired=bool(wegler_protein),
        compartment="liver microsomes, homogenate and hepatocytes",
        note=("Individual protein rows are retrievable (PRIDE PXD014131, CC0). "
              "No RNA was measured in these donors, and the intrinsic-clearance "
              "values are not in the deposit, so the outcome rows are not "
              "acquired either."),
    ))

    if eatris:
        rna = eatris.get("mRNA")
        protein = eatris.get("proteomics-ipp")
        shared = sorted(set(rna or []) & set(protein or []))
        cases.append(case(
            registry_key="eatris_plus_v4",
            source_id="eatris_plus_v4",
            rna_donors=rna,
            protein_donors=protein,
            outcome_donors=None,
            outcome_kind=NO_OUTCOME,
            source_rows_acquired=True,
            compartment=("BLOOD. Leukocyte transcriptome and plasma proteome. "
                         "Not hepatic, so even a donor-structure pass would say "
                         "nothing about drug metabolism."),
            note=(f"The ONLY registered source whose RNA and protein are "
                  f"same-donor: {len(shared)} individuals carry both. It fails "
                  f"on the outcome, not on the donors. The gate does not test "
                  f"compartment, so this verdict must be read with the "
                  f"compartment field."),
        ))

    # --- the fusion the module exists to prevent ----------------------
    if hlc_rna and wegler_protein:
        cases.append(case(
            registry_key=None,
            source_id="INVALID_FUSION_hlc_rna_plus_wegler_protein",
            rna_donors=hlc_rna,
            protein_donors=wegler_protein,
            outcome_donors=hlc_outcome,
            outcome_kind=ACTIVITY,
            source_rows_acquired=True,
            compartment="two unrelated cohorts",
            note=("Run deliberately. HLC's 466 liver donors and Wegler's 15 "
                  "donors are different people in different identifier spaces. "
                  "This is the exact error the module's docstring names, and it "
                  "MUST be refused. A fit_allowed of true here would mean the "
                  "gate is broken."),
        ))

    # --- the remaining registered sources -----------------------------
    # `acquired` separates two different refusals that must not be conflated:
    # a source whose individual rows this project holds but which simply has no
    # RNA or protein layer (ABSTAIN_NO_MATCHED_RNA, a property of the cohort),
    # and a source whose rows are not in hand at all
    # (ABSTAIN_SOURCE_ROWS_NOT_ACQUIRED, a property of our access). Reporting
    # both as "not acquired" would hide which cohorts could never help even if
    # fully downloaded.
    for key, source_id, kind, acquired, note in (
        ("ariana_metoprolol_47", "ariana_metoprolol_47", CLINICAL_PK, True,
         "Individual rows ARE held (the public workbook). The cohort has "
         "targeted metabolomics, genotype and clinical PK, and no RNA and no "
         "protein, so the layer is what is missing, not the data."),
        ("jarmusch_midazolam", "jarmusch_midazolam", CLINICAL_PK, True,
         "Individual rows ARE held through the GNPS processed output. "
         "Untargeted metabolomics and genotype only; no RNA, no protein."),
        ("epa_liver_100", "epa_liver_100", ACTIVITY, False,
         "RNA and activity, no protein layer. Individual rows not on disk in "
         "this session, so the refusal is about access as well as layers."),
        ("latham_tki_donors", "latham_tki_donors", CLINT, False,
         "Protein and activity, no RNA layer. Individual rows not on disk in "
         "this session."),
        ("rotroff_metformin_33", "rotroff_metformin_33", CLINICAL_PK, False,
         "Metabolomics and clinical PK; no RNA, no protein. Individual rows "
         "not on disk in this session."),
        ("grangeon_2025_cyp_protein", "grangeon_2025_cyp_protein", ACTIVITY,
         False,
         "Protein AND activity were measured in the same 16 donors, which "
         "would be half of what the gate needs -- but NO per-donor value is "
         "published, so the rows cannot be acquired at all without digitising "
         "scatter figures. This is the one source whose refusal is caused by "
         "the publication rather than by the biology."),
    ):
        cases.append(case(
            registry_key=key, source_id=source_id,
            rna_donors=None, protein_donors=None, outcome_donors=None,
            outcome_kind=kind, source_rows_acquired=acquired, note=note,
        ))

    allowed = [c for c in cases if c["fit_allowed"]]
    fusion = [c for c in cases
              if c["source_id"].startswith("INVALID_FUSION")]
    if any(c["fit_allowed"] for c in fusion):
        raise SystemExit(
            "the gate allowed a cross-cohort fusion; refusing to write a report")

    report = {
        "report": ("same-donor hepatic multiview readiness, executed against "
                   "every registered source"),
        "label": "GATE_EXECUTION",
        "built_on": date.today().isoformat(),
        "module": "src/omics2pk/hepatic_multiview_readiness.py",
        "what_this_is": (
            "The first execution of a gate that has existed with unit tests but "
            "had never been run on a registered source. It fits no model and "
            "produces no performance number; it records which cohorts can carry "
            "a same-donor RNA-plus-protein-plus-outcome comparison."),
        "verdict": {
            "sources_assessed": len(cases),
            "fit_allowed": len(allowed),
            "fit_allowed_sources": [c["source_id"] for c in allowed],
            "cross_cohort_fusion_refused": bool(fusion) and not any(
                c["fit_allowed"] for c in fusion),
        },
        "headline": (
            "No registered source permits a same-donor RNA-versus-protein "
            "comparison against a drug-specific outcome. EATRIS-Plus is the "
            "only source with RNA and protein in the same individuals and it "
            "has no pharmacological outcome and is blood rather than liver; "
            "HLC has the transcriptome and no protein; Wegler has the protein "
            "and no transcriptome. The gate refuses the HLC-plus-Wegler fusion, "
            "which is the error it was written to prevent."
        ),
        "cases": cases,
        "how_to_read": (
            "fit_allowed is the only field that licenses an analysis. A status "
            "beginning ABSTAIN_ is a recorded scientific result, not a failure "
            "of execution. The gate tests donor structure and outcome presence; "
            "it does NOT test tissue compartment, so the compartment field must "
            "be read beside every verdict."),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    print(f"written {out}\n")
    print(f"{'source':46s} {'status':44s} fit  rna prot out")
    for row in cases:
        print(f"{row['source_id']:46s} {row['status']:44s} "
              f"{str(row['fit_allowed']):5s} {row['n_rna']:4d} "
              f"{row['n_protein']:4d} {row['n_outcome']:4d}")
    print(f"\nfit allowed on {len(allowed)} of {len(cases)} sources")
    print(f"cross-cohort fusion refused: "
          f"{report['verdict']['cross_cohort_fusion_refused']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
