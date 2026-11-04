"""Published Topletz vinorelbine n20 HLM genotype/protein/outcome research gate.

Original article: DOI 10.1124/dmd.113.051094, Table 2.
This accepts only original source-qualified rows supplied by the caller.
There is NO bundled individual donor table, no internet fetch, no clinical PK,
and no claim that a manually transcribed table was machine-replayed.

The observed endpoint is the substrate-specific 1-uM vinorelbine metabolism
rate. Source-reported CLint and CLh are mathematically derived from rate and
physiological assumptions, not independent experimentally observed endpoints.
"""
from __future__ import annotations

from math import isfinite, log
from statistics import mean
import re
from typing import Mapping, Sequence

from omics2pk.functional_activity_benchmark import _ridge_predict
from omics2pk.functional_knowledge import assess_cyp3a5_evidence


class TopletzSourceError(ValueError):
    pass


DOI = "10.1124/dmd.113.051094"
PMC = "PMC3876810"


def _number(item: object, field: str, *, allow_zero: bool = False) -> float:
    if isinstance(item, bool):
        raise TopletzSourceError(field + " cannot be boolean")
    try:
        n = float(item)
    except (ValueError, TypeError, OverflowError):
        raise TopletzSourceError(field + " must be numeric") from None
    if not isfinite(n) or n < 0 or (not allow_zero and n == 0):
        raise TopletzSourceError(field + " must be finite and positive/zero as declared")
    return n


def qualify_topletz_table2(
    rows: Sequence[Mapping[str, object]],
) -> list[dict]:
    """Validate 20 independent donor records, keeping donor IDs private."""
    if len(rows) != 20:
        raise TopletzSourceError("Expected exactly 20 individual HLM donor rows")
    out=[]
    donors=set()
    for row in rows:
        if not isinstance(row, Mapping):
            raise TopletzSourceError("Original donor row must be a mapping")
        if row.get("source_doi") != DOI:
            raise TopletzSourceError("Original Topletz DOI mismatch")
        donor=row.get("donor_id")
        if not isinstance(donor,str) or not re.fullmatch(r"IUL-\d{1,3}",donor):
            raise TopletzSourceError("Unexpected original HLM donor ID format")
        if donor in donors:
            raise TopletzSourceError("Duplicate original HLM donor")
        donors.add(donor)
        dip=row.get("cyp3a5_diplotype")
        try:
            call=assess_cyp3a5_evidence(diplotype=dip,call_status="complete")
        except ValueError as err:
            raise TopletzSourceError("Invalid CYP3A5 source diplotype") from err
        if call["diplotype_classification"] is None:
            raise TopletzSourceError("Unresolved CYP3A5 functional allele")
        group=row.get("original_measured_CYP3A5_expression_group")
        if group not in {"high","low"}:
            raise TopletzSourceError("Original measured CYP3A5 group required")
        out.append({
            "donor_id": donor,
            "genotype": float(dip.split("/").count("*1")),
            "cyp3a4_protein": _number(row.get("cyp3a4_protein_pmol_mg"),
                                      "measured CYP3A4 protein"),
            # Original Table 2 includes a real measured CYP3A5 **zero**.
            "cyp3a5_protein": _number(row.get("cyp3a5_protein_pmol_mg"),
                                      "measured CYP3A5 protein",allow_zero=True),
            "observed_rate": _number(row.get("vinorelbine_rate_pmol_min_mg"),
                                    "original measured 1uM vinorelbine rate"),
            "group": group,
            "diplotype": dip,
        })
    if sum(x["group"]=="high" for x in out)!=8 or sum(x["group"]=="low" for x in out)!=12:
        raise TopletzSourceError("Original 8 high + 12 low source protein groups mismatch")
    if sum(x["diplotype"]=="*1/*6" for x in out)!=1 or sum(
            x["diplotype"]=="*1/*7" for x in out)!=1:
        raise TopletzSourceError("Original *1/*6 and *1/*7 genotype controls absent")
    if not any(x["diplotype"]=="*1/*3" and x["group"]=="low" for x in out):
        raise TopletzSourceError("Original low measured CYP3A5 in *1/*3 donor missing")
    if not any(x["cyp3a5_protein"] == 0 for x in out):
        raise TopletzSourceError("Original zero measured CYP3A5 protein source control absent")
    return sorted(out,key=lambda x:x["donor_id"])


ARMS = (
    "unconditional_donor_mean",
    "cyp3a5_genotype_only",
    "cyp3a5_protein_only",
    "cyp3a4_protein_only",
    "cyp3a4_protein_plus_cyp3a5_genotype",
    "cyp3a4_plus_cyp3a5_proteins",
    "cyp3a4_plus_cyp3a5_proteins_and_genotype",
)


def benchmark_topletz_table2(rows: Sequence[Mapping[str, object]]) -> dict:
    """Fixed LODO log-activity Ridge(alpha=1), training-only feature scaling."""
    selected=qualify_topletz_table2(rows)
    n=len(selected)
    outcomes={k:[] for k in ARMS}
    def vector(r,arm):
        pgx=r["genotype"]
        p4=r["cyp3a4_protein"]
        p5=r["cyp3a5_protein"]
        return {
            "cyp3a5_genotype_only":[pgx],
            "cyp3a5_protein_only":[p5],
            "cyp3a4_protein_only":[p4],
            "cyp3a4_protein_plus_cyp3a5_genotype":[p4,pgx],
            "cyp3a4_plus_cyp3a5_proteins":[p4,p5],
            "cyp3a4_plus_cyp3a5_proteins_and_genotype":[p4,p5,pgx],
        }[arm]
    for i,held in enumerate(selected):
        train=[row for j,row in enumerate(selected) if j!=i]
        for arm in ARMS:
            if arm=="unconditional_donor_mean":
                prediction=mean(row["observed_rate"] for row in train)
            else:
                prediction=_ridge_predict(
                    [vector(row,arm) for row in train],
                    [log(row["observed_rate"]) for row in train],
                    [vector(held,arm)])[0]
            outcomes[arm].append(abs(prediction-held["observed_rate"]))
    mae={arm:round(mean(errs),6) for arm,errs in outcomes.items()}
    return {
        "status":"ORIGINAL_SOURCE_DECLARED_NOT_AUTOMATICALLY_RETRIEVED_E2",
        "doi":DOI,"pmcid":PMC,
        "n_genuine_biological_donors":n,
        "same_LODO_donors_for_all_arms":True,
        "source_protein_group_counts":{"high":8,"low":12},
        "CYP3A5_functional_allele_expected_vs_protein_measured_separate":True,
        "has_measured_CYP3A5_zero_protein_donor":True,
        "endpoint":"observed_in_vitro_vinorelbine_depletion_at_1uM",
        "endpoint_unit":"pmol_vinorelbine_min_mg_microsomal_protein",
        "not_derived_Clint_or_CLh":True,
        "heldout_MAE_by_arm":mae,
        "genotype_only_worse_than_CYP3A4_protein_heldout":\n            mae["cyp3a5_genotype_only"] > mae["cyp3a4_protein_only"],
        "no_automatic_source_byte_verification":True,
        "no_individual_rows_donor_ids_or_predictions_exported":True,
        "not_clinical_PopPK_or_tacrolimus":True,
    }
