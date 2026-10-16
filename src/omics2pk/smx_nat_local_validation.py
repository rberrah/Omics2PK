"""Local-only source-gated hepatocyte NAT phenotype -> NA-SMX benchmark.

Primary source: Lanphier et al., Drug Metab Dispos 2025,
DOI 10.1016/j.dmd.2025.100195, PMC12861138, PRIDE PXD066786.

Outcome is actual in-vitro metabolite *formation*, NOT clinical PK,
plasma clearance, DILI incidence or donor-independent transport.
No publisher donor observations or peptide-level raw data bundled here.
"""
from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
from statistics import mean

SOURCE_DOI = "10.1016/j.dmd.2025.100195"
SOURCE_PXD = "PXD066786"
REQUIRED = (
    "donor_id", "row_kind", "source_doi", "na_smx_pmol_min_million_cells",
    "na_paba_pmol_min_million_cells", "na_smz_pmol_min_million_cells",
    "nat2_peptide_eflnshllpk_pmol_mg", "nat2_peptide_dntdlvefk_pmol_mg",
    "nat2_genotype_phenotype", "age_years", "sex", "bmi",
)
PHENOTYPES = {"slow", "intermediate", "rapid"}
EXPECTED_GENOTYPE_COUNTS = {"slow": 9, "intermediate": 6, "rapid": 5}


class NATSourceError(ValueError):
    """Insufficient or inconsistent original same-donor metabolite evidence."""


def _value(value, name, *, required=False):
    raw = "" if value is None else str(value).strip()
    if raw in ("", "-", "NA", "N/A", "ND", "NR"):
        if required:
            raise NATSourceError("Missing original measurement: " + name)
        return None
    try:
        x = float(raw)
    except ValueError as exc:
        raise NATSourceError("Non-numerical measurement: " + name) from exc
    if not math.isfinite(x) or x <= 0:
        raise NATSourceError("Original measurement must be finite and positive: " + name)
    return x


def qualify_original(rows):
    if len(rows) != 26:
        raise NATSourceError("Source specifies 26 original single-donor PHH preparations")
    people = []
    seen = set()
    for row in rows:
        if not set(REQUIRED).issubset(row):
            raise NATSourceError("Missing explicit original NA-SMX/NAT assay columns")
        if str(row["row_kind"]).strip().lower() != "individual":
            raise NATSourceError("Pooled PHH / repeat wells cannot be independent donors")
        if str(row["source_doi"]).strip() != SOURCE_DOI:
            raise NATSourceError("Unverified NAT study DOI")
        donor = str(row["donor_id"]).strip()
        if not donor or donor in seen:
            raise NATSourceError("Duplicate or absent biological donor ID")
        seen.add(donor)
        phenotype = str(row["nat2_genotype_phenotype"]).strip().lower()
        if phenotype in ("", "-", "na", "n/a", "unknown"):
            phenotype = None
        elif phenotype not in PHENOTYPES:
            raise NATSourceError("Unsupported original NAT2 acetylator classification")
        sex = str(row["sex"]).strip().lower()
        if sex not in ("male", "female"):
            raise NATSourceError("Sex field must match original study demographic coding")
        people.append({
            "id": donor,
            "y": _value(row["na_smx_pmol_min_million_cells"], "NA-SMX", required=True),
            "paba": _value(row["na_paba_pmol_min_million_cells"], "NAT1/NA-PABA"),
            "smz": _value(row["na_smz_pmol_min_million_cells"], "NAT2/NA-SMZ"),
            "peptide1": _value(row["nat2_peptide_eflnshllpk_pmol_mg"], "NAT2 peptide1"),
            "peptide2": _value(row["nat2_peptide_dntdlvefk_pmol_mg"], "NAT2 peptide2"),
            "phenotype": phenotype,
            "age": _value(row["age_years"], "age", required=True),
            "female": 1.0 if sex == "female" else 0.0,
            "bmi": _value(row["bmi"], "BMI"),
        })
    genotype = {name: sum(p["phenotype"] == name for p in people) for name in PHENOTYPES}
    if genotype != EXPECTED_GENOTYPE_COUNTS:
        raise NATSourceError("Original 20 genotyped donors must be 9 slow/6 intermediate/5 rapid")
    if sum(p["paba"] is not None and p["smz"] is not None for p in people) != 21:
        raise NATSourceError("Source reports 21 original donors with NAT probe assays")
    return people


def _model_data(row, arm):
    if arm == "clinical":
        vals = (row["age"], row["female"], row["bmi"])
    elif arm == "nat_probes":
        vals = (row["paba"], row["smz"])
    elif arm == "nat2_proteins":
        vals = (row["peptide1"], row["peptide2"])
    elif arm == "nat2_genotype":
        if row["phenotype"] is None:
            return None
        vals = (float(row["phenotype"] == "intermediate"),
                float(row["phenotype"] == "rapid"))
    else:
        raise NATSourceError("Model not prespecified")
    return vals if all(v is not None for v in vals) else None


def _linear_solve(matrix, vector):
    """Simple full pivot-safe Gaussian elimination, only tiny ridge matrices."""
    n = len(vector)
    a = [list(map(float, matrix[i])) + [float(vector[i])] for i in range(n)]
    for i in range(n):
        pivot = max(range(i, n), key=lambda j: abs(a[j][i]))
        if abs(a[pivot][i]) < 1e-12:
            raise NATSourceError("Degenerate ridge design")
        a[i], a[pivot] = a[pivot], a[i]
        d = a[i][i]
        for k in range(i, n + 1):
            a[i][k] /= d
        for j in range(n):
            if j == i:
                continue
            multiplier = a[j][i]
            for k in range(i, n + 1):
                a[j][k] -= multiplier * a[i][k]
    return [a[i][-1] for i in range(n)]


def _ridge(train_x, train_y, held_x, alpha=1.0):
    if not train_x:
        raise NATSourceError("No training subjects")
    d = len(held_x)
    mu = [mean(x[k] for x in train_x) for k in range(d)]
    sd = [math.sqrt(mean((x[k]-mu[k])**2 for x in train_x)) or 1.0 for k in range(d)]
    z = [[(x[k]-mu[k])/sd[k] for k in range(d)] for x in train_x]
    ym = mean(train_y)
    gram = [[sum(row[i]*row[j] for row in z)+(alpha if i==j else 0.0)
             for j in range(d)] for i in range(d)]
    cov = [sum(row[j]*(y-ym) for row,y in zip(z,train_y)) for j in range(d)]
    coefficients = _linear_solve(gram,cov)
    prediction = ym + sum(coefficients[k]*(held_x[k]-mu[k])/sd[k] for k in range(d))
    return max(0.0, prediction)


def evaluate(original_rows):
    people = qualify_original(original_rows)
    output = {
        "status": "SOURCE_CONTRACT_LOCAL_IN_VITRO_NAT_TO_NA_SMX_FORMATION_NOT_CLINICAL_PK",
        "original_biological_donors": 26,
        "n_genotyped": 20,
        "n_assayed_with_both_NAT_probes": 21,
        "outcome": "30min NA-SMX formation, pmol/min/million cells",
        "source_doi": SOURCE_DOI,
        "proteomics_source": SOURCE_PXD,
        "model": "donor-held-out fixed-alpha=1 Ridge with training-only scaling",
        "outcomes_cannot_be_treated_as_PK_CL_or_DILI": True,
        "arms": {},
    }
    for arm in ("clinical", "nat_probes", "nat2_proteins", "nat2_genotype"):
        valid = [(p,_model_data(p,arm)) for p in people]
        valid = [(p,x) for p,x in valid if x is not None]
        if len(valid) < 10:
            output["arms"][arm]={"status":"ABSTAIN_INSUFFICIENT_COMPLETE_CASES",
                                 "complete_case_donors":len(valid)}
            continue
        base_errors, model_errors = [], []
        for index, (person, x) in enumerate(valid):
            train = [(p,xx) for j,(p,xx) in enumerate(valid) if j != index]
            y = person["y"]
            baseline = mean(p["y"] for p,_ in train)
            prediction = _ridge([xx for _,xx in train],[p["y"] for p,_ in train], x)
            base_errors.append(abs(y-baseline))
            model_errors.append(abs(y-prediction))
        output["arms"][arm] = {
            "status":"DONOR_LODO_EXPLORATORY",
            "n_original_complete_case_donors":len(valid),
            "baseline_MAE_matched_cases":round(mean(base_errors),6),
            "model_MAE_matched_cases":round(mean(model_errors),6),
            "delta_MAE_model_minus_baseline":round(mean(model_errors)-mean(base_errors),6),
        }
    # Distinct assays and genotyping are not observed on the same 26 people.
    # Comparisons between arms require identical complete-case biological donors.
    arms=("clinical","nat_probes","nat2_proteins","nat2_genotype")
    common=[p for p in people if all(_model_data(p,arm) is not None for arm in arms)]
    if len(common)<10:
        output["strict_four_arm_shared_donor_comparison"]={
            "status":"ABSTAIN_INSUFFICIENT_IDENTICAL_COMPLETE_CASE_DONORS",
            "n_identical_complete_cases":len(common),
        }
    else:
        losses={name:[] for name in ("baseline",*arms)}
        for index,heldout in enumerate(common):
            training=[p for i,p in enumerate(common) if i!=index]
            ys=[p["y"] for p in training]
            y=heldout["y"]
            losses["baseline"].append(abs(y-mean(ys)))
            for arm in arms:
                prediction=_ridge([_model_data(p,arm) for p in training],ys,
                                  _model_data(heldout,arm))
                losses[arm].append(abs(prediction-y))
        output["strict_four_arm_shared_donor_comparison"]={
            "status":"IDENTICAL_DONOR_LODO_EXPLORATORY",
            "n_identical_complete_cases":len(common),
            "MAE":{name:round(mean(errors),6) for name,errors in losses.items()},
            "interpretation":"Identical donor subset; no independent population validation",
        }
    output["original_rows_redistributed"]=False
    return output


def run_local_csv(path):
    payload = Path(path).read_bytes()
    if len(payload) > 300_000:
        raise NATSourceError("Unexpected source size; original table is only 26 donors")
    rows = list(csv.DictReader(payload.decode("utf-8-sig").splitlines()))
    return {**evaluate(rows), "local_curation_sha256":hashlib.sha256(payload).hexdigest()}
