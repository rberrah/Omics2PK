"""First original public human tacrolimus C0/D source baseline (NO OMICS FIT).

PLOS ONE 2026 DOI 10.1371/journal.pone.0343817 supplemental table
provides 36 REAL individual whole-blood tacrolimus concentrations and dose,
not a longitudinal concentration-time PK history. The original 2026
clinical metabolomic/metagenomic sample crosswalk is NOT yet verified.
All donor rows are handled ephemerally and NEVER exported.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from statistics import mean, median
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from scripts.audit_public_tacrolimus_paired_omics_2026 import fetch_official, M

EXPECTED_SHA = "32c8a2d10c7c0448aebcb3ce12c32085098e24139b70c5b4c9fcd27ed3b37e49"
EXPECTED = (
    "No.","Age","Gender","Transplant Time Stage",
    "Daily Dose (mg)","Height (cm)","Weight (kg)","BMI (kg/m²)",
    "Daily Dose per kg (mg/kg)","Concentration (ng/ml)","C/D (ng·kg/ml·mg)",
    "Antibacterial Drugs","Immunosuppressive Drugs",
    "Type of Biliary Reconstruction","Portal Hypertension",
    "Conventional or Piggyback Technique",
)


class TacrolimusClinicalSourceError(ValueError):
    pass


def _index(ref):
    import re
    m=re.match(r"^[A-Z]+",ref)
    if not m:raise TacrolimusClinicalSourceError("Malformed source cell reference")
    n=0
    for letter in m.group(0):n=n*26+ord(letter)-64
    return n-1


def _xlsx_first_sheet(blob):
    with ZipFile(io.BytesIO(blob)) as z:
        stringtable=[]
        if "xl/sharedStrings.xml" in z.namelist():
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            stringtable=["".join(v.itertext()) for v in root.findall(M+"si")]
        sheet=ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        records=[]
        for r in sheet.findall(".//"+M+"sheetData/"+M+"row"):
            values={}
            for c in r.findall(M+"c"):
                col=_index(c.attrib["r"])
                v=c.find(M+"v")
                if v is None:
                    i=c.find(M+"is")
                    val="".join(i.itertext()) if i is not None else None
                elif c.attrib.get("t")=="s":
                    val=stringtable[int(v.text)]
                elif c.attrib.get("t") in ("str","inlineStr"):
                    val=v.text
                else:
                    try:val=float(v.text)
                    except (ValueError, TypeError):val=None
                values[col]=val
            if values:
                records.append([values.get(i) for i in range(max(values)+1)])
        return records


def _num(v, field):
    try:x=float(v)
    except (TypeError,ValueError):
        raise TacrolimusClinicalSourceError("Missing numeric original "+field) from None
    if not math.isfinite(x) or x<=0:
        raise TacrolimusClinicalSourceError("Nonpositive or invalid original "+field)
    return x


def qualify_clinical_rows(raw):
    if len(raw)!=37 or tuple(raw[0])!=EXPECTED:
        raise TacrolimusClinicalSourceError("Official 36-person Tac paper source layout changed")
    cases=[]
    seen=set()
    for r in raw[1:]:
        sid=str(r[0]).strip()
        if not sid or sid in seen:
            raise TacrolimusClinicalSourceError("Original rows must have 36 distinct patient identifiers")
        seen.add(sid)
        age=_num(r[1],"age")
        dose=_num(r[4],"daily dose")
        height=_num(r[5],"height")
        weight=_num(r[6],"weight")
        bmi=_num(r[7],"BMI")
        mg_kg=_num(r[8],"daily dose mg/kg")
        trough=_num(r[9],"tacrolimus C0")
        cd=_num(r[10],"dose-adjusted concentration")
        if not 0<age<=110:
            raise TacrolimusClinicalSourceError("Source age nonpositive or outside plausible human range")
        if not 5<=weight<=180:
            raise TacrolimusClinicalSourceError("Source weight outside expected adult range")
        if not 55<=height<=220:
            raise TacrolimusClinicalSourceError("Source height outside expected adult range")
        if not 0.1<=trough<=100:
            raise TacrolimusClinicalSourceError("Source tacrolimus C0 outside expected ng/mL range")
        if abs(dose/weight-mg_kg)>max(0.025*mg_kg,0.01):
            raise TacrolimusClinicalSourceError("Source dose per kg and measured weight inconsistent")
        computed=trough/mg_kg
        if abs(computed-cd)>max(0.05*cd,0.25):
            raise TacrolimusClinicalSourceError("Source C/D is not concentration divided by mg/kg daily dose")
        sex=str(r[2]).strip().lower()
        if sex not in {"male","female","m","f"}:
            raise TacrolimusClinicalSourceError("Source gender coding not qualified")
        cases.append({
            "source_id":sid,
            "age":age, "female":float(sex in ("f","female")),
            "BMI":bmi, "C0_ng_ml":trough, "daily_mg":dose,
            "dose_mg_kg":mg_kg, "CD_ng_kg_ml_mg":cd
        })
    if len(cases)!=36:
        raise TacrolimusClinicalSourceError("Need all 36 real source rows")
    return cases


def baseline_loocv(cases):
    """Intercept-only log C/D vs fixed age/sex/BMI ridge: original source only."""
    from omics2pk.epa_liver_cyp_validation import _ridge
    errors0=[]
    errors_clin=[]
    for i, case in enumerate(cases):
        tr=[r for j,r in enumerate(cases) if j!=i]
        ytrain=[math.log(x["CD_ng_kg_ml_mg"]) for x in tr]
        truth=case["CD_ng_kg_ml_mg"]
        p0=math.exp(mean(ytrain))
        xtrain=[[x["age"],x["female"],x["BMI"]] for x in tr]
        test=[case["age"],case["female"],case["BMI"]]
        p1=math.exp(max(-20,min(20,_ridge(xtrain,ytrain,test,alpha=1.0))))
        errors0.append(abs(p0-truth))
        errors_clin.append(abs(p1-truth))
    return {
        "n_patient_disjoint_holdout_folds":len(cases),
        "target":"observed tacrolimus dose-adjusted trough (C0 divided by daily mg/kg)",
        "baseline_reference":"training-only geometric mean C0/D",
        "clinical_reference":"fixed age+sex+BMI 3-covariate ridge log C0/D, alpha 1",
        "baseline_MAE_original_CD_units":round(mean(errors0),6),
        "clinical_MAE_original_CD_units":round(mean(errors_clin),6),
        "clinical_minus_intercept_MAE_delta":round(mean(errors_clin)-mean(errors0),6),
        "no_omics_model_fitted":True,
        "no_longitudinal_Concentration_time_history":True,
        "no_popPK_or_patient_dosing_claim":True,
    }


def run():
    original=fetch_official()
    if hashlib.sha256(original).hexdigest()!=EXPECTED_SHA:
        raise TacrolimusClinicalSourceError("PLOS DOI original clinical supplement SHA changed")
    cases=qualify_clinical_rows(_xlsx_first_sheet(original))
    return {
        "status":"REAL_HUMAN_TACROLIMUS_TROUGH_CLINICAL_BASELINE_NOT_OMICS_PREDICTION",
        "source_DOI":"10.1371/journal.pone.0343817",
        "source_sha256":EXPECTED_SHA,
        "n_original_patients":len(cases),
        "n_original_single_trough_concentrations":len(cases),
        "n_original_recipients_younger_than_18":sum(x["age"]<18 for x in cases),
        "units":"C0 ng/mL, daily dose mg and mg/kg, C/D ng*kg/(mL*mg)",
        "biological_unit":"one human liver transplant recipient, cross-sectional",
        "median_observed_tacrolimus_trough_ng_ml":round(median(x["C0_ng_ml"] for x in cases),5),
        "median_observed_dose_adjusted_trough":round(median(x["CD_ng_kg_ml_mg"] for x in cases),5),
        "original_clinical_C0_over_dose_formula_checked":True,
        "baseline":baseline_loocv(cases),
        "paired_stool_microbiome_RNA_or_metabolite_measurements_verified":False,
        "sequencing_original_bioproject":"PRJNA1416570",
        "original_individual_source_rows_or_ids_printed":False,
    }


if __name__=="__main__":
    print(json.dumps(run(),sort_keys=True,indent=2))
