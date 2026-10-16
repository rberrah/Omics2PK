"""Source-agnostic, noncompartmental descriptive PK NCA for ONE observed dose.

This module computes only partial AUC(0,last), Cmax and Tmax from true
recorded concentration-time pairs. It does not estimate systemic clearance,
AUC(0,infinity), accumulation or a population PK model. Concentration at t=0
may be a nonzero pre-dose reported measurement and must not be overwritten.
"""
from __future__ import annotations

from math import fsum, isfinite
from statistics import mean, median
from typing import Sequence


class LongitudinalNCAError(ValueError):
    pass


def single_dose_descriptive_nca(points: Sequence[tuple[float,float]]) -> dict:
    if len(points)<3:
        raise LongitudinalNCAError("At least three genuine original concentration-time observations")
    data=sorted(points)
    previous=None
    for time,conc in data:
        if not isfinite(time) or not isfinite(conc) or time<0 or conc<0:
            raise LongitudinalNCAError("Invalid observation value/time or censored value passed to NCA")
        if previous is not None and time<=previous:
            raise LongitudinalNCAError("Duplicate or invalid PK concentration time")
        previous=time
    if data[0][0] != 0.0:
        raise LongitudinalNCAError("Original profile must include a genuine t=0 baseline sample")
    auc=fsum((x2-x1)*(y1+y2)/2.0 for (x1,y1),(x2,y2) in zip(data,data[1:]))
    cmax=max(v for _,v in data)
    tmax=next(t for t,v in data if v==cmax)
    return {
        "n_original_observations":len(data),
        "last_original_sample_h":data[-1][0],
        "Cmax_original_mg_l":cmax,
        "Tmax_h":tmax,
        "AUC0_last_original_mg_h_l_linear_trapezoid":auc,
        "does_not_impute_terminal_AUC":True,
    }


def aggregate_original_theoph(rows) -> dict:
    groups={}
    for row in rows:
        sub=row["Subject"]
        groups.setdefault(sub,[]).append((float(row["Time"]),float(row["conc"])))
    if len(groups)!=12:raise LongitudinalNCAError("Expected 12 original Theoph subjects")
    summaries=[single_dose_descriptive_nca(v) for v in groups.values()]
    return {
        "status":"REAL_ORIGINAL_SINGLE_DOSE_DESCRIPTIVE_NCA_NOT_POPPK",
        "n_original_human_subjects":len(summaries),
        "n_observations_total":sum(x["n_original_observations"] for x in summaries),
        "median_Cmax_mg_l":round(median(x["Cmax_original_mg_l"] for x in summaries),5),
        "median_Tmax_h":round(median(x["Tmax_h"] for x in summaries),5),
        "median_AUC0_last_mg_h_l":round(median(x["AUC0_last_original_mg_h_l_linear_trapezoid"] for x in summaries),5),
        "n_subjects_with_nonzero_original_predose_concentration":sum(next(y for t,y in groups[k] if t==0)>0 for k in groups),
        "n_subjects_followed_up_to_last_time_between_23_and_26h":sum(23<=x["last_original_sample_h"]<=26 for x in summaries),
        "no_patient_subject_level_results_exported":True,
        "limitation":"AUC0-last varying follow-up lengths, nonzero observed predose uncorrected, not AUCinf, CL/F, dose recommendation or population NLME",
    }
