#!/usr/bin/env python3
"""Publishable aggregate-only figures: omics information budget and SHAP heatmap.

PK is a continuous target, so ROC/AUC is NOT appropriate. Show the actual
heldout MAE-vs-feature-count Pareto and descriptive fold-training SHAP heatmap.
No row-level omics, patient labels, concentrations or test predictions stored.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def plots(report,outdir):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    endpoints=list(report["endpoints"].items())
    if len(endpoints)!=2:
        raise ValueError("Expected 2 genuine metformin PK endpoints")
    fig,axes=plt.subplots(1,2,figsize=(13,5),layout="constrained")
    labels={"metformin_AUClast_0_24_original_units":"Metformin AUC 0–24",
            "metformin_CLr_renal_clearance_original_units":"Metformin renal CLr"}
    for ax,(label,value) in zip(axes,endpoints):
        series=value["fixed_k_full_nested_outer_heldout_curve"]
        ks=sorted(int(k) for k in series)
        y=[series[str(k)]["MAE_original_PK_table_units"] for k in ks]
        ax.plot(ks,y,marker="o",linewidth=1.4,markersize=4,label="Fixed k, heldout n=32")
        base=series["0"]["MAE_original_PK_table_units"]
        ax.axhline(base,linestyle=":",linewidth=1.3,label="Clinical only")
        pareto={int(row["k"]) for row in value["posthoc_descriptive_heldout_pareto_do_not_use_for_selection"]["non_dominated_models"]}
        for k in pareto:
            ax.scatter([k],[series[str(k)]["MAE_original_PK_table_units"]],
                       marker="s",s=55,zorder=3)
        for fraction,plotmarker in (("80","^"),("90","D")):
            model=value["adaptive_gain_retention_budgets"][fraction]
            ax.scatter([model["mean_selected_omic_features"]],
                       [model["observed_outer_patient_heldout_performance"]["MAE_original_PK_table_units"]],
                       marker=plotmarker,s=80,label=f"Nested {fraction}% inner gain")
        ax.set_xscale("symlog",linthresh=2)
        ax.set_xticks([0,1,2,5,10,20,40,80,160,384])
        ax.set_xticklabels(["0","1","2","5","10","20","40","80","160","384"],rotation=45)
        ax.set_xlabel("Number of omic features retained")
        ax.set_ylabel("Heldout MAE (original source numeric units)")
        ax.set_title(labels.get(label,label))
        ax.grid(alpha=.25)
        ax.legend(fontsize=7,loc="best")
    fig.suptitle("Original human Rotroff PK: molecular information vs complexity\n"
                 "Nested outer LOOCV; descriptive, no external validation",fontsize=11)
    fig.savefig(outdir/"pk_nested_information_pareto.png",dpi=175)
    plt.close(fig)

    # Only publicly named source feature metadata + aggregated absolute SHAP
    # from source original training folds, not a heatmap of patient profiles.
    top={}
    for endpoint,d in endpoints:
        top[endpoint]=d["top25_training_fitted_SHAP_ranking_frequency"]
    rows=sorted(set().union(*(set(v) for v in top.values())),
                key=lambda k:-sum(x.get(k,{}).get("average_train_SHAP_abs_log1p_PK",0)
                                  for x in top.values()))[:20]
    if not rows:
        return ["pk_nested_information_pareto.png"]
    mat=np.asarray([[top[key].get(name,{}).get("average_train_SHAP_abs_log1p_PK",0)
                     for key,_ in endpoints] for name in rows],dtype=float)
    denom=np.maximum(np.max(mat,axis=0),1e-12)
    scaled=mat/denom
    fig,ax=plt.subplots(figsize=(8,max(5,len(rows)*.35+1.7)),layout="constrained")
    im=ax.imshow(scaled,aspect="auto",vmin=0,vmax=1,cmap="viridis")
    ax.set_xticks(range(len(endpoints)))
    ax.set_xticklabels([labels.get(k,k) for k,_ in endpoints],rotation=15,ha="right")
    ax.set_yticks(range(len(rows)));ax.set_yticklabels(rows,fontsize=8)
    ax.set_title("Fold-training SHAP importance by original feature\n"
                 "(normalised within each endpoint, not predictive gain)")
    fig.colorbar(im,ax=ax,label="Relative mean absolute training Tree SHAP (0–1)")
    fig.savefig(outdir/"pk_shap_endpoint_feature_heatmap.png",dpi=175)
    plt.close(fig)
    return ["pk_nested_information_pareto.png","pk_shap_endpoint_feature_heatmap.png"]

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--aggregate",required=True)
    p.add_argument("--outdir",required=True)
    arg=p.parse_args()
    r=json.loads(Path(arg.aggregate).read_text(encoding="utf8"))
    if not r.get("status","").startswith("ORIGINAL_"):
        raise SystemExit("Original source benchmark not qualified")
    print(json.dumps({"figures":plots(r,arg.outdir),
                      "no_patient_omic_level_values_or_PK_results_rendered":True},indent=2))
