#!/usr/bin/env python3
"""Functional STRING overlay for original NCI60 source ranked genes, no PK claim.

Visual categories: original RNA high SHAP, original protein high SHAP, and
predefined doxorubicin pharmacology gene priors. Only STRING v12.5 returned
functional edges, not invented edges. No patient data or raw omics in figure.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scripts.audit_nci60_string_shap_mechanistic_network import (
    TOP_SHAP_SOURCE_RNA,TOP_SHAP_SOURCE_PROTEIN,LITERATURE_DOX_DRUG_MECHANISM
)

def plot_source_network(source,path):
    result=source
    if result.get("status") not in (
            "SOURCE_STRING_FUNCTIONAL_ASSOCIATIONS_NOT_CAUSAL",
            "NO_INPUT_GENE_INTERACTIONS_QUALIFIED",
            "NO_INTERACTIONS_AT_FIXED_SCORE_CUTOFF"):
        return {"status":"ABSTAIN_EXTERNAL_NETWORK_NOT_AVAILABLE","file_created":False}
    groups=[
       ("RNA / SHAP",TOP_SHAP_SOURCE_RNA,0.0),
       ("Protein / SHAP",TOP_SHAP_SOURCE_PROTEIN,2.4),
       ("Literature drug genes",LITERATURE_DOX_DRUG_MECHANISM,4.8),
    ]
    positions={}
    fig,ax=plt.subplots(figsize=(13.3,7),layout="constrained")
    for name,genes,x in groups:
        color={"RNA / SHAP":"#D3E7F5","Protein / SHAP":"#E3DDF3",
               "Literature drug genes":"#F2E5CD"}[name]
        for index,gene in enumerate(genes):
            y=(len(genes)-1)/2-index
            positions[gene]=(x,y)
            ax.text(x,y,gene,ha="center",va="center",fontsize=9,weight="bold",
                    bbox=dict(boxstyle="round,pad=0.45",facecolor=color,edgecolor="#66717D",linewidth=.8))
        ax.text(x,max(3.6,len(genes)/2+.5),name,ha="center",va="bottom",
                fontsize=12,weight="bold")
    for relation in result.get("edges",[]):
        a,b=relation["gene_A"],relation["gene_B"]
        if a not in positions or b not in positions:
            continue
        xa,ya=positions[a];xb,yb=positions[b]
        score=relation["confidence_functional_association"]
        ax.plot([xa,xb],[ya,yb],linestyle="--",alpha=.15+.75*score,
                linewidth=.8+1.8*score,color="#466879",zorder=0)
    ax.set_xlim(-.8,5.6);ax.set_ylim(-4.1,5.2)
    ax.set_xticks([]);ax.set_yticks([])
    for spine in ax.spines.values():spine.set_visible(False)
    ax.set_title("NCI-60 model-ranked genes and prior doxorubicin proteins\n"
                 "STRING 12.5 functional associations only; not a causal or PK network",fontsize=12)
    ax.text(0,-3.7,"Dashed edges = original STRING functional associations ≥0.400; "
                 "unlinked genes intentionally retained.\n"
                 "RNA/protein source experiments not contemporaneous. GI50 is cell-line PD, not clinical PK.",
            ha="left",va="center",fontsize=8,color="#404E5D")
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=170,facecolor="white")
    plt.close(fig)
    return {"status":"ORIGINAL_STRING_FUNCTIONAL_NETWORK_FIGURE_ONLY",
            "file_created":True,
            "n_input_genes":sum(len(g) for _,g,_ in groups),
            "n_qualified_source_functional_edges_plotted":result.get("n_edges",0),
            "physical_binding_or_PK_function_not_claimed":True,
            "no_original_assay_patient_profiles_plotted":True}

if __name__=="__main__":
    arg=argparse.ArgumentParser()
    arg.add_argument("--input",required=True)
    arg.add_argument("--output",required=True)
    o=arg.parse_args()
    result=json.loads(Path(o.input).read_text(encoding="utf8"))
    print(json.dumps(plot_source_network(result,o.output),sort_keys=True))
