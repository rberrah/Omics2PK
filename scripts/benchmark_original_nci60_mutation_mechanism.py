#!/usr/bin/env python3
"""Real NCI-60 somatic mutation-conditional mechanistic per-view nulls (GI50 only)."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from benchmark_public_nci60_twoomics import load_dataset
from benchmark_original_nci60_threeomics_swath import load_official_proteome,distinct_map
from omics2pk.mutation_mechanism_controls import compare_source_priors_by_mutation


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    ids,classes,y,rna_genes,met_names,rna,met=load_dataset(args.data)
    with (args.data/"mutation.csv").open(newline="",encoding="utf8") as f:
        records=list(csv.reader(f))
    labels=records[0][1:]
    if (records[0][0]!="subject_id" or labels[0]!="mut_TP53" or
            labels[-1]!="mutation_burden" or len(records)!=48 or
            tuple(row[0] for row in records[1:])!=tuple(ids)):
        raise ValueError("ABSTAIN original CellMiner mutation source no exact 47-line match")
    mutation=np.asarray([[float(x) for x in row[1:]] for row in records[1:]],dtype=float)
    names,tissues,proteins,sha,genes=load_official_proteome(return_gene_symbols=True)
    base=distinct_map(ids,"original GI50 cross-study line")
    pro=distinct_map(names,"original protein source line")
    ix=[pro[k] for k in base if k in pro]
    if len(ix)!=47 or any(tissues[ix[i]]!=classes[i] for i in range(47)):
        raise ValueError("ABSTAIN original source cross-lab protein lineage mismatch")
    out=compare_source_priors_by_mutation(
        ids,classes,y,mutation,labels,
        {"rna":rna_genes,"prot":genes,"met":met_names},
        {"rna":rna,"prot":proteins[ix],"met":met},
        repeats=3,controls=32)
    out["source_original_SWATH_SHA256"]=sha
    out["source_CellMiner_binary_mutation_release"]="Bioconductor rcellminerData 2.34.0"
    out["n_unrestricted_measured_features_per_layer"]={
        "rna":rna.shape[1],"prot":proteins.shape[1],"met":met.shape[1]}
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(out,indent=2,sort_keys=True))


if __name__=="__main__":main()
