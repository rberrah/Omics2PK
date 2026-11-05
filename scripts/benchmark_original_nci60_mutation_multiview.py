#!/usr/bin/env python3
"""Original cross-study 47-line somatic mutation + RNA/protein/metabolite GI50.
Individual source files only in temporary CI workspace; aggregate output only.
"""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from benchmark_public_nci60_twoomics import load_dataset
from benchmark_original_nci60_threeomics_swath import (
    load_official_proteome,distinct_map,PIN
)
from omics2pk.mutation_multiview import evaluate_mutation_modality_budget


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--data",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    ids,cls,y,_,_,rna,met=load_dataset(args.data)
    with (args.data/"mutation.csv").open(newline="",encoding="utf8") as f:
        table=csv.reader(f)
        header=next(table)
        records=list(table)
    if (not header or header[0]!="subject_id" or header[1]!="mut_TP53" or
            header[-1]!="mutation_burden" or
            len(set(header))!=len(header) or len(records)!=47 or
            tuple(row[0] for row in records)!=tuple(ids) or
            any(len(row)!=len(header) for row in records)):
        raise ValueError("ABSTAIN: original CellMiner mutation identities and calls do not align exact NCI60 original 47 lines")
    mut=np.asarray([[float(v) for v in row[1:]] for row in records],dtype=float)
    if not np.isfinite(mut).all():
        raise ValueError("ABSTAIN: missing or non-finite CellMiner mutation status")
    names, classes, protein, sha=load_official_proteome()
    ref=distinct_map(ids,"original GI50/multiomics line")
    pindex=distinct_map(names,"original SWATH line")
    idx=[pindex[k] for k in ref if k in pindex]
    if len(idx)!=47 or len(set(idx))!=47 or any(classes[idx[i]]!=cls[i] for i in range(47)):
        raise ValueError("Original NCI60 SWATH identity histology discordance")
    results=evaluate_mutation_modality_budget(
        ids,cls,y,header[1:],mut,
        {"rna":np.asarray(rna),"prot":np.asarray(protein)[idx],"met":np.asarray(met)},
        repeats=3)
    results["source_omic_features"]={"rna":rna.shape[1],"prot":protein.shape[1],"met":met.shape[1]}
    results["original_SWATH_SHA256"]=sha
    results["rcellminerData_binary_mutation_gene_source"] = "Bioconductor rcellminerData release 2.34.0, PMID 23856246"
    results["original_47_line_join_requires_exact_canonical_names"]=True
    results["individual_original_source_rows_committed_or_output"]=False
    args.out.write_text(json.dumps(results,indent=2,sort_keys=True)+"\n",encoding="utf8")
    print(json.dumps(results,indent=2,sort_keys=True))


if __name__=="__main__":main()
