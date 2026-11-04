#!/usr/bin/env python3
"""Experimental STRING v12.5 network overlay of original NCI60 XGB SHAP genes.

This is NOT an enrichment/causal mechanism test: top candidates were chosen
using the GI50 observations in a prior exploratory CV analysis. Uses no
clinical data, patient IDs, expression values or new pharmacokinetics.

STRING functional association != direct physical binding/catalysis or drug
transport; no new neighbor nodes are added and species is fixed human 9606.
"""
from __future__ import annotations
from collections import Counter
import csv,io,json
from urllib.parse import urlencode
from urllib.request import urlopen,Request
from urllib.error import HTTPError,URLError

URL="https://version-12-5.string-db.org/api/tsv/network"
TOP_SHAP_SOURCE_RNA=("TMEM139","TNFRSF21","GOLT1A","ATP8B1","CLEC2B")
TOP_SHAP_SOURCE_PROTEIN=("IL18","PLXNB2","LGALS3","BAIAP2L1","JUP")
LITERATURE_DOX_DRUG_MECHANISM=("TOP2A","TOP2B","ABCB1","ABCC1","ABCG2","CBR1","CBR3")
ALL=tuple(dict.fromkeys(TOP_SHAP_SOURCE_RNA+TOP_SHAP_SOURCE_PROTEIN+LITERATURE_DOX_DRUG_MECHANISM))
SUBSET=set(TOP_SHAP_SOURCE_RNA+TOP_SHAP_SOURCE_PROTEIN)
PRIOR=set(LITERATURE_DOX_DRUG_MECHANISM)

def parse_string_output(source_tsv):
    if not isinstance(source_tsv,str) or len(source_tsv)>250_000:
        raise ValueError("Unqualified STRING response size")
    if not source_tsv.strip():
        return {"status":"NO_INTERACTIONS_AT_FIXED_SCORE_CUTOFF",
                "n_edges":0,"edges":[],"n_high_shap_genes_connected_to_drug_priors":0}
    rows=list(csv.DictReader(io.StringIO(source_tsv),delimiter="\t"))
    if not rows:
        return {"status":"NO_INTERACTIONS_AT_FIXED_SCORE_CUTOFF",
                "n_edges":0,"edges":[],"n_high_shap_genes_connected_to_drug_priors":0}
    mandatory={"preferredName_A","preferredName_B","score"}
    if not mandatory.issubset(rows[0]):
        raise ValueError("STRING source field names changed")
    result=[]
    seen=set()
    bridges=set()
    types=Counter()
    for row in rows:
        a=row["preferredName_A"].strip();b=row["preferredName_B"].strip()
        if not a or not b or a==b or a not in ALL or b not in ALL:
            # Protein synonym/other species auto-mapping must not create
            # a made-up original measured gene or phenotype target.
            continue
        try:confidence=float(row["score"])
        except (TypeError,ValueError):
            raise ValueError("STRING confidence score malformed") from None
        if not 0.4<=confidence<=1:
            raise ValueError("STRING confidence score inconsistent with query")
        pair=tuple(sorted((a,b)))
        if pair in seen:continue
        seen.add(pair)
        relation=("SHAP_to_SHAP" if a in SUBSET and b in SUBSET else
                  "DRUG_PRIOR_to_DRUG_PRIOR" if a in PRIOR and b in PRIOR
                  else "SHAP_to_DRUG_PRIOR")
        types[relation]+=1
        if relation=="SHAP_to_DRUG_PRIOR":
            bridges.add(a if a in SUBSET else b)
        result.append({
            "gene_A":pair[0],"gene_B":pair[1],
            "confidence_functional_association":round(confidence,4),
            "relationship_class":relation,
            "is_direct_physical_interaction_confirmed":False,
            "is_drug_PK_mechanism_validated":False,
            "is_independently_GI50_predictive":False,
        })
    result.sort(key=lambda edge:-edge["confidence_functional_association"])
    return {
        "status":"SOURCE_STRING_FUNCTIONAL_ASSOCIATIONS_NOT_CAUSAL" if result
                  else "NO_INPUT_GENE_INTERACTIONS_QUALIFIED",
        "n_edges":len(result),"edges":result[:30],
        "n_high_shap_genes_connected_to_drug_priors":len(bridges),
        "original_high_shap_to_drug_prior_bridge_gene_symbols":sorted(bridges),
        "edge_class_counts":dict(types),
        "total_edges_not_limit_truncated":len(result)<=30,
    }


def main():
    # Exactly 17 source-measured/independent-lab gene symbols from the executed
    # original n47 GI50 XGB SHAP benchmark plus prespecified doxorubicin genes.
    url=URL+"?"+urlencode({
      "identifiers":"\r".join(ALL),
      "species":9606,
      "network_type":"functional",
      "required_score":400,
      "add_nodes":0,
      "caller_identity":"Omics2PK-source-research",
    })
    try:
        with urlopen(Request(url,headers={"User-Agent":"Omics2PK literature-linked human functional network research"}),timeout=25) as response:
            raw=response.read(250_001)
        network=parse_string_output(raw.decode("utf-8"))
    except (HTTPError,URLError,TimeoutError,ValueError,UnicodeError) as exc:
        return {"status":"ABSTAIN_STRING_SOURCE_API_UNAVAILABLE",
                "error_class":type(exc).__name__,"n_edges":None,
                "biological_network_function_confirmed":False,
                "molecular_fusion_gain_confirmed":False}
    return {
       **network,
       "database":"STRING 12.5 functional network",
       "reference_api":URL,"species":9606,
       "n_actual_SHAP_source_gene_symbols":len(SUBSET),
       "n_prespecified_doxorubicin_mechanism_genes":len(PRIOR),
       "n_source_genes_queried":len(ALL),
       "network_score_threshold":0.4,"add_external_neighbor_nodes":False,
       "output_is_functional_association_not_a_direct_physical_interaction":True,
       "no_unknown_metabolite_bin_assigned_to_a_gene":True,
       "patient_PK_or_drug_transport_rate_measured":False,
    }

if __name__=="__main__":
    print(json.dumps(main(),indent=2,sort_keys=True))
