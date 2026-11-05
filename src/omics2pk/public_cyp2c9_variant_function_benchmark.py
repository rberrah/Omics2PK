"""Original publication-anchored CYP2C9 variant-level protein-abundance/function test.

Source: Amorosi et al. 2021, AJHG DOI 10.1016/j.ajhg.2021.07.001,
dunhamlab/CYP2C9 data/CYP2C9_activity_abundance_scores.csv,
Git blob SHA 1a3041077f399dba6b4115877d31bd9def115787.

Outcome: click-seq CYP2C9 activity score measured in engineered yeast.
Abundance: VAMP-seq steady-state protein signal in human cells. The
measurements match variant identity but NOT a human-donor biospecimen.
No transcriptomic data, human plasma PK/PD, warfarin dose, matched RNA,
or original patient phenotype. Abundance does NOT equal specific activity.

Predeclared deterministic five-fold holdout by AMINO-ACID POSITION, not
random synonymous/missense variants sharing the same residue. All
standardization/ridge fits use the training positions only. The sequence
arm is intentionally SIMPLE and not a trained deep variant effect model.
No individual source variant rows, predictions, or assay data are exported.
"""
from __future__ import annotations
import csv
import io
from math import sqrt, isfinite
from random import Random
from statistics import mean
from collections import defaultdict

SOURCE_BLOB_SHA="1a3041077f399dba6b4115877d31bd9def115787"
SOURCE_DOI="10.1016/j.ajhg.2021.07.001"
SOURCE_REPOSITORY="dunhamlab/CYP2C9"
EXPECTED_ORIGINAL_ALL_VARIANTS=8701
EXPECTED_ORIGINAL_PAIRED_MISSENSE=4421
EXPECTED_ORIGINAL_PAIRED_POSITIONS=464
AA="ACDEFGHIKLMNPQRSTVWY"
HYDRO={"A":1.8,"C":2.5,"D":-3.5,"E":-3.5,"F":2.8,"G":-.4,
       "H":-3.2,"I":4.5,"K":-3.9,"L":3.8,"M":1.9,"N":-3.5,
       "P":-1.6,"Q":-3.5,"R":-4.5,"S":-.8,"T":-.7,"V":4.2,
       "W":-.9,"Y":-1.3}
VOLUME={"A":88,"C":108,"D":111,"E":138,"F":189,"G":60,
        "H":153,"I":167,"K":168,"L":167,"M":163,"N":114,
        "P":112,"Q":143,"R":173,"S":89,"T":116,"V":140,
        "W":227,"Y":193}
REQUIRED={"variant","class","position","start","end","activity_class",
          "activity_score","abundance_class","abundance_score"}
KNOWN_CLASSES={
    "wt-like","possibly_wt-like","possibly_decreased",
    "decreased","possibly_nonsense-like","nonsense-like","increased",
}
class PublicVariantDMSUnavailable(ValueError):
    pass


def _finite(s):
    if s is None or s=="NA" or s.strip()=="":
        return None
    try:
        x=float(s)
    except (ValueError, OverflowError):
        raise PublicVariantDMSUnavailable("Non-numeric original assay value") from None
    if not isfinite(x):
        raise PublicVariantDMSUnavailable("Nonfinite original assay")
    return x


def _charge(aa):
    return 1 if aa in "KR" else (-1 if aa in "DE" else 0)


def _features(seq):
    a,b=seq["start"],seq["end"]
    return [
        HYDRO[a],HYDRO[b],abs(HYDRO[a]-HYDRO[b]),
        _charge(b)-_charge(a),abs(VOLUME[a]-VOLUME[b])/100,
        float(b=="P")-float(a=="P"),
        float(b=="G")-float(a=="G"),
        float(b=="C")-float(a=="C"),
        float(b in "FWY")-float(a in "FWY"),
    ]


def parse_public_original(csv_bytes, *, strict_original=True):
    if not isinstance(csv_bytes,(bytes,bytearray)):
        raise PublicVariantDMSUnavailable("Original public CSV bytes required")
    if len(csv_bytes)>5_000_000:
        raise PublicVariantDMSUnavailable("Unexpected original source file size")
    try:
        text=bytes(csv_bytes).decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise PublicVariantDMSUnavailable("Original public CSV not UTF-8") from exc
    reader=csv.DictReader(io.StringIO(text))
    if reader.fieldnames is None or not REQUIRED.issubset(reader.fieldnames):
        raise PublicVariantDMSUnavailable("Unqualified original CYP2C9 activity/abundance headers")
    if len(reader.fieldnames)!=len(set(reader.fieldnames)):
        raise PublicVariantDMSUnavailable("Duplicate original source CSV header")
    out=[]
    count=0
    seen=set()
    for rec in reader:
        count+=1
        if None in rec:
            raise PublicVariantDMSUnavailable("Original CSV row has unexpected extra columns")
        variant=rec["variant"]
        if not isinstance(variant,str) or variant in seen or not variant:
            raise PublicVariantDMSUnavailable("Duplicate or empty original protein variant")
        seen.add(variant)
        if rec["class"] not in {"missense","synonymous","nonsense","wt"}:
            raise PublicVariantDMSUnavailable("Original genetic variant annotation unknown")
        if rec["class"]!="missense":
            continue
        if any(rec[k] not in AA for k in ("start","end")):
            raise PublicVariantDMSUnavailable("Missense variant amino acids not standard nuclear protein residues")
        try:
            pos=int(rec["position"])
        except (ValueError,TypeError) as exc:
            raise PublicVariantDMSUnavailable("Unresolved original protein position") from exc
        if pos<=0 or pos>1000 or variant!=rec["start"]+str(pos)+rec["end"] or rec["start"]==rec["end"]:
            raise PublicVariantDMSUnavailable("Original variant identity and amino acid coordinates conflict")
        yy=_finite(rec["activity_score"])
        ab=_finite(rec["abundance_score"])
        if (yy is None or ab is None or rec["activity_class"]=="NA"
                or rec["abundance_class"]=="NA"):
            continue
        if rec["activity_class"] not in KNOWN_CLASSES or rec["abundance_class"] not in KNOWN_CLASSES:
            raise PublicVariantDMSUnavailable("Original functional category not publication-qualified")
        out.append({
            "variant":variant,"position":pos,"start":rec["start"],"end":rec["end"],
            "activity":yy,"abundance":ab,
            "activity_class":rec["activity_class"],
            "abundance_class":rec["abundance_class"],
        })
    if strict_original and (count!=EXPECTED_ORIGINAL_ALL_VARIANTS
            or len(out)!=EXPECTED_ORIGINAL_PAIRED_MISSENSE
            or len({r["position"] for r in out})!=EXPECTED_ORIGINAL_PAIRED_POSITIONS):
        raise PublicVariantDMSUnavailable(
            "Original-public-source qualified row / paired-missense / position count changed")
    if len(out)<40 or len({r["position"] for r in out})<20:
        raise PublicVariantDMSUnavailable("Insufficient independently held-out variant positions")
    return out


def describe_known_class_controls(csv_bytes, *, strict_original=True):
    """Measured syn/missense/stop class controls on same original assay scales.

    This is *not* inferred human mRNA NMD: nonsense variants were assayed
    in engineered expression systems (yeast click-seq and human-cell VAMP-seq).
    No individual variant or original measurement row is released.
    """
    if not isinstance(csv_bytes,(bytes,bytearray)) or len(csv_bytes)>5_000_000:
        raise PublicVariantDMSUnavailable("Original variant class controls require source-qualified CSV bytes")
    try:
        stream=io.StringIO(bytes(csv_bytes).decode("utf-8-sig"))
    except UnicodeDecodeError as exc:
        raise PublicVariantDMSUnavailable("Source controls must be UTF-8") from exc
    reader=csv.DictReader(stream)
    if reader.fieldnames is None or not REQUIRED.issubset(reader.fieldnames):
        raise PublicVariantDMSUnavailable("Original genetic-code category source headers missing")
    classified=defaultdict(list)
    total=0
    for row in reader:
        total+=1
        kind=row["class"]
        if kind not in {"missense","synonymous","nonsense","wt"}:
            raise PublicVariantDMSUnavailable("Unexpected original variant genetic code class")
        activity=_finite(row["activity_score"])
        abundance=_finite(row["abundance_score"])
        if (activity is not None and abundance is not None and
                row["activity_class"]!="NA" and row["abundance_class"]!="NA"):
            classified[kind].append((activity,abundance))
    expected={"missense":4421,"synonymous":164,"nonsense":57,"wt":1}
    if strict_original and (total!=EXPECTED_ORIGINAL_ALL_VARIANTS or
           {key:len(classified[key]) for key in expected}!=expected):
        raise PublicVariantDMSUnavailable("Original source genetic-code control counts changed")
    def med(values):
        ordered=sorted(values)
        n=len(ordered)
        return ((ordered[(n-1)//2]+ordered[n//2])/2 if n else None)
    results={}
    for kind in ("missense","synonymous","nonsense","wt"):
        measurements=classified[kind]
        if not measurements:
            results[kind]={"status":"ABSTAIN_NO_QUALIFIED_ORIGINAL_CLASS_CONTROLS"}
        else:
            results[kind]={
                "n_source_variant_substitutions":len(measurements),
                "source_assay_median_activity_score":round(med([a for a,b in measurements]),6),
                "source_assay_median_abundance_score":round(med([b for a,b in measurements]),6),
                "clinical_or_native_tissue_variant_effect_proven":False,
            }
    return {
        "source_genetic_code_class_controls_not_human_mRNA_NMD":results,
        "stop_codon_construct_abundance_not_extrapolated_to_real_patient_NMD":True,
        "source_only_aggregate_no_original_variant_rows":True,
        "distinct_physical_assay_systems_yeast_activity_vs_human_cells_abundance":True,
    }


def _solve(a,b):
    n=len(b)
    a=[r[:] for r in a]
    b=b[:]
    for i in range(n):
        q=max(range(i,n),key=lambda k:abs(a[k][i]))
        a[i],a[q]=a[q],a[i];b[i],b[q]=b[q],b[i]
        pivot=a[i][i]
        if abs(pivot)<1.e-12:
            raise PublicVariantDMSUnavailable("Singular original ridge normal equations")
        for j in range(i,n):
            a[i][j]/=pivot
        b[i]/=pivot
        for k in range(n):
            if k==i:continue
            v=a[k][i]
            for j in range(i,n):
                a[k][j]-=v*a[i][j]
            b[k]-=v*b[i]
    return b


def _train_ridge(train, feature, *, alpha=30.0):
    ym=mean(r["activity"] for r in train)
    z=[feature(t) for t in train]
    p=len(z[0])
    if not p:
        return lambda _:ym
    means=[mean(row[j] for row in z) for j in range(p)]
    sds=[sqrt(mean((row[j]-means[j])**2 for row in z)) or 1.0
         for j in range(p)]
    standardized=[[(v-m)/sd for v,m,sd in zip(row,means,sds)]
                  for row in z]
    x=[r["activity"]-ym for r in train]
    gram=[[sum(row[i]*row[j] for row in standardized)
           +(alpha if i==j else 0.0)
           for j in range(p)] for i in range(p)]
    rhs=[sum(row[j]*v for row,v in zip(standardized,x))
         for j in range(p)]
    coefs=_solve(gram,rhs)
    def estimate(record):
        feat=feature(record)
        return ym+sum((v-m)/sd*b for v,m,sd,b
                      in zip(feat,means,sds,coefs))
    return estimate


def _ranks(x):
    indices=sorted(range(len(x)),key=lambda i:x[i])
    out=[0.0]*len(x)
    a=0
    while a<len(x):
        b=a+1
        while b<len(x) and x[indices[b]]==x[indices[a]]:
            b+=1
        for i in indices[a:b]:
            out[i]=(a+b+1)/2
        a=b
    return out


def _correlation(x,y):
    mx,my=mean(x),mean(y)
    covariance=sum((a-mx)*(b-my) for a,b in zip(x,y))
    varx=sum((a-mx)**2 for a in x)
    vary=sum((b-my)**2 for b in y)
    return covariance/sqrt(varx*vary) if varx>0 and vary>0 else None


def _auroc(predictions,case):
    if len(predictions)!=len(case) or len(case)<20:
        raise PublicVariantDMSUnavailable("Original heldout AUROC requires matched classes")
    positive=sum(case);negative=len(case)-positive
    if min(positive,negative)<10:
        raise PublicVariantDMSUnavailable("Original functional activity category counts too small")
    ranks=_ranks(predictions)
    return (sum(rank for rank,c in zip(ranks,case) if c)
            -positive*(positive+1)/2)/(positive*negative)


def benchmark_public_original(rows):
    if len(rows)<40 or len({r["position"] for r in rows})<20:
        raise PublicVariantDMSUnavailable("Insufficient unique protein residues")
    if len({r["variant"] for r in rows})!=len(rows):
        raise PublicVariantDMSUnavailable("Duplicate protein substitution in source")
    functions={
        "fold_training_mean":lambda _:[],
        "AA_physicochemical_change_only":_features,
        "protein_abundance_only":lambda r:[r["abundance"]],
        "AA_change_plus_protein_abundance":lambda r:_features(r)+[r["abundance"]],
    }
    pred={name:[None]*len(rows) for name in functions}
    heldout_positions_by_fold=[]
    for fold in range(5):
        train=[r for r in rows if ((r["position"]*2654435761)&0xFFFFFFFF)%5!=fold]
        test=[(i,r) for i,r in enumerate(rows) if
              ((r["position"]*2654435761)&0xFFFFFFFF)%5==fold]
        if len(train)<30 or len(test)<4:
            raise PublicVariantDMSUnavailable("Fold lacks independently held-out amino acid positions")
        if set(r["position"] for r in train)&{r["position"] for _,r in test}:
            raise PublicVariantDMSUnavailable("Position leakage into original test folds")
        heldout_positions_by_fold.append(len({r["position"] for _,r in test}))
        for name,feature in functions.items():
            model=_train_ridge(train,feature)
            for i,r in test:
                pred[name][i]=model(r)
    if any(not all(v is not None and isfinite(v) for v in ps)
           for ps in pred.values()):
        raise PublicVariantDMSUnavailable("Unpredicted original leave-position-out row")
    strict=[i for i,r in enumerate(rows) if
            r["activity_class"] in {"wt-like","decreased","nonsense-like"}]
    strict_y=[int(rows[i]["activity_class"]=="wt-like") for i in strict]
    outcomes=[r["activity"] for r in rows]
    metrics={}
    for name,ps in pred.items():
        errors=[a-b for a,b in zip(ps,outcomes)]
        metrics[name]={
            "held_position_out_MAE_activity_score":round(mean(abs(x) for x in errors),6),
            "held_position_out_RMSE_activity_score":round(sqrt(mean(x*x for x in errors)),6),
            "held_position_out_Spearman":round(
                _correlation(_ranks(ps),_ranks(outcomes)),6),
            "held_position_out_binary_activity_class_AUROC_not_PK_AUC":round(
                _auroc([ps[i] for i in strict],strict_y),6),
        }
    # Position-cluster paired uncertainty. There are 464 distinct residues,
    # NOT 4,421 independent gene locations. No refit or selection correction.
    groups=defaultdict(list)
    for i,r in enumerate(rows):
        groups[r["position"]].append(i)
    positions=sorted(groups)
    random=Random(20261105)
    contrasts={}
    for candidate,comparator in (
        ("protein_abundance_only","fold_training_mean"),
        ("AA_change_plus_protein_abundance","protein_abundance_only"),
        ("AA_physicochemical_change_only","fold_training_mean"),
    ):
        a=pred[candidate];b=pred[comparator]
        observed=mean(abs(a[i]-outcomes[i])-abs(b[i]-outcomes[i])
                      for i in range(len(rows)))
        simulated=[]
        for _ in range(999):
            sample=[random.choice(positions) for _ in positions]
            ia=[i for pos in sample for i in groups[pos]]
            simulated.append(mean(abs(a[i]-outcomes[i])-abs(b[i]-outcomes[i])
                                  for i in ia))
        simulated.sort()
        contrasts[candidate+"_minus_"+comparator]={
            "delta_MAE_negative_is_better":round(observed,6),
            "exploratory_95pct_resampled_by_independent_protein_position":[
                round(simulated[24],6),round(simulated[974],6)],
            "bootstrap_refitted_or_selection_adjusted":False,
        }
    strict_discordant=sum(r["activity_class"]=="decreased" and
                          r["abundance_class"]=="wt-like" for r in rows)
    true_damaging=sum(r["activity_class"] in {"decreased","nonsense-like"} and
                      r["abundance_class"] in {"wt-like","possibly_wt-like"}
                      for r in rows)
    return {
        "status":"ORIGINAL_PUBLIC_VARIANT_DMS_OUT_OF_FOLD_FUNCTION_NOT_HUMAN_PK",
        "source_doi":SOURCE_DOI,
        "upstream_repo":SOURCE_REPOSITORY,
        "original_git_blob_sha":SOURCE_BLOB_SHA,
        "original_source_variant_and_activity_rows_released":False,
        "n_individual_cyp2c9_missense_variants_with_both_scores":len(rows),
        "n_independent_protein_positions_grouped_for_holdout":len(positions),
        "heldout_unique_positions_per_fold":heldout_positions_by_fold,
        "n_original_public_missense_decreased_activity_but_wt_like_abundance":strict_discordant,
        "fraction_decreased_activity_wt_like_abundance":round(strict_discordant/len(rows),6),
        "n_lower_function_but_potentially_preserved_abundance_broader_nonconfirmatory":true_damaging,
        "activity_score_source_yeast_clickseq_probe":True,
        "protein_abundance_source_human_cell_VAMPseq_not_protein_activity":True,
        "scores_share_variant_identity_not_same_patient_or_tissue_sample":True,
        "not_transcriptome_metabolome_or_human_drug_PK":True,
        "no_specific_drug_metabolism_or_substrate_activity_measured":True,
        "consequence_missense_chemical_descriptors_no_direct_function_truth":True,
        "models_selected_without_outer_performance_comparison":True,
        "heldout_by_amino_acid_position_not_variant_random_split":True,
        "n_strict_source_function_classes_for_AUROC":len(strict),
        "heldout_metrics":metrics,
        "descriptive_paired_position_bootstraps":contrasts,
        "no_individual_allele_activity_or_dose_recommendation":True,
    }
