"""Real MultiOmics V1 individual matrices → fold-safe mono/multiomic PK arms.

All patient-specific PK/PGx/time sidecars are caller-supplied in memory and
stay outside GitHub and the CLI. No guesswork linking specimen IDs, dates
or anonymized GEO identifiers. This is NOT an automatic PopPK fit.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

from .exchange import inspect_exchange, _csv, _file, ExchangeError
from .multimodal_pk_ablation import (
    compare_modalities_pk, MultimodalPKBenchmarkError
)


class MultimodalExchangePKError(ValueError):
    pass


def _require(ok, explanation):
    if not ok:
        raise MultimodalExchangePKError(explanation)


def run_multimodal_exchange_pk(
    exchange_dir, clinical_pk_sidecar, *,
    selected_omics, pk_target_name, pk_target_unit,
    source_kind, source_sha256, fixed_scores=None,
    allow_complete_case=False, include_xgboost=False,
    n_outer_splits=5, validation_mode="participant_kfold",
):
    """Read same-person MultiOmics source V1, with strict source clinical join.

    clinical_pk_sidecar: list of authorized patient dictionaries, one record
    per original unique human, each with:
      subject_id: source original pseudonym key
      assay_ids_by_omic: {omic: original assay_id}
      sample_ids_by_omic: {omic: original source specimen ID}
      collected_at_by_omic: {omic: timezone-aware original ISO timestamp}
      index_dose_at: ISO timestamp of predictive index administration
      measured_pk: actually measured patient PK outcome (e.g. AUC, CL/F)
      clinical_values: fully measured prespecified numeric covariate vector
      genotype_dosages: all authentic 0/1/2 variants or absent for ALL people

    Source V1 manifest SHA attests the files, not the clinical sidecar.
    An original-data label is an *attestation by the data steward*, not
    independent scientific authentication.
    """
    root=Path(exchange_dir).resolve()
    report=inspect_exchange(root)
    mf=json.loads((root/"manifest.json").read_text(encoding="utf-8"))
    _require(mf["study"]["biological_unit"]=="participant",
             "Cross-individual human PK requires original participant identity")
    _require(isinstance(selected_omics,(list,tuple)) and
             1<=len(selected_omics)<=3 and
             len(selected_omics)==len(set(selected_omics)) and
             set(selected_omics)<={"transcriptomics","proteomics","metabolomics"} and
             set(selected_omics)<=set(report["matrices"]),
             "Select explicit existing omic layers without duplicate labels")
    if source_kind=="original_source_declared":
        _require(report["validation"]=="qualified_for_declared_scope",
                 "Original clinical association requires upstream source qualification")
    _require(isinstance(clinical_pk_sidecar,(list,tuple)) and len(clinical_pk_sidecar)>=12,
             "Original private clinical sidecar with >=12 people required")
    by_person={}
    for rec in clinical_pk_sidecar:
        _require(isinstance(rec,dict) and
                 isinstance(rec.get("subject_id"),str) and rec["subject_id"].strip()
                 and rec["subject_id"] not in by_person,
                 "No duplicate/missing original biological people in clinical sidecar")
        by_person[rec["subject_id"]]=rec
    head,samples=_csv(_file(root,mf["files"]["samples"]),
                      ["assay_id","sample_id","subject_id","omic"])
    original_samples=defaultdict(dict)
    source_population=set()
    for row in samples:
        person=row["subject_id"]
        source_population.add(person)
        omic=row["omic"]
        if omic not in selected_omics:
            continue
        _require(person not in original_samples[omic],
                 "Repeated specimens require a frozen visit/design before PK evaluation")
        _require(row.get("tissue","").strip(),
                 "Actual omics specimen tissue is mandatory")
        original_samples[omic][person]=row
    _require(set(by_person)<=source_population,
             "Clinical sidecar contains source-unknown biological patient")
    overlap=set(by_person)
    for omic in selected_omics:
        overlap &= set(original_samples[omic])
        compartments={r["tissue"] for r in original_samples[omic].values()}
        _require(len(compartments)==1,
                 "One molecular assay modality cannot pool multiple tissues")
    excluded=len(by_person)-len(overlap)
    _require(not excluded or allow_complete_case,
             "Multiomics complete-case attrition requires explicit researcher opt-in")
    _require(len(overlap)>=12, "Insufficient matched complete original individuals")
    ids=sorted(overlap)
    for p in ids:
        rec=by_person[p]
        assays=rec.get("assay_ids_by_omic")
        specimens=rec.get("sample_ids_by_omic")
        collected=rec.get("collected_at_by_omic")
        _require(isinstance(assays,dict) and isinstance(specimens,dict)
                 and isinstance(collected,dict),
                 "Original assay/specimen/time maps required; no inferred T0")
        for omic in selected_omics:
            row=original_samples[omic][p]
            _require(assays.get(omic)==row["assay_id"],
                     "Clinical and omic assay-to-patient identity mismatch")
            _require(specimens.get(omic)==row["sample_id"],
                     "Clinical and omic specimen-to-patient identity mismatch")
            _require(isinstance(collected.get(omic),str) and collected[omic],
                     "Original absolute collection time required for every modality")
            # V1 may preserve true absolute source collection timestamps.
            # A numeric/relative label is NOT inferred to be that timestamp.
            source_time=row.get("collection_time","").strip()
            if source_kind=="original_source_declared":
                _require(source_time and "T" in source_time and
                         "-" in source_time and ":" in source_time,
                         "Original source absolute specimen collection clock missing")
            if source_time and ("T" in source_time) and (
                    "-" in source_time and ":" in source_time):
                from datetime import datetime
                def dt(s):
                    try: return datetime.fromisoformat(s.replace("Z","+00:00"))
                    except ValueError: return None
                src=dt(source_time)
                private=dt(collected[omic])
                _require(src is not None and private is not None and
                         src.tzinfo is not None and private.tzinfo is not None and
                         src==private,
                         "Original MultiOmics and clinical molecular sampling clock disagree")
    _,feature_rows=_csv(_file(root,mf["files"]["features"]),
                        ["omic","feature_id","id_namespace"])
    matrices={}
    for omic in selected_omics:
        entries=[m for m in mf["matrices"] if m["omic"]==omic]
        _require(len(entries)==1,"Ambiguous V1 omics matrix")
        entry=entries[0]
        header,rows=_csv(_file(root,entry["file"]),["feature_id"])
        original_features=[f["feature_id"] for f in feature_rows if f["omic"]==omic]
        _require(len(rows)==len(original_features) and
                 {r["feature_id"] for r in rows}==set(original_features),
                 "Original measured feature inventory differs from manifest")
        assays=[original_samples[omic][p]["assay_id"] for p in ids]
        _require(len(assays)==len(set(assays)) and set(assays)<=set(header[1:]),
                 "Original source assay columns are inconsistent")
        # features by source table identity, patient mapping by authentic
        # assay key; row sorting is NOT used as an identity join.
        features=[r["feature_id"] for r in rows]
        values=[[float(r[assay]) for r in rows] for assay in assays]
        matrices[omic]={
            "feature_ids":features,"values":values,
            "collected_at":[by_person[p]["collected_at_by_omic"][omic] for p in ids],
            "tissue":original_samples[omic][ids[0]]["tissue"],
            "value_type":entry["value_type"],
        }
    _require(validation_mode in ("participant_kfold","leave_one_site_out"),
             "Explicit participant or whole-site validation required")
    original_sites=[by_person[p].get("site_id") for p in ids]
    if validation_mode=="leave_one_site_out":
        _require(all(isinstance(x,str) and x.strip() for x in original_sites),
                 "Every PK person needs a source-attested site for leave-one-site-out")
    else:
        _require(not any(x is not None for x in original_sites),
                 "Original site metadata detected: do not silently use mixed-site folds")
        original_sites=None
    pgx=[by_person[p].get("genotype_dosages") for p in ids]
    if all(g is None for g in pgx):
        pgx=None
    else:
        _require(all(isinstance(g,(list,tuple)) for g in pgx),
                 "Partial/or imputed PGx would invalidate reference comparator")
    try:
        benchmark=compare_modalities_pk(
            subject_ids=ids, modalities=matrices,
            clinical=[by_person[p]["clinical_values"] for p in ids],
            genotype_dosages=pgx,
            measured_pk=[by_person[p]["measured_pk"] for p in ids],
            index_dose_at=[by_person[p]["index_dose_at"] for p in ids],
            pk_target_name=pk_target_name,pk_target_unit=pk_target_unit,
            source_kind=source_kind,source_sha256=source_sha256,
            fixed_scores=fixed_scores,n_outer_splits=n_outer_splits,
            include_xgboost=include_xgboost,
            validation_mode=validation_mode,site_ids=original_sites)
    except (ValueError,KeyError,TypeError) as exc:
        if isinstance(exc,MultimodalPKBenchmarkError):
            raise
        raise MultimodalExchangePKError("Invalid original source clinical or PK measurements") from None
    benchmark["source_exchange_validated_sha256"]=True
    benchmark["source_exchange_qualified_original_declared_by_producer"]=(
        report["validation"]=="qualified_for_declared_scope")
    benchmark["n_private_clinical_people_offered"]=len(by_person)
    benchmark["n_complete_case_people_used"]=len(overlap)
    benchmark["n_people_excluded_for_missing_modalities"]=excluded
    benchmark["all_test_arms_share_identical_complete_case_people"]=True
    benchmark["clinical_sidecar_linkage_originality_caller_attested"]=True
    return benchmark
