# Fold-frozen RNA score vs graph vs direct-expression PK benchmark (2026-10-19)

Research-only **technical comparator** added in `src/omics2pk/predose_representation_benchmark.py`. Nothing in this software check constitutes a new successful empirical prediction of PK. It is intended for the next **original same-human RNA + original measured clinical PK** cohort once legal, temporal, identity and molecular QC audits are complete.

## Study data and physical meaning

Call `compare_predose_rna_pk` with a matched *private* patient cohort, one biological person per column: original gene-by-person **normalized individual transcriptome**; clinical age/weight/other numeric baseline; optional genuine called 0/1/2 genomic allele dosages; measured individual positive AUC/CL or other declared PK outcome **with original physical unit**; tissue; exact source predose sampling and first dosing timestamps (with timezones); original source digest. Gene list and gene-score weights must be **prespecified before outcome analysis**, not selected from this cohort's outcomes. Identifiers, individual measurements, fitted subject-specific parameters, source file and predictions are NEVER exported from this function. For full original research, use institutional/private execution, not public workflow data artifacts.

The benchmark rejects repeat biological people, post-dose RNA, invalid times, duplicate genes and fake non-dosage genotype codes. It does **not** independently prove that upstream MultiOmics normalization was inductive, that a proposed gene panel is causally valid, that clinical data/PGx is comprehensive, that the source unit is physically true or that a claimed source SHA came from the legitimate supplier. Those source-level qualifications are separate mandatory gates.

## Frozen leave-one-person-out model arms

For every identical outer test person: the exact same training people, fixed Ridge(alpha=1) of log(PK outcome), train-only clinical and expression scaling, train-only RNA PCA(2), and train-only SWEET-inspired inductive reference with 5-fold crossfitted training graphs. The comparator reports source-unit MAE/RMSE/bias and delta-MAE to **clinical+PGx** if the original genotype vector is available, or clinical baseline otherwise. No hyperparameter selection using the test donor.

- Clinical only
- Clinical + measured PGx variants, *only if truly observed*
- Reference + unsupervised two RNA PCs (global transcript state)
- Reference + direct selected RNA gene expression
- Reference + gene strength from the **experimental SWEET-inspired** graph
- Reference + direct RNA and experimental graph
- Reference + externally prespecified weighted gene activity score, when supplied

**Important:** SWEET-inspired graph is NOT the upstream original SWEET algorithm, and its reference cohorts, marginal feature distributions, gene weights and tissue biological legitimacy still require independent validation. Expression, score and graph do not estimate enzyme flux, in-vivo CL/F, bioavailability or a PopPK random effect in isolation. It is not M3/M4 NLME. A true genotype-vs-omics superiority claim additionally needs full PGx alleles, same subjects, validation set independent of model discovery, uncertainty and prospective performance.

## Execution (synthetic tests only in public CI)

```bash
python -m pip install -e .
python -m pip install 'numpy==2.3.5' 'scikit-learn==1.7.2'
python -m unittest tests.test_predose_representation_benchmark -v
```

No original Rotroff/Jarmusch/busulfan patients are exported or joined by these synthetic tests. The already executed public **Rotroff n32 predose METABOLOMICS** benchmark has no authentic transcriptome, so DO NOT apply RNA-only SWEET or claim that the graph improves that PK benchmark. Its real observed clinical+metabolomics result remains negative/inconclusive. A future separate metabolomics-specific score/transport/excretion benchmark must not be mislabeled SWEET.

**Next clinical gate:** recover original Jarmusch pre-treatment same-human feature intensities (GNPS quantitative file), or obtain a donor-matched full RNA clinical PK study; predeclare comparisons and use the authentic specimen-to-subject crosswalk. For hepatic donor RNA/CLint datasets, report ex-vivo metabolic activity only, never measured systemic PK.

## Direct native MultiOmics V1 input (private, no intermediate patient exports)

`src/omics2pk/matched_multiomics_pk.py` now reads the **actual** `exchange/` directory from MultiOmics, verifies its hashes with `inspect_exchange`, and joins one unique transcriptome assay to one actual biological person through both `subject_id` and `assay_id`. It rejects repeated visits/assays, mixed tissues/ID namespaces, raw counts and bad checksum; **no silent sample-order matching**. The original source clinical sidecar must be held solely in an authorized private process, with one per-person record containing these exact fields:

```python
{
    "subject_id": "private-original-id",
    "assay_id": "private-exact-MultiOmics-assay-id",
    "collected_at": "2026-01-01T09:00:00+01:00",
    "first_dose_at": "2026-01-02T09:00:00+01:00",
    "measured_pk": 1.0,  # placeholder only; real ORIGINAL measured PK required
    "clinical_values": [1.0, 1.0],  # prespecified numeric baseline example
    "genotype_dosages": [0, 1],  # optional *actual measured* allele dosages
}
```

**The above is a fabricated schema example, NOT original measurements.** Internal call:

```python
from omics2pk.matched_multiomics_pk import run_matched_rna_exchange_pk
report = run_matched_rna_exchange_pk(
    "authorized_exchange_directory", original_clinical_sidecar,
    target_name="study_specific_PK_target",
    target_unit="original_verified_physical_unit",
    source_kind="original_source_declared",
    source_sha256="verified_64_character_lowercase_original_source_digest",
    prespecified_gene_weights=None,
)
```

For `original_source_declared`, the upstream exchange validation status must already be `qualified_for_declared_scope`; software cannot confer that status or independently verify a sidecar source. **No real private clinical source is part of public GitHub CI**. Fictional integration tests verify real hash checking, sample/person crosswalk, shuffled sidecar row order, pre-dose windows and negative controls. The adapter returns aggregate metrics only.

This V1 strictly handles a single eligible RNA specimen per person. Multiple timepoints, technical repeats, multicentre batch adjustments and complete PopPK dose-to-C(t) likelihood must be addressed explicitly in a future validated contract, not silently collapsed. DNA genotypes here are *external measured PGx sidecars*, not falsely inferred from RNA. A joint proteome/metabolome score or M3 Pillar-A/PopPK is a distinct planned extension.
