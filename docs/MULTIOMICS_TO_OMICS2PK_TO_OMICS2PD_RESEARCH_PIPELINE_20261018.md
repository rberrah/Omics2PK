# MultiOmics → Omics2PK → optional Omics2PD: research platform architecture

**State 18 Oct 2026:** MultiOmics and Omics2PK are **separate repos**, and a true native JavaScript producer→Omics2PK V1 ZIP contract is already integration-tested on fictitious 8-person data in `docs/NATIVE_MULTIOMICS_INTEGRATION.md`. Omics2PD is a **planned independent scientific stage**, NOT a production-deployed software repository or validated full PK/PD prediction engine. This architecture makes the eventual staged pipeline explicit without touching or destabilizing the MultiOmics/PMx websites.

## Data-flow architecture

```mermaid
flowchart TD
    A["Original individual omics data & documented context"] --> B["MultiOmics: QC, source-aware normalisation, batch/tissue metadata, pathway analysis"]
    B --> C["DOWNSTREAM_EXCHANGE_V1: SHA256 manifest, samples, features, assay-level matrices"]
    C --> D{"Same biological person + timestamps + protocol + source permissions"}
    D -->|No verified join| STOP["ABSTAIN: no PK or PD claim"]
    D -->|Matched source| E["Omics2PK: fold-frozen molecular scores, genotype vs phenotype comparisons"]
    F["Clinical sidecar: actual dose events, concentration times, PD observations and units"] --> D
    E --> G{"Actual outcome available?"}
    G -->|Measured C(t) / AUC| PK["Omics2PK: nested CV, PopPK covariates, external heldout PK"]
    G -->|Longitudinal clinical PD endpoint| PD["Omics2PD later: risk/response, true exposure-effect if PK-linked"]
    PK --> PD
    G -->|Only in-vitro CYP turnover| E2["E2 biochemical mechanism evidence ONLY"]
```

## Interface contracts and division of responsibility

| Component | Owns/exports | Must **NOT** silently infer |
|---|---|---|
| **MultiOmics** | assay/biological sample/subject metadata; actual RNA/protein/metabolite values after modality-appropriate QC; gene metabolite namespace, units, batch/organ/tissue, provenance, upstream raw normalization; qualified individual signatures only when fitted inductively | PK AUC, clinical dose schedules, clinical outcomes, causal pathway effect from differential expression p-values. GROUP log2FC ≠ person covariate |
| **Omics2PK** | same-person, predose time gates; clinical+PGx baseline; trained-on-train molecular covariate scores; source-specific CYP phenotyping; E2/real E3/PopPK E4; trained network version/hash and likelihood-based model comparison; PK predictions | tissue-independent extrapolation, parent drug feature measured postdose as PREdose risk factor, derived substrate turnover that uses future test drug AUC, transductive SWEET fit on heldout patients |
| **Future Omics2PD** | same clinical subject post-therapy observation, explicit treatment allocation, PD response/tolerability/survival, time-to-event and competing-risk handling; drug exposure-effect when measured PK available | claiming efficacy from prognosis-only classifier, time-leakage of post-outcome omics, treatment heterogeneity from unrandomized group cross-validation, clinical dosing decisions |

### Native producer (already exists)

The existing actual MultiOmics `exchange/` V1 files:

- `manifest.json` with `schema_version, study/data_access, software.commit, validation.status, files` + SHA256 and per-matrix `value_type`;
- `samples.csv` with `subject_id, sample_id, assay_id, omic` plus biologically attested specimen tissue, modality-specific dates/units and original assay replicates;
- `features.csv` with `omic,feature_id,id_namespace` and annotation/namespace provenance;
- `matrices/transcriptomics.csv`, `proteomics.csv`, `metabolomics.csv` with **measured individual assay** values.

Real native importer: `src/omics2pk/exchange.py` (`inspect_exchange`, `prepare_sweet`), and `src/omics2pk/multiomics_pk_readiness.py`. Existing current producer may omit relative `timepoint` categorical T0/T12 from downstream normalized V1 and say `software.commit=unknown`. **Do not change MultiOmics silently from this Omics2PK PR**. Source-specific sidecars must document a shared dose/sampling time origin, not invent hours or genotype.

### Two externally supplied clinical sidecars

**PK study records (user-controlled private/authorized source):**
```text
doses.csv
  subject_id,dose_id,drug,administered_at,amount,dose_unit,route
plasma_concentrations.csv
  subject_id,observation_id,dose_id,drug,sampled_at,concentration,
  concentration_unit,assay_method,BLQ_flag
```

PK identity/time validation is already implemented in `src/omics2pk/pk_link_contract.py`. Plasma, serum and whole-blood analyte compartments MUST be distinguished; with repeated doses and multiple-dose steady-state histories, a true PopPK event table needs complete dosing history and event codes, not "nearest dose" heuristic. BLQ requires prespecified rules, not silent substitution.

**PD study records (no PopPK required just to study prognosis):**
```text
original_treatment_events.csv
  subject_id,intervention,allocation,started_at
original_outcomes.csv
  outcome_event_id,subject_id,endpoint,assessed_at,observed_value
```

The added `src/omics2pk/staged_pkpd_contract.py` checks original ID matching and absolute, timezone-aware pre-treatment omic sample dates against subsequent outcome date. It rejects non-unique molecular assay selection, missing molecular modalities or post-treatment assay leakage. This is only a SOURCE GATE; whether outcome is well-defined, time-to-event censoring handled or treatment randomized is a further scientific requirement. RACE/CaPS original historical files use documented **T0 categorical baseline**, so a dataset-specific source-protocol adapter can affirm pre-treatment without inventing absolute ISO dates; never fabricate timestamps to pass the generic clinical event gate.

Privacy/provenance: any file with original human rows, source ID crosswalk, genotypes or molecular expression values stays in an **authorized research workspace**, not GitHub, PRs, logs, external public CI or website analytics. Public code and synthetic fixtures only; original-data scientific scripts may operate on openly licensed public sources ephemerally in-memory without redistributing original study rows.

## Prespecified genotype versus molecular covariate experiments

For **a clinically measured PK outcome on the same people**, use these arms with identical OUTER patient-grouped CV and, when possible, external donor/center validation:

1. `C0`: reference clinical covariates + true dosing/occasion metadata;
2. `C1`: C0 + full authentic PGx genotype/star-allele/haplotype/diplotype (e.g. CYP3A4/5, CYP2C19, CYP2D6, pharmacological mechanism as relevant), independently QC'ed call orientation and CNVs;
3. `M1`: C1 + prespecified **targeted** CYP transcript/pathway scores;
4. `M2`: C1 + target CYP proteomic abundance (if same-human original source exists);
5. `M3`: C1 + measured baseline metabolomic state and/or full omic tissue-specific signatures, if independently validated on train;
6. `M4`: C1 + RNA+protein and/or paired RNA+metabolome, with predeclared multiomic fusion, and matched complete-case participants;
7. `F`: prior/independent nonclinical probe CYP functional phenotyping, only if measured before clinical target outcome and same subject (a negative control for outcome leakage).

**Scientific primary contrast:** `(M1/M2/M3/M4) minus C1` on heldout CL/F, PK AUC0–24, original concentration–time RMSE, predicted drug exposure bias and calibration, all in the same patient subset. Always compare genotype alone vs genotype+phenotype on exactly the same patients; population omics n vs separate PGx n is NOT a valid comparison. Do not tune fold-specific features on heldout test outcomes. A prespecified target-gene qPCR count is **targeted RNA**, not whole transcriptome; an independently measured CYP probe activity is functional phenotype, NOT proteomics.

For PopPK covariate formulation, classical base `CL_i=CL_pop * exp(eta_i + beta_PGX genotype_i)`; then prospective biologically interpretable `exp(beta_RNA score_i)` etc, with true assay-time/tissue and interindividual distribution, shrinkage/reliability checks. (V,F,K_a) alternatives should be tested only if information-supported. E2 microsomal turnover should not be labeled true patient CL/F.

**For PD**, use measured outcomes (rejection/ABMR/TCMR, toxicities, blood cell response, graft loss, tumor size, longitudinal measurements, survival with death/censor). Add drug exposure as a dynamic covariate **only where dose+C(t) are truly available and time order is valid**. Baseline phenotype prognosis and patient-specific *causal* treatment benefit are different inferential targets: the latter needs trial assignment/treatment×biomarker designs and external robustness.

## Qualification stages and launch bar

- **T0** MultiOmics file import, original source namespace, assay QC, train-only features and temporal provenance.
- **T1** original clinical genotype / dose / concentration/PD patient ID matching; source SHA and legal access; hard abstention if no true original matching.
- **T2** genotype baseline and each molecular modality ablation against *same* measured endpoints, leak-controlled donor CV, repeated patient grouping, prespecified endpoints, negative controls (shuffle molecular rows within nuisance strata).
- **T3** external cohort/laboratory/population validation, reliability/cali, clinically meaningful absolute and relative errors.
- **T4** only thereafter PopPK mechanistic integration and joint drug concentration–effect/PD model; NEVER infer prospective clinical dosing safety from research E2/E3 internal CV alone.

## Worked real study source examples (evidence must remain accurately labeled)

- [Rotroff 2016](E3_REAL_METFORMIN33_PREDOSE_384_OMICS_AUC_CLR_HELDOUT_20261018.md): n32 measured metformin plasma AUC0–24 and CLr, source-paired predose GC-TOF; near-zero predictive gain (not liver metabolism) — **actual E3 source validation**.
- [HLC / Sage original](HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS_20261018.md): n≈400 CYP transcriptome→observed microsomal CYP activity E2; genuine RNA+PGx 171–172 subset, no measured CYP proteins or patient systemic PK.
- [RACE / CaPS observed drug response](CAPS_RACE_PUBLIC_ORIGINAL_PD_BENCHMARK_20261018.md): real n225 trial human 90d mortality clinical benchmark vs baseline NMR/LCMS (same metabolome); ~0.724 to 0.764 original held-person AUROC for NMR arm, but no clear statistically significant delta Brier/causal treatment effect. **Actual prospective PD research proof**, not true clinical drug PK.
- [Jarmusch 2020](DEEP_PUBLIC_HUMAN_OMICS_PK_PD_COHORT_CATALOGUE_20261018.md): strongest genotype+16S+metabolomics→midazolam PK cohort for further original source-author data request; unverified original numeric patient PK key.
