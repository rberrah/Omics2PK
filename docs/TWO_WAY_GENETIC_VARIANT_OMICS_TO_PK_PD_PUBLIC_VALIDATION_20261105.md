# A unified, bidirectional evidence pipeline: MultiOmics ↔ Omics2PK

**5 November 2026 | public functional variants | original-source research audit.**

The unified central scientific question is not “Which omic is best?” but:

1. **MultiOmics direction:** Given an original experimentally qualified genetic
   variant, what does each measured molecular layer actually reveal about its
   *functional effect*? Can an omic predict a **separately assayed functional
   phenotype**, not merely the variant label, with honest held-out validation?
2. **Omics2PK direction:** Given an original drug, tissue, genotype (known,
   partial or unmeasured) and assays, **which layer, if any, improves the
   prediction of actual PK/PD** beyond clinical+genuine PGx/PopPK, and is that
   gain consistent with source-qualified biology?

The two directions have different targets and evidence thresholds. A measured
protein or omic is not in itself independent confirmation of enzyme function.
A purported “function score” cannot be built and scored using the *same* assay
readout; functional activity used as an outcome is not an independent predictor.
RNA and protein expression may be unchanged when activity is lost, and
a synonymous variant may disrupt splice usage.

## The exact causal-versus-statistical boundary

```
Qualified reference & genome assembly
  → authentic gene variant(s), transcript consequences, splice/codon/LOF knownness
  → per-assay RNA isoforms / protein coverage and localization /
    chemically qualified metabolites / independent function readout
  → MultiOmics: QC, features, functional hypotheses, uncertainty
  → optional source-qualified variant-effect evidence sidecar (not a new
    categorical "functional normal" truth inferred solely from omics)
  → Omics2PK: clinical + qualified PGx ± selected measured omics
  → independent held-patient PK/PD + biological consistency auditing
```

An unlabelled feature is not thrown away from a black-box experiment: it
may be predictive with no molecular identity. SHAP/STRING may suggest
experiments but cannot establish transporter specificity, a loss-of-function
mutation or causal effects. The choice to include an omic must be based on
**inner-fold predictive value**; the biological ranking is a prior/hypothesis,
not a fixed bonus to clinical activity.

## E1 actually achieved: original public CYP2C9 activity versus protein abundance

**Official source:** Amorosi et al., AJHG 2021,
DOI [10.1016/j.ajhg.2021.07.001](https://doi.org/10.1016/j.ajhg.2021.07.001);
the authors' original open [dunhamlab/CYP2C9](https://github.com/dunhamlab/CYP2C9),
`data/CYP2C9_activity_abundance_scores.csv`, Git blob SHA
`1a3041077f399dba6b4115877d31bd9def115787`.
The source was retrieved live from the immutable official GitHub Git blob API,
git-content SHA-verified, and benchmarked in
[GitHub Actions #41991782251](https://github.com/rberrah/Omics2PK/actions/runs/41991782251)
(**SUCCESS**, 2026-11-05). All source rows stay transient in the job; no
large third-party measurements, individual variant predictions or loss vectors
are committed to the repository.

**What is actually measured?** Activity comes from Click-seq with an
activity-based probe in engineered *yeast*. Protein abundance is an
independent VAMP-seq fluorescence proxy measured in human cultured cells.
These scores are matched by original **protein variant identity**, not
same human biospecimen; **no RNA, no endogenous metabolome, no direct clinical
warfarin/tacrolimus PK, no dose**. Assays and substrate dependence matter.

### Known genetic-code classes on the same original functional assay scales

| Original class | n with both measurements | Median activity score | Median protein abundance score |
|---|---:|---:|---:|
| Synonymous | 164 | 0.999905 | 0.991771 |
| Missense | 4,421 | 0.566604 | 0.782232 |
| Nonsense / stop | 57 | 0.010945 | 0.003991 |
| Reference control | 1 | 0.985036 | 0.993072 |

These are **experimental construct classes**, not a demonstration
of patient mRNA NMD or unique substrate-specific hepatic clearance.
Nonsense represents a strong *assay-class* control, not a declaration
that every human late-stop triggers NMD. The source has other variants
with only one measured assay: restrict paired comparisons to the
**same 4,421 missense variants**.

**Key counterexample:** **400/4,421 (9.05%)** original missense
variants are explicitly `decreased` for activity yet `wt-like`
for abundance according to the source's own measurement categories.
This directly falsifies any universal rule
“measured protein abundance ≈ active pharmacological protein.”
The reference-to-variant translation alone also cannot assign
a numeric percent of activity loss. Broadening to uncertain source
labels increases the discrepancy count but weakens its interpretation.

### True held-out experiment, at independent amino acid position

Unlike a random variant split, variants sharing an original protein
**residue** are kept together. Fixed five-fold position-disjoint
training/evaluation over **464 independent protein positions**,
no outcome-guided feature selection; fixed `Ridge(alpha=30)`,
standardization on training positions only.

| Inputs to predict independently assayed CYP2C9 activity score | Original held-position-out MAE | RMSE | Spearman | Source strict activity-class AUROC (not drug AUC) |
|---|---:|---:|---:|---:|
| Fold training mean | 0.328278 | 0.365678 | −0.038465 | 0.490867 |
| Simple amino-acid physicochemical difference | 0.311349 | 0.355013 | 0.245051 | 0.621551 |
| **Measured protein abundance** | **0.190947** | **0.242590** | 0.748420 | 0.846676 |
| Abundance + simple amino-acid properties | 0.190694 | 0.242442 | 0.749876 | 0.852587 |

AUROC uses **3,078** source variants with strict functional classes
`wt-like` versus `decreased`/`nonsense-like`, excluding
ambiguous `possibly_*` classifications.
The simple sequence model uses fixed hydropathy and side-chain
descriptors; it is not a full sequence-to-function foundation model,
nor a demonstration that all variant-sequence predictors are inferior.

**Paired original-protein-position bootstrap** (999 draws; interval
descriptive, no model refit / source-selection correction):
- abundance vs fold training mean: ΔMAE **−0.137331**,
  95% interval **[−0.149254; −0.125732]**;
- abundance plus AA features vs abundance only:
  ΔMAE **−0.000253**, **[−0.001793; +0.001357]**.

The abundance signal is informative, but **not universally functional**;
additional simple amino-acid features provide **no established incremental
gain** here. The strongly correlated abundance/activity relationship
was already reported in the original paper; the new contribution is
a source-replayed, residue-held-out **method qualification**, not a
novel biological discovery or clinical validation.

**Implementation:** `src/omics2pk/public_cyp2c9_variant_function_benchmark.py`,
`scripts/benchmark_public_cyp2c9_variant_activity_abundance.py`,
`tests/test_public_cyp2c9_variant_function_benchmark.py`,
`.github/workflows/cyp2c9-public-known-variants-abundance-vs-activity.yml`.

## E2 already partially measured: real human hepatocyte genotype → proteome/probe → drug depletion

The existing [Latham source replay](LATHAM_PHH_SOURCE_REPLAY_REAL_20261104.md)
contains **14 primary human hepatocyte donors, 13 complete protein records**,
CYP3A5 *3/*6/*7 genotype classes, CYP2C8/3A4/3A5 targeted peptide
abundances, CYP2C8 and **shared CYP3A4/5** probe activities, and
real *in vitro* imatinib apparent-clearance outcome.

Original **13 same-donor LODO** comparisons: PGx only MAE **2.906576**;
PGx+targeted proteins **3.003704**; PGx+two shared/probe activities
**2.762944**; PGx+both **3.262812** (source printed numerical scale).
In nine CPIC nonexpresser-class donors with complete protein values,
**9/9** have a positive CYP3A5-targeted peptide signal; all **10/10**
nonexpresser-class donors have positive **shared CYP3A** probe activity.
This is **not** evidence of residual CYP3A5-specific activity; it is an
assay-specific non-equivalence example.

This is n13, in vitro, non-selective probe and **without RNA**;
no genotype-stratum-specific omic selection, human tacrolimus plasma
PK or full biological mediation has been proven.

## Orthogonal public known-effect sources: different questions

| Real source | Known effect and independent observable | Valid inference | Limit |
|---|---|---|---|
| CYP3A5 *3 CRISPR HuH-7, Dorr et al. 2017 [PMCID PMC5518718](https://pmc.ncbi.nlm.nih.gov/articles/PMC5518718/) | Repair/alter cryptic splice allele restores CYP3A5 RNA/protein and increases **tacrolimus and midazolam metabolism** | Perturbation strengthens allele→functional metabolism evidence | Engineered clones, not patient-level PK; no generalization to other alleles |
| CYP3A5 *3/*6 fetal liver, [PMCID PMC6067386](https://pmc.ncbi.nlm.nih.gov/articles/PMC6067386/) | Total, normally spliced and abnormal SV1 RNA, protein, CYP3A-shared midazolam metabolites; 382 original RNA samples, subset functional assays | Test whether RNA isoform specificity improves **biological interpretation** | Figure/article-level summaries not a verified original 382× joint individual matrix; no individual heldout CV should be invented |
| CYP3A5 rs776746 eExon-QTL in primary hepatocytes, [PMCID PMC7192504](https://pmc.ncbi.nlm.nih.gov/articles/PMC7192504/) | n60 original donor genotype + RNA exon effects, source GEO GSE124076/GSE147628 | Exon-resolution can reveal CYP3A5 splice effect invisible in total-gene eQTL | Need exact source genotype/RNA join; no CYP3A5 selective function or PK confirmed |
| NUDT15, Suiter et al. 2020 [PMCID PMC7071893](https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/) | Public MaveDB score sets for ~2,844 variant-matched stability and **thiopurine cell sensitivity**; clinical toxicity relevance | Strong next **PD** replication: abundance ≠ drug susceptibility | Cytotoxicity is not isolated hydrolase activity; variant-labeled cultured-cell screens are not 2,844 patients |
| SLCO1B1 rs4149056 p.Val174Ala, [PMCID PMC10087178](https://pmc.ncbi.nlm.nih.gov/articles/PMC10087178/) | Known reduced OATP1B1 function and observed rosuvastatin PK association | Clinical PK check on transporter **missense** | Joint subject-level transcriptome/proteome/PK not yet source-verified |
| SLC22A1 OCT1 p.Arg61Cys / M420del, [JCI original paper](https://www.jci.org/articles/view/30558) | Published variants with similar measured RNA but reduced **substrate-dependent metformin transport** | Counterexample to gene-expression-only rule; substrate-specific function matters | Cross-variant experiments and clinical genotype response are not a matched multiomic PK cohort |

**Prioritization:** (1) CYP2C9 original variant-level replay DONE; (2)
NUDT15 original MaveDB reanalysis and variant-position-held-out
thiopurine PD benchmark NEXT; (3) source-row qualification of the
CYP3A5 *3 exonic/eQTL or splice-specific assay and direct CYP3A
function; (4) direct real same-patient tacrolimus/rosuvastatin PK
with complete pharmacogenetic panel. Do not manufacture a combined
genotype–RNA–protein–PK dataset by matching on a gene/variant only.

## Experimental design for a joint MultiOmics → Omics2PK scientific pathway

### Lane 1: MultiOmics (variant → molecular readout, and observed omics → functional phenotype)

- Preserve the genome build, HGVS/mapped allele, transcript, strand,
  zygosity/phase, authentic source; **no-variant-seen ≠ WT/normal function**.
- Include synonymous, missense, splice/NMD, nonsense and reference
  comparators, and **variants without observed function** as unlabelled.
- Compare total versus isoform-specific RNA, targeted peptide/whole
  protein, function-specific probes, identified metabolite/flux, in
  **matched original biological systems** only. Document tissue,
  mass-spectrometry coverage, timing and allele specificity.
- Predict **independently measured enzyme/transporter activity**
  or drug-cell PD, NOT a training label derived from the same abundance
  measurement. Enforce heldout protein residue / allele and, for
  clinical human samples, donor/hospital split.
- Use MAE/RMSE/ranks for continuous function; AUROC/PR-AUC,
  calibration and balanced error for **independently defined binary
  functional or toxicity groups**. Do not report drug-exposure AUC
  from a classification AUROC.
- Evidence grades: measured function > validated isoform-specific
  protein and substrate transport > indirect abundance/omics proxy >
  sequence consequence hypothesis > biological network association.

### Lane 2: Omics2PK (variant effect/assay relevance → clinical PK/PD prediction)

- Begin with PK process-specific physiological knowledge (CYP3A5/4
  hepatic and intestinal tacrolimus CL/F, OCT2–MATE metformin CLr,
  OATP1B1 hepatic rosuvastatin uptake).
- Baseline *actual* clinical covariates and complete original PGx,
  including allele status and original sample timing; allow **zero**
  omics. Compare original gene-type, codon/protein, protein-function
  and pathway omic hypotheses with full black-box all-feature arms.
- Only feed **pre-landmark** omics, with independent source-marked
  confidence. A drug-exposure-derived metabolite or postdose activity
  cannot serve as a “predose omic” without leakage.
- Train view selection **inside** donor-disjoint CV. For tacrolimus
  evaluate actual C0/AUC/CL/F as identifiable from sampling and dose;
  include log ratio/fold bias, MAE/RMSE, extreme errors, target
  attainment, and dosing harm. For PD focus on independently observed
  clinical toxicity/efficacy with calibration and decision utility.
- Assess incremental benefit over baseline and **against** black-box
  alternatives, including model selection instability across genotype,
  assay, tissue and substrate strata.
- A valid negative result (more omics worsens finite-sample estimates)
  is informative. A mechanistic inconsistency triggers a
  *hypothesis*, not automatic causality or elimination of a useful
  unassigned feature.

### Minimal compatible exchange (optional sidecar, not a schema-breaking change)

Keep **MultiOmics Research v1 and its existing `DOWNSTREAM_EXCHANGE_V1`
contract frozen**. Add a *versioned optional* annotation sidecar
for prospective research, conceptually `variant_omic_evidence_v1`:
`source_variant_id`, `assembly`, `transcript`, `consequence`,
`assay_feature_id`, `tissue`, `sample_landmark`,
`assay_method`, `isoform_or_peptide_specificity`,
`molecular_identity_confidence`, `source_evidence_url`,
`measured_vs_inferred`, `functional_endpoint`,
`evidence_level`, `conflicts_and_abstentions`.
Patient identifiers, raw data and patient-level predictions must remain
outside the public GitHub repository.

This adapter is a **design**, not a shipped MultiOmics implementation.
Biological records may be imported from Pillar A evidence-curated
sources. Pillar C compares predicted mechanism vs observed
generalization; Pillar B supplies model building and Pillar D
contextualized clinical state.
