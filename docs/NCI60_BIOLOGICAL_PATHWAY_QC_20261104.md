# Omics2PK — biological pathway and chemical-identity QC: original NCI-60 data
**Source-verified executed analysis:** 4 November 2026 · [successful workflow #41850397076](https://github.com/rberrah/Omics2PK/actions/runs/41850397076) · artifact original aggregate-only `13334327146` · PR #137 **research/draft**.

## Central interpretation for Pierre Marquet

The previous NCI-60 **RNA > automatic fusion** result is not by itself evidence that the transcriptome is intrinsically the best layer for a drug. A fair experiment must evaluate the **source measurement, chemical identity confidence, covered substrate/reaction biology, directly observed molecular assays and missing parts of the biological cascade**.

This audit now distinguishes:
- DNA genotype (not measured in the original NCI-60 source package) versus qPCR/gene-expression (transcript) versus SWATH peptide-based protein abundance (may miss isoforms, truncations or active-site function) versus metabolite pools.
- **direct doxorubicin pharmacological mechanisms** (TOP2A/TOP2B, ABC transport, anthracycline reduction) versus *indirect* redox/central-carbon/nucleotide metabolism.
- source-labeled, named candidate metabolites versus source-ambiguous/isobar/artefact features.
- an **experimentally complete cascade** (not yet available) versus a partial mapping of selected measured entries in a named pathway.

A positive/negative per-pathway predictive result **does not validate the causal pathway**, confirm GSH/GSSG redox flux, predict transporter function, or measure human drug PK.

### Original source and audit contract

**Metabolite names:** `Mathelab/NCI60_GeneMetabolite_Data@4c140708c18c41910458414a626b1a51d3dd2a85`, original `fData.metab.csv`, required to match exactly the source `metabData.csv` measured feature columns. The original IntLIM version filtered source Metabolon technical replicates to median CV below 0.3, retaining 280 log2-transformed features. The CSV `fData.metab.csv` has only `metabnames,id`; the two values are the **same human-readable text**, not standardized InChIKey/ChEBI/HMDB accession or MS/MS/standard-level proof.

**Measured response:** prespecified *in-vitro* doxorubicin GI50 from `mixOmics-org/mixOmics@9f9fd66b69105fdb618b54588b638b567e85c0d6`. Source cell lines rather than same-time aliquots: **47 lineage-concordant original NCI-60 identities**.

**Measured proteome:** original 2019 Guo et al. Table S1E, [10.1016/j.isci.2019.10.059](https://doi.org/10.1016/j.isci.2019.10.059), publisher XLSX SHA256 `e9741865cfb12cced0a471f534966e2bcd2f9c6979f591f454daf9c0206779cb`, 3,171 original protein accessions with the original gene symbol field.

**Standards:** Salek et al. 2013, [10.1186/2047-217X-2-13](https://pmc.ncbi.nlm.nih.gov/articles/PMC3853013/) (MSI levels 1–4). No source-supplied authentic-standard retention + MS/MS matching is present in the reviewed `fData.metab.csv`, so **none of its 111 named entities is claimed MSI level 1**. Original metabolomics provider may have additional metadata not released in this two-field file; *absence of evidence in the reviewed file is not proof that the providers failed to identify their metabolites*.

### What the 280 'metabolites' actually are

| Source annotation status (mutually exclusive, **source metadata only**) | Count | Meaning |
|---|---:|---|
| Chemical **name in source** (MSI identity confidence **unverified**) | **111** | Candidate chemical identity; mapping to a reaction or KEGG/HMDB/ChEBI is not independently verified from this source file |
| Original `X-...` feature | **161** | Unidentified original numbered MS feature, including tentative suffixes; cannot assign pathway via substring speculation |
| Source `Isobar ... includes ...` | **7** | Several candidate compounds share signal; cannot pick a favorable one to assign to a pathway |
| Explicit possible artifact not already `X-...` | **1** | Must not be treated as a biological molecule |

Thus **169/280 (60.4%)** features cannot be assigned a single biochemically evidenced pathway from the public labels alone. The `X` group includes possible/retired signals; the four source labels containing the words artifact/retired overlap the numbered group except the one counted separately.

### Source-checked candidate pathway coverage

*These six compact reaction-neighborhood panels are prespecified test hypotheses, **not** a curated complete Reactome/KEGG reaction graph. Source gene symbols and protein gene symbols are exact-string checks, not confirmed enzymatic function. One metabolite may map to more than one panel.*

| Mechanistic neighborhood | Named candidate metabolites used | RNA genes present / nominated | SWATH protein gene symbols present / nominated | Major missing intermediates |
|---|---:|---:|---:|---|
| Glutathione biosynthesis / redox | **6** | **10/10** | **6/10** | Oxidized glutathione (GSSG), NADPH and NADP+ |
| Glycolysis / pentose phosphates | **6** | **8/9** | **9/9** | 6-phosphogluconate, NADPH/NADP+ |
| TCA / mitochondrial carbon | **5** | **8/9** | **9/9** | Acetyl-CoA, isocitrate |
| Purine / pyrimidine pools | **8** | **7/8** | **7/8** | dATP, dGTP, dCTP; TOP2A abundance absent from SWATH source protein table |
| Methionine / one-carbon | **6** | **7/7** | **4/7** | Tetrahydrofolate and 5-methyl-THF |
| Polyamines | **4** | **5/5** | **2/5** | Spermine |

**High-value direct drug check:** exact source gene symbols were measured for **TOP2A, TOP2B, ABCB1, ABCC1, ABCG2, CBR1 and CBR3** in the gene array. Original SWATH protein table reports matching gene symbols only for **TOP2B, ABCB1, CBR1** among those seven. There is **no quantified doxorubicin or doxorubicinol feature** in the original 280-metabolite panel. No functional TOP2, ABCB1 transport or CBR activity measured. A transcript gene symbol does not establish variant function, full active protein, or the activity of the original drug pathway. Biological pharmacology reference: Thorn et al., *Pharmacogenetics and Genomics* 2011 [Doxorubicin pathways: pharmacodynamics and adverse effects](https://pmc.ncbi.nlm.nih.gov/articles/PMC3116111/).

**Glutathione as concrete counterexample:** `glutathione, reduced`, cysteine, glutamate, glycine, gamma-glutamyl-cysteine and 5-oxoproline are present; **oxidized glutathione and NADPH are not annotated**. A reliable actual redox ratio, NADPH recycling rate, GST detoxification flux, reactive oxygen species level or drug-specific metabolite formation **cannot** be derived merely from these six abundance profiles.

### Original n47 predictive experiment after chemistry QC

All **23 arms** tested on **exactly the same 47 cell-line IDs**, original same GI50, 3 × 5 grouped heldout splits (3 correlated repetitions), fixed Ridge(alpha=10), **3 PCA components learned within each train split** for each omic panel, train-only feature variance ranking and scaling. Original full three-omics and RNA/other comparator MAEs were required to replay the previous source benchmark to within 0.000005 before results were accepted.

| Response model (all cancer-lineage adjusted) | Actual heldout MAE (lower is better) |
|---|---:|
| Cancer lineage only | **0.342488** |
| + RNA only | **0.334666** |
| + original all 280 metabolite signals | 0.355840 |
| + only 111 named candidate metabolites | **0.345659** |
| + 169 ambiguous/isobar/artefact candidate signals | 0.347366 |
| + RNA + all 280 metabolite signals | 0.351537 |
| + RNA + 111 named candidate metabolites | **0.344727** |
| + RNA + all protein signals | 0.341053 |
| + RNA + proteins + all metabolites | 0.358124 |
| + RNA + proteins + 111 named metabolites | **0.350081** |

Per-pathway metabolite ablation:

| **Pathway panel only + tumor lineage** | MAE | **Tumor lineage + RNA + pathway panel** | MAE |
|---|---:|---|---:|
| Glutathione/redox | 0.361753 | RNA + glutathione | 0.351715 |
| Glycolysis/PPP | 0.354961 | RNA + glycolysis/PPP | 0.350484 |
| TCA | 0.354170 | RNA + TCA | 0.345168 |
| Nucleotides | **0.344997** | RNA + nucleotides | **0.338575** |
| One-carbon | 0.359534 | RNA + one-carbon | 0.346685 |
| Polyamines | 0.352598 | RNA + polyamines | 0.345153 |

**All six biochemical candidate panels failed to improve upon tumor-lineage + original RNA-alone MAE 0.334666 in this fixed protocol**. No result should be turned into a blanket statement about these biochemical pathways: each panel lacks several reaction-relevant molecules, source annotations do not certify standards, samples come from separate studies, and molecular source batch effects differ.

### Paired uncertainty — no claim of statistical proof

Exploratory **1,999 resamples of the 47 unique cell-line heldout losses** (repeat predictions averaged within cell-line before resampling; not full model bootstrap; no multiplicity adjustment):

| New model minus comparator | Delta MAE (positive worse) | Descriptive paired 95% CI |
|---|---:|---|
| Named 111 metabolite-only minus all 280 metabolite-only (each tumor-adjusted) | **−0.010181** | **[−0.018241, −0.001674]** |
| RNA + named 111 minus RNA-alone (tumor-adjusted) | **+0.010061** | [−0.002059, +0.023116] |
| RNA + all 280 minus RNA-alone (tumor-adjusted) | **+0.016871** | [+0.005573, +0.027710] |
| RNA + proteins + named 111 minus RNA + proteins | **+0.009028** | [−0.004426, +0.023120] |
| RNA + proteins + all 280 minus RNA + proteins | **+0.017071** | [+0.004287, +0.029630] |
| RNA + proteins + named 111 minus RNA + proteins + all 280 | **−0.008043** | [−0.016471, +0.000833] |
| RNA + glutathione panel minus RNA alone | **+0.017049** | [+0.007707, +0.025845] |
| RNA + nucleotide panel minus RNA alone | **+0.003908** | [−0.020645, +0.027964] |

**No unconditional statistical conclusion** because many arms were examined, bootstrap is descriptive, and the annotation-based and pathway analyses were created *after* the previous overall multiomic benchmark results had been seen. The improved named-only result may be dimensionality, non-specific noise rejection, source/platform artefacts, chemical coverage or chance; not evidence of a particular molecular pathway.

### What this changes scientifically

1. **Genetics versus RNA remains a different, source-qualified HLC question**: same 171 original HLC livers, age/sex baseline MAE 706.771; + one observed CYP3A5 rs776746 SNP 697.803; + observed CYP3A4/5 gene-level RNA 641.842; + RNA+SNP 641.836. **RNA is numerically better than that one polymorphism for shared CYP3A midazolam microsomal activity**, not better than complete CYP3A4/5 germline diplotypes, variants, expression/splicing or tacrolimus plasma PK. [Original HLC data report](HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS_20261018.md).
2. **Layer comparability first:** before fitting an omic, audit specimen/source date, target tissue, primer or array exon, alternative isoforms, proteomic tryptic peptide position/species, substrate, reaction, MS annotation identity grade, specimen relative drug timing, denominator and measured PK outcome.
3. **Do not mistake context for a complete metabolic cascade:** for glutathione, GSH without GSSG/NADPH is an incomplete redox measurement. For doxorubicin, TOP2A and ABC transport are far closer to the drug's direct mechanism than static citrate/lactate abundance, but in the source even those gene/protein abundances are incomplete functional proxies.
4. **No budget needed to qualify more source datasets**, but a meaningful RNA-versus-actual complete PGx comparison or direct same-person CYP3A5 activity/tacrolimus PK **cannot be fabricated** from currently missing genotyping/paired-tacrolimus assays.

## Implementation and reproducibility

- Pure source-label evidence gate: `src/omics2pk/nci60_biochemical_coverage.py`.
- Source-locked real matched benchmark: `scripts/benchmark_original_nci60_pathway_annotation.py`.
- Source original SWATH parser enriched with gene-symbol return (no backward-compatibility change for existing 3-omics benchmark): `scripts/benchmark_original_nci60_threeomics_swath.py`.
- Synthetic-only label guards: `tests/test_nci60_biochemical_coverage.py`.
- [Live original-data replay workflow, SUCCESS #41850397076](https://github.com/rberrah/Omics2PK/actions/runs/41850397076); only JSON aggregate artifact, **no cell-line-specific assay values or predictions committed or uploaded**.
- A full **reaction-edge** graph requires authoritative compound-ID mapping with original MS evidence and enzyme/reaction citations, not lexical match. Existing MultiOmics pathway annotation outputs may be used **only after** preserving evidence/provenance and verifying keys, never substituted for confirmed NCI60 source identities.

**Final interpretation:** the negative NCI-60 threeomic result is scientifically explainable *in part* by source and annotation limitations, but this audit **does not identify the causal reason** for the lack of gain. It demonstrates why one must assess the measured substrate-specific cascade before concluding that a complex biological signal is redundant.
