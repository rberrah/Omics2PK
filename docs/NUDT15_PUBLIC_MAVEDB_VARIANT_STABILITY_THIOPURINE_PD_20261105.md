# NUDT15: measured variant protein stability versus independent thiopurine cellular PD

**2026-11-05 — original-source reproduction, research E1/PD-cell only, NOT human plasma PK or clinical toxicity.**

## Authoritative source and data qualification

- Suiter et al., *PNAS* 2020, DOI [10.1073/pnas.1915680117](https://doi.org/10.1073/pnas.1915680117), [PMCID PMC7071893](https://pmc.ncbi.nlm.nih.gov/articles/PMC7071893/).
- MaveDB original meta-analysis [`urn:mavedb:00000055-0-1`](https://www.mavedb.org/score-sets/urn:mavedb:00000055-0-1), merging protein abundance/stability VAMP-seq and thiopurine sensitivity **separately measured** in engineered HEK293T cellular constructs.
- Official [MaveDB scores endpoint](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000055-0-1/scores), pinned source SHA-256 `2e248907c0dcf54372d6ebf79c4d8a91b9ec2f2cd50d6c8aa365abec39d85a59`.
- Original source contains **3,015** rows. Strict original sample qualification retains **2,844** positions-parseable, complete, unique single-AA **missense** constructs with both measured scores, at **163 unique protein positions**; excludes **170** rows missing one of the scores and **1** remaining non-missense protein identity (`p.=` reference). No ambiguous row is converted to a fabricated missense. The earlier Omics2PK pilot erroneously truncated NUDT15 at residue 163 and excluded all **19 valid missense substitutions at residue 164**. The official protein is 164 amino acids long. This was corrected and rerun on the original source; the recovered **2,844** paired variants now match the published independent benchmark. Do not equate the 3,015 total rows with complete paired missense variants.
- Original numerical rows, individual variant identities, individual predictions and bootstrap loss arrays are transient in CI; none are committed or exported.

Original workflow: [successful corrected live paired MaveDB replay #42016076614](https://github.com/rberrah/Omics2PK/actions/runs/42016076614).

**Prior-work positioning (29 August 2026):** Kırboğa, *G3*, DOI [10.1093/g3journal/jkag248](https://doi.org/10.1093/g3journal/jkag248), already reproduced the 2,844 paired NUDT15 screens, reported 48% variant-level discordance, and tested ESM2/ESM1v/AlphaMissense/structure/contact-graph prediction of biochemical assay scores. **Our replication is a contract and feasibility control, not a new discovery about NUDT15.** Omics2PK's distinct future contribution must be *drug-specific, patient-heldout* assay choice conditioned on genotype, analytical methodology and independently observed PK/PD.

## Experimental biological question

**Does a genetically altered protein's measured abundance predict its independent pharmacological cellular function?**

The outcomes are different:
1. VAMP-seq fluorescence **protein stability/abundance** of a NUDT15 variant; this is not mRNA measurement or measured proteolysis rate.
2. Thioguanine-exposed HEK293T cell survival/resistance from a **second screen** of the same *variant identity*. Lower `sensitivity_score` indicates the drug is more toxic to the cell; **higher score means tolerance**. This is **cell-level thiopurine PD**, **not** direct isolated NUDT15 nucleotide hydrolase activity, not in-vivo measured TG nucleotide pool, not measured human patient marrow toxicity.

Constructs are matched by amino acid variant, not same original human tissue/clinical donor. A discordance does **not** mean an assay is invalid or prove protein degradation versus substrate dysfunction without an orthogonal experiment.

## Matched measurements: real original numerical results

**Fixed method:** Ridge(alpha=30) with three-view representations, same **5 outer folds separated by 162 protein positions** (so no other mutation at a heldout residue can appear in training), normalization fitted **train-only**. No hyperparameters learned on heldout values; phenotype is never used as an input feature.

| Inputs to predict thiopurine cell tolerance | Heldout MAE ↓ | RMSE ↓ | Spearman ↑ | AUROC* ↑ |
|---|---:|---:|---:|---:|
| Fold training mean | 0.305701 | 0.372867 | −0.142310 | 0.429567 |
| Simple physicochemical amino-acid descriptors | 0.300793 | 0.371984 | 0.107423 | 0.571698 |
| **Measured protein stability/abundance** | **0.271991** | **0.345336** | 0.308120 | 0.684932 |
| Stability + simple AA descriptors | 0.274854 | 0.348773 | 0.319197 | 0.689862 |

\*AUROC uses an **illustrative score ≥0.8**, not independent clinical toxicity or a validated therapeutic threshold. It is a classification AUROC, **never a pharmacokinetic AUC**.

**Descriptive paired bootstrap by protein position**, 999 re-samples, no model refitting/multiplicity correction:

- Protein stability versus training mean: ΔMAE **−0.032695**, 95% exploratory interval **[−0.043037; −0.022045]**.
- Stability + naive amino-acid descriptors versus stability alone: ΔMAE **+0.002022**, interval **[−0.002845; +0.006646]**. **No convincing additional predictive improvement** from those naive sequence features.

The modest heldout gain from protein stability is important, but this **does not mean protein stability is sufficient** to infer drug sensitivity from a missense variant. Nor do these results exclude stronger sequence/structure-based models; the sequence arm is deliberately minimal.

### Original variant-level contradictory molecular profiles

Prespecified *exploratory* thresholds on the source's normalized scoring scales:

- **131/2,844 variants**: measured protein stability ≥0.8 while measured thiopurine tolerance <0.3. **Protein apparently stable, cell highly sensitive to thiopurine.** Candidate loss of function through catalysis, nucleotide binding, localization or other process: **not directly proven from these two screens alone**.
- **222/2,844 variants**: measured protein stability <0.3 while measured thiopurine tolerance ≥0.8. **Low protein abundance despite preserved resistance**; important falsification of an automatic `low protein → ineffective drug protection` rule. Requires clone-level functional and assay validation.

These thresholds are *descriptive*, not a CPIC NUDT15 clinical allele function classifier. No inference that all missense variants reduce function or that NUDT15 protein is degraded by the proteasome in the patient.

Original Suiter 2020 reports **17 independent residues** at which thiopurine response shifts without a measured destabilization effect; this replay complements, **does not claim novelty over**, that peer-reviewed finding.

## Consequences for a unified MultiOmics and Omics2PK method

**MultiOmics E1:** source-qualified mutation → measured protein abundance → **separately observed pharmacological cell response**. A discordance can prioritize a functional assay versus additional protein quantitation. An optional evidence sidecar should distinguish:
- `MEASURED_PROTEIN_ABUNDANCE` (observed)
- `MEASURED_CELL_THIOPURINE_TOLERANCE` (observed and separate)
- `HYPOTHESIS_SUBSTRATE_SPECIFIC_FUNCTION_DEFECT` (hypothesis)
- `LOSS_OF_FUNCTION_PROVEN_IN_PATIENT=False`
- `PROTEIN_DECAY_RATE_MEASURED=False`.

**Omics2PK E2/E3:** the direct consequence to test is whether the functional signal is more useful than genotyping alone for *independent, real clinical dosing/exposure or toxicity*. This public screen does **NOT** have human NUDT15 genotype/omics/plasma thiopurine PK/observed patient toxicity in the same people. Do **not** train real clinical dosing or claim Omics2PK superiority from the cell-line scores. For a real cohort compare clinical+true complete NUDT15/TPMT diplotypes vs addition of source-timed RNA, protein and independently measured TG nucleotide pools on *matched patients*; predeclare marrow toxicity endpoints/calibration and external validation.

**Code:** `src/omics2pk/public_nudt15_variant_pd_benchmark.py`,
`scripts/benchmark_public_nudt15_cell_PD.py`, synthetic tests `tests/test_public_nudt15_variant_pd_benchmark.py`, source schema audit `scripts/audit_nudt15_mavedb_source.py` and two source CI workflows.

**Status:** Qualifies the E1 **biological function versus omic-proxy hypothesis** using authentic public constructs, with an independently measured *cellular PD* response. Full patient PK/PD validation remains open.
