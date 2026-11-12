# Qualification of three candidate public sources, and what each can carry

**Date:** 2026-11-09
**Status: source audit. No model was fitted except one descriptive correlation (§4), which is labelled.**
Every number is a read of an authoritative endpoint — a repository API, a
deposited file, or a published table — not of an abstract.

**Summary of verdicts:**

| Source | Individual-level data public? | Verdict | Which axis it can carry |
|---|---|---|---|
| **EATRIS-Plus v4** (Zenodo 19136013) | **yes**, 12 modalities | **USABLE** | A (multi-omic method), **not** B or C |
| **Wegler 2021** (PRIDE PXD014131) | **yes**, better than the registry recorded | **USABLE** | B (functional) |
| **Grangeon 2025** (Pharmaceuticals 18(12):1789) | **no** — genotype only | **COMPARATOR ONLY** | B, as a published benchmark |

---

## 1. EATRIS-Plus v4 — usable, and the first multi-layer cohort in this project's reach

**Zenodo record 19136013**, DOI `10.5281/zenodo.19136013`, published 2026-03-20,
title *"EATRIS-Plus multi-omics data of a human reference cohort"*.
**Licence MIT. `access_right: open`, files unrestricted.** No account, no
application, no terms click-through.

Five files are deposited. The one that carries the data:

| File | Size | SHA-256 (first 16) |
|---|---|---|
| `EATRIS_PLUS_mudata.h5mu` | 141,704,166 B | `a764a6cd7246edbd` |
| `mae_experiments.h5` | 116.33 MB | (not downloaded) |
| `eatris_isa_filled.xlsx` | 211,761 B | ISA metadata |
| `eatris_isa_filled_JH.ttl` | 3.47 MB | RDF serialisation |
| `mae_mae.rds` | 2.22 MB | R MultiAssayExperiment |

### 1.1 What is actually in the container

Enumerated from the HDF5 directly — **twelve modalities**, 127,620 features in
all:

| Modality | individuals | features |
|---|---|---|
| `acylcarnitines` | 125 | 33 |
| `aminoacids` | 125 | 46 |
| `fattyacids` (VLCFA) | 125 | 5 |
| `lipidomics-neg` | 125 | 164 |
| `lipidomics-pos` | 125 | 196 |
| `proteins` (targeted) | 125 | 1,393 |
| `lipidomics-spline` | 124 | 257 |
| `proteomics-ipp` | 123 | 6,831 |
| `EM-seq` (methylome) | 117 | 100,000 |
| **`mRNA`** | **90** | 16,386 |
| `miRNA-mature` | 89 | 2,130 |
| `miRNA-qRT-PCR` | 86 | 179 |

Person-level variables (24): `Age`, `Sex`, `BMI`, `Height`, `Weight`,
`Smoking.Status`, `ABO.Blood.Group`, `Rh.Blood.Group`, eleven haematology
indices (`Leukocytes`, `Neutrophils`, `Lymphocytes`, `Hematocrit`,
`Erythrocytes`, `Hemoglobin`, `MCV`, `MCH`, `MCHC`, `MXD`, `Platelets`),
`Genomic_Abberation`, `PhenotypicFeature`, `subject.id`. In the ISA workbook
these are complete for 127/127 on demographics and 122–126/127 on haematology.

### 1.2 Two facts the metadata gets wrong, and the number that matters

**The ISA workbook claims 128 rows (127 individuals) for every assay sheet.
The deposited data does not agree.** The container's union is **125**
individuals, and `mRNA` covers **90**, `miRNA-mature` 89, `miRNA-qRT-PCR` 86.
Anyone sizing an analysis from the ISA sheet would over-count by up to 41
individuals on the transcriptome.

The complete-case intersection, measured:

| Layers retained | Complete-case n |
|---|---|
| all **12** | **55** |
| 11 (drop `miRNA-qRT-PCR`) | **79** |
| 10 (drop both miRNA layers) | 79 |
| 9 (drop both miRNA **and mRNA**) | **114** |

Coverage is clean: 55 individuals have all twelve layers, 32 have eleven, 25
have ten, 13 have nine, and **none has fewer than nine**. The 12-versus-11
layer choice costs or buys 24 individuals, and the 9-layer option buys 59 more
at the price of the transcriptome — which is the layer any gene-level claim
needs. That trade-off is a design decision to make **before** running
anything, not after.

### 1.3 What it cannot carry

**Blood only, healthy reference cohort, single biobank** (all 125 identifiers
carry the `CZC` prefix — one site, so no centre confound and no internal
replication either).

**There is no drug, no dose, no concentration, no time, no enzyme activity and
no clinical endpoint of any kind.** The only outcome-like variables are
demographics and haematology. The user's own statement that this dataset has
no confirmed clinical PK endpoint is correct and in fact understates it: it
has no *pharmacological* endpoint at all.

Two consequences, both binding:

- EATRIS-Plus can support **Axis A** — complementarity, redundancy,
  reconstruction of one modality from the others, and parsimony — because
  reconstruction needs no endpoint and redundancy is a property of the layers
  themselves. Where a supervised target is needed it must be a demographic or
  haematological variable, declared as such.
- It cannot support **Axis B or C**. A leukocyte transcriptome and a plasma
  proteome do not measure hepatic drug-metabolising activity, and no analysis
  on this cohort may be presented as evidence about drug metabolism or
  pharmacokinetics. This is the compartment caveat the commissioning note
  raised, and it is not a hedge: it is the reason this dataset is in Axis A and
  nowhere else.

---

## 2. Wegler 2021 (PXD014131) — usable, and the registry understated it

**PRIDE PXD014131**, *"Human liver homogenate, liver microsomes, low-speed
centrifugation pellet and hepatocytes LC-MS/MS from 15 matched donors"*,
submitted 2019-06-05, public 2021-03-29. **Licence CC0 (public domain).**
Reference: Wegler C, Matsson P, Krogstad V, Urdzik J, Christensen H, Andersson
TB, Artursson P. *Mol Pharm* 2021, PMID 33739838.

The submission is flagged `submissionType: PARTIAL`, which suggests raw spectra
only. **It is not.** Of 100 deposited files, 83 are `.raw` and 17 are
processed MaxQuant output:

| File | Size | What it is |
|---|---|---|
| `Hom_proteinGroups.txt` | 24.14 MB | homogenate, **per-sample intensities** |
| `Mic_proteinGroups.txt` | 24.99 MB | microsomes, per-sample intensities |
| `Hep_proteinGroups.txt` | 22.75 MB | hepatocytes, per-sample intensities |
| `Pellet_proteinGroups.txt` | 9.22 MB | low-speed pellet |
| `*_peptides.txt` ×4 | 25–50 MB | peptide-level evidence |
| `*_experimentaldesigntemplate.txt` ×4 | <1 KB | **raw file → donor mapping** |
| `*_summary.txt`, `*_parameters.txt` | small | run QC and search settings |

The experimental design templates were downloaded and parsed. The donor
structure is published and unambiguous:

| Preparation | donors | runs |
|---|---|---|
| Homogenate | **15** (`1_Homogenate` … `15_Homogenate`) | 30 (two proteases each) |
| Microsomes | **15** | 30 |
| Hepatocytes | **15** | 30 |
| Low-speed pellet | **5** (donors 1, 7, 11, 12, 15) | 15 |

**Verdict: individual-level quantitative proteomics is publicly accessible,
CC0, for 15 matched donors in three complete liver preparations**, with the
donor mapping deposited rather than inferred. The registry's `wegler_15` entry
records the cohort but not that the per-sample quantitative matrices are
retrievable; that is an understatement to correct.

What still needs separate retrieval: the measured intrinsic clearance for the
five drugs, which lives in the article rather than in PRIDE. The registry
records it as present; it is not in this deposit.

---

## 3. Grangeon 2025 — a comparator, not a dataset

**Pharmaceuticals 2025, 18(12), 1789**, PMID 41471278, PMC12735975, **CC BY
4.0**. *"Protein Levels of 16 Cytochrome P450s and 2 Carboxyl Esterases Using
Absolute Quantitative Proteomics."*

MDPI blocks direct supplement download (HTTP 403) and the PMC download route
returns an interstitial page; the file was retrieved through the Europe PMC
`supplementaryFiles` endpoint, the same route this project already uses for
Jarmusch.

### 3.1 The supplement is one PDF, and it was read

`pharmaceuticals-18-01789-s001.zip` → a single 2,253,487-byte PDF, 9 pages:

| Page | Content |
|---|---|
| 1 | **Table S1 — per-donor genotypes, 16 HLM subjects** |
| 2 | Table S2 — proteotypic peptides, MRM transitions |
| 3 | Table S3 — intra/inter-day precision and accuracy |
| 4 | Table S4 — matrix effect |
| 5 | Table S5 — autosampler stability |
| 6–9 | images only (Figures S1–S4) |

Table S1 is the **only** per-donor table in the entire supplement, and it is
thinner than its caption implies:

- `CYP2C9`, `CYP2C19`, `CYP2D6`: called for **HLM01–HLM10 only**. HLM11–16 are
  `Unknown`. So **n = 10**, not 16.
- `CYP3A5`: called for 11 donors — and **every one of them is `3A5*3/*3`**.
  **Zero variance.** As a predictor it carries no information whatsoever.
- The five donors with 2C9/2C19/2D6 calls and the eleven with 3A5 calls overlap
  only partly, so no donor has a complete four-gene genotype plus a usable
  comparison group.

### 3.2 What is not published

- **No per-donor protein abundance.** Table 2 gives N, mean ± SD and min–max.
  Figure 1 plots individual donors as points without listing values.
- **No per-donor activity.** Table 3 is not an activity table at all: it is a
  **correlation** table, with columns `CYP450 | Probe | r | p | N`.
- The mass-spectrometry data is in PeptideAtlas `PASS05929`, not PRIDE, and
  would be raw spectra rather than a quantitative matrix.
- Activity values "were provided by the supplier", not measured by the authors.

**Verdict: Grangeon is not an exploitable individual-level dataset.** Reading
abundances or activities off its scatter figures would be exactly the
"individual tables artificially reconstructed from published means" the mission
constraints forbid, and it is not done here.

What it *is* good for is a published benchmark, which §4 uses.

---

## 4. The consequence: how much of CYP activity does mRNA leave unexplained?

> **SUPERSEDED IN PART, 2026-11-10.** The published protein figure below comes
> from a single source that publishes **no per-donor value**, so it has no
> interval and cannot be recomputed. An independent cohort that *can* be
> recomputed — Wegler 2021 Table S3, 15 donors, verified against twenty of the
> article's own published correlations — gives a **microsomal median Spearman
> of 0.412**, not 0.944, and the HLC transcript correlation falls **inside**
> that cohort's interval on **5 of 5** shared enzymes. The "9 of 9" statement
> below therefore compares a figure with an interval against a figure without
> one, and must not be quoted on its own. See
> `docs/OMICS2PK_WEGLER_PER_DONOR_CLINT_QUALIFICATION_20261110.md` §4. The
> table below is left unaltered because it is what was measured on
> 2026-11-09; only its interpretation is withdrawn.

**Label: `DESCRIPTIVE_CROSS_COHORT_COMPARISON`.** Artefact
`reports/data/rna_versus_protein_activity_correlation.json`, script
`scripts/compare_rna_versus_published_protein_activity_correlation.py`.

Grangeon publishes the Spearman correlation between absolutely-quantified
protein content and probe-substrate activity **in the same donors**. The Human
Liver Cohort allows the same descriptive statistic for the **transcript**. Both
are whole-cohort correlations with no model, no selection and no tuning, so
they are like-for-like in form.

| Endpoint | mRNA ρ (HLC) | 95 % CI | protein r (published, HLM) | gap | protein r above the CI? |
|---|---|---|---|---|---|
| CYP3A4 midazolam | 0.586 | [0.514, 0.652] | 0.8168 | 0.230 | yes |
| CYP2D6 | 0.548 | [0.467, 0.621] | 0.8291 | 0.281 | yes |
| CYP2C8 | 0.542 | [0.462, 0.613] | 0.9176 | 0.375 | yes |
| CYP2B6 | 0.539 | [0.460, 0.611] | 0.9794 | 0.440 | yes |
| CYP2A6 | 0.508 | [0.430, 0.586] | 0.9471 | 0.439 | yes |
| CYP3A4 testosterone | 0.498 | [0.417, 0.574] | 0.9559 | 0.458 | yes |
| CYP2C19 | 0.439 | [0.355, 0.518] | 0.9500 | 0.511 | yes |
| CYP2C9 | 0.431 | [0.338, 0.513] | 0.6500 | 0.219 | yes |
| CYP2E1 | 0.313 | [0.216, 0.407] | 0.9441 | 0.631 | yes |
| CYP1A2 | `NOT_EVALUABLE` | — | 0.6814 | — | — |

- **Median mRNA ρ = 0.508** (0.313 to 0.586), n ≈ 390–398 per enzyme.
- **Median published protein r = 0.944** (0.650 to 0.979), n = 16.
- **On 9 of 9 evaluable enzymes the published protein correlation lies entirely
  above the bootstrap interval of this cohort's transcript correlation.**

> ~~In rank-variance terms, the median enzyme: mRNA accounts for about 26 %,
> protein for about 89 %.~~
>
> **WITHDRAWN, 2026-11-10 — this sentence was not a valid statistic.** It
> squares a **Spearman** coefficient and reads the result as a fraction of
> variance explained. The square of a Pearson coefficient is a variance
> fraction; the square of a rank correlation is not, because the ranks are not
> the quantity whose variance is being apportioned. No percentage of explained
> variance may be derived from any ρ in this project, and none appears
> anywhere else. The correlations themselves (the table above) stand as
> measured; only this derived reading is struck.

`CYP1A2_activity` is `NOT_EVALUABLE` here for the same reason as in the
F-series: none of its four probes is admissible under the frozen identity rule.
Grangeon's protein r for it (0.6814) is reported but has nothing to be compared
against.

### 4.1 Three explanations, and why the project cannot choose between them

The 9-of-9 consistency makes a cohort accident unlikely, but **at least three
mechanisms would produce this gap and the public record cannot separate them**:

1. **Biology.** Protein abundance is closer to catalytic capacity than mRNA is;
   translation and turnover decouple the two. This is the interesting
   explanation and the one a funding case would rest on.
2. **Assay precision.** Grangeon used targeted MRM absolute quantification;
   HLC used a 2008-era Agilent microarray. Measurement error attenuates a
   correlation, and the two platforms are not comparable in precision. **Part
   of this gap is certainly instrumental.**
3. **Shared preparation.** Grangeon measured protein *and* activity in the same
   microsomal preparation, so both inherit that preparation's yield and
   membrane integrity — a shared nuisance that inflates their correlation. HLC's
   transcript comes from tissue RNA while activity comes from microsomes, so
   they share no preparation step and no such inflation.

Explanations 2 and 3 both push in the same direction as 1, which is precisely
why **this comparison motivates a measurement rather than settling a question**.
Separating them requires mRNA and protein measured on the **same** donors with
a common activity assay. No public dataset has that: HLC has the transcriptome
in 466 donors, Wegler the proteome in 15, Grangeon the proteome in 16 — all
different people.

That is the experimental justification for paired proteomics, stated as a
number rather than an intuition: **the layer this project cannot measure
accounts for roughly three times the activity variance of the layer it can,
and three mechanisms compete to explain it.**

---

## 5. What this changes in the registry

Four amendments, all additive:

1. Register **EATRIS-Plus v4** with its measured modality table, the
   complete-case curve (55 / 79 / 114), the `h5mu` checksum, MIT licence, and
   `objectives_supported: []` — it supports no current objective, because none
   of H1–H5 admits a blood-only cohort with no pharmacological endpoint.
2. Record the **ISA-versus-data discrepancy** (127 claimed per assay; 125 union,
   90 for mRNA) as a caveat on that entry.
3. Correct **`wegler_15`** to record that per-sample quantitative matrices are
   retrievable under CC0, with the four-preparation donor counts (15/15/15/5),
   and that the intrinsic-clearance values are *not* in the deposit.
4. Register **Grangeon 2025** as a `published_comparator` and not a dataset,
   with the measured reasons: per-donor genotype for 10 donors, zero-variance
   CYP3A5, no per-donor abundance or activity.

---

## 6. What should happen next, in order

The commissioning note's own order is the right one and this audit does not
change it, only its inputs.

1. **Execute the readiness and ablation modules on the real accessible data**
   — `hepatic_multiview_readiness` and `multimodal_pk_ablation` exist with unit
   tests but have never been run against a registered source, so the project
   has no recorded verdict. This is the cheapest remaining deliverable and it
   needs no new data and no new method.
2. **Axis A on EATRIS-Plus**, with the layer/n trade-off chosen and frozen
   first, and every target declared non-pharmacological.
3. **Axis B on Wegler**, which now has retrievable individual proteomics; the
   intrinsic-clearance values must be extracted from the article and qualified
   separately before any protein→CLint analysis.
4. **Axis C unchanged** — ARIANA, Jarmusch, Rotroff, with the frozen protocols
   and the negative results kept as they stand.

What should **not** happen: no analysis joining individuals across cohorts, no
Grangeon figure digitisation, and no statement about drug metabolism derived
from EATRIS-Plus blood layers.
