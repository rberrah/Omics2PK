# Wegler 2021 — the per-donor intrinsic clearances are accessible, and they correct our strongest number

**Version 1.0.0 — 2026-11-10.** Evidence types: **source qualification** and
**independent replication attempt**. Axis B.
Artefacts: `reports/data/wegler_table_s3_per_donor_kinetics_audit.json`,
`reports/data/wegler_protein_clint_replication.json`.
Scripts: `scripts/audit_wegler_table_s3_per_donor_kinetics.py`,
`scripts/compare_wegler_protein_clint_against_published_and_transcript.py`.

---

## 1. The answer to the question that was asked

> *Search the original sources, the supplements and the associated
> repositories for the individual CLint values. If they are really accessible
> and matchable, this analysis can provide an important additional functional
> proof. Otherwise keep the blocked status and document exactly what is
> missing.*

**They are accessible and they are matchable.** The registry's `PARTIAL` status
for Wegler is wrong and is corrected.

**Table S3 of the supporting information**, "Kinetics and proteomics data for
the 15 donors", publishes for **each of the 15 donors** (`UU1`–`UU15`) and each
of **five probe substrates**:

| Column | Meaning |
|---|---|
| `cyp_pmol_per_mg_HLM` / `_HH` / `_HL` | the matching CYP's concentration in microsomes, hepatocytes and homogenate (pmol/mg protein) |
| `cyp_pmol_in_incubation_HLM` / `_HH` | the CYP amount actually in the incubation |
| `kpuu_HH` | unbound intracellular-to-extracellular ratio |
| `clint_mic_eq7` / `_eq8` / `_eq9` | **intrinsic clearance in microsomes**, three scaling equations (mL/min/kg bw) |
| `clint_hep_eq7` / `_eq8` | **intrinsic clearance in hepatocytes**, two equations |

Five substrates × 15 donors × 11 columns = **825 per-donor values recovered**,
with `ND` recorded explicitly where the article reports none (bufuralol
`CLint,mic` for two donors).

| Substrate | CYP | `CLint,mic` eq. 7 median | fold range | `CLint,hep` eq. 7 median | fold range |
|---|---|---|---|---|---|
| Midazolam | CYP3A4 | 374.7 | 4.0 | 134.0 | 9.1 |
| Diclofenac | CYP2C9 | 494.1 | 6.1 | 392.3 | 2.4 |
| Omeprazole | CYP2C19 | 39.8 | **64.2** | 16.6 | 32.1 |
| Bufuralol | CYP2D6 | 37.0 (n = 13) | 3.5 | 65.7 | 6.6 |
| Bupropion | CYP2B6 | 24.5 | 8.3 | 68.9 | 2.9 |

And **Data S1** (`mp1c00053_si_002.xlsx`, 5.1 MB) publishes, for the **same 15
donors**, absolute protein concentrations by the Total Protein Approach
(fmol/µg protein, which is the same unit as pmol/mg):

| Preparation | donors | proteins |
|---|---|---|
| Microsomes | 15 | 5,172 |
| Homogenate | 15 | 4,991 |
| Hepatocytes | 15 | 4,784 |
| Low-speed ("discard") pellet | 5 of 15 | 4,682 |

This is better than the PRIDE deposit the registry recorded: PXD014131 carries
MaxQuant output, while Data S1 carries the authors' own absolute
quantification, already donor-labelled `1`–`15`, matching the `Sample_ID`
numbering of the PRIDE experimental-design files.

**So Wegler is, as far as this project has been able to establish, the only
public source holding same-donor proteome and drug-specific intrinsic
clearance.** It is also very small: **15 donors.**

### 1.1 What is *not* available, exactly

- **No per-donor demographics.** Table S1 gives characteristics as aggregates
  only: age mean 64 (39–79), BMI mean 25.7 (20.1–32.9), 10 male and 5 female,
  13 colorectal cancer, 1 breast cancer, 1 GIST. **No per-donor age, sex, BMI
  or diagnosis.** No demographic baseline block can be built for these donors,
  which changes what an Axis B contrast can look like (§5).
- **No transcriptome.** These donors have no RNA, so Wegler cannot answer
  "does protein explain variability the transcriptome does not capture" within
  one cohort. That question stays cross-cohort and therefore descriptive.
- **No clinical pharmacokinetics.** The endpoint is an in-vitro intrinsic
  clearance, not an exposure in a patient.
- **No per-donor genotype.**
- Data S2 (`mp1c00053_si_003.xlsx`) is DAVID enrichment output, not data.

---

## 2. How it was recovered, and why the recovery needed verifying

The article is open access (PMC8041379). The Europe PMC supplementary endpoint
throttled to a stop after delivering all three supporting-information files but
before the trailing figures and the archive's end-of-central-directory record,
so the archive would not open. The three files were recovered from the
truncated stream by reading each member's local header and inflating it;
`mp1c00053_si_001.pdf` has SHA-256
`5beaab1117a32bb281734504dc19237bd264b5db0857837dfd794a9319892bab`.

**A text dump of the PDF is not usable, and this matters.** `pdftotext
-layout` breaks the row correspondence of Table S3: the clearance cells sit
slightly above the batch labels, so a text-order reading assigns one donor's
clearance to another, silently, and leaves eight donors with no clearance at
all. Anyone transcribing this table from a text dump would produce a plausible
file with mis-assigned donors.

The table was therefore reconstructed from **word coordinates** — rows
clustered by the vertical centre of each word, columns assigned by eleven
left-edge bands read off the page geometry — with a hard failure if any row
does not yield exactly eleven values or if two values land in one column.

### 2.1 The verification, which is the part that makes this trustworthy

The article publishes **Table S4**: Pearson (on log-transformed values) and
Spearman correlations between each substrate's clearance and each CYP's
concentration, in both systems. The five pairings of a substrate with its own
CYP give **twenty published anchors**, every one of which would move far more
than any plausible tolerance if a single donor row were mis-assigned.

**All twenty reproduce.**

| Substrate | system | coefficient | published | recomputed | difference |
|---|---|---|---|---|---|
| Midazolam | HH | r | 0.94 | 0.9396 | 0.0004 |
| Midazolam | HH | rs | 0.89 | 0.8893 | 0.0007 |
| Midazolam | HLM | r | 0.84 | 0.8362 | 0.0038 |
| Midazolam | HLM | rs | 0.76 | 0.7643 | 0.0043 |
| Omeprazole | HH | r | 0.88 | 0.8936 | 0.0136 |
| Omeprazole | HH | rs | 0.90 | 0.8964 | 0.0036 |
| Omeprazole | HLM | r | 0.32 | 0.3251 | 0.0051 |
| Omeprazole | HLM | rs | 0.47 | 0.4714 | 0.0014 |
| Diclofenac | HH | r | 0.73 | 0.7492 | 0.0192 |
| Diclofenac | HH | rs | 0.75 | 0.7679 | 0.0179 |
| Diclofenac | HLM | r | 0.18 | 0.1772 | 0.0028 |
| Diclofenac | HLM | rs | 0.11 | 0.1107 | 0.0007 |
| Bupropion | HH | r | 0.40 | 0.4039 | 0.0039 |
| Bupropion | HH | rs | 0.37 | 0.3345 | 0.0355 |
| Bupropion | HLM | r | 0.33 | 0.3310 | 0.0010 |
| Bupropion | HLM | rs | 0.15 | 0.1464 | 0.0036 |
| Bufuralol | HH | r | 0.13 | 0.1271 | 0.0029 |
| Bufuralol | HH | rs | 0.36 | 0.3728 | 0.0128 |
| Bufuralol | HLM | r | 0.14 | 0.1279 | 0.0121 |
| Bufuralol | HLM | rs | 0.41 | 0.4121 | 0.0021 |

Eighteen of twenty agree to within 0.02; the worst is 0.0355. Independently,
the Results quote a median midazolam `CLint,mic` of **375** mL/min/kg bw; the
recovered table gives **374.7** across all 15 donors.

The audit script **refuses to emit the table** if any anchor disagrees by more
than 0.05 or if the median check fails. This is the project's own standard
applied to a transcription: a value that cannot be verified against its source
is not a value.

---

## 3. Licence and handling

Europe PMC records the licence as **CC BY**. The supporting-information PDF
carries no licence statement of its own, so that is the repository's assertion
and not a term read from the document — and **silence about terms is not
permission**. The per-donor rows are therefore kept **out of this repository**
regardless of the licence: the audit script refuses to write them to any path
under the repository root, writes them to the local data directory instead, and
commits only aggregates plus the SHA-256 of the canonical per-donor table
(`750a43095e223bff…`), which is enough to prove a later run recovered the same
825 values without publishing them.

---

## 4. What this does to our strongest number

**This is a correction, and it goes against the funding case's previous
headline.**

On 2026-11-09 this project recorded, as the number a funding request would rest
on: median transcript→activity Spearman **0.508** in the Human Liver Cohort
(n ≈ 390–398 per enzyme) against a median published protein→activity Spearman
of **0.944** (Grangeon 2025, n = 16 microsomal donors), with the published
protein value above the HLC bootstrap interval on **9 of 9** evaluable enzymes.
Three explanations were recorded as inseparable from the public record:
biology, assay precision, and a shared microsomal preparation.

Wegler supplies what Grangeon could not: **per-donor values that can be
recomputed.** Grangeon publishes no per-donor number at all, so its
correlations cannot be checked, cannot be given an interval, and — a detail
that matters more now — its activity values were *provided by the commercial
supplier*, not measured by its authors.

Recomputed from the verified Wegler table, with a 2,000-draw donor bootstrap:

| Substrate | CYP | Wegler HLM ρ [95 %] | Wegler HH ρ [95 %] | Grangeon HLM | HLC transcript |
|---|---|---|---|---|---|
| Midazolam | CYP3A4 | **0.764** [+0.330, +0.953] | 0.889 [+0.620, +0.996] | 0.817 | 0.586 |
| Omeprazole | CYP2C19 | 0.471 [−0.144, +0.830] | 0.896 [+0.657, +0.989] | 0.950 | 0.439 |
| Bufuralol | CYP2D6 | 0.412 [−0.277, +0.904] | 0.373 [−0.192, +0.819] | 0.829 | 0.548 |
| Bupropion | CYP2B6 | 0.146 [−0.440, +0.631] | 0.335 [−0.300, +0.775] | 0.979 | 0.539 |
| Diclofenac | CYP2C9 | 0.111 [−0.489, +0.696] | **0.768** [+0.439, +0.931] | 0.650 | 0.431 |

- **Microsomal median ρ = 0.412** (0.111 to 0.764), against Grangeon's
  published microsomal median of **0.944**.
- **Hepatocyte median ρ = 0.768** (0.335 to 0.896).
- Grangeon's published value falls **outside** the interval these 15 donors
  support on **2 of 5** enzymes — CYP2C19 (0.950 against [−0.144, +0.830]) and
  CYP2B6 (0.979 against [−0.440, +0.631]).
- **The HLC transcript correlation falls *inside* the interval these 15 donors
  support on 5 of 5.**

That last line is the one that changes the funding narrative. In the
**microsomal** compartment — the like-for-like one, since Grangeon's figure is
microsomal — **these 15 donors cannot tell protein apart from transcript as an
explanation of enzyme function, for any of the five enzymes.** The previous
"9 of 9" statement rested on a comparison in which only one side had an
interval.

**"Cannot tell apart" is not "the same".** An interval wide enough to contain
both values is a statement about this cohort's resolution, not about the two
layers being equivalent; at n = 15 that interval contains almost everything, so
it would be wrong to read non-exclusion as evidence of equality. No sentence in
this project may conclude RNA/protein equivalence from a wide interval. What
the data supports is narrower and is the whole point: **the size of the
difference is not established, and the sample sizes available are why.**

**What is *not* being claimed:**

- Not that Grangeon is wrong. Its arithmetic cannot be checked because it
  publishes no per-donor value; and 3 of 5 of its figures do fall inside
  Wegler's interval.
- Not a like-for-like replication. The endpoints differ — probe **activity**
  in Grangeon and the HLC, intrinsic **clearance** by substrate depletion in
  Wegler — the probes differ for some enzymes, and the preparations differ.
- Not that protein is no better than transcript. **It is that the public record
  does not establish the size of the difference**, and that n = 15 and n = 16
  are why.

**What is claimed, narrowly:** CYP3A4 is the one enzyme where both independent
sources agree on a strong protein→function relationship (Grangeon 0.817,
Wegler microsomal 0.764 [+0.330, +0.953], Wegler hepatocyte 0.889). For
CYP2B6, CYP2C9, CYP2C19 and CYP2D6 the two sources disagree or the intervals
are uninformative. **A single-enzyme number is not a median, and the median is
0.412 in microsomes.**

### 4.1 A compartment effect that is its own finding

Hepatocyte correlations exceed microsomal ones on 3 of 5 enzymes, and
decisively so for CYP2C9 (0.768 [+0.439, +0.931] against 0.111 [−0.489,
+0.696]) and CYP2C19 (0.896 against 0.471). The article's own subject is
exactly this microsome-versus-hepatocyte discordance. For a funding case it
means a measurement choice is not neutral: **which preparation the proteome and
the clearance are measured in changes the apparent protein→function
relationship more than the choice between protein and transcript does, in these
data.**

### 4.2 One internal oddity, recorded not explained

CYP2D6 protein spans **587-fold** across the 15 microsomal donors (0.10 to
58.70 pmol/mg — poor metabolisers are present), while bufuralol `CLint,mic`
spans only **3.5-fold** and is `ND` for two donors. A 587-fold range in the
enzyme with a 3.5-fold range in its substrate's clearance is not consistent
with the enzyme being rate-limiting in that assay, and is the most likely
reason the CYP2D6 correlation is near zero in microsomes. This project has not
investigated it further and does not claim an explanation.

---

## 5. What Axis B can now be — and what it still cannot

**Unblocked.** The design that the data now permits, stated so that its
weakness is visible before anyone runs it:

- **Endpoint:** per-donor `CLint,mic` (eq. 7), one model per substrate.
  `drug_specific_clint` under the readiness gate's vocabulary.
- **Baseline:** the target enzyme's **own** absolute protein concentration,
  per donor, from Table S3. This is the honest baseline and the only one
  available, since there are no per-donor demographics (§1.1).
- **Incremental arm:** the baseline plus the rest of the microsomal proteome
  (5,171 further proteins from Data S1), reduced inside training folds.
- **The question:** does the proteome beyond the target enzyme add to the
  target enzyme's own abundance in predicting its substrate's intrinsic
  clearance? This is the hepatic analogue of the F-series question and the
  Axis B question the commissioning note named.
- **Replication within the source:** five substrates, and a second preparation
  (hepatocytes) in the same donors — but **five substrates in 15 donors are not
  five independent validations**, and the two CYP3A4-adjacent readings are not
  two cohorts.
- **Controls required:** donor permutation, and — given §1.2 of the Axis A
  results — a permutation null for **every** reported contrast, because the
  incremental arm is wider than the baseline and differential shrinkage alone
  produces a positive apparent gain.

**What it cannot become.** n = 15 with leave-one-donor-out gives 14 training
donors against 5,172 features. Any interval will be very wide; the Axis A
experience says the honest prior is that nothing will be detected and that the
permutation null will be the informative output. The result would be
**`exploratory`**, a **functional** proof at best, and **never** a clinical-PK
proof: there is no patient, no dose and no exposure in this source.

**It must not be fused with anything.** The readiness gate already refuses
HLC RNA plus Wegler protein at real counts (466 / 15 / 488, same-donor
intersection 0). Nothing here changes that.

---

## 6. Registry consequences

| Field | Before | After |
|---|---|---|
| `status` | `PARTIAL` | `ACQUIRED` |
| `outcome_kind` | not recorded | `drug_specific_clint` |
| per-donor outcome | "must be extracted from the article and qualified" | **held, verified against 20 published anchors** |
| per-donor protein | PRIDE MaxQuant matrices | **Data S1 absolute quantification, 15 donors, 4,682–5,172 proteins** |
| per-donor demographics | not recorded | **absent — aggregates only** |
| `objectives_supported` | — | Axis B; none of H1–H5 |

---

## 7. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First write-up. Establishes that the per-donor intrinsic clearances are published in Table S3 and recovers all 825 values; records that a `pdftotext` reading mis-assigns donors and that the coordinate reconstruction was verified against twenty of the article's own published correlations and the quoted median; records that no per-donor demographics exist; and **corrects the project's previous headline number** by showing that an independent cohort with recomputable data gives a microsomal protein→clearance median of 0.412 against the single published 0.944, and that the HLC transcript correlation lies inside this cohort's interval on 5 of 5 enzymes. |
