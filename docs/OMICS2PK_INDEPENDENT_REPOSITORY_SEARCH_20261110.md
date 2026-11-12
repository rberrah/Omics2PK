# Searching past the publishers — three usable datasets, and a withdrawn conclusion

**Version 1.1.0 — 2026-11-10.** Evidence type: **source search**.
Script: `scripts/audit_independent_repository_search.py` (fits nothing).
Artefact: `reports/data/independent_repository_search.json`.

---

## 0. The finding, and the claim it withdraws

The candidate qualification of 2026-11-10 searched seven named citations
publisher-first and concluded that step 2 of the commissioning order — *execute
only the comparisons that are genuinely identifiable* — had **"nothing to
run."**

**That conclusion is withdrawn.** A search over eight publisher-independent
repositories, and over the Europe PMC open-access-plus-supplement filter,
found:

| | |
|---|---|
| **A usable three-layer cohort** | per-donor hepatic protein, per-donor plasma metabolite and per-donor measured activity in the **same fifteen donors** |
| **A usable clinical-PK cohort** | **twelve patients** with per-patient venetoclax dose, plasma concentration **and** a plasma CYP3A metabolite — the project's first clinical-PK source outside ARIANA |
| **A usable two-compartment proteome** | 9,027 proteins × 37 jejunum donors and 6,982 × 37 matched liver donors, per donor, with **no** pharmacological endpoint |

The earlier search was too narrow. The public record was not too thin.

That matters more than the datasets themselves, because the funding case rests
on a claim about what public data does **not** contain. A claim of that shape
is only as good as the search behind it, and the previous search stopped one
channel too early.

---

## 1. What was searched

Eight repositories and indexes, with the query strings recorded in the script
so that "we searched the repositories" is checkable rather than asserted:

| Repository | Outcome |
|---|---|
| **PRIDE** | one deposition of the cohort already held; eight depositions of per-donor liver proteomics **with no deposited outcome** |
| **ProteomeXchange** | nothing beyond PRIDE |
| **GEO** | nothing relevant — cell lines, organoids, iPSC-derived hepatocytes, unrelated disease series |
| **BioStudies / ArrayExpress** | **the productive one**, because its `S-EPMC` records wrap PMC supplementary files and serve them over the EBI FTP mirror |
| **MetaboLights** | nothing pairing metabolomics with a drug concentration in the same individuals beyond MTBLS3128, already held |
| **Zenodo** | nothing — reviews and figures |
| **Figshare** | nothing |
| **Dryad** | nothing |

---

## 2. PRIDE: the layer is deposited, the outcome never is

| Accession | Verdict | What it is |
|---|---|---|
| `PXD014131` | `ALREADY_HELD` | **the Wegler cohort itself** — "Human liver homogenate, liver microsomes, low-speed centrifugation pellet and hepatocytes LC-MS/MS from 15 matched donors", submitter Christine Wegler, lab head Per Artursson |
| `PXD020910` | `LAYER_WITHOUT_OUTCOME` | per-donor liver microsome proteomics, Barber group |
| `PXD021025` | `LAYER_WITHOUT_OUTCOME` | adult human livers, Achour / Barber |
| `PXD038776` | `LAYER_WITHOUT_OUTCOME` | the proteomics behind *Clin Pharmacol Ther* 2022, which publishes aggregates |

Finding Wegler in PRIDE confirms the cohort and adds nothing analysable: the
per-donor clearances this project needed came from the article supplement, not
from the repository.

The Manchester group's eight depositions are the structural point. Raw mass
spectrometry for individual donors is deposited; **no per-donor activity or
clearance is deposited anywhere**. A layer without an outcome cannot answer a
protein-to-function question, however many donors it holds.

---

## 3. The access ladder, and the lesson it actually teaches

The brief's instruction was right: an HTTP 403 does not establish that
individual data is unavailable through another public channel. Tested on the
hardest case:

**Bhatt 2019 (PMC6222000) — UGT abundance *and* activity in 171 liver donors,
with genotype.** The largest paired cohort in this literature.

| Channel | Result |
|---|---|
| publisher (*Clin Pharmacol Ther*, Wiley) | HTTP 403 to every automated client |
| Europe PMC `supplementaryFiles` | refused: *"Article with id PMC6222000 is not open access one"* |
| NCBI author-manuscript download endpoint | **proof-of-work challenge** (`cloudpmc-viewer-pow`, difficulty 4) |
| NCBI author-manuscript bulk collection | requires registration |
| **BioStudies → EBI FTP mirror** | **served without restriction**, 2,547,704 bytes |

So the file was obtainable, and Europe PMC's own open-access flag was not the
last word.

**And then the lesson inverts.** Opened and read, the supplement gives:

- **Table S2** — abundance and activity of six UGTs as **mean ± SD by age
  category** (neonatal, infancy, early childhood, …);
- **Table S3** — a genuine per-donor row for each of the 171 samples, whose
  `Activity` and `Proteomics` columns are **availability ticks, not values**.

The pairing was measured. The article publishes a per-donor table that records
*which* donors have each measurement and *no measurement*.

**The verdict is therefore `ABSTAIN_NO_PER_DONOR_VALUES`, not
`ABSTAIN_ACCESS`.** Climbing the ladder changed a wrong reason into a right
one. It did not change the outcome. Had the search stopped at the 403, the
project would have recorded an access barrier where the real barrier is a
publication convention — and the funding argument would have been aimed at the
wrong target.

### 3.1 What was not done

NCBI's proof-of-work challenge was **not** defeated. It is an anti-automation
control, and working around it is not a channel this project will use. The
author-manuscript bulk collection requires registration, which the standing
constraint on authorisation-free sources excludes.

---

## 4. The two paywalled sources, re-tested and unchanged

| Source | Result |
|---|---|
| Drozdzik 2018, paired liver and intestine, n = 9 | `ABSTAIN_ACCESS` stands |
| Wenzel 2024, intra-subject jejunum and liver, n = 11 | `ABSTAIN_ACCESS` stands |

Both have **no PMCID at all**, so there is no PMC record for BioStudies to
mirror and no file for the EBI FTP channel to serve. Four channels were tried
for each and are recorded in the artefact.

The audit re-queries both on every run and **fails** if either acquires a
PMCID, because a PMCID would mean the supplement may now be mirrored and the
verdict would have to be revisited. The verdict is not a conclusion; it is a
standing check.

---

## 5. The two usable datasets

### 5.1 Lee / Jackson CYP3A panel — `PMC9926076`, n = 15, CC BY-NC-ND

Per donor, for the **same** fifteen hepatocyte donors:

- hepatic CYP3A4 and CYP3A5 absolute abundance (supplement Table S2);
- plasma 4β-hydroxycholesterol, cholesterol and their ratio (**main-text
  Table 1**);
- midazolam 1′-hydroxylation and ibrutinib M37 / M34 / M25 formation
  (supplement Table S3).

Qualified in full in
`docs/OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_20261110.md`, including
three correlations the article publishes and the extraction reproduces.

**The plasma values are in the main text, not the supplement.** A search that
only harvests supplementary files would have missed the layer that makes this
cohort interesting.

**One cohort, not two.** Lee et al. 2024 (`PMC11645475`) reports venetoclax in
the **same** fifteen donors and the same twenty microsome lots. Two drugs in
the same donors are one cohort with two substrates — the error the brief
explicitly forbids. It is recorded as `ALREADY_HELD`.

### 5.2 Lee 2024 venetoclax patients — `PMC11645475`, n = 12, CC BY-NC

**Two cohorts live in that one article and must be separated.**

Its **in-vitro** half is the same fifteen hepatocyte donors and twenty
microsome lots as §5.1, with venetoclax as a further substrate. One cohort,
two substrates — counting it as a second cohort is the error the brief
forbids, and it is recorded as such.

Its **in-vivo** half is a **different set of people**: twelve patients with
chronic lymphocytic leukaemia carrying, per patient and per visit,

- venetoclax daily dose and plasma concentration (main-text table), and
- plasma 4β-HC, cholesterol and their ratio (**supplement Table S3**).

That is a clinical pharmacokinetic endpoint paired with a plasma metabolite in
the same individuals — the first such source this project has found outside
ARIANA. Executed in
`docs/OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_CONTRAST_20261110.md`, where it
returns a negative result with controls.

**This is the dataset §7 of the first version of this document said did not
exist.** The Europe PMC supplementary-file endpoint delivers only images for
this article; the BioStudies channel delivers the whole 2 MB supplement. The
error was to write a conclusion from the first channel after having just
documented the second.

### 5.3 Matched jejunum and liver proteome — `PMC9310776`, n = 37, CC BY-NC

The supplementary workbook arrives as **11,163,005 bytes over the EBI FTP
mirror** against under 1 MB and a corrupt archive from the Europe PMC endpoint
— the same channel difference as in §3. Inside:

| Sheet | Shape |
|---|---|
| `Proteins_Jejunum_Obese` | 9,027 proteins × **37 per-donor columns**, fmol/µg |
| `Proteins_Liver_Obese` | 6,982 proteins × **37 per-donor columns**, fmol/µg |
| `Proteins_Liver_NonObese`, `ADMEProteins_Liver_NonObese` | no per-donor columns |

The article's own tables publish **medians and ranges** for 37 donors. The
per-donor matrix is in the workbook.

**It has no pharmacological endpoint.** No activity, no clearance, no drug, no
concentration. It can support a cross-compartment redundancy question — does
one tissue's proteome predict another's in the same individual — and it can
support **nothing** about drug metabolism. It is recorded as available and is
**not** executed here, because the standing rule is that each experiment must
answer a new question, and a new question needs a declared target first.

---

## 6. The pattern, now on nine sources

| | |
|---|---|
| sources that **measured** a per-donor layer paired with a per-donor outcome | **9** |
| sources that **published the per-donor rows** | **3** |

The three are **Wegler 2021**, **Lee 2023** and **Lee 2024**, and all three
published rows by accident of format rather than by policy: Wegler's Table S3
happens to be a per-donor table, Lee 2023's plasma values happen to sit in a
main-text table where a word limit did not squeeze them into an aggregate, and
Lee 2024's patient biomarkers sit in a supplement that one retrieval channel
delivers as images and another delivers whole.

Nine studies measured the pairing. Three published it. That is the sentence the
funding case should argue from, and it is stronger than the version this
project held yesterday — because the barrier is a publication convention, which
can be addressed, rather than an absence of measurement, which cannot.

---

## 7. What this search did not settle

- ~~**Twelve real-world patients** in `PMC11645475` have per-patient venetoclax
  dose and plasma concentration in a main-text table, and the article reports
  their plasma 4β-HC/cholesterol ratio **only in a figure**. A clinical-PK
  contrast is **not** unblocked.~~

  **Struck 2026-11-10, within hours of being written. It was an inference, not
  a check, and it was wrong.** Retrieving that article's supplement through the
  BioStudies channel described in §3 — which this document had just finished
  recommending — gives `CTS-17-e70106-s001.docx`, whose **Table S3 publishes
  per-patient plasma 4β-HC, cholesterol and their ratio for all twelve
  patients at both visits**. Together with the main text's per-patient dose and
  venetoclax plasma concentration, that is a **clinical pharmacokinetic
  endpoint paired with a plasma metabolite in the same individuals** — the
  first such source this project has found outside ARIANA. It is qualified and
  executed in
  `docs/OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_CONTRAST_20261110.md`.

  The lesson is the same one as §3, turned on this document: a channel that
  returns images is not evidence that a file holds no table. The failure was
  to write "therefore not available" from a truncated download instead of
  trying the channel that had already been shown to work.
- The **171-donor** Bhatt cohort remains the single most valuable unpublished
  per-donor table identified anywhere in this project. It is the obvious target
  of a data request, which the standing constraint excludes and which the
  commissioning note has deferred.
- The eight PRIDE depositions could in principle be reprocessed into per-donor
  abundances. That would produce a layer and still no outcome, so it would not
  answer any question this project has.

---

## 8. Change log

| Version | Date | Change |
|---|---|---|
| 1.1.0 | 2026-11-10 | Second issue, hours after the first. Adds §5.2: the twelve venetoclax patients, which §7 of version 1.0.0 had wrongly concluded were unusable. Three usable datasets, not two; three sources publish per-donor rows, not two. |
| 1.0.0 | 2026-11-10 | First issue. Withdraws the "nothing to run" conclusion of `OMICS2PK_CANDIDATE_PAIRED_LIVER_SOURCE_QUALIFICATION_20261110.md` §4 and records the access ladder together with the fact that climbing it did not change the 171-donor verdict. |
