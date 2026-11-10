# The Lee / Jackson CYP3A panel — a three-layer per-donor public dataset

**Version 1.0.0 — 2026-11-10.** Evidence type: **source qualification**.
Script: `scripts/audit_lee_jackson_cyp3a_panel.py` (fits nothing).
Artefact: `reports/data/lee_jackson_cyp3a_panel.json`.
Per-donor table: **outside this repository**.

Source: Lee J, Fallon JK, Smith PC, Jackson KD. "Formation of CYP3A-specific
metabolites of ibrutinib in vitro is correlated with hepatic CYP3A activity and
4β-hydroxycholesterol/cholesterol ratio." *Clin Transl Sci*
2023;16(2):279–293. PMCID PMC9926076, PMID 36350327, **CC BY-NC-ND**.

---

## 0. Why this source is different from the nine before it

For the **same fifteen** hepatocyte donors it publishes, each at the level of
the individual donor:

| | Where | Unit | Range across donors |
|---|---|---|---|
| CYP3A4 abundance | supplement Table S2 | pmol/mg protein | 4.05 – 545.23 (**135-fold**) |
| CYP3A5 abundance | supplement Table S2 | pmol/mg protein | 0.31 – 106.8 (345-fold) |
| plasma 4β-HC / cholesterol | **main-text Table 1** | ×10⁴ | 0.05 – 0.73 (15-fold) |
| midazolam 1′-hydroxylation | supplement Table S3 | pmol/10⁶ cells/min | 1.32 – 138.23 (**105-fold**) |
| ibrutinib M37 formation | supplement Table S3 | pmol/10⁶ cells/min | 2.46 – 105.51 (43-fold) |
| age, sex, race, BMI | supplement Table S2 | — | 1 – 63 years |

**Two molecular layers and a measured outcome, in the same individuals.** Five
of the seven candidates qualified earlier that day measured a pairing and
published only aggregates. This one publishes the rows.

**The second layer is a plasma measurement.** That is what makes the cohort
worth more than its size: 4β-hydroxycholesterol is obtainable from a blood
sample in a living patient, where CYP3A4 abundance requires liver tissue. The
cohort therefore supports a question no other public source has let this
project ask — whether an accessible layer adds to, or could stand in for, an
inaccessible one.

---

## 1. The extraction is verified against the article's own numbers

An extraction that cannot reproduce the source's published figures is not
evidence of anything. Three anchors, all reproduced:

| Anchor | n | Published | Recomputed | Gap |
|---|---|---|---|---|
| ibrutinib M37 vs plasma ratio, excluding the 1-year-old | 14 | 0.62 | **+0.6176** | 0.0024 |
| ibrutinib M25 vs plasma ratio, excluding the 1-year-old | 14 | 0.67 | **+0.6703** | 0.0003 |
| CYP3A4 abundance vs M37 | 15 | 0.97 | **+0.9679** | 0.0021 |

Tolerance 0.02, which leaves room for the two-decimal rounding of the
published values and for nothing else. `verify_published_anchors` raises if any
gap exceeds it, so nothing downstream can run on a mis-parsed table.

The first two anchors require excluding donor `VYU`, aged 1, because the
article excludes that donor from the correlation it reports. Reproducing a
published number means reproducing the published exclusion, not a convenient
one.

---

## 2. The circularity check, done before any model was fitted

Wegler cost this project two of three candidate endpoints: equation 8 divides
the clearance by the probed CYP and multiplies back by the homogenate CYP
concentration, and equation 9 adds K<sub>puu</sub>. Both would have made a
CYP-abundance predictor predict its own input.

**This source has no such quantity.** The outcome is normalised **per million
cells**; the predictor is an abundance **per mg of protein**. The outcome is
divided by no CYP quantity.

**And that difference has a cost, which is recorded rather than hidden.** The
two are expressed per *different denominators*, so each donor's
protein-per-cell ratio enters any relationship between them as an uncontrolled
term. It weakens the comparison. It does not make it circular.

---

## 3. What the source permits, as rank correlations

Computed on all fifteen donors. **These are in-sample correlations and are not
performance.**

| Relationship | ρ |
|---|---|
| CYP3A4 abundance → ibrutinib M37 | **+0.9679** |
| CYP3A4 abundance → midazolam 1′-hydroxylation | **+0.7286** |
| plasma ratio → midazolam 1′-hydroxylation | +0.5286 |
| plasma ratio → ibrutinib M37 | +0.4357 |
| **plasma ratio → CYP3A4 abundance** | **+0.4071** |
| CYP3A5 abundance → midazolam 1′-hydroxylation | +0.3857 |
| age → CYP3A4 abundance | +0.0751 |
| age → midazolam 1′-hydroxylation | +0.1662 |

Two of these decide whether the cohort is interesting:

**The plasma metabolite and the hepatic protein correlate at only +0.41.** They
are therefore *not* two measurements of one thing. Each may carry information
the other lacks, which is why an incremental contrast between them is a real
question rather than an arithmetic certainty.

**Age explains almost nothing here** (+0.075 against abundance). The two
donors under four are a genuine confounder in principle; in this cohort the
rank relationship with age is close to nil, which weakens — without removing —
the ontogeny objection.

---

## 4. Independence, stated in both directions

**Independent of Wegler.** Different donors (BioIVT letter codes against
Wegler's `UU` codes), different laboratory, different year, different
quantification, and an activity rate rather than a fitted intrinsic clearance.

**Not independent of the venetoclax paper.** Lee et al., *Clin Transl Sci* 2024
(PMC11645475), reports venetoclax metabolism in the **same** fifteen hepatocyte
donors and the same twenty microsome lots. Two drugs measured in the same
donors are **one cohort with two substrates**. The brief forbids counting them
as two, and the artefact records this under `independence`.

---

## 5. What this source cannot do

- **No intrinsic clearance.** The outcome is a formation rate at a single
  substrate concentration — not V<sub>max</sub>/K<sub>m</sub>, not a fitted
  CL<sub>int</sub>, and not on Wegler's scale.
- **No pharmacokinetics.** No dose, no plasma drug concentration, no exposure.
- **No microsomal outcome.** Supplement Table S1 gives abundance for twenty
  single-donor microsome lots and **no activity** for them. The microsomal arm
  Wegler could test is not testable here — which matters, because Wegler's
  microsomal result was the negative one.
- **Procurement plasma.** The 4β-HC comes from deceased organ donors at
  procurement. It is not interchangeable with plasma from a living patient:
  agonal state, intensive-care drug exposure and cholesterol handling all
  differ.
- **Ontogeny is present and cannot be removed.** Ages 1 to 63, two donors
  under four, fifteen individuals.
- **CC BY-NC-ND.** Analysis and the reporting of aggregates are unaffected.
  Redistributing the per-donor table or a derivative of it is not permitted, so
  it is written outside this repository and never committed — which is this
  project's practice regardless of licence.

---

## 6. Status

| Field | Value |
|---|---|
| status | `ACQUIRED` |
| access tier | `open_public` |
| n analysed | 15 |
| outcome kind | `drug_specific_metabolic_activity` |
| axis | **B — functional** |
| ceiling | exploratory; an independent replication at n = 15 |

The executed analysis is reported in
`docs/OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_RESULTS_20261110.md` under
the protocol frozen in
`docs/OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_PRESPECIFICATION_20261110.md`.

---

## 7. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Found through the repository and open-access-supplement search recorded in `OMICS2PK_INDEPENDENT_REPOSITORY_SEARCH_20261110.md`. |
