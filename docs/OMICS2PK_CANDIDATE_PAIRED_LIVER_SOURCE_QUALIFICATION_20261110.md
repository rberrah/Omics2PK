# Seven candidate paired liver sources — none qualifies, and the reason is not what we assumed

**Version 1.0.1 — 2026-11-10.** Evidence type: **source qualification**.
Script: `scripts/audit_candidate_paired_liver_sources.py` (fits nothing).
Artefact: `reports/data/candidate_paired_liver_sources.json`.

Three candidates were put forward in the commissioning note of 2026-11-10;
four more were found while checking them. One rule, declared before use, was
applied to all seven.

---

## 0. The finding, which is not the one that was expected

**Of seven candidates, none qualifies as a dataset.** But the reasons separate
into two very different kinds, and the split is the result:

| Reason | Candidates |
|---|---|
| `ABSTAIN_ACCESS` — paywalled, publisher blocks automated clients | 2 |
| `ABSTAIN_NO_PER_DONOR_VALUES` — no per-donor layer paired with a per-donor outcome | **5** |

**Six of the seven measured both a molecular layer and either a second layer
or an enzyme activity in the *same* donors.** The measurements exist. What is
missing is the individual-level sharing: four open-access, CC BY articles with
supplements publish medians, ranges, group means and correlation coefficients,
and no per-donor row.

The seventh is the opposite case and must be separated from the others.
**Al-Majdoub publishes per-donor *activity* and measured abundance on pooled
commercial lots**, identified by lot number — so in that source the
donor-level pairing was never created at all, rather than created and
withheld. Its donor count (four, three after the article's own exclusion) is a
second and independent disqualification.

That is a materially different obstacle from "this data has never been
collected", and it changes what the funding case should argue. It is also the
reason Wegler was usable at all: its Table S3 happens to be a per-donor table,
which among these eight studies makes it the exception rather than the norm.

---

## 1. The rule, declared before it was applied

A source becomes a **dataset** only if **all four** hold:

1. retrievable **without** subscription, authorisation request or
   collaboration;
2. at least one measured layer **and** the outcome published or deposited **at
   the level of the individual donor** — a per-donor outcome paired with a
   **pooled** layer does not satisfy this, because no individual pairing
   exists to analyse;
3. at least **10 donors** carrying the required layers;
4. an admissible outcome kind — `drug_specific_clint`,
   `drug_specific_metabolic_activity`, `clinical_pk` — or, for Axis A, a
   declared non-pharmacological target.

**Why the floor is 10.** Leave-one-donor-out with an inner leave-one-out for
the penalty leaves n − 2 donors in the inner training set. Below ten donors
that is under eight, and neither the penalty choice nor the interval is
interpretable. The floor is a property of the design, not of the biology, and
it is stated so that a small cohort is refused by a rule rather than by taste.

A source failing **only** on (2) may still be a `published_comparator`. A
source failing on (1) is `citation_only`, because nothing about its contents
can be asserted at all.

**Two kinds of fact, handled differently.** Access and licence are re-queried
from Europe PMC **at every run**, and the audit **fails** if the fields moved;
an access verdict that is never re-checked is one that has quietly expired.
Whether per-donor values exist cannot be read from an index, so it is recorded
as evidence naming the file retrieved, its SHA-256 and its contents.

---

## 2. The seven, in the order that matters

### 2.1 Paywalled: the two that would have closed the biggest gap

**Drozdzik 2018** — *Clin Pharmacol Ther* 104(3):515-524, PMID 29205295.
Nine organ donors, paired liver and duodenum/jejunum/ileum/colon, **mRNA
(qRT-PCR) and protein (LC-MS/MS) in the same donors**.

This is the most conceptually valuable of the seven. Both layers in the same
livers is exactly what would turn this project's RNA-versus-protein comparison
from cross-cohort into within-cohort — the requirement named in the Axis B
results document. **It fails on access alone.** Europe PMC records
`isOpenAccess N`, `inPMC N`, `hasSuppl N`; the publisher answers **HTTP 403**
to every automated client, both for the article page and for the
`downloadSupplement` endpoint, behind a Cloudflare interstitial.

A web summary asserted that Tables S1 and S2 carry the per-donor protein and
mRNA values for the nine donors. **That assertion could not be checked against
the file, so it is not recorded as a fact.** If it is true, the dataset this
project needs for one of its three requirements exists and is one subscription
away; that is worth saying plainly, and it is not something this audit can
establish.

**Wenzel 2024** — *Clin Pharmacol Ther* 115(2):221-230, PMID 37739780.
Eleven donors, jejunum and liver, mRNA and protein in the same subjects, no
activity. Same group, same journal, same barrier: `isOpenAccess N`,
`hasSuppl N`, and the licence as indexed is **CC BY-NC-ND**, which would
forbid redistributing a derivative even if the file were reachable.

### 2.2 Pooled abundance, and far too few donors: the paper closest to our own question

**Al-Majdoub 2026** — *Drug Metab Dispos* 54(4), PMID 41855594, PMC13197946,
CC BY, open. Protein abundance **and** k<sub>cat</sub> / specificity constant
in matched liver and intestinal microsomes, seven enzymes.

The supplement was retrieved and opened (`mmc1.docx`, 934,810 bytes, SHA-256
prefix `80b321c6f197214c9c741464`), and it splits in a way the abstract does
not reveal:

- **Supplemental Tables 8 and 9 give intrinsic clearance and V<sub>max</sub>
  "from matched individual tissue donors"** — the *outcome* is per donor.
- **Supplemental Tables 11, 12, 13 and 14 give abundance for "pooled donor"
  microsomes identified by lot number**, for liver, intestine and kidney — the
  *predictor* is pooled.

**So the abundance-to-activity pairing per donor does not exist in this source,
independent of its size.** That is a different failure from the other four: not
a per-donor measurement summarised away, but a measurement never made per
donor. With four hundred donors it would still not support the contrast.

Separately, and sufficient on its own: **four paired donors, and the article's
own main comparison excludes one for low intestinal activity, leaving three.**
The commissioning note described this as "only three donors for some individual
comparisons"; the correct reading is that three to four *is* the whole cohort,
not a subset of a larger one.

Scientifically it remains the closest published work to this project's question
— it measures activity **per unit of enzyme** and concludes that scaling
activity through abundance needs care, which is the same concern Wegler's
microsome-versus-hepatocyte discordance raised. It is a citation for the
argument and cannot be a dataset.

### 2.3 Open, read, and aggregate-only: the four where the pairing exists but is not shared

**Zhang 2015** — *Sci Rep* 5:17671, PMC4669488, CC BY. **128 livers**, MPPGL
and ten CYP activities measured in **both** microsomes **and** liver tissue.

By far the largest candidate. The supplement was retrieved from the Springer
CDN (`41598_2015_BFsrep17671_MOESM1_ESM.pdf`, 1,004,095 bytes, SHA-256
`23c1ceba…bf7b`) and read in full: **41 lines of text**, holding one table of
aggregates — ranges, means, SD and 95 % prediction intervals for five scaling
parameters on the n = 78 prediction subset — and three figure captions. **No
per-donor value for MPPGL or for any CYP activity.** Recovering the individual
points would mean digitising scatter figures, which this project forbids.

**Its published direction concords with this project's own Axis B finding**,
and that is worth recording carefully. Zhang reports that MPPGL correlates
strongly with CYP activity measured in liver **tissue** (r ≥ 0.6, P < 0.001 for
most CYPs, moderate for CYP2A6, CYP3A4/5 and CYP2C8) and shows **no
significant correlation** with activity measured in **microsomes**. The Axis B
pilot found, in different donors with a different measurement, that the
protein-to-clearance relationship holds in hepatocytes and not in microsomes.
Two cohorts pointing the same way on the preparation is **concordance in
direction**; it is **not a replication**, because no per-donor value is
available to recompute anything, and nothing was fitted here.

**Drozdzik 2023** — *Int J Mol Sci* 24(5):4543, PMC10002520, CC BY.
Fifty-eight hepatitis C livers by Child-Pugh class, mRNA and protein. On paper
the candidate that would have closed the gap. The supplement was retrieved and
opened: the Europe PMC archive carries `ijms-24-04543-s001.zip` (133,511 bytes,
SHA-256 prefix `537c94b7d60348d9ca689dce`) containing exactly one file,
`ijms-2146685-supplementary.pdf` (141,723 bytes). Supplemental Table S1 reports
**median, min, max, mean, SD and CV % per enzyme per group** — group-level
aggregates only, no per-donor row for either layer.

**Drozdzik 2021** — *Pharmaceutics* 13(9):1334, PMC8471929, CC BY. About 97
livers across five pathologies (controls 20, HCV 21, ALD 20, AIH 20, PBC 10,
PSC 6 — the paper never states a total in one sentence). Read from the full
text: its **Data Availability Statement is "Not applicable"**, and the
supplementary tables it names are subject characteristics and **between-group
differences** in abundance and mRNA. It publishes per-enzyme
mRNA-to-protein correlations in liver at that scale, which is comparator
content for the layer question — the same role Grangeon plays for
protein-to-activity.

**Song 2018** — *Data Brief* 18:831-834, PMC5996730, CC BY. Table 1 is headed
"Individual human liver microsomes (n = 16)" and reports only the maximum, the
minimum, the max/min ratio and group means by age and sex. The supplementary
files are Mascot-processed MS identifications in `.doc` form, one per organ —
protein identifications, not per-donor quantification. Sixteen individual
donors with both proteomics and activity would have been a second Wegler; the
donor-level pairing does not survive the summarising.

---

## 3. What this changes, and what it does not

**Changes.** The project's data gap was previously stated as "no public cohort
has paired multi-omic and PK measurements". That is true but too coarse. The
sharper version, now evidenced:

> **The measurements exist in several cohorts. The individual-level values are
> almost never published.** Four open-access, CC BY articles covering 16, 58,
> 97 and 128 donors measured both layers, or a layer and an activity, in the
> same donors and published summaries. Two more, with the exact two layers the
> project needs in the same donors, are behind a subscription. A seventh
> measured the activity per donor and the abundance on pooled lots, so its
> pairing never existed.

For a funding case that is a stronger and more precise argument. It is also an
argument that can be acted on in a second way the project has not tried: five
of these seven have reachable corresponding authors, and a data-sharing
request is a different instrument from a grant. **This audit does not make
such a request** — the standing constraint is to use only data accessible
without a request or collaboration — but it names the option and the five
sources it would apply to, so the decision is the commissioner's.

**Does not change.** No conclusion from any earlier work moves. Nothing was
fitted here. Wegler remains the only qualified Axis B dataset, the readiness
gate's verdict of `fit_allowed` on 0 of 10 registered sources stands, and the
Axis B pilot's results are untouched.

**Does not claim.** Not that these studies are flawed: publishing group-level
summaries is a normal and legitimate choice. Not that no further source
exists: this is seven candidates examined, not a census. Not that the two
paywalled supplements lack per-donor values — their contents could not be
checked, and those two verdicts are about access and nothing else.

---

## 4. What should happen next, in the commissioning note's own order

The note's order is right and this audit does not change it, only its inputs.

1. **Step 1 is closed for these seven**, with evidence, and no new comparison
   became available. Step 2 — "run only the genuinely identifiable new
   comparisons" — therefore has **nothing to run**: there is no new per-donor
   dataset to run anything on. The single qualitative concordance with Zhang
   is recorded in §2.3 and is deliberately **not** turned into a computation,
   because a concordance of published directions is not a contrast and
   dressing it as one would be the error this project exists to avoid.
2. **Step 3, the evidence table**, is the next thing with content: it needs no
   new data and it is the deliverable a reviewer would read first.
3. **Step 4, the sizing**, should follow the evidence table and must be built
   per question rather than from the Axis B pilot's single number, which sizes
   one enzyme-abundance contrast on an in-vitro endpoint and nothing else.

---

## 5. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.1 | 2026-11-10 | **One verdict corrected by opening the file rather than trusting the abstract.** Al-Majdoub was first recorded as `ABSTAIN_TOO_FEW_DONORS`, on the assumption that its abundance and activity were both per donor. The supplement shows otherwise: Supplemental Tables 8 and 9 give per-donor clearance and V<sub>max</sub>, but Supplemental Tables 11 to 14 give abundance for **pooled donor** microsomes identified by lot number. The reason is therefore `ABSTAIN_NO_PER_DONOR_VALUES` — the pairing was never created — with the donor count recorded as a second, independent disqualification. The audit's own consistency check did not catch this, because the recorded field and the recorded reason were derived from the same wrong belief; what caught it was reading the file. The rule's second condition now states explicitly that a per-donor outcome paired with a pooled layer does not satisfy it, and the field is renamed `per_donor_layer_and_outcome` so the same confusion cannot recur. Counts move from 2 / 4 / 1 to **2 access and 5 no-per-donor-pairing**; no other verdict changed. |
| 1.0.0 | 2026-11-10 | First write-up. Qualifies seven candidates — the commissioning note's three plus four found while checking them — against a rule declared before use, with access and licence re-verified live against Europe PMC and per-donor availability recorded as checksummed file evidence. All seven abstain: two on access, four on the absence of per-donor values, one on donor count. Records three corrections to the note's framing: Drozdzik 2018 is 2018 rather than 2017 and is paywalled rather than merely "to be verified"; Al-Majdoub's cohort is four donors and three after exclusion rather than three "for some comparisons"; and Zhang 2015 publishes no individual values at all, so its 128 livers cannot be used. Records the concordance in direction between Zhang's microsome-versus-tissue finding and the Axis B pilot's, explicitly as concordance and not replication. |
