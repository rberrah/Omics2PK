# H6 feasibility — which omic layer deserves to be measured?

**Date:** 2026-11-09
**Status: feasibility study. No code was written and no model was fitted.**
Every number below is a read of a registered source or of the repository.

**Verdict in one line: H6 as formulated is not executable on the data this
project has registered. Two of its three questions substantially duplicate H3,
H4 and H5, and the third cannot be posed because no registered cohort carries
three omic layers in the same individuals.** One narrow pilot survives, and it
is reported with the reason it still cannot answer the headline question.

---

## 0. One claim in the commissioning note that does not hold

The note states that *"la demande précédente concernant H6 n'a pas encore été
traitée"*.

Checked: in the whole session transcript, `H6` appears as a hypothesis label
**only in that note itself**. Every other occurrence is a gene symbol (`ADH6`,
`ALDH6A1`, `MSH6`, `H6PD`, …). The protocol defines §4 H1, §5 H2, §6 H3, §7 H4
and §8 H5, and **no H6**. The repository contains no H6 document, section or
test.

So nothing was left undone. This changes nothing operationally — H6 is being
asked for now and is answered below — but the record should be accurate about
whether a request was dropped.

---

## 1. What is actually paired, measured from the registry

The three H6 questions all presuppose that two or more omic layers exist **in
the same individuals**. Here is every registered dataset with an omic layer,
with what it actually carries.

| Dataset | Access | Individual data | n analysed | Layer A | Layer B | Endpoint class |
|---|---|---|---|---|---|---|
| `ariana_metoprolol_47` | open | yes | 47 | targeted metabolomics (3 analytes) | CYP2D6 genotype + copy number | **`pk_clinical`** |
| `jarmusch_midazolam` | open | yes | 14 with PGx | untargeted metabolomics (GNPS) | genotype | **`pk_clinical`** |
| `muller_potato_meal_20` | open | yes | 20 | targeted metabolomics (3 analytes) | CYP2D6 genotype | `biomarker_only` |
| `hlc_syn4499` | registered | yes | 466 expression / 195 genotyped | transcriptomics (40,638 probes) | genotype, **2 usable rsIDs** | `functional_in_vitro` |
| `wegler_15` | open | — | 15 | proteomics | **none** | `functional_in_vitro` |
| `latham_tki_donors` | open | — | — | protein + activity | **none** | `functional_in_vitro` |
| `rotroff_metformin_33` | open | — | 33 | metabolomics | **none** | `pk_clinical` |
| `epa_liver_100` | open | — | — | — | **none** | `functional_in_vitro` |
| `nci60` | open | — | — | — | **none** | `pd_cell_line` |

**The structural fact that decides H6:**

> **No registered cohort pairs two *measured* omic layers in the same
> individuals.** The maximum is two, and one of the two is always **genotype**.
> Transcriptome + proteome: nowhere. Transcriptome + metabolome: nowhere.
> Proteome + metabolome: nowhere.

Proteomics exists exactly once, in 15 donors, with no second layer and against
an in vitro endpoint. That is the whole of the project's protein evidence.

---

## 2. Question by question

### Q1 — which layer carries the most information for a given endpoint?

**`NOT_EVALUABLE`.** A within-subject ranking of layers requires ≥ 2 measured
layers in the same people, which no cohort has. A between-cohort ranking would
compare, for instance, HLC transcriptomics (n ≈ 397, in vitro, liver
microsomes, microarray) against ARIANA metabolomics (n = 47, clinical trough,
plasma, LC-MS) — confounding the **layer** with the drug, the endpoint class,
the matrix, the assay, the cohort and the sample size simultaneously. No
adjustment available here separates those.

Reporting such a ranking would be the between-cohort fusion that
`src/omics2pk/hepatic_multiview_readiness.py` was written to refuse. Its
docstring names the exact error: *"to prevent invalid cross-study fusion such
as combining EPA RNA donors with unrelated TKI/Wegler protein donors."*

### Q2 — which combinations are genuinely complementary?

**Answerable for exactly one pair, and already answered — negatively.**

The only layer pair with same-subject data and a clinical endpoint is
**genotype + targeted metabolomics**, in ARIANA (n = 47). That is precisely
what the M4 / M5 / M5c arms measured, and the answer is in
`docs/OMICS2PK_H1_ARIANA_MECHANISTIC_HYBRID_RESULTS_20261109.md` §3:

| Arm | MAE | layers used |
|---|---|---|
| M4 | 39.44 | genotype only |
| M5 | 32.83 | genotype + metabolomics |
| **M5c** | **32.47** | **metabolomics only** |

M5c, using **no genotype at all**, matches M5. In the one cohort where the
question is posable the two layers are **redundant, not complementary**, and
the metabolomic layer is the one carrying the signal.

The second instance of this pair, `muller_potato_meal_20` (n = 20), has
endpoint class `biomarker_only` — no pharmacokinetics — so it cannot extend the
finding, and its registry entry records that the cohort's only poor
metaboliser is below the limit of quantification at all nine timepoints for
both metabolites.

### Q3 — can the most parsimonious multimodal panel be chosen automatically, accounting for tissue, sampling time and analytical cost?

**Cannot be posed.** Selecting *among* modalities requires at least three in
the same subjects to have anything to select from; the maximum is two, one of
them genotype. With two layers the "panel choice" has four possible answers and
one of them is already known from Q2.

Separately, most of this question is **already protocol**, and the overlap is
not marginal:

| Component of Q3 | Already specified in |
|---|---|
| performance against panel size, `k80` / `k90`, selector comparison | **§6 H3** |
| **analytical cost of the panel** | **§6 H3**, reported explicitly |
| redundancy across transcript / protein / metabolite | **§7 H4**, a listed condition |
| robustness to removing one modality | **§7 H4**, a reported quantity |
| tissue mismatch, sampling-time mismatch | **§7 H4** conditions and **§8 H5** |

What is left of Q3 after removing what H3, H4 and H5 already cover is a
**cost-weighted objective function over modalities** — a genuinely new element,
and one with no data on which to calibrate a cost or validate a choice.

---

## 3. What already exists and must not be rebuilt

The commissioning note asks to reuse existing MultiOmics capability. Measured:
the capability is substantially present, and the note's status line
*"Méthodologie à concevoir"* is therefore too pessimistic about the code and too
optimistic about the data.

| Module | What it already does | Run under the protocol? |
|---|---|---|
| `hepatic_multiview_readiness.py` | fail-closed gate: are RNA, protein and a drug-specific outcome measured on the **same** donors from one authenticated source? Returns explicit codes `ABSTAIN_NO_MATCHED_RNA`, `ABSTAIN_NO_MATCHED_PROTEIN`, `ABSTAIN_NO_MATCHED_DRUG_OUTCOME`, `ABSTAIN_NO_SAME_DONOR_INTERSECTION`, `ABSTAIN_INCOMPLETE_CASE_MISMATCH`, `READY_SAME_DONOR_RNA_PROTEIN_DRUG_OUTCOME` | **no** |
| `multimodal_pk_ablation.py` | fold-safe mono- versus multi-omic ablations against observed pharmacokinetics, one biological individual per row per modality, transforms fitted inside training partitions only | **no** |
| `multiomics_pk_readiness.py`, `multimodal_exchange_pk.py`, `matched_multiomics_pk.py` | source-hash-verified exchange audit and the individual-matrix → PK arm plumbing | **no** |
| `hepatic_protein_activity_fusion.py`, `multicyp_rna_protein_activity_source_gate.py` | protein/activity source gates | **no** |

All have unit tests (`tests/test_hepatic_multiview_readiness.py`,
`tests/test_multimodal_pk_ablation.py`, and six more) and those tests exercise
the **logic on synthetic inputs**. None has been run against the registered
sources, and `reports/data/` contains no artefact from any of them.

So the honest gap is not the method. It is that **the gate has never been
executed, so the project has no recorded verdict** — only my reading of the
registry above.

---

## 4. The one pilot that survives, and why it still does not answer Q1

HLC is the only cohort with genotype **and** transcriptome **and** an endpoint
in the same donors. I measured the three-layer intersection rather than
assuming it:

| | rs1799853 (`CYP2C9*2`) | rs776746 (`CYP3A5*3`) |
|---|---|---|
| genotype donors called | 178 of 195 | 192 of 195 |
| …also having expression | 161 | 175 |
| …also having the matching activity | **155** (CYP2C9) | **168** (CYP3A4 midazolam) |
| genotype distribution | CC 137, CT 41, **TT 0** | CC 155, CT 35, TT 2 |

(Expression donors 466; genotyped 195; intersection **178** — the October
source audit's figure, reproduced.)

**A three-arm comparison on CYP2C9 activity is therefore executable: genotype
alone, transcript alone, both, n = 155, same donors, same folds.** It would
reuse the F-series machinery directly.

**And it cannot answer Q1.** The "PGx layer" here is **one SNP with two
observed levels and 41 carriers, no homozygous variants at all**. `CYP2C9*2` is
a modestly reduced-function allele. A null result would be a statement about
*that single SNP*, not about pharmacogenetics as a layer, and reporting it as
"transcriptomics beats PGx" would be a category error. The CYP3A5 arm is worse:
its genotype predicts **CYP3A5** activity, which this cohort does not measure —
only CYP3A4 — and the registry already records that the array's reference
allele orientation is not independently confirmed.

So the pilot is worth running as a **bounded, honestly labelled** experiment —
"one SNP versus one transcript on one in vitro endpoint" — and worth nothing as
a layer-ranking.

---

## 5. The data requirement, stated so it can be checked

For H6 to become answerable rather than rephrased, one of these must exist,
publicly downloadable without authorisation:

1. **≥ 2 measured omic layers** (not counting genotype) in the same
   individuals, with a pharmacokinetic endpoint — the minimum for Q1 and Q2 at
   a clinical endpoint;
2. **≥ 3 layers** in the same individuals — the minimum for Q3 to be a
   selection problem at all;
3. failing both, **≥ 2 measured layers in the same individuals against an in
   vitro endpoint**, which would at least make Q1 posable within-subject.

The project's registry currently satisfies **none** of the three. The nearest
miss is HLC plus Wegler — transcriptome in 466 donors, proteome in 15 donors,
**different people** — which is exactly the fusion the readiness gate refuses.

`estonian_biobank_pk_recall` remains excluded: its access tier is
`governed_application`, and the mission constraints forbid cohorts requiring
authorisation.

---

## 6. Recommendation

The commissioning note's own recommendation — stop developing new
transcriptomic architectures — is supported by the F-series result and I agree
with it. I would go one step further on H6 itself.

**Do not open H6 as a sixth protocol objective.** As formulated it is one new
question (Q1) that the data cannot answer, one question already answered
negatively in the only cohort that can pose it (Q2), and one question that is
mostly H3 + H4 + H5 with a cost term added (Q3). A new objective would create
the appearance of a new line of evidence without a new line of data.

**Do these three things instead, in order:**

1. **Execute the readiness gate and commit its verdict.** One run of
   `hepatic_multiview_readiness` across every registered source, producing the
   `ABSTAIN_*` codes as a dated artefact. This converts §1 of this document
   from my reading into the project's recorded, reproducible position, and it is
   the cheapest honest deliverable available. It needs no new data and no new
   method.
2. **Add the layer inventory of §1 to the registry as a derived view**, so that
   "which layers are paired in whom" is a queryable registered fact rather than
   something re-derived by reading nine entries.
3. **Treat the CYP2C9 pilot as optional and bounded.** If run, it is labelled
   "one SNP versus one transcript", capped `exploratory`, and prespecified
   before execution like the F-series. It is not a layer ranking and must not be
   written up as one.

What I would **not** do: build a cost-aware multimodal selector. There is no
dataset on which to calibrate a cost, no dataset with three layers on which to
exercise a selection, and H3 already reports panel cost for the single-layer
case.

**The direction the note draws is right and this study sharpens it.** The
question worth asking is indeed which omic information deserves to be measured
rather than which representation of an omic is best. But on the data registered
today that question has one answerable instance, ARIANA's genotype-versus-
metabolite pair, and its answer is already recorded: the cheaper layer won and
the other was redundant.
