# Automating the pharmacological reasoning that chooses which omic to measure

**Version 1.1.0 — 2026-11-11.** Evidence type: **architecture and protocol;
the first demonstrator is executed, and its diagnosis has been tested.** Answers the six items commissioned on
2026-11-11. Companion to
`OMICS2PK_MAVEDB_VARIANT_MECHANISM_RESULTS_20261111.md`.

**H1–H5 and every frozen result are untouched.** The Fusing branch is
suspended, not deleted: `nested_fusion_evaluation.py`, its protocol and its
negative results stay in the repository as the record that *fusing predictions*
does not pay on HLC. This document is about a different operation — *selecting
information* — and that negative result does not transfer.

---

## 0. The one-paragraph answer

The commissioned priority was to show that *pharmacological knowledge helps
choose which biological information is useful, rather than making prediction
models more complex without measurable gain.* **Executed on public data, it is
confirmed — once the knowledge layer is adequate.**

A variant within three residues of an annotated functional site is **6.74×
[4.82, 9.42]** more likely to keep protein abundance while losing activity, so
the mechanistic rule is right. With that **sequence** annotation it bought no
prediction, and the reason was measurable: the annotation reached only **48.8 %**
of the discordant cases on the better-annotated protein and **12.5 %** on the
other. The diagnosis — *the binding constraint is annotation coverage, not the
algorithm* — was a prediction, and it was tested by adding two crystal
structures. **It holds on both halves:** coverage rises to **95.9 %** and
**87.5 %**, enrichment to **55.35× [26.95, 251.41]** and **29.70× [10.78,
105.08]**, and the contrast that matters — *does knowledge-derived information
add to the best single measurement?* — goes from **+5.23 % / −2.98 %** (a sign
flip between the two validation splits) to **+15.83 % / +9.47 %, reaching the
declared margin under both.**

**So the commissioned priority is met**, on the primary gene, under the harder
validation, against a measurement that was already the best single predictor.
One part still fails: **conditional weighting** (`M5`) reaches the margin
nowhere, with either knowledge layer. Knowledge helps **select** information
and does not help **weight** it — the same split this project found for
prediction fusion.

---

## 1. Deliverable 1 — audit: what is reusable

### 1.1 Omics2PK, 67 modules

| What the new direction needs | Where it already exists | Verdict |
|---|---|---|
| random walk with restart on a protein graph | `mechanistic_graph_selection.py`: STRING v12 ≥700, restart 0.30, column-stochastic transition, **degree-stratified** permutation FDR | **reuse directly** |
| knowledge proposes / data disposes | `hybrid_mechanistic_selection.py`: shortlist 3k, retention by partial rank correlation **residualised on the baseline**, abstention allowed | **reuse as a comparator only** — it is detectably worse than a univariate filter in 21 of 72 cells |
| the arm to beat | univariate Spearman filter, in `scripts/run_h2a_hlc_selection_comparison.py` and `qualify_selected_probes.py` | **reuse; it is still the arm to beat** |
| LASSO, Elastic Net, XGBoost comparators | `run_h2a_hlc_selection_comparison.py`, `run_h1_ariana_incremental_omic_value.py`, `benchmark_public_*.py` | reuse |
| nested evaluation, paired intervals, the six-rung verdict ladder | `nested_fusion_evaluation.py`, `fusion_headroom_diagnostics.py` | **reused already**: the new module imports `ridge_predict` and `paired_gap` and shares the ladder and the 5 % margin |
| ceilings — is anything recoverable at all | `fusion_headroom_diagnostics.py` | reuse before claiming any gain |
| abstention | `hybrid_mechanistic_selection.py`, `hepatic_multiview_readiness.py`, `epa_eight_cyp_sensitivity.py` | reuse — abstention is a required output, §6 |
| admissibility gating of public sources | `public_data_evidence_registry.py` and the `*_source_gate.py` family | reuse unchanged |
| Reactome and CTD ingestion | `mechanistic_graph_selection.py` | reuse |

**Absent and needed:**

| Missing | Status |
|---|---|
| variant-level annotation and mechanism rules | **built: `variant_mechanism_omic_selection.py`** |
| **3D structural proximity to catalytic machinery** | **built: `structural_annotation.py`** — a minimal PDB reader, declared catalytic centres, and a numbering guard that fails unless a fact requiring correct numbering holds. This is what converted the rule from true-but-useless to useful |
| pharmacogenetic knowledge ingestion (CPIC, PharmVar) | **absent.** §3 shows why it cannot yet be a feature |
| CollecTRI / OmniPath regulatory ingestion | data is downloaded; **no module reads it** |
| multi-layer supervised selection (sPLS, block sPLS) | absent from Omics2PK; exists in MultiOmics as DIABLO |
| choosing the **selector itself** inside training folds | **absent.** This is the §6 success criterion and nothing implements it |

### 1.2 MultiOmics, 117 R files

`MOFA2` and `DIABLO` representation fusion with fold-safe tuning and selection
stability; DESeq2 references with null calibration; nutrimouse DIABLO/MOFA
validation. **Conceptually a comparator, not a dependency**: it is R, and the
commissioned requirement is that Omics2PK stay independent. Nothing is
imported. `DIABLO` remains the mandatory comparator for any multi-layer
representation work.

### 1.3 What the audit changes about the plan

The pieces for *statistical* selection and for *honest evaluation* are all
present and reusable. **The missing half is the knowledge side**, and
specifically its *coverage*, not its machinery. A new algorithm is not what is
short.

---

## 2. Deliverable 2 — an autonomous architecture

Five layers, no dependency on Thesis, NetPK, MetaPK or a GraphRAG, and no step
that requires a model that does not exist yet.

```
  L1  KNOWLEDGE          versioned, checksummed, offline-replayable
      ClinPGx/CPIC · PharmVar · UniProt · STRING · Reactome · CollecTRI · MaveDB
                            |
  L2  VARIANT ANNOTATION    what the variant is and where it sits
      position · substitution physicochemistry · distance to annotated
      machinery · motif membership · terminus proximity
                            |
  L3  MECHANISM HYPOTHESES  a RANKED LIST with prior weights, never a label
      nonsense · site-proximal · motif-internal · core-destabilising ·
      surface-tolerant · terminal-anchor
                            |
  L4  MEASUREMENT PROPOSAL  which omic TYPE should carry the information,
      and then which specific biomarkers, with an explicit ABSTAIN
                            |
  L5  SUPERVISED SELECTION  everything above is a PRIOR entering inside the
      and EVALUATION        training folds only; comparators, controls, ladder
```

### 2.1 The three rules the architecture is built to respect

1. **A ranked hypothesis is not a fact.** L3 returns candidates with declared
   prior weights and an argument in words. Nothing downstream may collapse
   them to a label. `mechanism_hypotheses` is tested to return a sorted list,
   and `M3` — mechanism features alone — is a comparator that is allowed to
   lose, which is how the rules' usefulness is measured rather than assumed.
2. **An interaction is not a cause, and a confidence score is not a
   probability of predictive usefulness.** A STRING combined score of 900 means
   the interaction is well supported, not that the partner is a useful
   predictor. This is the measured π_m failure mode of MetaPK, transposed:
   *network proximity measures how well studied a gene is.* The
   **degree-stratified** permutation null in `mechanistic_graph_selection.py` is
   therefore mandatory on every network-derived claim, and a UniProt
   binding-site annotation enters only as a hypothesis whose predictive value
   `M3` has to earn.
3. **Every supervised choice happens inside the training folds**, including the
   choice of method. A procedure that picks its own selector is only honest if
   that pick is itself nested.

### 2.2 Why the knowledge engine is files, not a service

L1 is a directory of versioned, checksummed downloads plus a reader module, not
a live API call inside an analysis. Reasons, all of them learned the hard way
in this project: an analysis that calls a network at run time is not
reproducible; a provider can change an annotation silently; and the admissible
public channel for a source is not always the obvious one. The new benchmark
therefore **embeds** the UniProt annotation it used, with the accession and the
date, and offers `--verify-uniprot` to re-fetch and fail loudly if the embedded
copy has gone stale.

---

## 3. Deliverable 3 — the public data, access verified today

Every row below was probed on **2026-11-11**, not inferred from a publication.

| Source | Channel | Status | What it gives |
|---|---|---|---|
| **MaveDB** | `api.mavedb.org`, CC0 | **ACCESSIBLE, used** | variant→function at scale; and, rarely, **two assays on one library** |
| **RCSB PDB** | `files.rcsb.org`, `search.rcsb.org` | **ACCESSIBLE, used** | experimental coordinates **with ligands**: 1R9O, 5LPG |
| **AlphaFold DB** | `/files/` path **403** → `/api/prediction/` **200** | accessible by the API path; **not used** | predicted coordinates, but **no ligands**, so it cannot define a catalytic centre |
| **UniProt** | `rest.uniprot.org` | **ACCESSIBLE, used** | sequence length, binding/active sites, motifs, domains |
| **CPIC** | `api.cpicpgx.org` | **ACCESSIBLE** | allele functional status: 94 CYP2C9 alleles, 23 NUDT15 |
| **PharmVar** | API **401, key requires an account** → public download file **200** | **ACCESSIBLE by the file channel only** | haplotype definitions: CYP2C9 v6.2.29, 345 variant rows; NUDT15, 57 |
| **STRING v12** | local, registered | accessible | protein interactions |
| **Reactome, CTD** | local, registered | accessible | pathways, chemical–gene |
| **CollecTRI / OmniPath** | local, downloaded | **accessible, unused** | regulatory relations |
| **Ensembl VEP REST** | timed out twice today | **UNVERIFIED** | would be needed to map PharmVar genomic coordinates to protein positions |
| HLC, EPA, Wegler, Lee | local, registered | accessible | enzyme activity, the §7 target |

### 3.1 PharmVar is another instance of the access ladder

The PharmVar **API** now answers `401 — API Key is invalid or missing … Keys
may be generated with a PharmVar Account`. An account registration is a
request for authorisation, so that channel is **inadmissible** under this
project's constraints. The **public download file for the same data returns
200.** As with BioStudies before it: *a 401 is a statement about the channel,
not about the data.*

### 3.2 The finding that reshapes the plan: the PGx layer cannot be a feature

PharmVar lists **345** variant rows for CYP2C9 and **57** for NUDT15. CPIC
lists **94** and **23** alleles. The functional screens measure **4 478** and
**2 844** variants. Even with a perfect coordinate-to-protein mapping — which
is not in hand, since VEP did not answer — **the pharmacogenetic knowledge
layer annotates a low single-digit percentage of the assayed variants.**

It therefore enters as an **evaluation stratum**, not a feature block: do the
predictions behave sensibly on the handful of clinically annotated alleles?
Building a feature that is 97 % missing and calling it pharmacological
knowledge would be the kind of move this project exists to avoid.

### 3.3 What separates the three questions

| Question | Data that can answer it | Status |
|---|---|---|
| is the **variant annotation** relevant? | MaveDB + UniProt | **executed**, §5 of the results |
| is the **omic type selection** right? | paired assays on one library (MaveDB), paired layers on one donor set (HLC, Wegler, Lee) | **executed on MaveDB**; HLC is §7 |
| does it improve **functional prediction**? | enzyme activity on the same donors | HLC, **not yet run under this architecture** |
| does it improve **clinical PK**? | omics + genotype + dose + sampling + concentration, same individuals | **`NOT_EVALUABLE`** on admissible public data; unchanged |

---

## 4. Deliverable 4 — comparators, criteria, controls

### 4.1 Comparators, in the order they must be beaten

| | Information used | Role |
|---|---|---|
| `M0` | nothing, training mean | anchors every relative gain |
| `M1` | the measurement closest to the gene, alone | **the arm to beat** |
| `M2` | variant physicochemistry only | is the variant itself enough? |
| `M3` | mechanism features only | **do the knowledge rules carry anything?** |
| `M4` | fixed union of all blocks | complexity without conditioning |
| `M5` | abundance used **conditionally** on mechanism (interaction form) | the new idea |
| `M6` | an interpretable per-variant **trust weight** on abundance | the new idea, in readable form |

`M5` collapses **exactly** to `M4` when its conditional coefficients vanish,
which is what makes `M5` vs `M4` a clean test of the conditional part and
nothing else. A test asserts it does not win where abundance is uniformly
informative — otherwise the extra columns would be buying flexibility rather
than knowledge.

### 4.2 Success criteria

| Criterion | How it is measured |
|---|---|
| out-of-sample performance | absolute error, position-grouped and region-grouped |
| **margin, not just significance** | the six-rung ladder: `GAIN_REACHES_MARGIN` needs the **lower** bound at 5 % |
| parsimony and stability | number and recurrence of selected variables across folds |
| biological relevance and traceability | every mechanism hypothesis carries a prior weight and an argument |
| **conditional** gain | always against the information already available, never against zero |
| transportability | a second gene, a harder split, and in §7 a second data type |
| **ability to abstain** | reported, and an abstention is a valid outcome |

### 4.3 Negative controls, declared before running

| Control | Must happen | Result |
|---|---|---|
| **permuted mechanism** | `M5`'s advantage disappears | holds: −0.06 %, −0.10 %, +0.03 % |
| **permuted abundance** | `M1` degrades | holds |
| **region-grouped split** | conclusions survive a harder split | **does not hold for `M4` vs `M1`** — see §5 of the results |
| **specificity**: predict abundance instead of activity | mechanism features help *less* | run and reported |
| **positive control** | the harness recovers a conditional rule that exists | holds: **+26.6 %**, lost under permutation |
| degree-stratified network null | mandatory on any future network claim | not yet applicable |

The positive control is not optional. Without it a null `M5` vs `M4` could
equally mean the data carry nothing or the implementation is broken — the
lesson the Fusing work paid for.

---

## 5. Deliverable 5 — the demonstrator, executed

Built, run and reported in
`OMICS2PK_MAVEDB_VARIANT_MECHANISM_RESULTS_20261111.md`. Headline: **the
mechanistic rule is confirmed at 6.74× and refuted as a predictive feature.**

---

## 6. Deliverable 6 — what is frozen, what is exploratory

| | Status |
|---|---|
| H1–H5 and every frozen contrast | **unchanged.** Nothing here touches them |
| the Fusing negative results and protocol v3.0.0 | **kept**, branch suspended |
| the MaveDB variant-mechanism benchmark | **exploratory.** First run, two genes, one data type, and the two genes are not independent evidence of each other |
| the 6.74× enrichment | **exploratory but robust**: interval excludes 1, and it is a descriptive statistic, not a trained model |
| the `M5`/`M6` nulls | **exploratory**, and dependent on sparse annotation — a null here is not a null for a better-annotated protein |
| anything about pharmacokinetics | **not addressed.** A functional-screen activity score is not a clearance, an AUC or a concentration |
| anything about multi-omic measurement in individuals | **not addressed.** These are engineered single variants in cell and yeast systems |

---

## 7. What to do next, in the order the measurements imply

1. ~~**Fix the knowledge layer, not the model.**~~ **DONE 2026-11-11 and it
   worked.** Structural proximity to the heme and substrate pocket was added
   from 1R9O and 5LPG; coverage went 48.8 %→95.9 % and 12.5 %→87.5 %,
   enrichment 6.74×→55.35× and 9.55×→29.70×, and `M4` vs `M1` on NUDT15 went
   from a sign flip to **+15.83 % / +9.47 %, reaching the margin under both
   splits**. CYP2C9's interval moved from `[0.00, 23.65]` to `[10.78, 105.08]`.
   **Still outstanding from this item:** cross-species conservation, and the
   CYP substrate-recognition-site definitions, which would be literature-derived
   and must be labelled as such rather than database-derived.
1a. **The structural layer's own limits now bound it.** It needs a crystallised
   ligand complex: a protein without one gets nothing. 327 CYP2C9 and 120
   NUDT15 assayed variants sit at **unresolved** positions — for a microsomal
   CYP that is the truncated membrane anchor — and are flagged, never silently
   imputed, because "absent from the structure" and "far from the machinery"
   are opposite claims.
2. **Map PharmVar to protein positions** so the PGx stratum can be used at all.
   Needs a coordinate→protein route; VEP did not answer today, so this is a
   data-access task before it is an analysis task.
3. **Then the HLC enzyme-activity benchmark**, with the same ladder, the same
   controls and donor-grouped folds: RNA only, protein only, functional proxies,
   genotype where genuinely available, and a supervised combination. **Do not
   attribute genotypes to HLC donors that are not actually published**, and do
   not join individuals across cohorts.
4. **Implement selector choice inside the folds** — the §6 criterion that
   nothing currently implements, and the only part of the brief that genuinely
   needs new machinery.
5. **Only then** consider block sPLS or graph-guided supervised selection, with
   DIABLO as the comparator. **No GNN**, per the brief and because nothing
   measured so far would justify one.

Boltz and molecular dynamics are not needed for any of the above, and item 1 is
the cheapest large gain available.

---

## 8. Change log

| Version | Date | Change |
|---|---|---|
| **1.1.0** | **2026-11-11** | §0 rewritten: the version 1.0.0 diagnosis was a prediction, it was tested with two crystal structures, and it holds on both halves — coverage 48.8 %→95.9 % and 12.5 %→87.5 %, enrichment to 55.35× and 29.70×, and `M4` vs `M1` from a sign flip to +15.83 %/+9.47 % reaching the margin under both splits. `structural_annotation.py` added to the audit. RCSB recorded accessible and used; AlphaFold recorded accessible by its API path but unusable here because its models carry no ligands. Roadmap item 1 marked done with its residue — conservation and SRS definitions — and item 1a added for the structural layer's own limits. Conditional weighting (`M5`) still reaches the margin nowhere. |
| 1.0.0 | 2026-11-11 | First issue. Answers the six commissioned items; audit of 67 Omics2PK modules and 117 MultiOmics R files; a five-layer architecture with no Thesis dependency; nine data sources probed today, of which PharmVar is admissible only by its file channel and VEP is unverified; seven comparators, seven success criteria, six negative controls including a mandatory positive control; the demonstrator executed. |
