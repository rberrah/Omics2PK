# Omics2PK — funding preparation dossier

**Version 1.0.5 — 2026-11-11.**
Consolidates every executed result, open question and sizing the project holds.
Nothing new is computed here; every figure is traceable to a named document and
a named artefact, and each is labelled with its level of proof.

**Written for:** reviewers and co-investigators assessing whether the missing
cohort is worth funding, and what it should look like.

---

## Summary in one page

**What the method does.** Given individual-level molecular measurements, it
selects a small mechanistically motivated predictor set inside training folds
and tests out-of-sample whether that set adds to what is already known, with a
permutation null and an interval that refits.

**What the public record now supports.** The pipeline runs end-to-end on real
individual-level multi-omic data under protocols frozen before the first fit.
A molecular measurement predicts function in the direction mechanism says it
should, consistently, across three cohorts and two layers — and **no public
cohort is large enough to establish its size.** One clinical cohort shows a
large exploratory gain from two metabolite ratios; a second, on a different
drug, shows none.

**What it does not support.** No multi-omic complementarity, in three attempts
on three pairs of layers. No clinical-PK proof: the single C-level contrast is
negative. No transfer between drugs, enzymes, tissues or cohorts.

**What is missing, and why no amount of searching will supply it.** Nine public
liver cohorts measured a per-donor molecular layer paired with a per-donor
outcome. **Three published the rows**, all three by accident of format rather
than policy, and the largest of them has fifteen donors. No public source
anywhere pairs two molecular layers with a drug exposure in the same
individuals at a size that can establish an effect.

**What to ask for.** A staged cohort: **≈120 donors** for the functional
question, which the literature sizes; and **≈400–1,000 patients** for the
clinical question, which only scenarios size. The honest version of the
request names which of those two numbers is measured and which is assumed.

---

## 1. Public results established, with their level of proof

Full detail and every interval: `docs/OMICS2PK_EVIDENCE_LEDGER_20261110.md`
(v1.2.0), whose figures are read back out of their artefacts by
`tests/test_evidence_ledger.py`, so the ledger cannot drift silently.

| Level | Count | The strongest row at that level |
|---|---|---|
| **R** reproduction | 3 | published results and artefact identities recomputed from their own sources |
| **V** independent validation | 1 positive, **1 failed** | the CYP transcript adds to age and sex for in-vitro CYP activity in **100 EPA donors**: 2 of 8 enzymes after Bonferroni, with one enzyme going the wrong way |
| **F** functional proof | 5 | one enzyme's absolute abundance cuts its substrate's hepatocyte clearance error from **×1.55 to ×1.18** (CYP3A4) and **×2.46 to ×1.55** (CYP2C19), Wegler n=15 |
| **C** clinical-PK | **1, negative** | a plasma CYP3A metabolite does **not** predict a dose-normalised venetoclax trough in 12 patients, against a baseline of no information: fold error **2.609 → 2.541**, paired effect **+0.0006** log₁₀, refitting interval [−0.111, +0.190] |
| **E** exploratory | 10 | two solanidine ratios cut an observed metoprolol trough's error by **38 %**, surviving Bonferroni over eleven contrasts, ARIANA n=47 |

### 1.1 The three statements a funding case can make, and their exact form

**(a) The pipeline handles real individual-level multi-omic data.** Executed on
EATRIS-Plus v4: 12 layers, 127,620 features, 125 individuals, three
prespecified configurations whose sizes and checksums match the commissioned
ones. Fold-safe preprocessing, dimension reduction and penalty selection
throughout.

**(b) A molecular measurement predicts function where mechanism says it should
— consistently in direction, and at no establishable size.** Three cohorts:

| Cohort | Layer | Result |
|---|---|---|
| Wegler n=15 hepatocytes | absolute protein | 2 of 5 substrates after Bonferroni |
| EPA n=100 | transcript | 2 of 8 enzymes after Bonferroni |
| **Lee / Jackson n=15 hepatocytes** | absolute protein | **0 of 6 contrasts**, every estimate improving, five of six clearing their null |

The third is an **attempted independent replication of the first and it did not
confirm it**. The defensible claim is therefore the weaker one: *the
relationship runs the expected way everywhere it has been looked at, and no
public cohort can size it.* That is a better argument for new data than the
stronger claim was.

**(c) Measurement beats representation where data is scarce.** In ARIANA, two
measured metabolite ratios cut an observed trough's error by 38 % while the
reimbursed CYP2D6 genotype adds **nothing** on top of them (+5.7 %, a
degradation; a ≥10 % PGx benefit is excluded by the interval).

### 1.2 Two methodological findings that constrain every number above

**A permutation null is positive whenever the compared arms differ in width.**
A wider arm shrinks harder under an inner-CV penalty and so scores better with
*no information at all*. Measured null medians run to +0.012 ΔR². Three of the
project's own bootstrap-significant contrasts were withdrawn on this ground.
**Consequence:** any incremental-value claim without a permutation null is
uninterpretable.

**A non-refitting paired bootstrap understates its interval by a measured
factor.** Computing both estimators on the same six contrasts: a refitting
out-of-bag interval is **1.78–2.09×** wider for a one-predictor arm and
**3.76–6.32×** wider for two. The project's own earlier positive findings were
declared on the narrow estimator and should be discounted accordingly.

### 1.3 Eleven claims this project has withdrawn on its own evidence

Including its strongest published comparison (a protein-versus-mRNA headline
withdrawn for comparing a figure with an interval against one without), an
invalid statistic (a squared Spearman read as a variance fraction), a source
verdict corrected after opening the file, a search conclusion withdrawn as too
narrow, and the downgrade of its own best functional result. The register is
§2 of the evidence ledger.

**This is offered as evidence about the project, not as an apology.** A ledger
that only grows is not a ledger.

---

## 2. The scientific questions still open

| # | Question | Why it is still open | Sized? |
|---|---|---|---|
| **A** | Does an enzyme's own abundance predict its activity well enough to be useful? | the direction is consistent across three cohorts; the largest has 100 donors and a transcript, the two with protein have 15 each | **yes, from literature** |
| **B** | Does any molecular layer add **beyond** clinical covariates and pharmacogenetics? | three attempts, three non-detections — and each was underpowered or on a surrogate endpoint | **scenarios only** |
| **C** | Does a molecular layer improve prediction of a **clinical** exposure? | one exploratory positive (n=47, one drug), one negative (n=12, another drug) | **scenarios only** |
| **D** | Can a **plasma** measurement substitute for a tissue one? | newly askable: in 15 donors the plasma metabolite ranked *ahead* of the hepatic protein on one endpoint and behind it on another | **not sized** |
| **E** | Which **marker class** works — substrate-derived or turnover-dependent? | the public record has one of each, with opposite results, and no cohort carrying both | **not sized** |

Questions **D** and **E** are new as of 2026-11-10 and neither was on the
project's list before the data search. Both are cheap to design for and
expensive to discover afterwards.

### 2.1 Why the three non-detections of complementarity do not settle question B

| Attempt | Why it is not evidence of absence |
|---|---|
| Axis A, EATRIS-Plus, 79 individuals × 127,620 features | the sizing shows a 20,000-feature block with a 10 % true gain has **3 % power at n = 800**. At 79 individuals this design could not have found it |
| Axis B, Wegler, two principal components of 4,037 proteins | two components of a standardised wide block; a mechanism carried by a few low-variance proteins is invisible in the dominant directions |
| Axis B2, Lee / Jackson, plasma metabolite + hepatic protein | n = 15, and the refitting interval spans [−0.33, +0.74] |

**None of the three had the power to detect the effect the funding case cares
about.** Quoting them as evidence that layers do not complement would be a
misuse of the project's own results.

---

## 3. The minimum data needed to resolve them

### 3.1 What exists and why it cannot be stretched further

Nine public liver cohorts measured a per-donor molecular layer paired with a
per-donor second layer or outcome. **Three published the per-donor rows**:
Wegler 2021 (15 donors), Lee 2023 (15), Lee 2024 (12 patients). Two are
paywalled with no PMC record at all. The remaining four publish medians,
ranges and correlation coefficients.

**The barrier is individual-level sharing, not measurement**, and that matters
because it is addressable:

| Instrument | What it could yield | Status |
|---|---|---|
| **a data-sharing request** to the four aggregate-only sources, above all the **171-donor** Bhatt UGT cohort with abundance, activity *and* genotype in the same livers | would answer question **A** outright, at no measurement cost | **named, not used.** Outside the standing authorisation-free constraint; deferred by the commissioning note of 2026-11-10 pending exhaustion of public channels — which this search has now done |
| **a grant** | the cohort no level of sharing can supply | the subject of §6 |

### 3.2 The minimum requirement per question

| Question | Minimum data | Could sharing supply it? |
|---|---|---|
| **A** | ≥120 donors with same-donor absolute enzyme abundance and a drug-specific activity in the **same preparation** | **yes** — the Bhatt cohort alone has 171 |
| **B** | two molecular layers **and** clinical covariates **and** genotype in the same individuals, with a pharmacological endpoint | **no.** No public source has this, shared or not |
| **C** | dose, sampling time, concentration, genotype, clinical covariates and ≥1 omic layer in the same patients | **no** |
| **D** | a plasma layer and a tissue layer in the same individuals, with a measured function | **partly** — Lee 2023 has it at n=15 |
| **E** | both marker classes measured in the same patients with an exposure endpoint | **no** |

**That table is the funding argument.** Questions A and D are bounded by data
that exists and is not shared. Questions **B, C and E — the ones the programme
is actually about — are bounded by data that has never been collected.**

---

## 4. Sample-size scenarios, with their assumptions and their uncertainty

Full detail: `docs/OMICS2PK_COHORT_SIZING_RESULTS_20261110.md`. Method: the
whole pipeline is simulated — selection inside folds, penalty by inner
cross-validation, out-of-sample scoring, multiplicity — because shrinkage and
selection change the distribution of the quantity actually reported in ways no
closed form captures.

> **Two levels, and the difference is load-bearing.** Question A's inputs are
> published effects of the same kind on the same sort of endpoint — **level 1**.
> Questions B and C run on **declared scenarios** — level 2 — because no public
> source gives this project an out-of-sample incremental effect for an omic
> layer on a pharmacokinetic endpoint. **No in-vitro effect is carried into
> either.**

### 4.1 Question A — level 1, from published correlations

Input: 28 published Spearman correlations between each CYP's abundance and its
own probe reaction in individual human preparations (*Pharmaceuticals*
2025;18:1789, CC BY). Used only to set a simulation's generating correlation,
through a declared Gaussian-copula mapping; never read as performance, and
never squared into a variance fraction.

| Rank correlation | Required n (80 % power) |
|---|---|
| 0.20 — the weakest published hepatocyte value's neighbourhood | **902** |
| 0.40 | 189 |
| **0.53** — published hepatocyte value for CYP3A4/midazolam | **95** |
| 0.70 | 45 |
| 0.90 | 14 |

**One assumption decides the answer.** The correlation moves the required size
by a factor of 64; tripling the analytical dispersion costs 11 extra donors and
20 % missingness costs 28. A reviewer should interrogate *which enzyme and which
preparation*, not assay precision.

**And ≈95 is a floor, not an estimate.** The confirmatory sizing computed from
the effect this project actually *measured* in the Lee / Jackson cohort gives
**406** donors for midazolam and **53** for ibrutinib M37 — the same fifteen
donors, two endpoints. In that cohort the in-sample correlation was **+0.73**,
above the published 0.53, and out-of-sample prediction was still off by a factor
of three.

> **Question A, stated for a reviewer: ≈100 donors if the published
> correlations transfer, up to ≈400 if the measured out-of-sample effect is the
> better guide. Ask for ≈120 and prespecify an interim analysis.**

### 4.2 Question B — level 2, declared scenarios

| Assumed reduction in mean absolute error | Required n |
|---|---|
| 5 % | **not reached below 1,200** |
| 10 % | **976** |
| 20 % | 529 |
| 30 % | 365 |

Reference scenario: baseline R² 0.30, 200 candidate features of which 10
informative, redundancy 0.30 with the baseline, analytical noise 0.30, 5 %
missing, 10 features selected per training fold. **Every one of those is an
assumption.** Two simplifications favour the layer — whole-column mean
imputation and a clean shared signal — so each figure is a **lower bound**.

**Dimensionality is the binding constraint, and this is the sharpest number in
the dossier:**

| Layer width (10 informative features, 10 % true gain) | Power at n = 800 |
|---|---|
| 200 features | **0.58** |
| 2,000 features | 0.20 |
| 20,000 features | **0.03** |

A transcriptome-scale block is **undetectable** at 800 individuals. This is the
quantitative case for mechanistic restriction of the candidate set: it is the
difference between 0.58 and 0.03, not a stylistic preference.

**A layer that is 85 % redundant with the baseline is unreachable at any size
inside the search range.** Recruitment cannot rescue a layer that repeats what
is already known.

### 4.3 Question C — level 2, and the three endpoints are not interchangeable

| Endpoint | 10 % reduction | 20 % reduction |
|---|---|---|
| single trough concentration | **1,051** | 538 |
| AUC from dense sampling | **939** | 501 |
| clearance from a population model | **929** | 501 |

**The three barely differ, which contradicts what the project assumed.**
Halving the measurement noise and raising the baseline R² from 0.20 to 0.50
moves the required size by about 12 %, because a cleaner endpoint reduces both
the error the layer could remove and the error it must beat, roughly in
proportion.

**Consequence for a protocol:** dense sampling must be justified by what the
endpoint *means* to a clinician, not by a cohort-size saving, because there is
almost none. **A single-concentration design at ~1,050 patients is not
meaningfully weaker than an AUC design at ~940**, and it is far cheaper per
patient.

**Observed anchor, not an assumption:** a dose-normalised venetoclax trough
carries a fold error of ≈2.5 with no covariate model at all (n=12), consistent
with the 0.20 baseline R² assumed above.

### 4.4 What the sizing cannot tell you

- Every question-B and question-C figure is a **scenario**. Quote the
  assumption with the number or do not quote the number.
- A required n from a simulation is a statement about the assumed generating
  process, not about biology.
- Questions **D** and **E** are **not sized**: the project has no effect
  estimate for either, and producing one would need the data being requested.

---

## 5. Modalities, compartments and endpoints — advantages and limits

### 5.1 Modality

| Modality | For | Against | Evidence here |
|---|---|---|---|
| **absolute proteomics** | closest to the mechanism; a CYP's abundance is the quantity in the clearance equation | expensive; needs tissue; targeted panels pre-commit to the enzymes | 2 of 5 substrates in Wegler; 0 of 6 in Lee / Jackson |
| **transcriptomics** | cheap, broad, widely deposited | one step further from function; 20,000 features is where power dies | 2 of 8 enzymes at n=100; **no public cohort pairs it with protein per donor** |
| **targeted metabolomics** | obtainable from **blood**; can be a direct phenotype readout | marker-specific; turnover-dependent markers are confounded by physiology | ARIANA **+38 %** on metoprolol; venetoclax **no gain** |
| **genotype** | cheap, reimbursed, stable | captures only coded variation | adds **nothing** on top of metabolite ratios (ARIANA) |

**The ordering is not fixed.** In Lee / Jackson the plasma metabolite alone
predicted hepatic activity slightly *better* than the hepatic protein alone on
the primary endpoint (fold error 2.796 against 2.976) and slightly worse on the
secondary (2.233 against 2.116). Neither ordering is established, and a blood
sample against a liver biopsy is the cheapest high-value contrast the project
has identified.

### 5.2 Compartment

Hepatocyte and microsomal preparations do not give the same answer, and the
public record is inconsistent about which is better:

| Source | Hepatocytes | Microsomes |
|---|---|---|
| Wegler n=15, out-of-sample | 2 of 5 substrates | **0 of 5** |
| *Pharmaceuticals* 2025, in-sample ρ for CYP3A4/midazolam | 0.54 | **0.82** |
| Zhang n=128, MPPGL against activity | strong in **tissue** | none in microsomes |

Wegler's out-of-sample result favours hepatocytes; an independent in-sample
correlation favours microsomes. **The discordance is unresolved and is recorded
rather than averaged away.** Part of it is certainly the difference between an
in-sample correlation and an out-of-sample contrast — the distinction the Axis
B2 run demonstrated at ρ = 0.73 and a factor-of-three error.

Intestinal mucosa shows essentially no abundance–activity relationship for
CYP3A4 (published ρ 0.13 midazolam, 0.23 testosterone, n=15). **A compartment
can abolish the relationship entirely**, so a protocol must declare its
compartment before it declares its hypothesis.

### 5.3 Endpoint

| Endpoint | For | Against |
|---|---|---|
| in-vitro CL<sub>int</sub> | mechanistically interpretable; recomputable | says nothing about a patient; two of three published definitions are **circular** against a CYP-abundance predictor |
| in-vitro activity rate | simple, widely measured | one substrate concentration; not V<sub>max</sub>/K<sub>m</sub>; not comparable across labs |
| single clinical concentration | cheap; what TDM actually uses | sampling-time and adherence error; lowest baseline R² |
| AUC | what exposure means | expensive; **buys almost no cohort-size saving** (§4.3) |
| model-derived clearance | least noisy | the baseline is partly the same model; the apparent gain depends on the structural model chosen |

**Two traps, both encountered and both expensive:** Wegler's equations 8 and 9
scale clearance by the very CYP abundance used as predictor — a circular
endpoint that would have produced a spuriously favourable result; and donor
identifiers in the two halves of that source are permuted for 6 of 15 donors,
so assuming identity mis-pairs them. **A new cohort should fix its endpoint
definition and its identifier mapping before the first measurement.**

---

## 6. A minimal cohort, then an ideal one

Staged deliberately: the exploratory phase must be able to **fail cheaply**,
and the validation phase must be designed before the exploratory phase reports.

### 6.1 Stage 1 — minimal, exploratory

| | |
|---|---|
| **n** | **120 donors** with liver tissue, or ≥150 if archival material allows |
| **layers** | absolute targeted proteomics (≤30 enzymes and transporters) **and** transcriptomics of the same panel, **same donors, same preparation** |
| **endpoint** | drug-specific activity for ≥5 probe substrates, **plus** a declared non-circular CL<sub>int</sub> definition fixed in advance |
| **also** | genotype; a plasma sample for **both** a substrate-derived and a turnover-dependent marker |
| **answers** | **A** outright; **D** and **E** for the first time |
| **sizing basis** | level 1, published correlations; ≈95 at ρ = 0.53 with a margin to ≈120 for the gap between correlation and measured effect (§4.1) |
| **cost control** | a 30-protein targeted panel, not a discovery proteome. §4.2 shows a 20,000-feature block is undetectable at this size, so a wide block would be money spent on power that does not exist |

**Prespecified failure criterion.** If the abundance arm does not beat the mean
out-of-sample for ≥2 of 5 substrates with a refitting interval excluding zero,
the programme does not proceed to stage 2. Writing that down now is what makes
stage 1 worth funding.

### 6.2 Stage 2 — the clinical cohort, independent validation

| | |
|---|---|
| **n** | **400 patients** minimum; **1,000** for a 10 % gain under the reference scenario (§4.2–4.3) |
| **design** | one drug with a known CYP route and a reimbursed genotype; dose, **recorded sampling time**, concentration |
| **layers** | the stage-1 panel, measured in **blood**, plus both plasma marker classes |
| **baseline** | clinical covariates **and** genotype — the comparator must be what clinical practice already has, which neither ARIANA nor the venetoclax cohort could fully provide |
| **endpoint** | a single timed concentration is **statistically adequate** (§4.3); add dense sampling only where exposure itself is the clinical question |
| **answers** | **B** and **C** |
| **split** | prespecified: exploratory half, then **locked** validation half. Three attempts at complementarity have now been reported by this project and none was confirmatory by design |

### 6.3 The protocol requirements this project has paid for

Each of these cost a withdrawn claim or a failed replication to learn:

1. **Select the layer's features on the residual after the baseline**, inside
   the training folds — not on the raw outcome. Selecting on the outcome
   **rewards redundancy**: selection accuracy rises to 1.000 exactly where the
   layer is least useful (§4 of the sizing results).
2. **Compute a permutation null for every reported contrast**, not only the
   significant ones. A wider arm beats a narrower one with no information.
3. **Report an interval that refits the model inside each draw.** The
   non-refitting version is 1.8–6.3× too narrow.
4. **Fix the endpoint definition in advance and check it is not circular**
   against its own predictor.
5. **Re-derive identifier mappings from the data**, never assume them.
6. **Declare the compartment and the preparation before the hypothesis.**
7. **Treat repeated visits as repeated measures** — the patient is the
   resampling unit, not the observation.

### 6.4 What a reviewer should press on

- **Question A can be answered without this grant** if the 171-donor Bhatt
  cohort were shared. That request should be made in parallel, and the dossier
  says so rather than hiding it.
- **The 38 % ARIANA gain is exploratory, n=47, one drug**, and a second drug
  showed no gain. Stage 2 is designed to test it, not to assume it.
- **Three non-detections of complementarity** are on the record. §2.1 shows
  none had the power to find the effect stage 2 is sized for — which is an
  argument for stage 2 and also the reason stage 2 must be sized honestly.
- **No in-vitro effect anywhere in this dossier is used to size a clinical
  study.**

---

## 7. Source documents

| Topic | Document |
|---|---|
| every claim by level of proof, and eleven withdrawals | `OMICS2PK_EVIDENCE_LEDGER_20261110.md` |
| the eight-repository search and the access ladder | `OMICS2PK_INDEPENDENT_REPOSITORY_SEARCH_20261110.md` |
| the publisher-first search it supersedes | `OMICS2PK_CANDIDATE_PAIRED_LIVER_SOURCE_QUALIFICATION_20261110.md` |
| the new three-layer cohort | `OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_20261110.md` |
| the failed replication and the interval measurement | `OMICS2PK_AXIS_B2_PLASMA_BIOMARKER_VALIDATION_RESULTS_20261110.md` |
| the first clinical-PK contrast | `OMICS2PK_AXIS_C2_VENETOCLAX_CLINICAL_PK_CONTRAST_20261110.md` |
| all sizing, and the selection-criterion finding | `OMICS2PK_COHORT_SIZING_RESULTS_20261110.md` |
| the original functional pilot | `OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_RESULTS_20261110.md` |
| the multi-omic feasibility run | `OMICS2PK_AXIS_A_EATRIS_PLUS_LAYER_COMPLEMENTARITY_RESULTS_20261110.md` |
| the clinical ablation | `OMICS2PK_ARIANA_STATISTICAL_VERSUS_MECHANISTIC_ABLATION_RESULTS_20261110.md` |
| every source and its status | `registry/PUBLIC_DATA_EVIDENCE_REGISTRY.json` v1.22.0 |

---

## 8. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-10 | First issue. Supersedes no earlier dossier. |
| 1.0.1 | 2026-11-10 | Registry citation moved to v1.20.0 after two further sources were qualified (both ABSTAIN). No figure in this dossier changes: neither new source is usable, so the data gap it describes is unchanged. Caught by `tests/test_funding_preparation_dossier.py`, which is why that test exists. |
| 1.0.2 | 2026-11-11 | Registry citation moved to 1.20.0, which added the `variant_functional_screen` endpoint class and two MaveDB paired-assay score sets. No sizing figure and no dossier conclusion changes: deep-mutational-scanning score sets measure engineered variants, not individuals, so they cannot enter a cohort-size calculation. Caught by `tests/test_funding_preparation_dossier.py`. |
| 1.0.3 | 2026-11-11 | Registry citation moved to 1.20.0, which registered two RCSB PDB structures as a structural knowledge layer. No sizing figure and no dossier conclusion changes: a crystal structure carries no individuals, so it cannot enter a cohort-size calculation. Caught by `tests/test_funding_preparation_dossier.py`. |
| 1.0.4 | 2026-11-11 | Registry citation moved to **1.21.0**, which registered a pre-registered independent test of the structural selection rule and its outcome — a refutation on BRAF, where the rule's loss-of-function premise does not hold. No sizing figure and no dossier conclusion changes: a deep-mutational-scanning score set measures engineered variants and a crystal structure carries no individuals, so neither can enter a cohort-size calculation. The dossier's **data gap** is unchanged, and the scientific case is now narrower and better evidenced: the mechanistic selection claim is scoped to enzymes whose variants disrupt catalysis at a ligand-binding site, which is a measured boundary rather than an assumed generality. Caught by `tests/test_funding_preparation_dossier.py`. |
| 1.0.5 | 2026-11-11 | Registry citation moved to **1.22.0**, which registered the ARAF generalisation test and RCSB 9AXM. No sizing figure and no dossier conclusion changes: a deep-mutational-scanning score set measures engineered variants and a crystal structure carries no individuals, so neither can enter a cohort-size calculation. **But one claim of this dossier is narrowed and one explanation withdrawn.** The scope of the mechanistic selection claim is no longer “enzymes whose variants cause loss of function”: ARAF is a loss-of-function library and the rule is refuted on it anyway, so the measured boundary is **metabolic enzyme versus signalling kinase**, tested on two kinases with opposite screen mechanisms. And the HSP90 explanation of the BRAF refutation is **withdrawn**: abundance under HSP90 inhibition is a detectably worse proxy than abundance under vehicle, so that refutation is unexplained again. Caught by `tests/test_funding_preparation_dossier.py`. |