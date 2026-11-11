# Q1 answered on the Human Liver Cohort: a contextual gate does not beat stacking

**Version 1.0.0 — 2026-11-11.** Evidence type: **exploratory, single cohort,
pre-declared comparators, four negative controls, one positive control.**
Artefact: `reports/data/nested_fusion_evaluation.json`.
Companion to `OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md` and
`OMICS2PK_FUSION_HEADROOM_RESULTS_20261111.md`.

**No frozen protocol or result is modified. No mechanistic prior is fitted, so
Q2 (`F` vs `X`) is untested.**

---

## 0. The headline, before the caveats

On 397 liver donors, ten activity endpoints, three context sets and two
temperature grids, a contextual gate learned without any mechanistic prior
**never beats constant-weight stacking**:

| Contrast | Cells | Reaching the declared 5 % margin |
|---|---|---|
| **`X` vs `K`** — Q1, primary | 30 | **0** (25 refuted, 5 underpowered) |
| `X` vs `K`, widened temperature grid | 10 | **0** (8 refuted, 2 underpowered) |
| `X` vs `S` | 30 | **0** (24 refuted, 6 underpowered) |
| `K` vs `S` | 10 | **0** (1 positive-but-below-margin, 1 undecided, 7 refuted, 1 underpowered) |
| `U` vs `S` | 10 | **0** (8 refuted, 2 underpowered) |
| context probes | 40 | **0** (28 refuted, 8 positive-below-margin, 4 undecided) |

**Nothing in this run reaches the margin it was measured against.** That
includes stacking, which is the comparator the gate had to beat, and it
includes the three endpoints the previous diagnostics had flagged as the
promising ones.

The most informative single number is not a loss. It is that the gate
**disagrees with simple selection about which arm should win on 34 % to 76 %
of donors** and is no better for it. The gate is not inert — it is making
hundreds of different choices per endpoint, and they are not better choices.

---

## 1. What the commissioning review asked for, and what it changed

Four defects were named in commit `5eda7c4`. All four were real. Their
consequences differ by two orders of magnitude, and saying so is the point of
this section.

| # | Defect | Status | Measured consequence |
|---|---|---|---|
| **1** | the probe chose its compared arm pair from the **whole cohort** before splitting | **fixed** — chosen inside each fold from training rows only | **negligible here**: mean Δρ **−0.0032**, mean \|Δρ\| **0.0047**, median 0.0000, verdict changed in 2 of 40 cells and in opposite directions |
| **2** | one **cached prediction layer**: stacking weights were fitted on rows whose base models had trained on the rows being scored | **fixed** — arms refit inside every outer fold, 210 fits per endpoint instead of 20 | **this is the one that mattered**: the stacking gain on CYP3A4/testosterone falls from **+167.5 [+29.6, +366]** to **+58.2 [−31.5, +151.5]**, and the verdict from a recovered gain to `REFUTED_GAIN_BELOW_MARGIN` |
| **3** | an adjusted interval excluding zero was named `CONTEXT_CARRIES_SIGNAL` at ρ as low as 0.17, below the declared 0.30 | **fixed** — the floor now applies to both sides, six rungs | **all four previously "robust" probes are reclassified**; none reaches the margin |
| **4** | the three endpoints were selected after seeing the probes on these same donors | **not repairable, only labelled** | they are marked `prespecified_but_selected_on_these_data`; all ten endpoints are reported |

### 1.1 Why defect 1 mattered so little, which the review could not have known

The four HLC arms have well-separated mean losses, so the two best arms are the
same in every fold on **7 of 10** endpoints (`arm_pair_stability` = 1.00). The
leak had almost nothing to leak. That is a property of this cohort, not a
defence of the old code: on a cohort with closely matched arms the same leak
would bite, which is why it is fixed and why the fix is tested.

### 1.2 A claim of mine that the measurement contradicted

I wrote in the first patch that the arm-selection leak "ran in the optimistic
direction". On synthetic null cohorts it runs the **other** way: the
globally-best arm is selected to have the lowest total loss, so the advantage
it defines is systematically negative, and over 12 null cohorts the
whole-cohort pair gives a **lower** ρ than the in-fold pair. The leak is a bias
with no reliable sign. The sentence is **struck, not deleted**, and the module
now says the direction has to be measured. On these HLC predictions the
isolated effect was +0.003 on average — neither optimistic nor pessimistic to
any useful degree.

### 1.3 The same error, found in my own new code

Defect 3 is a verdict ladder that awards a name implying the margin was met
when only zero was excluded. While writing the contrast ladder for this run I
reproduced it: `K` vs `S` on CYP1A2 came out `GAIN_PRESENT` on a **+0.62 %**
gain against a declared **5 %** margin. The contrast ladder now mirrors the
probe ladder rung for rung — `GAIN_REACHES_MARGIN`,
`GAIN_PRESENT_MARGIN_UNDECIDED`, `GAIN_PRESENT_BUT_BELOW_MARGIN` — and a test
fails if a bare `GAIN_PRESENT` is ever reintroduced.

---

## 2. The design, and why a null from it is readable

### 2.1 Nothing crosses the outer split

For each of 5 outer folds: inner 5-fold out-of-fold predictions are produced
on the **outer-training rows only**, from arms refit on subsets of those rows;
every comparator is fitted on those predictions and nothing else; the arms are
then refit on all the outer-training rows to predict the outer-test rows.
`fit_predict` is called 30 times per endpoint per configuration, 210 in total.

A test asserts that no fit ever receives a row it will predict. A second test
adds an arm that returns the exact truth for any row it was trained on: under
a correct design that arm is never trained on a row it predicts, so it is
useless everywhere and attracts under 1 % of the constant weight. Under the
previous cached-layer design it would have looked perfect to the weight-fitter.

### 2.2 The positive control, without which a null means nothing

`nested_fusion_evaluation.self_test` builds a cohort where two arms bracket the
truth on opposite halves and the context says which half. Constant weights
cannot help there by construction. The harness recovers a **91.7 % loss
reduction** for `X` over `K`, in all three context modes, and **loses it
entirely** when the context is permuted.

So a null `X` vs `K` on HLC is a statement about HLC, not about this
implementation. That is the only reason the rest of this document can be read
as evidence.

### 2.3 The gate was given every reasonable chance

| Given | Detail |
|---|---|
| three context sets | clinical (age, sex); the arms' own predictions, spread, mean and sd; both |
| a non-linear option | the probe runs a quadratic expansion as well as a linear one |
| a tuned sharpness | temperature chosen on a held-out split of the training rows, from a scale-free grid |
| a widened grid | the frozen grid chose its own edges in 70 of 150 folds, so the run was repeated with a grid two octaves wider on each side |
| room to degenerate | the softmax reduces to a constant mixture when the ridge coefficients vanish, and to hard per-individual selection as the temperature falls |

Under the widened grid the chosen temperatures scatter instead of pinning to
the edges, and the null is unchanged. **The null is not a grid artefact.**

One honest consequence of the grid being free: CYP2E1, the negative control,
produces the run's largest point estimate (+6.38 %) and its widest interval
(**[−30.6, +383.6]**), because on the two folds where the gate did sharpen,
nothing it could sharpen towards was informative. A nominally positive point
estimate on the negative control is a reminder that this design's power is
limited, not a hint of a signal.

---

## 3. Results

### 3.1 Q1: `X` vs `K`, context mode `both`, mean absolute error in source units

| Endpoint | Role | n | `S` | `U` | `K` | `X` | `X` vs `K` | Verdict |
|---|---|---|---|---|---|---|---|---|
| CYP1A2 | sensitivity | 397 | 394.8 | 396.2 | 392.3 | 400.2 | −2.00 % | refuted |
| CYP2A6 | sensitivity | 396 | 549.4 | 561.1 | 549.4 | 556.1 | −1.23 % | refuted |
| CYP2B6 | sensitivity | 397 | 207.3 | 209.0 | 208.2 | 209.4 | −0.58 % | refuted |
| CYP2C8 | sensitivity | 387 | 53.6 | 51.0 | 51.2 | 49.5 | +3.23 % | underpowered |
| **CYP2C9** | **prespecified** | 397 | 243.7 | 244.0 | 241.8 | 241.9 | −0.04 % | **refuted** |
| CYP2C19 | sensitivity | 397 | 42.0 | 42.9 | 42.4 | 43.0 | −1.33 % | refuted |
| **CYP2D6** | **prespecified** | 396 | 97.1 | 97.8 | 97.3 | 97.8 | −0.52 % | **refuted** |
| **CYP2E1** | **negative control** | 397 | 1214.3 | 1313.4 | 1316.6 | 1232.5 | +6.38 % | underpowered |
| **CYP3A4/midazolam** | **specificity control** | 396 | 705.5 | 701.0 | 683.7 | 686.8 | −0.45 % | refuted |
| **CYP3A4/testosterone** | **prespecified** | 396 | 3129.8 | 3132.1 | 3071.6 | 3080.0 | −0.27 % | **refuted** |

The three endpoints selected on the strength of the earlier probes are
**refuted on all three**, in all three context modes.

### 3.2 The negative controls all behave

| Control | What must happen | What happened |
|---|---|---|
| **permuted context** | `X` loses its advantage | 28 refuted, 2 underpowered of 30. No gain anywhere. |
| **decoy arm** (20 noise features) | the decoy attracts no weight | `K` gives it **0.000–0.050**. **`X` gives it 0.081–0.173** — see §4.1 |
| **CYP2E1**, all arms at the age/sex baseline | nothing to arbitrate | with clinical context the gate picked the **flattest temperature in the grid on all 5 folds**, i.e. it declined to gate — correct behaviour. With the arms' own predictions as context it did so on 3 of 5, and the two remaining folds are what makes this endpoint's interval the widest in the run, [−30.6, +383.6]. |
| **CYP3A4/midazolam** vs **/testosterone** | the two substrates are one endpoint | neither shows a gain, so they agree. The earlier contradiction disappears because the earlier signal does. |

### 3.3 There is headroom, and the gate recovers none of it

The per-individual convex ceiling on these nested predictions sits at
**64–79 %** of the selection loss, with bracketing rates 0.23–0.33. So 21–36 %
of the loss is theoretically reachable by *some* member of the fusion family.

`X` recovers **none** of it. A permissive ceiling licenses nothing, which the
headroom document already said; this run is what that statement looks like when
executed.

---

## 4. What is **not** established

### 4.1 The minimal gate has a defect of its own

A softmax over predicted losses cannot drive a weight to zero. With a
pure-noise decoy arm present, `X` still assigns it **8–17 %** of its weight
where `K` assigns **0–5 %**, and `X` becomes detectably worse than `K` on two
endpoints. **This is a limitation of the gate I built, not a measurement of
the hypothesis.** A gate with a sparsity mechanism or an explicit abstention
rule would not share it. It does not explain the primary null — the no-decoy
runs contain no useless arm, and the temperature was free to approach hard
selection — but it is a real reason not to read §3.1 as the last word on all
possible gates.

### 4.2 The null covers what was tested, and no more

| Not refuted | Why |
|---|---|
| gating on **context this cohort does not publish** | HLC publishes age and sex. That is a thin context set and the null is partly a statement about it. |
| **non-linear** gating beyond a quadratic expansion | the probe tests quadratic; the gate's loss model is linear in the expanded basis |
| gating on **other arms** | four arms from one model family on one cohort |
| **Q2**, the mechanistic prior | no prior was fitted. `F` vs `X` is untested, and the roadmap puts it after a positive Q1 that did not arrive. |
| contextual fusion **on clinical PK** | every endpoint here is in-vitro microsomal activity. Not a clearance, not an AUC, not a concentration. |
| contextual fusion at **larger n** | 5 of 30 primary cells are `NO_DETECTABLE_GAIN_UNDERPOWERED`, not refuted |

### 4.3 The provenance problem that cannot be fixed

CYP3A4/testosterone, CYP2D6 and CYP2C9 were named **after** the context probe
ran on these same 397 donors. This run governs what happens next; it does not
convert that choice into independent evidence. The fact that all three were
refuted makes the selection effect moot in this instance — but had one shown a
gain, that gain would not have been admissible as confirmation.

---

## 5. What this means for the programme

1. **Q1 is answered negatively on the only cohort where it could be asked.**
   Context-dependent weighting does not beat a fixed mixture here, and a fixed
   mixture barely beats picking one arm.
2. **Q2 should not be opened on this cohort.** The protocol makes `F` vs `X`
   conditional on context being exploitable at all. It is not, here. Fitting a
   mechanistic prior on weights that have nothing to weight would measure the
   prior's flexibility, not its knowledge — the π_m failure mode, transposed.
3. **The honest comparator recommendation stands and strengthens.** The
   commissioning note's own instruction — *if the concept shows no identifiable
   methodological advantage over an existing method, say so explicitly and
   prefer that method* — applies. On this cohort the method to prefer is
   **`S`, selecting one arm on training folds**. `U`, `K` and `X` all fail to
   beat it by anything worth having.
4. **The NetPK philosophy is untouched by this.** What is refuted is
   *contextual weighting of prediction outputs* on in-vitro activity with
   age and sex as context. Conditional **selection** of covariates and omic
   information — the broader ambition — is a different operation, and the
   commissioning note's own warning applies: the optimal algorithm may differ
   for each of them.

**The absence of gain is a valid scientific result**, which this project's
constraints state explicitly. This is one.

---

## 6. Reproduction

```
python scripts/run_nested_fusion_evaluation.py \
    --hlc-expression <dir>/curatedExpressionLiver.zip \
    --hlc-phenotype  <dir>/curatedPhenotypesynapse.zip \
    --out reports/data/nested_fusion_evaluation.json
```

Frozen before execution: seed 20261112, 5 outer × 5 inner folds, 2000 bootstrap
resamples over individuals, α 0.05, Bonferroni family 9 for the primary
contrasts and 40 for the probes, minimum interesting relative gain 5 %,
minimum interesting rank association 0.30, absolute-error loss. The arms,
designs, model and transform are those of `hlc_real_cyp_validation`, unchanged;
only the moment of fitting changed.

The artefact contains no individual-level values: per-individual losses and
predictions are held in memory and excluded from the written JSON.

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-11 | First issue. Q1 refuted on 25 of 30 cells and reaching the margin on none; four negative controls and one positive control; two leaks closed and their consequences measured separately. |
