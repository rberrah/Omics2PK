# Axis B results — the enzyme predicts its own clearance in hepatocytes, and two principal components of the rest of the proteome add nothing

**Version 1.1.0 — 2026-11-10.** Protocol:
`docs/OMICS2PK_AXIS_B_WEGLER_PROTEIN_CLINT_PILOT_PRESPECIFICATION_20261110.md`
(v1.0.1; v1.0.0 frozen before any model was fitted).
Artefact: `reports/data/axis_b_wegler_protein_clint_pilot.json`.
Status: **`exploratory`**. Evidence type: **functional proof attempt on
independent public data, in vitro**.

**Specific to in-vitro intrinsic clearance.** There is no patient, no dose and
no exposure in this source. Nothing here is a clinical-pharmacokinetic result.
Five substrates in 15 donors are **not** five independent cohorts, and the
microsomal and hepatocyte readings of the same donors are **not** two cohorts.

Reported in the three separated sections the protocol requires.

---

## 0. The three answers

| Question | Answer |
|---|---|
| **Q1. Does the target CYP's own abundance predict its substrate's intrinsic clearance?** | **In hepatocytes, yes, for 2 of 5 substrates after Bonferroni** (midazolam/CYP3A4, omeprazole/CYP2C19) and a third unadjusted (diclofenac/CYP2C9). **In microsomes, for none of the five.** |
| **Q2. Do *two prespecified principal components* of the rest of the microsomal proteome add?** | **No, in any of the ten cells.** Not one arm improved with an interval excluding zero; the only interval that excludes zero is a **degradation**. **This does not show the other 4,037 proteins are uninformative** — see §4. |
| **Q3. Do these relationships differ between the two systems?** | **Yes, paired on the same donors, for 3 of 5 substrates** — the enzyme's abundance is more informative in hepatocytes. **Not a causal compartment effect.** |

---

## 1. Predictive performance

### 1.1 Q1 — the direct measurement of the enzyme (`W1` versus `W0`)

Mean absolute error on `log10` CLint; a negative difference is an improvement
over predicting the training-fold mean.

| System | Substrate | CYP | n | Δ MAE (log₁₀) | unadjusted 95 % | Bonferroni 99.00 % |
|---|---|---|---|---|---|---|
| **HH** | Midazolam | CYP3A4 | 15 | **−0.1168** | **[−0.2024, −0.0455]** | **[−0.2276, −0.0303]** |
| **HH** | Omeprazole | CYP2C19 | 15 | **−0.2012** | **[−0.3283, −0.0679]** | **[−0.3646, −0.0204]** |
| HH | Diclofenac | CYP2C9 | 15 | −0.0393 | **[−0.0712, −0.0059]** | [−0.0800, +0.0059] |
| HH | Bupropion | CYP2B6 | 15 | +0.0009 | [−0.0097, +0.0121] | [−0.0131, +0.0160] |
| HH | Bufuralol | CYP2D6 | 15 | +0.0012 | [−0.0008, +0.0039] | [−0.0011, +0.0048] |
| HLM | Midazolam | CYP3A4 | 15 | −0.0518 | [−0.1146, +0.0054] | [−0.1352, +0.0226] |
| HLM | Omeprazole | CYP2C19 | 15 | +0.0084 | [−0.0359, +0.0785] | [−0.0427, +0.1045] |
| HLM | Diclofenac | CYP2C9 | 15 | +0.0004 | [−0.0020, +0.0034] | [−0.0026, +0.0042] |
| HLM | Bupropion | CYP2B6 | 15 | +0.0132 | [−0.0019, +0.0331] | [−0.0039, +0.0427] |
| HLM | Bufuralol | CYP2D6 | **13** | +0.0038 | [−0.0004, +0.0108] | [−0.0007, +0.0138] |

**Two substrates carry a claim after Bonferroni over five**, both in
hepatocytes. Read in interpretable units, the same result:

| Cell | `W0` fold error | `W1` fold error | `W0` MAE (mL/min/kg) | `W1` MAE |
|---|---|---|---|---|
| HH Midazolam | ×1.55 | **×1.18** | 63.2 | **24.8** |
| HH Omeprazole | ×2.46 | **×1.55** | 22.9 | **14.1** |
| HH Diclofenac | ×1.26 | ×1.15 | 100.7 | 62.9 |
| HLM Midazolam | ×1.40 | ×1.25 | 147.6 | 95.7 |
| HLM Omeprazole | ×2.21 | ×2.25 | 27.0 | 26.6 |

Knowing one protein's absolute concentration takes midazolam's hepatocyte
clearance from "wrong by a factor of 1.55" to "wrong by a factor of 1.18", and
omeprazole's from ×2.46 to ×1.55. That is the strongest single statement this
pilot supports.

### 1.2 Q2 — two components of the rest of the proteome (`W2` versus `W1`, k = 2) — the primary contrast

| System | Substrate | Δ MAE (log₁₀) | unadjusted 95 % | Bonferroni 99.00 % | permutation null (median, range) | rank p |
|---|---|---|---|---|---|---|
| HLM | Midazolam | +0.0161 | [−0.0002, +0.0357] | [−0.0047, +0.0463] | +0.0027 [−0.0249, +0.0258] | 0.900 |
| HLM | Omeprazole | +0.0128 | [−0.0002, +0.0280] | [−0.0020, +0.0331] | +0.0016 [−0.1056, +0.0459] | 0.896 |
| HLM | Diclofenac | +0.0077 | [−0.0009, +0.0203] | [−0.0019, +0.0242] | +0.0014 [−0.0420, +0.0237] | 0.876 |
| HLM | Bupropion | −0.0099 | [−0.0269, +0.0047] | [−0.0371, +0.0076] | −0.0101 [−0.0421, +0.0363] | 0.507 |
| HLM | Bufuralol | −0.0021 | [−0.0044, +0.0004] | [−0.0048, +0.0012] | +0.0008 [−0.0228, +0.0172] | 0.045 |
| HH | Midazolam | −0.0042 | [−0.0148, +0.0076] | [−0.0167, +0.0118] | +0.0027 [−0.0225, +0.0184] | 0.239 |
| HH | Omeprazole | −0.0129 | [−0.0379, +0.0117] | [−0.0436, +0.0184] | +0.0049 [−0.0444, +0.1014] | 0.119 |
| **HH** | **Diclofenac** | **+0.0196** | **[+0.0094, +0.0302]** | **[+0.0069, +0.0326]** | +0.0105 [−0.0195, +0.0300] | 0.930 |
| HH | Bupropion | +0.0076 | [−0.0024, +0.0174] | [−0.0045, +0.0201] | +0.0008 [−0.0273, +0.0113] | 0.945 |
| HH | Bufuralol | −0.0142 | [−0.0448, +0.0189] | [−0.0552, +0.0267] | +0.0014 [−0.0233, +0.0259] | 0.030 |

**No cell shows an improvement whose interval excludes zero.** The one interval
that does exclude zero, and survives Bonferroni, is **HH diclofenac at
+0.0196 — a degradation**, and its permutation rank p of 0.930 says the
degradation is larger than most of what permuted proteomic information
produces. Adding two components there actively hurts.

Two cells sit in the improving tail of their null — HH bufuralol
(rank p = 0.030) and HLM bufuralol (0.045) — but **their bootstrap intervals
span zero comfortably**, and the permutation p-values are themselves
**unadjusted across the ten cells**: at Bonferroni over ten, neither survives.
Bufuralol is also the one substrate where `W1` itself carries nothing, so what
the proteome would be adding to is a baseline with no signal.

The null's median is **positive in eight of the ten cells**, which is the
shrinkage artefact the protocol predicted: `W2` has three columns where `W1`
has one, and under an inner-CV penalty the wider arm shrinks harder and scores
slightly better with no information. Any study reporting a `log10` MAE
improvement of 0.002–0.005 between arms of different width, without this
control, is at risk of reporting the artefact.

### 1.3 Both declared sensitivities point the same way

**k = 1** (one component instead of two). Five cells now show a **small
degradation with an interval excluding zero** — HLM midazolam +0.0043
[+0.0006, +0.0080], HLM omeprazole +0.0124 [+0.0003, +0.0264], HLM diclofenac
+0.0011 [+0.0001, +0.0024], HH diclofenac +0.0040 [+0.0001, +0.0088], HH
bupropion +0.0055 [+0.0008, +0.0119]. So the conclusion does not depend on the
component count, and at k = 1 the proteome component reliably costs a little.

**Cytochrome-free block** (all 32 cytochrome and redox genes removed from the
microsomal block, 28 from the hepatocyte one). Results are essentially
unchanged — HH diclofenac +0.0196 against +0.0196, HLM midazolam +0.0158
against +0.0161. **`W2`'s absence of contribution is not about other
cytochromes being in the block.**

### 1.4 Q3 — the paired system comparison

A difference of differences over donors: the arm's benefit in microsomes minus
the same arm's benefit in hepatocytes, on the donors present in both. A
**positive** value means the benefit is larger in hepatocytes.

| Substrate | n | `W1`−`W0` in HLM | in HH | paired difference | 95 % |
|---|---|---|---|---|---|
| **Omeprazole** | 15 | +0.0084 | **−0.2012** | **+0.2095** | **[+0.0637, +0.3533]** |
| **Midazolam** | 15 | −0.0518 | **−0.1168** | **+0.0650** | **[+0.0014, +0.1279]** |
| **Diclofenac** | 15 | +0.0004 | **−0.0393** | **+0.0397** | **[+0.0060, +0.0717]** |
| Bupropion | 15 | +0.0132 | +0.0009 | +0.0123 | [−0.0079, +0.0363] |
| Bufuralol | 13 | +0.0038 | +0.0003 | +0.0035 | [−0.0002, +0.0103] |

**For three of five substrates the enzyme's abundance is informative in
hepatocytes and not in microsomes, and the paired difference excludes zero.**

For the primary contrast the system difference is small and mixed: midazolam
+0.0202 [+0.0009, +0.0439] and diclofenac −0.0119 [−0.0225, −0.0006] exclude
zero, both at magnitudes inside their cells' permutation nulls.

**This is not a causal statement about the compartment.** The two preparations
differ simultaneously in enzyme content, membrane environment, the presence of
cytosolic cofactors, and the scaling constant used to reach mL/min/kg. This
design separates none of them, and the protocol forbids attributing the
difference to any one.

---

## 2. Biological reading

Each statement carries its caveat, and none is pharmacokinetic.

- **One protein's absolute concentration is a genuinely useful predictor of its
  own substrate's clearance — in hepatocytes, for CYP3A4 and CYP2C19.** This is
  the mechanistic expectation holding where it should, measured
  out-of-sample rather than as a correlation.
- **It fails in microsomes, for every substrate including CYP3A4.** The same
  donors, the same enzyme, the same drug: only the preparation differs. In a
  project whose whole aim is predicting metabolism from molecular
  measurements, **which preparation the measurement comes from changed the
  answer more than the choice of molecular layer did.**
- **It fails for CYP2B6 and CYP2D6 in both systems.** For CYP2D6 the reason is
  visible: abundance spans 587-fold across the 15 microsomal donors while
  bufuralol `CLint,mic` spans 3.5-fold and is undetermined for two donors —
  the enzyme cannot be rate-limiting in that assay. This pilot records the
  oddity and does not explain it.
- **A two-component summary of the rest of the proteome adds nothing, and that
  much is a real negative rather than an underpowered shrug:** the direction is
  *against* a benefit in six of ten cells, the one significant interval is a
  degradation, and the k = 1 sensitivity shows five significant degradations.
  A design with no power tends to produce intervals around zero, not a
  consistent small cost. **What it does not show is that the other 4,037
  proteins carry no mechanism** (§4).
- **What this does not license.** Nothing about enzyme induction, about drug–
  drug interaction, about transporters, or about any in-vivo quantity. The
  endpoint is an in-vitro depletion clearance scaled with literature constants.

---

## 3. Analytical validity

- **The circularity check fired, and it changed the experiment.** The
  commissioning note required that the clearance definition must not contain
  the predictor. The article's Methods state that `eq. 8` divides clearance by
  the probed CYP amount **and** multiplies it back by the **homogenate** CYP
  concentration — two of Table S3's three CYP columns — and that `eq. 9`
  carries that normalisation and adds Kpuu. **Both are circular against a
  CYP-abundance predictor and are forbidden.** Only `eq. 7`, scaled by
  hepatocellularity or literature MPPGL, was used. The runner refuses each
  forbidden name and additionally asserts the guard still rejects them, so a
  guard that stopped working would fail the run rather than pass silently.
- **The donor mapping is not the identity, and assuming it would have
  mis-paired six of fifteen donors.** Table S3's `UU1`–`UU15` and Data S1's
  `1`–`15` agree for donors 1–9 and are permuted for 10–15
  (10→UU14, 11→UU15, 12→UU10, 13→UU11, 14→UU12, 15→UU13). The bijection is
  re-derived at every run from the CYP concentrations both files publish, and
  the run refuses to start unless every index resolves unanimously. In this
  execution all fifteen resolved on **8 of 8** voting vectors.
- **The target CYP is out of every proteome block.** Verified per cell in the
  artefact. For bufuralol the question does not arise: CYP2D6 is absent from
  Data S1 in both preparations.
- **The two undetermined values are excluded, never imputed.** Bufuralol's
  microsomal analysis runs on 13 donors and 13 folds; every other cell on 15.
- **Fold safety.** Column means and SDs, component loadings, predictor
  standardisation, the target mean and the ridge penalty are all estimated
  inside the training donors, including inside the inner leave-one-out that
  selects the penalty. Tests assert that replacing the held-out row leaves the
  training scores bit-for-bit identical and that the first component of a
  two-component extraction equals the one-component extraction.
- **Proteome coverage, no imputation.** Of 5,172 deposited microsomal rows,
  4,039 are complete across all 15 donors; one has a blank gene name; the
  block is 4,037 proteins after removing the target CYP (4,038 for bufuralol).
  Hepatocytes: 4,784 deposited, 3,786 complete, 3,784 in the block.
  Concentrations are strictly positive and span 4.06 × 10⁶, so they are
  `log10`-transformed.
- **No donor identifier, value or prediction is in the artefact**, and no
  numeric array of ten or more entries: the 200 permutation replicates are
  summarised by their deciles.
- **Runtime** ≈ 325 s for 60 fitted cells and 10 nulls of 200 replicates.

### 3.1 The sample-size projection, and a defect of the freeze

**The prespecified projection targets the primary contrast**, which is the one
with no effect. Its most favourable in-interval value is therefore near zero,
and the donors required to detect it run to **hundreds of thousands**
(midazolam: 610,719 at 80 % power with the five-substrate correction). That is
not an absurdity to be hidden; **it is the finding, stated as a design
requirement: the proteome-beyond-the-enzyme question is not worth sizing a
study for at this effect size.**

Taking instead the observed point estimates — which are degradations for three
of the five substrates — the projection says **51 to 131 donors** would
establish them (bufuralol 51, omeprazole 61, midazolam 72, diclofenac 100,
bupropion 131, at 80 % power with the five-substrate correction). A study
sized to confirm that extra proteomic components make the prediction slightly
worse is not a study anyone should run.

**This is a defect of the freeze, recorded in the protocol's §10 rather than
repaired.** §7.3 pointed the projection at the primary contrast; the contrast
that actually shows an effect is the secondary `W1` versus `W0`. A projection
for it was added to the artefact under
`post_hoc_q1_projection_not_prespecified`, with
`status: POST_HOC_ADDED_AFTER_SEEING_THE_RESULTS`, and is labelled as a
reporting addition everywhere it appears. The prespecified projection is
unaltered and no result was recomputed.

**The post-hoc projection, which is the one a future cohort would be sized
from.** Donors required for `W1` versus `W0`, at 80 % power, two-sided
α = 0.05 with the five-substrate Bonferroni correction. The third column is
the **least** favourable value inside this pilot's interval, which is the
conservative question when the effect is real.

| Cell | SD of paired differences | at the observed effect | at half of it | at the weakest plausible effect |
|---|---|---|---|---|
| **HH Midazolam** (CYP3A4) | 0.1624 | **23** | 91 | **149** |
| **HH Omeprazole** (CYP2C19) | 0.2773 | **23** | 89 | **196** |
| HH Diclofenac (CYP2C9) | 0.0676 | 35 | 139 | 1,517 |
| HLM Midazolam | 0.1245 | 68 | 271 | — (interval reaches a degradation) |
| HH Bupropion, HH Bufuralol, HLM Omeprazole, HLM Diclofenac, HLM Bupropion, HLM Bufuralol | — | — | — | not design targets: the observed effect is a degradation or indistinguishable from zero |

**Read as an order of magnitude: a cohort of roughly 150 to 200 donors would
establish the protein-to-clearance relationship for CYP3A4 and CYP2C19 in
hepatocytes even under the weakest version of the effect this pilot's interval
allows** — against 15 donors today. Twenty-three would suffice if the effect
is the size observed here, which is why the honest figure is the conservative
one.

For the cells whose observed effect is a degradation, "donors required to
detect the observed effect" is arithmetic without a purpose and is not quoted
as a target; the artefact carries the numbers with their signs.

All projections assume the future cohort's dispersion equals this pilot's —
estimated from 15 donors and therefore itself imprecise — the same assay, the
same preparation, the same scaling equation and the same arms, with
independent donors and a constant effect. They assume a paired two-sided test
on the per-donor difference in `log10` absolute error, and they size **one**
substrate at a time with the correction for five.

---

## 4. What may not be concluded

### 4.0 Two limits named by the reviewer of 2026-11-10, both correct

**(i) The null on `W2` is a statement about two principal components, not about
the other 4,037 proteins.** A mechanism carried by a handful of proteins with
little variance across these 15 donors would not appear in the leading
directions of a standardised 4,037-column block, and this design would not see
it. Nothing here licenses "the rest of the proteome is uninformative"; what it
licenses is "**the leading one or two directions of the rest of the proteome do
not improve on the target enzyme's own abundance, and adding them costs a
little.**" The restriction was deliberate — the commissioning note asked for a
parsimonious pilot rather than a free search over 5,172 proteins, and a free
search at n = 15 would have produced an uninterpretable winner — but the price
of that choice is exactly this limit on the conclusion. **Testing a
few-protein mechanism needs a hypothesis fixed in advance from outside these
data, or more donors; it is not available from this pilot and was never
claimed.**

**(ii) The intervals do not include the uncertainty of learning on 15 donors.**
`paired_bootstrap` resamples the **stored out-of-fold per-donor errors** and
does not refit the model inside each draw. It therefore propagates the
variability of the donors used for *scoring* and not the variability of the
components, the standardisation constants and the penalty that were estimated
on 14 donors. **Every interval in this document is therefore narrower than a
fully nested one would be.** It is reported as a descriptive interval, its
resampling unit is recorded in the artefact, and the conclusion it supports is
a negative one — a fully nested interval would be wider still and would not
turn a non-detection into a detection. Where the intervals carry a *positive*
claim (Q1, hepatocyte CYP3A4 and CYP2C19) the understatement matters more, and
those two claims should be read as **"detected under a descriptive interval
that does not price in the learning variance"**. Closing that gap means a
refitting bootstrap, which is a new computation and was not run: the
commissioning note's instruction for this stage was to stop optimising this
cohort.

- **An interval containing zero means no difference was detected.** It is never
  a demonstration that there is none, and at n = 15 the power is low by
  construction. For Q2 the evidence against a benefit rests on the *direction
  and consistency* of the estimates and on the significant degradations, not
  on non-detection alone.
- **No variance fraction from a rank correlation.** None is computed here, and
  a sentence that derived one in
  `docs/OMICS2PK_PUBLIC_SOURCE_QUALIFICATION_20261109.md` §4 was withdrawn on
  2026-11-10.
- **No RNA/protein equivalence from a wide interval.** This pilot measures no
  transcript. The separate finding that the HLC transcript correlation falls
  inside this cohort's protein interval on 5 of 5 enzymes is a statement about
  resolution, not about the layers being the same.
- **No clinical-pharmacokinetic claim.**
- **No causal compartment effect** (§1.4).
- **No biomarker panel.** Individual protein coefficients are not reported and
  no selection among the 4,037 proteins was performed: the extra information
  entered as two prespecified components and nothing else.

---

## 5. What this contributes to the funding case

**Evidence type: functional proof attempt on independent public data.**

**Supports:**

1. **A mechanistic prediction works where it should, out-of-sample.** One
   protein's absolute concentration cuts the hepatocyte clearance error of
   midazolam from ×1.55 to ×1.18 and of omeprazole from ×2.46 to ×1.55, with
   Bonferroni-adjusted intervals excluding zero at n = 15.
2. **It is enzyme-specific and preparation-specific.** It holds for CYP3A4 and
   CYP2C19 in hepatocytes, not for CYP2B6 or CYP2D6, and not in microsomes for
   any substrate. A single-enzyme demonstration would have been misleading;
   five substrates in two systems is what makes the limit visible.
3. **A parsimonious mechanistic representation is not improved by a
   two-component reduced-rank summary of the rest of the proteome** — in ten
   cells, with two declared sensitivities agreeing and a permutation null that
   explains what a naive reading would have called a small gain. The claim is
   about that representation, not about the proteome (§4).
4. **The public record's limit is now quantified rather than asserted.** The
   only public source with same-donor proteome and drug-specific intrinsic
   clearance has 15 donors, no demographics, no transcriptome and no patient.

**Argues for measuring, with a number.** The question this pilot cannot answer
is whether the proteome adds to the target enzyme *where the target enzyme
itself is informative* — which, in these data, is hepatocytes for two enzymes.
Testing that needs a cohort with (a) more donors, (b) a transcriptome on the
same donors so the layer comparison stops being cross-cohort, and (c) a
donor-matched in-vivo endpoint so the functional result can be connected to
exposure. **None of the three exists in any public dataset this project has
qualified**, and that is the gap the request should be written around.

The size of (a) is now quantified rather than asserted: **roughly 150 to 200
donors** would establish the CYP3A4 and CYP2C19 protein-to-clearance
relationship in hepatocytes under the weakest effect this pilot's interval
allows, at 80 % power with the five-substrate correction — against the 15 that
the entire public record provides (§3.1). That projection is post-hoc and
labelled as such, and it rests on a dispersion estimated from 15 donors, so it
is an order of magnitude and not a protocol.

**That number sizes one thing and one thing only.** It is the donor count for
`W1` versus `W0` — the abundance of **one** enzyme against the mean, for
**two** enzymes, on an **in-vitro** clearance, under a descriptive interval
(§4.0 ii). It is **not** a sample size for a multi-omic complementarity study,
**not** one for a clinical-pharmacokinetic study, and **not** one for a
few-protein mechanism (§4.0 i). Those designs have different effects,
different dispersions and different numbers of comparisons, and none of them
can be sized from this pilot. Quoting 150–200 for any of them would be a
misuse of it.

---

## 6. Change log

| Version | Date | Change |
|---|---|---|
| 1.1.0 | 2026-11-10 | **Three claims narrowed after review, no number recomputed.** (1) The Q2 conclusion said "the rest of the proteome adds nothing"; what was tested is **two prespecified principal components** of it, and a mechanism carried by a few low-variance proteins would be invisible in the leading directions of a standardised 4,037-column block. The title, §0, §1.2, §2 and §5 are rescoped and §4.0 (i) states the limit. (2) The intervals come from a bootstrap of **stored out-of-fold errors without refitting**, so they price the scoring donors' variability and not the variability of the components, the standardisation constants and the penalty estimated on 14 donors. **Every interval here is narrower than a fully nested one**, which matters most for the two positive Q1 claims; §4.0 (ii) says so and records that closing the gap needs a refitting bootstrap that was not run. (3) The 150-200 donor figure sizes `W1` versus `W0` for two enzymes on an in-vitro endpoint and **cannot** size a multi-omic complementarity study, a clinical-PK study, or a few-protein mechanism; §5 says so explicitly. No arm, endpoint, fold, loss, interval or projection changed. |
| 1.0.0 | 2026-11-10 | First write-up. Reports all ten primary cells with unadjusted and Bonferroni-adjusted intervals and their permutation nulls, both declared sensitivities, the paired system comparison, the three losses per arm, and the sample-size projection with the defect of the freeze that made a post-hoc projection necessary. Records that the circularity check forbade two of the article's three clearance equations and that the donor mapping is not the identity. |
