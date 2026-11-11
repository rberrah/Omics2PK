# The mechanistic rule is right, and it becomes useful when the ruler is three-dimensional

**Version 3.0.0 — 2026-11-11.** Evidence level:
**EXPLORATORY — EXPERIMENTALLY SUPPORTED ON THESE DATA, PENDING INDEPENDENT
VALIDATION.** Two genes, two knowledge layers, position- and region-grouped
validation, position-clustered intervals, five negative controls and one
positive control. Artefact: `reports/data/variant_mechanism_benchmark.json`.
Protocol and architecture:
`OMICS2PK_VARIANT_MECHANISM_SELECTION_ARCHITECTURE_20261111.md`.

The two structures used here were chosen **after** these proteins' discordant
residues had been examined. That does not make the analysis leak — no
functional score entered the choice of a crystal — but it does make it
**exploratory**, and no amount of internal validation converts it. Only a
protein declared before its functional data is read can do that; §0.3 records
the first such declaration.

---

## 0. What version 2.0.0 claimed too strongly

A reading of commit `ccda244` raised three points. The first two were measured
and both bite; the third is accepted as stated.

### 0.1 The resampling unit was wrong, and it widened every interval by 2.3–2.9×

The folds were grouped by protein position. **The bootstrap was not.** Up to
nineteen substitutions share one residue — same structural environment, same
assay artefacts — and a row bootstrap counted them as that many independent
observations. Every verdict in this project is decided by whether an interval's
**lower** bound clears a margin, so a too-narrow interval errs in the one
direction that manufactures a verdict.

Measured, by resampling **positions** with replacement instead of rows:

| Headline contrast | row bootstrap, as published | **position-clustered** | width |
|---|---|---|---|
| NUDT15, position, `M4` vs `M1`, sequence layer | [+0.0077, +0.0206] `GAIN_PRESENT_MARGIN_UNDECIDED` | **[−0.0041, +0.0334] `NO_DETECTABLE_GAIN_UNDERPOWERED`** | ×2.89 |
| NUDT15, region, `M4` vs `M1`, sequence layer | [−0.0160, −0.0007] `DETECTABLY_WORSE` | **[−0.0294, +0.0129] `REFUTED_GAIN_BELOW_MARGIN`** | ×2.76 |
| CYP2C9, position, `M4` vs `M1`, sequence layer | [+0.0020, +0.0073] `GAIN_PRESENT_BUT_BELOW_MARGIN` | **[−0.0016, +0.0104] `NO_DETECTABLE_GAIN_UNDERPOWERED`** | ×2.26 |
| CYP2C9, region, `M4` vs `M1`, sequence layer | [−0.0051, +0.0011] `REFUTED_GAIN_BELOW_MARGIN` | [−0.0101, +0.0057] `REFUTED_GAIN_BELOW_MARGIN` | ×2.55 |

**Withdrawn, and the withdrawal cuts both ways.** Version 1.0.0's statement
that the sequence layer made `M4` **detectably worse** under the region split
goes with it: the gain was never detectable in either direction. What survives
is weaker and more boring — with a sequence ruler, adding variant and mechanism
information to abundance **changes nothing measurable**.

**What survives unchanged.** The enrichment of discordance near catalytic
machinery is robust to the unit: with the 3D ruler, NUDT15 **55.35×
[21.09, 277.53]** position-clustered against [26.95, 251.41] row-wise, and
CYP2C9 **29.70× [9.31, 113.56]** against [10.78, 105.08]. The sequence ruler's
own enrichment widens a lot — 6.74× from [4.82, 9.42] to **[2.15, 17.83]** —
and still excludes 1. The mechanistic rule itself is not at stake here; the
published interval around it was simply too narrow.

A cluster bootstrap does **not** always widen an interval: it widens one whose
statistic is positively correlated inside a cluster, and can narrow one that is
not. The positive control is a case of the latter, and an assertion written the
other way round fired on it. The direction is measured, not assumed — the same
discipline this project applied to the arm-selection leak.

### 0.2 There was no direct contrast of "with structure" against "without"

`M4` beating `M1` says the union beats the best single measurement. It cannot
say what the **six structural variables** contributed, because `M4` already
held the sequence mechanism block. Version 2.0.0 reported only a difference of
two mean losses, with no interval and no pairing. The paired contrast now
exists, on identical rows and identical folds — see §1.4. It is the stronger
result, because it isolates the knowledge layer instead of the model.

### 0.3 The structures were chosen after the discordant residues were seen

Accepted. Two mitigations were run and one is declared for later.

* **The thresholds are not carrying the result.** Removing all three distance
  thresholds and keeping only the two continuous log-distances leaves the gain
  intact (§1.5). Enrichment excludes 1 at **every** cutoff from 4 Å to 15 Å.
* **8 Å is nonetheless the largest point estimate** in the NUDT15 sweep below
  10 Å, which is exactly the signature of a value chosen after seeing data. It
  is kept because it was declared in code before this run, and the conclusion
  does not depend on it.
* **A third protein is pre-registered**, in
  `OMICS2PK_BRAF_PREREGISTRATION_20261111.md`, written and committed before any
  BRAF score was downloaded. Until it runs, the claim here reads
  **experimentally supported on these data, pending independent validation**,
  and the word *confirmed* has been removed from every structural claim below.

---

**Not a pharmacokinetic result.** The target is an activity score from a
functional screen — not a clearance, not an AUC, not a concentration. **Not a
multi-omic measurement in individuals** either: these are engineered single
variants in cell and yeast systems.

---

## 1. Three findings with a sequence ruler, and what a structural ruler changed

| | Result |
|---|---|
| **The premise is true.** The measurement closest to the gene is not reliably the one that reports its function | ρ(abundance, activity) = **+0.340** for NUDT15 but **+0.757** for CYP2C9. On NUDT15, **4.3 %** of variants keep abundance in the top quartile and lose activity to the bottom |
| **The mechanistic rule holds** | a variant within 3 residues of an annotated functional site is **6.74× [2.15, 17.83]** more likely to be discordant. The interval excludes 1, and this is the position-clustered one — the published [4.82, 9.42] was 2.9× too narrow, §0.1 |
| **With a sequence ruler it was not a useful feature** | conditional use of abundance (`M5` vs `M4`) reaches the declared 5 % margin in **none** of the four cells, and predicting *where* the assays disagree gives ρ **+0.083 [+0.045, +0.120]** — real, and far too small to act on |

**And the reason the second and third coexist is measurable:** the annotation
reaches only **48.8 %** of NUDT15's discordant variants and **12.5 %** of
CYP2C9's. The rule is correct about the cases it can see and blind to the rest.

### 1.1 That diagnosis was refutable, and it was tested — it holds

Version 1.0.0 concluded that **annotation coverage, not the algorithm**, was
the binding constraint. That is a prediction: a ruler that measures distance in
**three dimensions** should reach the discordant residues a **sequence** ruler
misses. Tested by adding two crystal structures — **1R9O** (CYP2C9 with heme
and flurbiprofen, 2.0 Å) and **5LPG** (NUDT15 with Mg and 6-thio-GMP, 1.7 Å):

| | sequence ruler, ≤3 residues | **3D ruler, ≤8 Å** |
|---|---|---|
| NUDT15 discordant variants reached | 48.8 % | **95.9 %** |
| NUDT15 enrichment | 6.74× [2.15, 17.83] | **55.35× [21.09, 277.53]** |
| CYP2C9 discordant variants reached | 12.5 % | **87.5 %** |
| CYP2C9 enrichment | 9.55× **[0.00, 31.67]** | **29.70× [9.31, 113.56]** |

**CYP2C9's interval excludes 1 only after the structural layer is added**, so
on that gene the rule moved from *untestable* to *supported*. All intervals in
this table resample **positions**, not substitutions (§0.1).

**And prediction improved too, which version 1.0.0 did not promise.** On
NUDT15, position-grouped, the mechanism arm `M3` improves **+10.4 %** and:

| Contrast, NUDT15 | sequence annotation | **+ structure** |
|---|---|---|
| `M3` vs `M0`, position-grouped | +6.92 % undecided | **+16.61 % [+0.0273, +0.0750] `GAIN_REACHES_MARGIN`** |
| `M4` vs `M1`, position-grouped | +5.23 %, not detectable | **+15.83 % [+0.0205, +0.0684] `GAIN_REACHES_MARGIN`** |
| `M3` vs `M0`, **region-grouped** | +2.49 %, not detectable | +9.89 % [+0.0035, +0.0596] `GAIN_PRESENT_MARGIN_UNDECIDED` |
| `M4` vs `M1`, **region-grouped** | −2.98 % `REFUTED_GAIN_BELOW_MARGIN` | +9.47 % [+0.0001, +0.0521] `GAIN_PRESENT_MARGIN_UNDECIDED` |

**One cell reaches the margin and the other does not, and the difference is the
split.** Position-grouped, knowledge-derived information beats the best single
measurement by a margin-reaching amount. Region-grouped — the harder split,
which removes the structural neighbourhood from the training set as well — the
gain is **positive but the margin is undecided**: the lower bound is +0.0001,
barely clear of zero.

**Version 2.0.0 said this reached the margin under both splits. That was the
row bootstrap talking** (§0.1), and it is withdrawn. The surviving claim is the
commissioned priority met on the primary gene under the easier of two
validations, and supported but not established under the harder one.

**What still fails:** `M5`, the *conditional* arm, reaches the margin nowhere
even with structure (+4.05 % undecided, +1.25 % refuted on NUDT15; +1.14 %
below margin, −3.10 % worse on CYP2C9). So knowledge helps **select**
information and does **not** help **weight** it conditionally — the same split
this project found for prediction fusion.

### 1.2 Why a sequence ruler was always going to fail

| Residue | sequence distance to the nearest annotation | 3D to cofactor | 3D to substrate |
|---|---|---|---|
| CYP2C9 **Thr301**, most discordant | **134 residues** | **3.15 Å** | 3.83 Å |
| CYP2C9 Val113 | 322 residues | 3.65 Å | 3.72 Å |
| CYP2C9 Ser365 | 70 residues | 2.88 Å | 6.87 Å |
| NUDT15 Gly137 | 70 residues | 14.04 Å | **3.18 Å** |
| NUDT15 Gly15 | 32 residues | 10.05 Å | 4.90 Å |

Thr301 sits in the I-helix, which packs directly against the heme. It is in
contact with the catalytic centre and 134 positions away in sequence. A
one-dimensional ruler cannot see a folded protein, and NUDT15's missed residues
are reached mostly through the **substrate** pocket rather than the metal.

### 1.4 What the six structural variables themselves contribute

The contrast above compares two *arms*. This one compares one arm against
itself, with and without the structural block, on identical rows and identical
folds, paired and position-clustered. Nothing else differs between the two
fits.

| Arm, NUDT15 | position-grouped | region-grouped |
|---|---|---|
| `M3` mechanism only | +10.41 % [+0.0118, +0.0496] undecided | +7.59 % [−0.0067, +0.0530] not detectable |
| **`M4` fixed union** | **+11.19 % [+0.0118, +0.0494] undecided** | **+12.09 % [+0.0068, +0.0622] undecided** |
| `M5` conditional interaction | +13.46 % [+0.0172, +0.0540] **`GAIN_REACHES_MARGIN`** | +13.79 % [+0.0109, +0.0686] undecided |

| Arm, CYP2C9 | position-grouped | region-grouped |
|---|---|---|
| `M3` mechanism only | +2.43 % below margin | +3.84 % undecided |
| `M4` fixed union | +1.95 % `REFUTED_GAIN_BELOW_MARGIN` | +3.81 % undecided |
| `M5` conditional interaction | +2.73 % not detectable | +3.67 % undecided |

Three things follow that the `M4`-vs-`M1` contrast could not show.

1. **The structural block is detectably positive on NUDT15 under both
   splits** — +11.19 % and +12.09 %, both lower bounds above zero. It is the
   *isolated* contribution of the knowledge layer, and it is more stable across
   splits than the composite `M4`-vs-`M1` gain is.
2. **Under the region split the structural block contributes more (+12.09 %)
   than the whole union gains over abundance (+9.47 %).** The sequence
   mechanism block was costing accuracy there; the structure more than pays for
   it. The net figure was hiding two effects of opposite sign.
3. **CYP2C9 gains 2–4 %, below the margin.** Consistent with its ρ = 0.757
   between abundance and activity and 0.5 % discordance: where abundance is
   already a good proxy there is little for a mechanism layer to add. That is
   the system working, not failing.
4. **The arm the structure helps most is `M5`** — and `M5` still does not beat
   `M4`. Better knowledge raises the conditional arm's absolute accuracy
   without making conditional *use* of abundance pay off. The split between
   selecting and weighting is not an artefact of a thin knowledge layer.

### 1.5 The distance thresholds are not carrying the result

The structural block holds two continuous log-distances, three threshold
indicators (contact at 5 Å, second shell at 8 Å) and an unresolved flag.
Removing **all three thresholds**:

| `M4` fixed union | all six columns | **continuous distances only** |
|---|---|---|
| NUDT15, position-grouped | +11.19 % [+0.0118, +0.0494] | **+10.46 % [+0.0102, +0.0470]** |
| NUDT15, region-grouped | +12.09 % [+0.0068, +0.0622] | **+13.61 % [+0.0120, +0.0655]** |
| CYP2C9, position-grouped | +1.95 % | +1.66 % |
| CYP2C9, region-grouped | +3.81 % | +4.22 % |

**No threshold appears in the threshold-free model, so no threshold can be
carrying its gain.** The continuous distance to the catalytic machinery is the
whole of the signal.

The enrichment sweep says the same for the descriptive statistic. Every cutoff
from 4 Å to 15 Å gives an interval excluding 1:

| cutoff | NUDT15 enrichment, position-clustered | share of discordant reached |
|---|---|---|
| 4 Å | 12.58× [5.01, 33.73] | 52.1 % |
| 5 Å | 21.41× [9.53, 80.08] | 76.9 % |
| 6 Å | 17.81× [7.05, 92.61] | 80.2 % |
| **8 Å** | **55.35× [21.09, 277.53]** | **95.9 %** |
| 10 Å | 27.28× [10.42, 132.12] | 95.9 % |
| 12 Å | 77.39× [18.53, 113.79] | 99.2 % |
| 15 Å | 30.74× [7.13, 48.57] | 99.2 % |

**8 Å is the largest point estimate below 10 Å, which is the signature of a
threshold chosen after seeing data.** It is kept because it was frozen in code
before this run and because the conclusion holds at every value in the sweep —
not because it is optimal. On CYP2C9 the sweep runs 36.28× at 4 Å to 54.26× at
10 Å, and becomes unbounded at 12 Å and beyond, where **every** discordant
variant falls inside the cutoff and there is no comparison group left.

### 1.3 The numbering was verified, not trusted

A crystal structure's residue numbers need not match the reference sequence,
and an offset would produce **plausible wrong distances** rather than an error.
The guard asserts a fact that is only true under correct numbering: UniProt
annotates **Cys435 as CYP2C9's axial heme ligand**, so its sulphur must lie
within a coordinate bond of the iron. Measured: **2.378 Å**. For NUDT15 the
nearest Nudix-box residue lies **2.11 Å** from the catalytic Mg. Tests confirm
the guard fires when the numbering is shifted, the residue identity is wrong,
or the metal is moved away.

---

## 2. The data

Four MaveDB score sets, all **CC0**, downloaded from `api.mavedb.org` on
2026-11-11. The pairing is what makes the question answerable and it is rare:
the **same variant library** measured by one abundance assay and one activity
assay.

| URN | n | Assay |
|---|---|---|
| `urn:mavedb:00000055-a-1` | 3100 | NUDT15 protein stability, VAMP-seq |
| `urn:mavedb:00000055-b-1` | 3100 | NUDT15 **activity**, thiopurine cytotoxicity |
| `urn:mavedb:00000095-b-1` | 6821 | CYP2C9 abundance, VAMP-seq, HEK293T |
| `urn:mavedb:00000095-a-1` | 6524 | CYP2C9 **activity**, Click-seq, humanised yeast |

After parsing protein-level changes and intersecting the two assays:
**NUDT15 n = 2 844 over 163 positions**; **CYP2C9 n = 4 478 over 464
positions**. Synonymous variants are dropped; nonsense is kept, because a
premature stop is one of the mechanism classes under test.

**NUDT15 is the primary gene.** Its two assays are both in human cells.
CYP2C9's abundance assay is in HEK293T and its activity assay is in humanised
yeast, so a CYP2C9 discordance may be a difference of *system* rather than of
biology. This design cannot separate the two, so every CYP2C9 figure here is a
sensitivity analysis with that confound attached.

### 2.1 The annotation, and how thin it is

| Gene | UniProt | Annotated functional sites | Motifs |
|---|---|---|---|
| NUDT15 | Q9NV35, 164 aa | 47, 49, 63, 67 | Nudix box 48–69 |
| CYP2C9 | P11712, 490 aa | **435** — one residue, the heme axial ligand | none |

**This asymmetry is a finding, not a configuration detail.** It is why the two
genes give enrichment intervals of very different widths, and it is the subject
of §6.

---

## 3. Why the folds are grouped, and why it changed a conclusion

Four different substitutions at NUDT15 Arg34 sit in the discordant set. A
random split would train on `p.Arg34Ala`, test on `p.Arg34Gly` and call that
out-of-sample. Every split here is grouped by **protein position**, and
repeated with contiguous **sequence region** blocks, because neighbouring
residues share structural context too.

**The instruction to validate on unseen positions or regions was load-bearing.**
`M4` vs `M1` — does adding variant and mechanism information to abundance help?
— answers differently under the two splits:

| Gene, sequence layer | position-grouped | region-grouped |
|---|---|---|
| NUDT15 | +5.23 % [−0.0041, +0.0334] `NO_DETECTABLE_GAIN_UNDERPOWERED` | −2.98 % [−0.0294, +0.0129] `REFUTED_GAIN_BELOW_MARGIN` |
| CYP2C9 | +2.40 % [−0.0016, +0.0104] not detectable | −1.00 % [−0.0101, +0.0057] refuted |

With a position-clustered interval the sequence layer adds **nothing
measurable, in either direction, under either split**. Version 1.0.0 read the
pair as a sign reversal and version 2.0.0 repeated it; both were reading a
row-wise interval about 2.8× too narrow (§0.1).

The split still matters, and the structural layer shows where: there,
position-grouped reaches the margin and region-grouped does not (§1.1). A
position-grouped split is too easy — it leaves the structural neighbourhood in
the training set — so the region-grouped figure is the one to believe.

---

## 4. Results

### 4.1 NUDT15, primary

| Arm | position-grouped | region-grouped |
|---|---|---|
| `M0` intercept | 0.304 | 0.307 |
| `M1` abundance only | 0.269 | 0.272 |
| `M2` substitution only | 0.291 | 0.295 |
| `M3` mechanism only | 0.283 | 0.300 |
| `M4` fixed union | **0.255** | 0.280 |
| `M5` conditional interaction | 0.251 | 0.282 |
| `M6` trust gate | 0.259 | **0.269** |

| Contrast | position-grouped | region-grouped |
|---|---|---|
| **`M5` vs `M4`** | +1.53 % `GAIN_PRESENT_BUT_BELOW_MARGIN` | −0.71 % `REFUTED_GAIN_BELOW_MARGIN` |
| `M6` vs `M4` | −1.55 % `DETECTABLY_WORSE` | +4.23 % `GAIN_PRESENT_MARGIN_UNDECIDED` |
| `M4` vs `M1` | +5.23 % `GAIN_PRESENT_MARGIN_UNDECIDED` | −2.98 % `DETECTABLY_WORSE` |
| `M1` vs `M0` | **+11.26 % `GAIN_REACHES_MARGIN`** | **+11.42 % `GAIN_REACHES_MARGIN`** |
| `M3` vs `M0` | +6.92 % `GAIN_PRESENT_MARGIN_UNDECIDED` | +2.49 % `REFUTED_GAIN_BELOW_MARGIN` |
| `M5` vs `M4`, permuted mechanism | −0.06 % refuted | +0.03 % refuted |

### 4.2 CYP2C9, sensitivity — assay-system confound applies

| Arm | position-grouped | region-grouped |
|---|---|---|
| `M0` | 0.331 | 0.333 |
| `M1` abundance only | 0.190 | **0.192** |
| `M4` fixed union | **0.185** | 0.194 |
| `M5` conditional | 0.185 | 0.199 |
| `M6` trust gate | 0.212 | 0.218 |

| Contrast | position-grouped | region-grouped |
|---|---|---|
| `M5` vs `M4` | +0.34 % refuted | −2.95 % `DETECTABLY_WORSE` |
| `M6` vs `M4` | −14.33 % `DETECTABLY_WORSE` | −12.67 % `DETECTABLY_WORSE` |
| `M1` vs `M0` | **+42.60 % `GAIN_REACHES_MARGIN`** | **+42.32 % `GAIN_REACHES_MARGIN`** |

**The only contrast that reaches the margin anywhere is `M1` vs `M0`:** the
abundance measurement itself. Everything built on top of it fails to add
anything worth having. On CYP2C9 abundance alone explains so much (+42.6 %)
that there is little left; on NUDT15 it explains +11.3 %, and the remainder is
*not* recovered by the variant or mechanism information available.

### 4.3 The two conditional forms disagree, and neither is reliable

`M5` (interaction) and `M6` (interpretable trust weight) swap places between
splits on NUDT15 and `M6` is badly worse on CYP2C9. **An unstable ranking
between two implementations of the same idea is itself evidence that the idea
is not carrying much here**: a real effect would not depend on which of two
equivalent parameterisations was used.

`M6`'s structural handicap is declared rather than hidden: it blends two
predictors that were each already pooled across regimes, so it cannot express
"trust abundance here, not there" as cleanly as an interaction. The first
version of `M5` had the same flaw and **failed its own positive control at
+3.2 %**, which is how the flaw was found.

### 4.4 The controls

| Control | Result |
|---|---|
| **positive control** | a synthetic protein where abundance is informative on one half only: `M5` recovers **+26.6 %**, `GAIN_REACHES_MARGIN`. **So a null below is about the data, not the harness** |
| permuted mechanism, synthetic | the +26.6 % collapses to −0.1 %, refuted |
| permuted mechanism, real data | −0.06 %, −0.10 %, +0.03 %: no residual advantage |
| permuted abundance | `M1` degrades, as it must |
| region-grouped split | **changes a conclusion**, §3 |
| specificity: predict abundance instead | run and reported in the artefact |

---

## 5. The finding that decides what to build next

| | NUDT15 | CYP2C9 |
|---|---|---|
| discordance rate **within 3 residues** of an annotated site | **0.168** | 0.045 |
| discordance rate beyond | 0.025 | 0.005 |
| **enrichment**, position-clustered | **6.74× [2.15, 17.83]** | 9.55× **[0.00, 31.67]** |
| variants near a site / far | 352 / 2 492 | 66 / 4 412 |
| **share of discordant variants the annotation reaches** | **48.8 %** | **12.5 %** |

Read in order:

1. **The rule is right.** On the primary gene, proximity to annotated
   machinery multiplies the odds of "protein present, activity lost" by 6.74,
   with an interval excluding 1. The commissioned table row — *missense near a
   catalytic site → protein present but inactive* — is **supported on public
   data**. The interval is [2.15, 17.83] clustered by position, not the
   [4.82, 9.42] published in version 1.0.0 (§0.1).
2. **CYP2C9 cannot test it.** Its interval is `[0.00, 31.67]` because only 66
   of 4 478 variants are within 3 residues of its single annotated residue. The
   point estimate is larger than NUDT15's and means almost nothing. **On
   CYP2C9 the rule is untestable here, not refuted** — that is a statement
   about UniProt's coverage of a cytochrome P450, not about CYPs.
3. **Coverage is the binding constraint.** Even on the better-annotated
   protein the annotation reaches **fewer than half** the discordant variants.
   A rule that is right about half the cases and silent on the rest cannot
   produce a large predictive gain, and did not.

The discordant NUDT15 positions the annotation misses — 113, 117, 137, 15 — are
32 to 70 residues from any annotated site. On CYP2C9 the most frequent
discordant position is **301**, 134 residues from the only annotation. These
are not random: they are the kind of residue that structural and
substrate-pocket knowledge would reach and a UniProt feature table does not.

---

## 6. What is not established

| | |
|---|---|
| that conditional selection cannot work | **no**, but it is now refuted with **good** annotation too: `M5` reaches the margin on none of eight cells across two knowledge layers. The remaining escape routes are a different conditional form or a different protein, not better annotation |
| that the structural gain transfers beyond these two proteins | **no.** Two proteins, both with a crystallised ligand complex. A protein with no structure, or none with a bound substrate, gets nothing from this layer |
| that a 3D distance is causal | **no.** It is a hypothesis that a variant affects catalysis rather than folding, and only the arms that use it earn it |
| that abundance is a poor proxy in general | **no.** It is the single most useful thing measured here, by a wide margin, on both genes |
| anything about CYP2C9 biology | the assay-system confound is unresolved and its annotation is one residue |
| that the two genes replicate each other | **no.** Two genes with different annotation density and different assay pairs. They are not independent evidence of one another |
| anything about pharmacokinetics | a functional-screen activity score is not a clearance, an AUC or a concentration |
| anything about measuring omics in individuals | engineered single variants in cell and yeast systems. The relevance of a measurement *type* and of a specific *biomarker* remain two questions |
| that the PGx layer was tested | **it was not.** PharmVar and CPIC annotate a low single-digit percentage of these variants; §3 of the architecture explains why they are an evaluation stratum, not a feature |

**Status ceiling: `exploratory`.**

---

## 7. Reproduction

```
python -I scripts/run_variant_mechanism_benchmark.py \
    --mavedb-dir <dir-with-the-four-csv-files> \
    --verify-uniprot \
    --out reports/data/variant_mechanism_benchmark.json
```

Frozen before execution: seed 20261112, 5 position folds, 8 region blocks,
2000 bootstrap resamples, α 0.05, Bonferroni family 8, minimum interesting
relative gain 5 %, absolute-error loss. The UniProt annotation is embedded with
its accession and `--verify-uniprot` re-fetches and fails if it has gone stale.
The score files stay outside the repository and the runner refuses to read them
from inside it. The artefact holds no per-variant values.

| Version | Date | Change |
|---|---|---|
| **3.0.0** | **2026-11-11** | **Three review points measured; two headline claims withdrawn, one strengthened.** (1) The bootstrap now resamples **protein positions**, not substitutions: intervals widen ×2.26–×2.89 and four verdicts move, including the withdrawal of `DETECTABLY_WORSE` for the sequence layer and of `GAIN_REACHES_MARGIN` for `M4` vs `M1` under the region split (now +9.47 %, margin undecided). (2) A paired **ablation** of the structural block itself, absent before, puts its isolated contribution at **+11.19 % / +12.09 %** on NUDT15 under the two splits — more stable than the composite contrast it replaces. (3) A **threshold-free** variant reproduces the gain (+10.46 % / +13.61 %), and enrichment excludes 1 at every cutoff from 4 Å to 15 Å, so no threshold is carrying the result; the word *confirmed* is withdrawn in favour of **experimentally supported, pending independent validation**, and BRAF is pre-registered as the independent test. |
| 2.0.0 | 2026-11-11 | **The version 1.0.0 diagnosis was tested and holds.** Adding a 3D structural layer (1R9O, 5LPG) raises discordant-variant coverage from 48.8 % to **95.9 %** (NUDT15) and 12.5 % to **87.5 %** (CYP2C9), enrichment from 6.74× to **55.35×** and from 9.55× **[0.00, 23.65]** to **29.70× [10.78, 105.08]** — so CYP2C9 moved from untestable to confirmed. Prediction improved as well: `M4` vs `M1` on NUDT15 goes from +5.23 %/−2.98 % (sign flip between splits) to **+15.83 %/+9.47 %, reaching the margin under both**. `M5`, the conditional arm, still reaches the margin nowhere. Residue numbering verified by the Cys435–Fe coordinate bond at 2.378 Å. |
| 1.0.0 | 2026-11-11 | First issue. The mechanistic rule confirmed at 6.74× [4.82, 9.42] and refuted as a feature — 0 of 4 cells reach the margin; annotation coverage identified as the binding constraint at 48.8 % and 12.5 %; the region-grouped split reverses the run's most promising contrast. |
