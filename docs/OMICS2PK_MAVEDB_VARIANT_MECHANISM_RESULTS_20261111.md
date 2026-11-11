# The mechanistic rule is right, and it becomes useful when the ruler is three-dimensional

**Version 2.0.0 — 2026-11-11.** Evidence type: **exploratory; two genes, two
knowledge layers, position- and region-grouped validation, five negative
controls and one positive control.** Artefact:
`reports/data/variant_mechanism_benchmark.json`. Protocol and architecture:
`OMICS2PK_VARIANT_MECHANISM_SELECTION_ARCHITECTURE_20261111.md`.

**Not a pharmacokinetic result.** The target is an activity score from a
functional screen — not a clearance, not an AUC, not a concentration. **Not a
multi-omic measurement in individuals** either: these are engineered single
variants in cell and yeast systems.

---

## 1. Three findings with a sequence ruler, and what a structural ruler changed

| | Result |
|---|---|
| **The premise is true.** The measurement closest to the gene is not reliably the one that reports its function | ρ(abundance, activity) = **+0.340** for NUDT15 but **+0.757** for CYP2C9. On NUDT15, **4.3 %** of variants keep abundance in the top quartile and lose activity to the bottom |
| **The mechanistic rule is confirmed** | a variant within 3 residues of an annotated functional site is **6.74× [4.82, 9.42]** more likely to be discordant. The interval excludes 1 comfortably |
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
| NUDT15 enrichment | 6.74× [4.82, 9.42] | **55.35× [26.95, 251.41]** |
| CYP2C9 discordant variants reached | 12.5 % | **87.5 %** |
| CYP2C9 enrichment | 9.55× **[0.00, 23.65]** | **29.70× [10.78, 105.08]** |

**CYP2C9's interval excludes 1 only after the structural layer is added**, so
on that gene the rule moved from *untestable* to *confirmed*.

**And prediction improved too, which version 1.0.0 did not promise.** On
NUDT15, position-grouped, the mechanism arm `M3` improves **+10.4 %** and:

| Contrast, NUDT15 | sequence annotation | **+ structure** |
|---|---|---|
| `M3` vs `M0`, position-grouped | +6.92 % undecided | **+16.61 % `GAIN_REACHES_MARGIN`** |
| `M4` vs `M1`, position-grouped | +5.23 % undecided | **+15.83 % `GAIN_REACHES_MARGIN`** |
| `M3` vs `M0`, **region-grouped** | +2.49 % refuted | **+9.89 % `GAIN_REACHES_MARGIN`** |
| `M4` vs `M1`, **region-grouped** | **−2.98 % `DETECTABLY_WORSE`** | **+9.47 % `GAIN_REACHES_MARGIN`** |

**The contrast that reversed sign under the harder split now reaches the margin
under both.** This is the commissioned priority — *that pharmacological
knowledge helps choose which biological information is useful* — met on the
primary gene, under the harder of two validations, against the measurement that
was already the best single predictor.

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

| Gene | position-grouped | region-grouped |
|---|---|---|
| NUDT15 | **+5.23 %**, `GAIN_PRESENT_MARGIN_UNDECIDED` | **−2.98 %**, `DETECTABLY_WORSE` |
| CYP2C9 | +2.40 %, below margin | −1.00 %, refuted |

The most promising cell in the whole run **reverses sign** under the harder
split. A position-grouped split is still too easy: it leaves the structural
neighbourhood in the training set. Had only the easier split been run, this
document would have reported a +5.2 % gain.

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
| **enrichment** | **6.74× [4.82, 9.42]** | 9.55× **[0.00, 23.65]** |
| variants near a site / far | 352 / 2 492 | 66 / 4 412 |
| **share of discordant variants the annotation reaches** | **48.8 %** | **12.5 %** |

Read in order:

1. **The rule is right.** On the primary gene, proximity to annotated
   machinery multiplies the odds of "protein present, activity lost" by 6.74,
   with an interval nowhere near 1. The commissioned table row — *missense near
   a catalytic site → protein present but inactive* — is **confirmed on public
   data**.
2. **CYP2C9 cannot test it.** Its interval is `[0.00, 23.65]` because only 66
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
| **2.0.0** | **2026-11-11** | **The version 1.0.0 diagnosis was tested and holds.** Adding a 3D structural layer (1R9O, 5LPG) raises discordant-variant coverage from 48.8 % to **95.9 %** (NUDT15) and 12.5 % to **87.5 %** (CYP2C9), enrichment from 6.74× to **55.35×** and from 9.55× **[0.00, 23.65]** to **29.70× [10.78, 105.08]** — so CYP2C9 moved from untestable to confirmed. Prediction improved as well: `M4` vs `M1` on NUDT15 goes from +5.23 %/−2.98 % (sign flip between splits) to **+15.83 %/+9.47 %, reaching the margin under both**. `M5`, the conditional arm, still reaches the margin nowhere. Residue numbering verified by the Cys435–Fe coordinate bond at 2.378 Å. |
| 1.0.0 | 2026-11-11 | First issue. The mechanistic rule confirmed at 6.74× [4.82, 9.42] and refuted as a feature — 0 of 4 cells reach the margin; annotation coverage identified as the binding constraint at 48.8 % and 12.5 %; the region-grouped split reverses the run's most promising contrast. |
