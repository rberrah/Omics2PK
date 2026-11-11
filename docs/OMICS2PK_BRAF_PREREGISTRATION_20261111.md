# Pre-registration: does the structural selection rule generalise to a third protein?

**Version 2.0.0 — 2026-11-11.** Status: **RUN. THE PRIMARY CONTRAST IS REFUTED.** This
document is committed **before any BRAF functional score has been
downloaded**. Nothing below may be changed after the data is read; a later
version may only *add* a results section and record deviations.

---

## 1. Why this document exists

`OMICS2PK_MAVEDB_VARIANT_MECHANISM_RESULTS_20261111.md` v3.0.0 reports that
three-dimensional distance to a protein's catalytic machinery adds predictive
information about variant activity beyond protein abundance, on NUDT15:
**+11.19 %** position-grouped and **+12.09 %** region-grouped for the isolated
structural block.

That result is **exploratory for a reason no internal validation can fix**: the
two crystal structures, 1R9O and 5LPG, were chosen *after* those proteins'
discordant residues had been examined. No functional score entered the choice
of a crystal, so nothing leaked — but the analyst had seen which residues
needed explaining. Two mitigations were run (a threshold-free model; an
enrichment sweep from 4 Å to 15 Å) and neither addresses the structure choice
itself.

**Only a protein declared before its functional data is read can.** This is
that declaration.

---

## 2. The protein, and why it is a real test rather than a safe one

**BRAF** (UniProt **P15056**), serine/threonine protein kinase.

| Criterion | BRAF |
|---|---|
| Abundance **and** activity on the *same* variant library | yes, MaveDB `urn:mavedb:00001315` |
| n | ≈ 7 400 abundance, ≈ 7 500 activity, one construct |
| Ligand-bound crystal structure with the genuine cofactor | yes, AMP-PNP + Mg |
| Same enzyme family as the two training proteins | **no** |

**The last row is the point.** NUDT15 is a Nudix hydrolase and CYP2C9 a
cytochrome P450; both are small-molecule-metabolising enzymes whose function is
a catalytic turnover at a buried active site. BRAF is a signalling kinase whose
substrate is **another protein**, and whose best-known variants act by
**relieving autoinhibition or forcing dimerisation**, not by touching the
nucleotide. Choosing a fourth hydrolase would have been a safer test and a less
informative one.

### 2.1 A reason this may fail that is not about the method

Declared in advance so it cannot become a post-hoc excuse: **if BRAF variants
change activity mainly through dimerisation and autoinhibition, distance to the
ATP site is a poor mechanism proxy, and the rule should fail here even though
it is sound for metabolic enzymes.** A null would then be evidence about the
*scope* of the rule — which proteins it applies to — and not about its truth.
That is a weaker, honest reading, and it is the one that will be reported if
the primary contrast comes back null.

---

## 3. The data, declared

| Role | URN | Declared n | Assay |
|---|---|---|---|
| **abundance** | `urn:mavedb:00001315-c-1` | 7 438 | BRAF (C-terminal construct) abundance, baseline |
| **activity** | `urn:mavedb:00001315-d-1` | 7 496 | BRAF (C-terminal construct) activity, baseline |

**One construct for both assays.** The pooled, standardised score sets
(`-0-3`, `-0-4`) mix the N- and C-terminal constructs, which would add a
construct effect to every contrast. The C-terminal pair is primary; the pooled
pair is a **declared sensitivity analysis only**, and will not be substituted
for the primary if the primary is null.

Licence: MaveDB is CC0. Score files are written **outside** the repository, in
their own empty directory, and no per-variant table enters any artefact. This
is a **functional screen**, registry endpoint class `variant_functional_screen`
— not a PK quantity, not a multi-omic measurement in individuals.

---

## 4. The structure, selected by a declared rule

**Rule, applied before reading any score:** among RCSB entries whose polymer
references UniProt P15056, keep those containing an ATP-analogue ligand
(`ANP`, `ATP`, `ADP`, `AGS`, `ACP`) **in the BRAF chain**, available in legacy
PDB format, and take the best resolution.

**Selected: `6PP9`, 2.59 Å** — BRAF:MEK1 complex. BRAF is **chain A**
(residues 447–722), with **ANP** (AMP-PNP, a non-hydrolysable ATP analogue) and
**MG** bound in that chain. Residue 600 is **VAL**, so this is wild-type BRAF,
not V600E.

| | CYP2C9 | NUDT15 | **BRAF** |
|---|---|---|---|
| cofactor | HEM, FE | MG | **ANP, MG** |
| substrate in the crystal | flurbiprofen | 6-thio-GMP | **MEK1, a protein chain** |
| resolution | 2.0 Å | 1.7 Å | **2.59 Å** |

### 4.1 Two declared weaknesses of this instrument

1. **Chain restriction is mandatory and is new code.** MEK1, chain B, carries
   its own ANP and Mg. Without restricting both residues and ligands to
   chain A, `to_cofactor` would silently measure the distance to the *other*
   kinase's nucleotide. This is the BRAF equivalent of a numbering offset:
   plausible wrong distances rather than an error.
2. **The substrate is a protein, so `to_substrate` changes meaning.** For the
   first two proteins it was the distance to a bound small molecule. Here it
   is the distance to the nearest MEK1 heavy atom — an *interface* distance,
   not a pocket distance. It is declared as a distinct quantity and reported
   separately; the primary analysis is run **with it** and the cofactor-only
   variant is a declared sensitivity.
3. **6PP9 is an inactive, αC-out conformation.** Measured before registration:
   the catalytic Lys483–Glu501 salt bridge is **absent, 12.69 Å**. Residues in
   the αC helix and activation segment therefore sit where an inactive kinase
   puts them. A crystal structure is one conformation; this one is the
   MEK1-bound inactive state.

### 4.2 The numbering guard, declared

The guard must assert a fact that is **true only under correct numbering**.
Verified against the coordinates before this registration — a structural check
that reads no functional score:

| Check | Required | Measured |
|---|---|---|
| residue 483 identity | LYS, the catalytic lysine | **LYS** |
| residue 501 identity | GLU, the αC glutamate | **GLU** |
| residue 529 identity | THR, the gatekeeper | **THR** |
| residues 594, 595, 596 | ASP, PHE, GLY — the DFG motif | **ASP, PHE, GLY** |
| residue 600 identity | VAL in wild type | **VAL** |
| Lys483 NZ to nearest ANP atom | ≤ 4.0 Å, a salt bridge to the nucleotide | **2.875 Å** |

**The guard is the nucleotide contact, not the Lys483–Glu501 salt bridge** —
because that salt bridge is broken in this conformation (§4.1), so declaring it
would have failed on correctly numbered coordinates. The failure mode a
numbering guard exists to catch is a plausible wrong answer; the failure mode it
must not have is a false alarm on a correct structure.

---

## 5. Everything else is unchanged and frozen

| Parameter | Value | Where from |
|---|---|---|
| contact / second-shell cutoffs | 5 Å / **8 Å** | `structural_annotation.FROZEN`, unchanged |
| margin, smallest interesting relative gain | **5 %** | project-wide |
| seed | 20261112 | project-wide |
| position folds / region blocks | 5 / 8 | unchanged |
| loss | absolute error | unchanged |
| **resampling unit** | **protein position**, cluster bootstrap | corrected in `593d96f` |
| family size | 8 | unchanged |
| comparators | `M0`–`M6`, unchanged | `COMPARATORS` |

**No parameter is tuned for BRAF.** If the cutoff that worked for a buried
hydrolase pocket is wrong for a kinase cleft, that is a result.

---

## 6. The primary contrast, declared singular

> **The paired ablation of the structural block on `M4_fixed_union`,
> contiguous-region-grouped, position-clustered.**

That is: the same arm, with and without the structural variables, on identical
rows and identical folds, under the **harder** of the two splits, with the
interval resampling positions.

Declared secondaries, reported but not decisive: the same ablation
position-grouped; `M3` vs `M0`; `M4` vs `M1`; `M5` vs `M4`; the cofactor-only
variant; the threshold-free variant; the pooled-construct sensitivity.

### 6.1 The prediction, with numbers

| | Predicted |
|---|---|
| discordance enrichment near the ATP site, 8 Å, position-clustered | **> 1, interval excluding 1** |
| share of discordant variants the 8 Å ruler reaches | **> 70 %** |
| **primary**: structural ablation on `M4`, region-grouped | **gain > 0 with a lower bound above 0** |
| `M5` vs `M4`, the conditional arm | **does not reach the margin**, as on both earlier proteins |
| ρ(abundance, activity) | **between +0.2 and +0.7** |

### 6.2 What refutes it

* **The primary contrast's interval contains zero** → the structural rule does
  not transfer to a signalling kinase. Report as such; do not substitute a
  secondary.
* **The primary contrast's upper bound falls below the 5 % margin** → a
  refutation proper, not an absence of power.
* **The ablation is positive position-grouped and not region-grouped** → the
  same fragility the sequence layer showed, and the rule is reported as
  split-dependent.
* **`M5` reaches the margin** → the selecting-versus-weighting split, which has
  held on five proteins' worth of cells, would be broken, and this document's
  §6.1 would be wrong in the interesting direction.
* **The enrichment interval contains 1 while coverage is high** → the ruler
  reaches the discordant variants and they are not enriched, which refutes the
  mechanism claim rather than the coverage claim.

### 6.3 What will be reported whatever happens

The primary contrast, its interval, its verdict, and the positive and negative
controls. **An absent gain is a valid scientific result** and is the outcome
this document exists to make reportable. No arm will be added, no cutoff moved
and no gene swapped after the scores are read; any deviation forced by the data
— a parse failure, an unusable construct — will be recorded in §7 as a
deviation, with the reason, before any contrast is quoted.

---

## 7. Deviations and additional declarations

Append-only. Nothing here was written after a BRAF score was read.

### 7.1 The sequence-layer annotation, declared at implementation time

Version 1.0.0 fixed the structural layer, the data, the parameters and the
primary contrast, but not BRAF's **sequence** annotation, which is new per-gene
data. It is declared here, before download, from UniProt P15056
(`rest.uniprot.org`, retrieved 2026-11-11), and mapped by the convention the
first two genes used — point sites to `functional_sites`, a **short**
functional range to `motifs`:

| UniProt class | Positions | Mapped to |
|---|---|---|
| Active site, proton acceptor | 576 | `functional_sites` |
| Binding site | 235, 248, 251, 261, 264, 269, 272, 280, 483 | `functional_sites` |
| Binding site, range | 463–471, the ATP glycine-rich loop | `motifs` |
| Domain | 155–227 RBD; 457–717 protein kinase | **recorded, not used as a feature** |

**Why the kinase domain is not a feature.** The assayed region lies inside
457–717, so an "in domain" indicator would be near-constant: no information,
dressed as knowledge. The RBD and kinase domain are recorded in the
annotation block as `domains_not_used_as_features` so the choice is visible
rather than silent. The eight cysteine-rich-region binding sites are **kept**
even though they are zinc sites in a different functional module, because
dropping them would be a hand-picked restriction of what "annotated functional
site" means.

### 7.2 One structural column changes meaning, and one arm removes it

`to_substrate` is an interface distance to the MEK1 chain here, not a pocket
distance to a small molecule (§4.1). The **cofactor-only** arm — log distance
to AMP-PNP/Mg, the cofactor contact indicator and the unresolved flag — is
therefore run as a declared sensitivity. The second-shell indicator is
**excluded** from that arm because it is computed from
`min(cofactor, substrate)` and would smuggle the substrate distance back in.

### 7.3 Noted before the run, from the coordinates alone

**Val600 is 9.61 Å from the nucleotide and 13.66 Å from the MEK1 interface** —
outside the 8 Å ruler. The best-known BRAF variant is one the structural rule
cannot see, which is §2.1's declared risk showing up in the geometry before
any score has been read.

### 7.4 One deviation: a control added after the scores were read

§6.3 says no arm will be added after the scores are read. One was, and the
reason matters: it was added to **attack an unexpected positive**, not to
rescue the null.

`M6`, the trust gate, came back positive on BRAF — the first time any
conditional-weighting arm has reached the margin on real data in this project.
Since `M4`, the union, is *worse* than abundance alone here, the obvious
alternative explanation is that `M6` wins by **regularisation**: it blends two
small predictors instead of fitting one large union. The control that separates
the two is a **constant-weight** blend of the same two components, with the
weight chosen inside the training folds. It was not pre-registered and it was
run after `M6`'s result was seen.

The brief this project works to says to seek refutation actively. Withholding
a control that could kill a surprise finding would be the opposite. It is
reported in §8.4 and **labelled post-hoc throughout**; it has no bearing on
the primary contrast, which was refuted before it was run.

*No other deviation from §§2–6 has occurred. No cutoff was moved, no gene
swapped, and the primary contrast was not substituted.*

---

## 8. Results

Run 2026-11-11 on `urn:mavedb:00001315-c-1` and `-d-1` as declared.
**n = 6 583 variants over 340 positions**, 1 237 at positions 6PP9 does not
resolve. The numbering guard passed: seven residue identities and Lys483 NZ at
**2.875 Å** from AMP-PNP. Artefact:
`reports/data/variant_mechanism_benchmark.json`.

### 8.1 Every quantitative prediction in §6.1 failed

| §6.1 predicted | Measured | |
|---|---|---|
| enrichment near the ATP site > 1, excluding 1 | **0.655, [0.42, 0.96]** | **excludes 1 from below** |
| share of discordant variants the 8 Å ruler reaches > 70 % | **21.2 %** | failed |
| **primary**: structural ablation on `M4`, region-grouped, lower bound > 0 | **−1.60 %, [−0.0333, +0.0040]** `REFUTED_GAIN_BELOW_MARGIN` | **refuted** |
| `M5` vs `M4` does not reach the margin | refuted or detectably worse in all four cells | held |
| ρ(abundance, activity) between +0.2 and +0.7 | **−0.290** (−0.376 excluding nonsense) | **wrong sign** |

**The primary contrast is a refutation proper, not an absence of power.** Its
upper bound, +0.0040, lies below the 5 % margin of 0.0436. Position-grouped it
is worse still: **−5.82 % [−0.0770, −0.0215]**, `DETECTABLY_WORSE`. Adding
three-dimensional distance to BRAF's catalytic machinery **costs** accuracy.

No cutoff rescues it. The sweep runs 0.55× at 4 Å to 0.92× at 15 Å and never
exceeds 1. Nor does the sequence ruler: 1.15× [0.46, 2.14], reaching 4.7 % of
discordant variants. The threshold-free and cofactor-only arms are negative or
null as well; the single exception is `M3` with the cofactor-only block
region-grouped, **+4.51 % [+0.0147, +0.0627]**, which says the **MEK1 interface
column was actively harmful** — a declared sensitivity doing its job.

### 8.2 Why it failed, and why that is a result rather than an excuse

§2.1 declared in advance that BRAF variants may act through dimerisation and
autoinhibition rather than at the nucleotide, in which case the rule should
fail here while remaining sound for metabolic enzymes. **The data says that is
what happened, and says it in three independent ways.**

1. **This is a gain-of-function library.** Mean activity score: nonsense
   **+0.729**, synonymous **+0.992** (wild type ≈ 1), missense **+1.420**. The
   average missense substitution makes BRAF **more** active than wild type.
   BRAF is autoinhibited; most substitutions relieve that autoinhibition.
2. **Abundance is an *inverted* proxy, not a weak one.** ρ = **−0.290**.
   Variants that destabilise tend to signal **more**. Both assays are oriented
   the same way — nonsense is lowest in both — so this is not a scale artefact.
3. **The rule under test presumes loss of function.** "Missense near catalytic
   machinery → protein present but inactive" is a statement about losing
   activity. On a library dominated by activation it describes a minority of
   variants, and the quadrant it names is **depleted** near the machinery
   rather than enriched.

**Val600 was the warning, visible before any score was read** (§7.3): the
best-known BRAF variant sits **9.61 Å** from the nucleotide and 13.66 Å from
the MEK1 interface, outside the 8 Å ruler.

### 8.3 What this costs the NUDT15 result, and what it does not

**It does not refute it.** A rule that holds for two small-molecule-metabolising
enzymes and fails for a signalling kinase whose variants act allosterically is
a rule with a **scope**, and the scope is now measured rather than assumed.

**It does remove the right to generalise.** The claim that survives is: *for an
enzyme whose variants act by disrupting catalysis at a ligand-binding site,
three-dimensional distance to that site adds information beyond protein
abundance.* Every clause is load-bearing, and the first one is what BRAF
violates. A fourth protein should be chosen to test the clause, not the
conclusion: an enzyme with a crystallised ligand complex **and** a
loss-of-function screen.

**The honest reading of the whole programme so far** is that the selection rule
needs the mechanism class to match, and that **nothing in this architecture
detects that mismatch automatically.** A system that chose measurements for
BRAF using this rule would have chosen badly and had no way to know. That is
the gap the next experiment should attack, and it is more interesting than a
fourth confirmation would have been.

### 8.4 One unexpected positive, post-hoc, and attacked

`M6`, the trust gate — which learns, **without** abundance, how much abundance
can be trusted and blends accordingly — is the best arm on BRAF under both
splits. It uses the **sequence** mechanism block only; no structure.

| Contrast, post-hoc | position-grouped | region-grouped |
|---|---|---|
| `M6` vs `M4` | +7.63 % [+0.0287, +0.0968] undecided | **+9.84 % [+0.0575, +0.1129] `GAIN_REACHES_MARGIN`** |
| `M6` vs `M1` | **+7.24 % [+0.0424, +0.0780] `GAIN_REACHES_MARGIN`** | +6.37 % [+0.0352, +0.0704] undecided |
| `M6` vs `M0` | **+10.13 % `GAIN_REACHES_MARGIN`** | **+8.66 % `GAIN_REACHES_MARGIN`** |

**The regularisation explanation was tested and does not account for it**
(§7.4). Against the best **constant** weight chosen inside the training folds:

| | position-grouped | region-grouped |
|---|---|---|
| gate vs constant weight | **+5.10 % [+0.0342, +0.0478]** undecided | **+3.62 % [+0.0236, +0.0350]** below margin |
| constant weight vs `M1` | +2.26 %, refuted below margin | +2.86 %, below margin |

The constant blend captures only a third to a half of the gate's gain over
abundance alone. **The conditionality is doing real work** — detectably, under
both splits, with every lower bound above zero. In-fold constant weights landed
at 0.50–0.80, so the gate is not merely reproducing a fixed compromise.

**Three reasons to hold this loosely.** It is **post-hoc**: `M6` was a declared
secondary, the constant-weight control was not declared at all, and it was run
after the result was seen. Region-grouped, the gate's advantage over a constant
weight is **below the 5 % margin** — present, and by this project's own
standard not yet interesting enough to build on. And it is **one protein**,
the one where abundance happens to be an inverted proxy, which is exactly the
situation a trust gate should exploit.

**It is nonetheless the first real-data evidence in this project that
conditional weighting beats constant weighting.** It stands against a
contextual gate failing on 0 of 30 HLC cells and `M5` failing on 0 of 8
MaveDB cells. The difference worth testing next: those were cases where the
best measurement was *already good*, and this is a case where it is
**actively misleading**. That is a pre-registrable hypothesis, and it needs
its own registration, not an extension of this one.

---

---

## Change log

| Version | Date | Change |
|---|---|---|
| **2.0.0** | **2026-11-11** | **Run. The primary contrast is REFUTED.** The structural ablation on `M4`, region-grouped and position-clustered, gives −1.60 % with an upper bound of +0.0040 — below the 5 % margin, so a refutation proper and not an absence of power; position-grouped it is `DETECTABLY_WORSE` at −5.82 %. Enrichment near the ATP site is **0.655 [0.42, 0.96]**, excluding 1 from below, and the 8 Å ruler reaches 21.2 % of discordant variants against the 70 % predicted. ρ(abundance, activity) is **−0.290**, the opposite sign to the prediction: this is a gain-of-function library (missense mean activity +1.420 against a wild type of ~1.0), so "protein present but inactive" is not the dominant mode of dysfunction and the rule under test does not apply. §2.1 declared this failure mode in advance. One unexpected positive, `M6` the trust gate, reported post-hoc in §8.4 with a constant-weight control that does not explain it away. |
| 1.0.1 | 2026-11-11 | Appended §7: BRAF's sequence annotation declared from UniProt P15056 with the kinase domain deliberately excluded as a feature; the cofactor-only sensitivity declared for a protein substrate; and the observation, from coordinates alone, that Val600 lies 9.61 Å from the nucleotide and so outside the ruler. Still no BRAF score downloaded. |
| 1.0.0 | 2026-11-11 | Declared. BRAF, `urn:mavedb:00001315-c-1` and `-d-1`, structure 6PP9 chain A selected by rule, numbering guard verified against coordinates, primary contrast and refutation conditions fixed. Committed before any BRAF score was downloaded. |
