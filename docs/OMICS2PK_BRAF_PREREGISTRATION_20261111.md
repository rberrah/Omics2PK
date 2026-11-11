# Pre-registration: does the structural selection rule generalise to a third protein?

**Version 1.0.1 — 2026-11-11.** Status: **DECLARED, NOT YET RUN.** This
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

*No deviation from §§2–6 has occurred.*

---

## Change log

| Version | Date | Change |
|---|---|---|
| 1.0.1 | 2026-11-11 | Appended §7: BRAF's sequence annotation declared from UniProt P15056 with the kinase domain deliberately excluded as a feature; the cofactor-only sensitivity declared for a protein substrate; and the observation, from coordinates alone, that Val600 lies 9.61 Å from the nucleotide and so outside the ruler. Still no BRAF score downloaded. |
| 1.0.0 | 2026-11-11 | Declared. BRAF, `urn:mavedb:00001315-c-1` and `-d-1`, structure 6PP9 chain A selected by rule, numbering guard verified against coordinates, primary contrast and refutation conditions fixed. Committed before any BRAF score was downloaded. |
