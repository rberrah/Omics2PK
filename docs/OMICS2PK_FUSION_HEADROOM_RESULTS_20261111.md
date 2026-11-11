# Fusion headroom — the corrected gate, executed on ARIANA, Lee and HLC

**Version 1.0.0 — 2026-11-11.** Evidence type: **executed diagnostic, exploratory**.
Protocol: `OMICS2PK_FUSING_REFUTATION_PROTOCOL_20261110.md` **v2.0.0**.
Module: `src/omics2pk/fusion_headroom_diagnostics.py`.
Runner: `scripts/run_fusion_headroom_diagnostics.py`.
Artefact: `reports/data/fusion_headroom_diagnostics.json`.

**No fusion method is implemented, no gate is fitted, and no frozen result is
modified.** Every number here is computed from arms that already exist.

---

## 0. The result, and the error it corrects

**The old stopping rule would have been wrong on this data.** Version 1.0.0 of
the protocol would have stopped the programme when the **constant**-weight
ceiling showed no detectable headroom. On the thirteen cells measured here that
ceiling is **refuted or undetectable in 9 of 13** — and the **per-individual**
ceiling, which is the one that bounds a gate, shows headroom in **13 of 13**.

| Ceiling | Verdict counts over 13 cells |
|---|---|
| **constant-weight** (bounds stacking only) | 5 `REFUTED_CEILING_BELOW_MARGIN`, 4 `NO_DETECTABLE_HEADROOM`, 4 `HEADROOM_PRESENT` |
| **per-individual convex** (bounds everything, gates included) | **13 `HEADROOM_PRESENT`** |

So a mathematically wrong ceiling would have discarded the hypothesis it was
meant to test, on real data, on five endpoints outright.

**And the informative test is not the ceiling.** Of 26 context probes,
**16 are `REFUTED_CONTEXT_UNINFORMATIVE`** and **4 carry signal after Bonferroni
adjustment over 20** — all four on three HLC endpoints. The headroom is real
and mostly unreachable from the context these cohorts publish.

**The one place both conditions hold** is **CYP3A4 / testosterone**: a
per-individual ceiling 32 % below the honest reference, constant-weight
stacking that actually recovers **+167.5 [+29.6, +366]** out of sample, and
**both** context probes carrying signal after adjustment (ρ +0.238 [+0.087,
+0.366] on clinical covariates and **+0.318 [+0.144, +0.476]** on the arms' own
disagreement).

---

## 1. The ceilings, all thirteen cells

Loss is mean absolute error in each cohort's own units. **Every ceiling is
computed knowing the answers and none may be quoted as a performance.**

| cohort / endpoint | n | honest best single | constant ceiling | **per-individual ceiling** | bracketing | nested stacking, recovered |
|---|---:|---:|---:|---:|---:|---:|
| ARIANA, repeat 0 | 47 | 38.318 | 35.614 | **22.063** | 0.34 | 37.677 |
| ARIANA, repeat 1 | 47 | 33.464 | 32.826 | **20.479** | 0.30 | 33.959 |
| Lee / Jackson | 15 | 0.429 | 0.396 | **0.287** | 0.20 | 0.451 |
| HLC CYP1A2 | 397 | 399.4 | 395.0 | **275.2** | 0.31 | 440.8 |
| HLC CYP2A6 | 396 | 550.2 | 544.2 | **377.0** | 0.27 | 545.6 |
| HLC CYP2B6 | 397 | 207.3 | 206.2 | **154.4** | 0.32 | 214.5 |
| HLC CYP2C8 | 387 | 52.25 | 49.54 | **35.21** | 0.29 | 51.08 |
| HLC CYP2C9 | 397 | 252.2 | 242.2 | **168.2** | 0.27 | 244.7 |
| HLC CYP2C19 | 397 | 41.94 | 41.75 | **33.03** | 0.24 | 42.01 |
| HLC CYP2D6 | 396 | 99.33 | 97.93 | **64.58** | 0.24 | 98.69 |
| HLC CYP2E1 | 397 | 1224.7 | 1224.7 | **921.8** | 0.23 | 1792.4 |
| HLC CYP3A4 midazolam | 396 | 741.9 | 689.8 | **457.9** | 0.33 | 695.1 |
| **HLC CYP3A4 testosterone** | 396 | 3234.6 | 3052.1 | **2188.3** | 0.26 | **3067.1** |

### 1.1 Three readings, all needed

**The per-individual ceiling sits 23 % to 42 % below the honest reference**, and
the bracketing rate is **0.20 to 0.34**: on a quarter to a third of individuals
the arms' predictions straddle the truth, so a gate that knew which way to lean
would be exactly right there. That is where the headroom comes from, and it is
not an artefact of many arms — HLC has only four.

**A permissive ceiling licenses nothing.** It is zero wherever the arms bracket
the truth, so it is a necessary condition and no more. It can refute. It cannot
promise. §2 is the test that discriminates.

**Constant-weight stacking mostly fails to recover anything, and twice makes
things worse.** Nested stacking is *worse* than honest selection on CYP1A2
(440.8 against 399.4) and badly so on CYP2E1 (**1792.4 against 1224.7**, +46 %)
— an endpoint where every arm sits at the clinical baseline and the fitted
weights are therefore noise. It helps on CYP3A4 testosterone (+167.5
[+29.6, +366]), CYP3A4 midazolam and CYP2C9. This reproduces, on omic arms,
MetaPK's measured finding that **averaging without good weights degrades**.

---

## 2. The context probe — the test that actually discriminates

Can any observable predict, out of sample, **which arm will win** for an
individual? Two context sets, because a refutation must not rest on having
looked in one place: the clinical covariates each cohort publishes, and the
arms' own predictions and disagreement.

Declared before the run: a minimum interesting rank association of **0.30**, and
a Bonferroni adjustment over a declared family of **20** probes.

| Probe verdict | Count of 26 |
|---|---|
| `REFUTED_CONTEXT_UNINFORMATIVE` | **16** |
| `CONTEXT_CARRIES_SIGNAL` (survives Bonferroni) | **4** |
| `CONTEXT_CARRIES_SIGNAL_UNADJUSTED_ONLY` | 2 |
| `NO_DETECTABLE_SIGNAL_UNDERPOWERED` | 4 |

### 2.1 Where the context carries signal

| endpoint | probe | ρ out of sample | 95 % | Bonferroni / 20 |
|---|---|---:|---|---|
| **CYP3A4 testosterone** | arms' disagreement | **+0.318** | [+0.200, +0.430] | **[+0.144, +0.476]** |
| **CYP3A4 testosterone** | clinical covariates | **+0.238** | [+0.133, +0.335] | **[+0.087, +0.366]** |
| CYP2D6 | clinical covariates | +0.191 | [+0.075, +0.304] | **[+0.011, +0.360]** |
| CYP2C9 | arms' disagreement | +0.174 | [+0.068, +0.273] | **[+0.009, +0.341]** |
| CYP2C9 | clinical covariates | +0.183 | [+0.070, +0.302] | [−0.015, +0.353] |
| CYP2D6 | arms' disagreement | +0.171 | [+0.046, +0.290] | [−0.023, +0.350] |

**Three endpoints of ten, and only one robustly.** CYP3A4 / testosterone is the
only endpoint where both probes survive adjustment. CYP2D6 and CYP2C9 each
survive on one probe and not the other, which is the pattern of a weak true
effect or of multiplicity, and these data cannot tell those apart.

**Four things must be said against over-reading this.**

1. **ρ ≈ 0.17 to 0.32 is a weak association.** The probe says there is
   *something* to condition on, not that a gate will convert it into a loss
   reduction. Those are different claims and only the first is tested here.
2. **Ten HLC endpoints are not ten experiments.** They are ten measurements on
   the same ~397 donors, with a spectral dependence index of about 3.55
   recorded in the H2a header. "Three of ten" is not three independent findings.
3. **The two CYP3A4 substrates disagree.** Testosterone carries signal on both
   probes; midazolam carries none on either. By this project's own rule the two
   CYP3A4 substrates count as **one** endpoint. A within-enzyme contradiction is
   either substrate-specific biology or noise, and this design cannot separate
   them — which materially weakens the one robust result.
4. **A null probe is not always a refutation.** On ARIANA (n = 47) and Lee
   (n = 15) four probes are `NO_DETECTABLE_SIGNAL_UNDERPOWERED`: the intervals
   reach past 0.30 and nothing is decided. The 16 refutations are almost all at
   n ≈ 390, which is the only place the probe has power.

### 2.2 The honest state of the two co-primary questions

| Question | Status |
|---|---|
| **Q1 — does context help choose between predictions?** | **Not refuted, and supported on one endpoint of ten.** Worth testing properly on CYP3A4 testosterone, CYP2D6 and CYP2C9, and refuted on the other seven HLC endpoints |
| **Q2 — does mechanistic knowledge help exploit that context?** | **Untested.** No gate exists, so no prior has been compared to no prior. Nothing here speaks to it either way |

---

## 3. What this changes about the roadmap

**Step 1 is complete and it did not stop the programme.** The corrected gate
returns headroom everywhere and reachable context on three endpoints, so the
condition the commissioning note set — *proceed only if the diagnostics
identify an exploitable gain and the contextual covariates are available* — is
**met on three HLC endpoints and on no other cohort**.

| Next | Why |
|---|---|
| **Step 3 on CYP3A4 testosterone, CYP2D6, CYP2C9 only** | the three endpoints where context carries out-of-sample signal. Fitting a gate where the probe is refuted would be fitting noise |
| `S`, `U`, `K` first, then `X` | Q1 before Q2, and `K` is already measured here: it recovers +167.5 on CYP3A4 testosterone and degrades on CYP1A2 and CYP2E1 |
| **not ARIANA, not Lee** | n = 47 and n = 15 cannot estimate a gate; their probes are underpowered, not informative. They supplied a clinical-endpoint headroom check and that is all they can supply |
| the mechanistic prior **last** | Q2 is only askable once a gate exists that works without it |

**And one negative control is now mandatory for a reason this run supplies:**
CYP2E1, where nested stacking degrades by 46 % because every arm sits at the
baseline, is the natural test that a gate abstains instead of fitting noise. It
was not designed as a control; the data made it one.

---

## 4. What this run does not establish

| Claim | Status |
|---|---|
| a fusion method would improve prediction | **not established.** Ceilings are hindsight bounds |
| contextual fusion is refuted | **no** — refuted on 7 of 10 HLC endpoints, not on 3 |
| constant-weight stacking is refuted | **on 5 of 13 cells**, and it demonstrably helps on 3 |
| the mechanistic prior adds anything | **untested** |
| anything about pharmacokinetics | **no.** HLC is in-vitro activity; ARIANA is one observed trough; Lee is an activity rate |
| these are independent replications | **no.** All three cohorts have been analysed before by this project; everything here is **exploratory** |

**Status ceiling:** `exploratory`. The HLC arms used are the four from
`hlc_real_cyp_validation`, not the eight from H2a — those store losses rather
than per-donor predictions, so their ceilings are not recoverable from the
frozen artefacts without a re-run.

---

## 5. Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-11 | First issue. Executes the corrected Gate 1 and the context probe on 13 cells across three cohorts. |
