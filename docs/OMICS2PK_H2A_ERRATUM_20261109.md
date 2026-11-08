# H2a — erratum: one endpoint withdrawn, and a platform defect that changes the reading

**Erratum version:** 1.0.0
**Date:** 2026-11-09
**Applies to:** the frozen H2a run, protocol v1.4.3, merged result `reports/data/h2a_hlc_selection.json`
**Status:** no H2a number is recomputed. One endpoint's mechanistic arms are withdrawn as uninterpretable, and one conclusion I reported is retracted.

---

## E1 — Withdrawn: the CYP1A2 result cannot be read as a result about CYP1A2

**What I reported.** That CYP1A2 and CYP2E1 were cases where the target enzyme's own transcript carried no more information than age and sex, and that NetPK failed on CYP1A2 at every panel size (non-inferiority 0 of 8).

**What is actually true.** Every probe on this array carrying the symbol `CYP1A2` is a **macaque** sequence. Verified by accession at NCBI Nucleotide on 2026-11-09:

| probe accession | organism | molecule |
|---|---|---|
| `AB185338` | *Macaca fuscata* | mRNA |
| `AF123056` | *Macaca fascicularis* | DNA |
| `AF123057` | *Macaca mulatta* | DNA |
| `D86474` | *Macaca fascicularis* | mRNA |

**There is no human CYP1A2 probe under that symbol on the platform.** The frozen one-feature-per-gene rule, which picks the probe with the highest training-fold mean, selected `AF123056` — a macaque *genomic* probe.

So the CYP1A2 "target transcript" arm never measured human CYP1A2 messenger RNA. Consequences:

1. **The claim that CYP1A2 mRNA is uninformative for CYP1A2 activity is withdrawn.** It was never tested. The measured equality with the baseline (492.6 against 492.4) reflects a non-human probe, not the biology of the human transcript.
2. **NetPK's 0-of-8 failure on CYP1A2 is withdrawn as evidence about mechanistic selection.** The ordering was seeded on a gene whose only available probes are from another species, so the arm was handicapped by the gene-to-probe mapping, not by the mechanism.
3. The statistical arms on CYP1A2 remain interpretable, because they chose their own probes. The endpoint is not void; the arms **seeded on CYP1A2** are.
4. Nothing can repair this by re-running: the platform does not carry the probe. The endpoint is reported with its mechanistic arms marked uninterpretable.

**This defect was mine, not the protocol's.** The label-resolution rule was designed to avoid hard-coding a gene list, and it does. It never checked that the probe behind the resolved symbol belongs to the organism under study, and nothing in the protocol required it to. That requirement is added below.

## E2 — CYP2E1 survives the same check

Three probes carry `CYP2E1`: two human (`NM_000773`) and one macaque (`AY635465`). The frozen highest-mean rule picks the **human** probe. The CYP2E1 result therefore stands as what I said it was: a genuine case in which the human transcript of an enzyme does not predict that enzyme's measured activity better than age and sex.

**One conclusion survives, one does not.** The admissible statement is: *on one of ten endpoints, the enzyme's own human transcript carried no more information than the clinical baseline.* Not two.

## E3 — The candidate universe contains probes no biological reading can use

Qualified by accession lookup, restricted to the probes some arm actually selected, so every verdict is auditable:

| verdict | distinct symbols selected |
|---|---:|
| human transcript | 374 |
| unverified (no retrievable record) | 74 |
| human genomic | 13 |
| **non-human** | **9** |
| **technical control** | **7** |
| mixed organism | 3 |

Inadmissible selections per arm, out of the symbols that arm chose:

| arm | inadmissible | which |
|---|---:|---|
| `univariate_spearman` | 10 of 152 | CYP2A23, CYP2A24, CYP2C74, CYP2C75, CYP2D17, CYP2D42, CYP3A64, CYP3A66 (all macaque), NegCont_210, NegCont_434 |
| `xgboost_gain` | 9 of 185 | CYP2A23, CYP2A24, CYP2D17, CYP3A64, CYP3A66, and four negative controls |
| `lasso_entry_order` | 5 of 182 | CYP2A24, CYP3A66, and three negative controls |
| every mechanistic arm | 1 of 20–37 | CYP1A2, by the defect of E1 |

**What this does and does not mean.** A non-human probe cannot carry independent human biology: on human liver RNA its signal is cross-hybridisation with the human orthologue, or noise. A technical control carries no transcript at all. But neither makes the arm that chose it *wrong about prediction* — averaging two cross-hybridising probes of one transcript is a legitimate way to reduce noise, and that is probably what the univariate arm was doing when it ranked CYP2D17 and CYP2D42 behind CYP2D6. What it means is that **those selections must not be reported as distinct biomarkers**, and that my earlier description of them as "a platform artefact" was too dismissive in one direction while assuming too much in the other.

**A human record is not dismissed because its symbol looks like a pseudogene.** CYP2A7P1 and CYP2B7P1 are human records and are left alone. Where no record was retrievable the verdict is `unverified`, not a guess — 74 of 480 selected symbols, which is itself a limit on how far any biological reading of these panels can go.

---

## What changes for the analyses still to run

**Required before H2d.** The gene-to-probe mapping must verify the organism of the probe behind a resolved symbol, and refuse a symbol whose probes are all non-human rather than silently selecting one. H2d must not inherit the technical controls or the cross-species probes.

**Three universes, kept separate and traced**, as an external review proposed:

1. **the frozen primary** — unchanged, including the withdrawn CYP1A2 arms, so the negative result stands as it was obtained;
2. **the common network pool** — every arm restricted to genes the public network contains, which separates ordering quality from pool size; already running when this erratum was written;
3. **the biologically qualified pool** — human probes only, technical controls removed, cross-hybridising paralogues flagged.

**The third must not be read as favouring either side.** Restricting to the network pool removes technical probes that only a statistical arm could reach, which helps NetPK; it also removes real human biomarkers absent from the network, which hurts it. Restricting to qualified human probes removes the cross-species probes a statistical arm used, which helps NetPK; it leaves the uncharacterised loci, which it can still use. Each restriction moves the comparison in both directions at once, and no single one of them is the honest answer on its own. All three are reported.

---

## What is unchanged

Every threshold, margin, comparator and fold of the frozen H2a protocol. The destroyed-signal diagnostic, which does not depend on the identity of any probe beyond the target gene, and whose finding — that an outcome-blind ordering keeps a variable whose information has been destroyed, on every endpoint, while a selector reading the endpoint was never once fooled in ten folds — is unaffected on the nine endpoints whose target probe is human.
