# H2a — erratum: one endpoint withdrawn, and a platform defect that changes the reading

**Erratum version:** 1.0.2
**Date:** 2026-11-09
**Applies to:** the frozen H2a run, protocol v1.4.3, merged result `reports/data/h2a_hlc_selection.json`
**Status:** no H2a number is recomputed. One endpoint's mechanistic arms are withdrawn as uninterpretable, and one conclusion I reported is retracted.

---

## E1 — Withdrawn: the CYP1A2 result cannot be read as a result about CYP1A2

**What I reported.** That CYP1A2 and CYP2E1 were cases where the target enzyme's own transcript carried no more information than age and sex, and that NetPK failed on CYP1A2 at every panel size (non-inferiority 0 of 8).

**What is actually true, stated as precisely as the evidence allows — and version 1.0.1 weakens what version 1.0.0 of this erratum claimed.** Two independent public annotations describe these four features, they **contradict each other**, and one of them is internally incoherent:

| feature | archive accession | organism of that accession | platform accession | platform gene assignment |
|---|---|---|---|---|
| `10033668809` | `AF123056` | *Macaca fascicularis* (DNA) | `DN745607` | dog `LOC479559` |
| `10033668817` | `AB185338` | *Macaca fuscata* (mRNA) | `CO597948` | dog `LOC479559` |
| `10033668906` | `AF123057` | *Macaca mulatta* (DNA) | `AM086385` | dog **`MYH9`, myosin heavy chain 9** |
| `10033668921` | `D86474` | *Macaca fascicularis* (mRNA) | — | dog `LOC479559` |

A cytochrome probe is not myosin, so the platform's gene assignment is wrong here; and the two sources cannot both be right about the accession. **Neither source distributes the oligonucleotide sequence**, so hybridisation specificity cannot be checked at all.

**Why version 1.0.0 was too strong.** It said the arm "never measured human CYP1A2 messenger RNA". That does not follow. A sixty-mer designed from a macaque cytochrome coding sequence, highly similar to its human orthologue, may well hybridise the human transcript, with unknown efficiency. The defensible claim is narrower and is the one that now stands: **the molecular identity of what these probes measure cannot be established from any available public source.** Not established as human, not established as non-human — unestablished.

That is a different reason from the one first given, and it is a sufficient one. The consequences are unchanged in substance:

1. **The claim that CYP1A2 mRNA is uninformative for CYP1A2 activity is withdrawn.** It was never tested: the measured equality with the baseline (492.6 against 492.4) came from a probe of unestablished identity.
2. **NetPK's 0-of-8 failure on CYP1A2 is withdrawn as evidence about mechanistic selection.** The ordering was seeded on a symbol behind which no probe has an establishable human identity, so the arm was handicapped by the gene-to-probe mapping rather than by the mechanism.
3. **The numbers are kept, not deleted.** CYP1A2's full result stays in the merged output exactly as obtained, marked uninterpretable *for this comparison*. The statistical arms on that endpoint remain interpretable, because they chose their own probes, and the enzyme activity itself remains a usable **secondary endpoint** for studying mechanisms other than its own transcript.
4. Nothing can repair this by re-running: the platform carries no probe for this symbol whose identity can be established.

**This defect was mine, not the protocol's.** The label-resolution rule was designed to avoid hard-coding a gene list, and it does. It never checked that the probe behind the resolved symbol belongs to the organism under study, and nothing in the protocol required it to. That requirement is added below.

## E2 — CYP2E1 survives the same check

Three probes carry `CYP2E1`: two human (`NM_000773`) and one macaque (`AY635465`). The frozen highest-mean rule picks the **human** probe. The CYP2E1 result therefore stands as what I said it was: a genuine case in which the human transcript of an enzyme does not predict that enzyme's measured activity better than age and sex.

**One conclusion survives, one does not.** The admissible statement is: *on one of ten endpoints, the enzyme's own human transcript carried no more information than the clinical baseline.* Not two.

## E3 — The candidate universe, qualified over every probe rather than every selected symbol

**Version 1.0.0 qualified only the symbols the first fold recorded as selected, and judged a symbol by pooling the accessions of all probes behind it.** Both limits were real: the ten folds select far more than the first fold records, and the analysis picks one *particular* probe per symbol per fold, so a verdict on a symbol certifies nothing about the probe a fold used. The qualification is redone **per probe, over all 40,638 features, before any selection is consulted**, tracing `feature → accession → organism → symbol` for each.

**Version 1.0.1's single `identified_human` verdict was itself too affirmative, and its conflict count was inflated.** Two further corrections, both from the same review:

* the verdict covered two unequal situations — two annotations agreeing on a human gene, and **one** annotation with nothing to agree with it. A lone source is not a concordance, and counting both as identified turned *annotation-qualified candidates* into a claim about *demonstrated molecular identities*;
* **two different accessions are not a conflict when they designate transcripts of the same gene.** Comparing accession strings manufactured conflicts out of an archive citing one transcript and a platform citing another. Accessions are now resolved to genes *before* anything is compared, using NCBI's own `gene2accession`, streamed over 291,554,102 rows and filtered to *Homo sapiens*, giving 2,107,223 accession-to-gene pairs in a local file with a recorded checksum.

The verdicts are graded, strongest first, and the difference the second correction makes is large:

| verdict | probes | what it establishes |
|---|---:|---|
| `concordant_human` | **23,261** | two annotations resolve to the same human gene |
| `single_source_human` | 4,215 | one source names a human gene; nothing agrees or disagrees |
| `conflicting` | **1,476** | the sources name different genes; identity unestablished |
| `undetermined` | 10,679 | no source resolves the feature to a human gene |
| `technical_control` | 1,007 | the array's own control feature |

**The conflict count falls from about 7,000 to 1,476** once accessions are resolved rather than compared as strings. Version 1.0.1's figure was five times too high, and the claim that "44 % of candidate symbols have no probe whose human identity can be established" is withdrawn with it: **22,237 of 31,862 symbols carry at least one admissible probe, so about 30 % do not.** The direction of the finding is unchanged; its size was overstated.

An example of the correction working, from the data: the CYP3A4 probe whose archive accession is `D00003` and whose platform accession is `NM_000776` was a "conflict" under the string rule. Both resolve to human gene 1576. It is concordant.

**The nine target genes, checked individually under the graded rule.** Eight have admissible probes. `CYP1A2` has **none**: all four of its features are `conflicting`, because neither accession resolves to a human gene and the platform assigns a non-human one. The withdrawal of E1 therefore stands under the stricter and fairer rule, not only under the first one.

**Per-fold traceability, which a symbol-level verdict cannot give.** For both mixed symbols the analysis actually used — `CYP2E1`, the target of one endpoint, and `CYP1A1`, which appears in the mechanistic expansion for another — the frozen highest-training-mean rule selected an `identified_human` probe in **10 of 10 folds**. The CYP2E1 result therefore rests on a human-identified probe in every fold, not merely on average.

**The earlier per-selection table is superseded; its reading stands.** For the record, among the symbols the first fold's selections reached:

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

**A human record is not dismissed because its symbol looks like a pseudogene.** Under the full qualification CYP2A7P1, CYP2B7P1, POR and FABP1 all come out `all_probes_identified_human` and are left alone. The symbols that fail are the ones that fail on evidence: CYP2D17, CYP2D42, CYP2A23, CYP2A24, CYP2C74, CYP2C75, CYP3A64 and CYP3A66 all have `no_identified_human_probe`. Where no record was retrievable the verdict is `unverified`, not a guess — 74 of 480 selected symbols, which is itself a limit on how far any biological reading of these panels can go.

---

## What changes for the analyses still to run

**Required before H2d.** The gene-to-probe mapping must verify the organism of the probe behind a resolved symbol, and refuse a symbol whose probes are all non-human rather than silently selecting one. H2d must not inherit the technical controls or the cross-species probes.

**Three universes, kept separate and traced**, as an external review proposed:

1. **the frozen primary** — unchanged, including the withdrawn CYP1A2 arms, so the negative result stands as it was obtained;
2. **the common network pool** — every arm restricted to genes the public network contains, which separates ordering quality from pool size; already running when this erratum was written;
3. **the biologically qualified pool** — human probes only, technical controls removed, cross-hybridising paralogues flagged.

**The third must not be read as favouring either side.** Restricting to the network pool removes technical probes that only a statistical arm could reach, which helps NetPK; it also removes real human biomarkers absent from the network, which hurts it. Restricting to qualified human probes removes the cross-species probes a statistical arm used, which helps NetPK; it leaves the uncharacterised loci, which it can still use. Each restriction moves the comparison in both directions at once, and no single one of them is the honest answer on its own. All three are reported.

## E3.1 — The pools H2d will use, built by probe and then grouped

Built at the **probe** level and only then grouped by the symbol the analysis
uses, so a gene with one admissible and one inadmissible probe keeps the
admissible one. Dropping such a gene whole would discard a real measurement
because a *different* probe of the same gene is badly annotated — `CYP2E1` is
exactly that case and keeps its two concordant probes.

| pool | admits | symbols | probes | symbols keeping only some probes |
|---|---|---:|---:|---:|
| **A** qualified human | concordant + single source | 22,237 | 27,245 | 2,075 |
| **A strict** | concordant only | 19,938 | 23,107 | 3,058 |
| **B** = A ∩ network | concordant + single source | 11,092 | 14,647 | 1,383 |

Every arm gets the **same** pool, the probe retained in each fold is traceable,
and the absence of technical controls is asserted rather than assumed: the
runner refuses to proceed if a control symbol survives into a pool.

**Provenance is frozen, not summarised.** The committed file carries the rule
version, the SHA-256 of the platform record, of the human gene reference and of
the accession-to-gene map, and the **complete list of admissible feature
identifiers** for each pool — not only counts. Every individual decision can be
reconstructed from it.

**Two sensitivities, with the undetermined stratum handled separately.** Pool A
answers whether the comparison holds among analytically admissible variables;
pool B answers it among variables both admissible *and* reachable by the
mechanistic ordering. Neither is the honest answer alone, for the reason already
recorded: each restriction moves the comparison in both directions at once.

## E4 — Status of this qualification, and the three levels it separates

This is a **secondary analysis, posterior to the H2a run**. It does not requalify the frozen H2a result as confirmatory, changes no threshold, and introduces no selection rule chosen in the light of an inspected endpoint. It establishes which candidate variables can carry a biological reading at all, which is a question about the platform and not about any hypothesis.

Three levels have to be kept apart, and the H2a experience is what separates them:

| level | question | what went wrong here |
|---|---|---|
| analytical quality | is the announced biomarker the molecule actually measured? | about 30 % of candidate symbols carry no probe whose human identity can be established |
| biological validity | does that molecule represent the supposed mechanism? | messenger RNA abundance need not represent enzyme activity |
| predictive usefulness | does it add information about the phenotype? | the only question H2a was built to answer |

A variable can pass the third and fail the first, which is exactly what the contested cytochrome probes did. A selection method cannot be judged on predictive usefulness alone when the first level is unsettled.

**And the converse, which matters as much.** A probe whose identity is undetermined is **not** a probe whose measurement is wrong. It is a measurement whose biological meaning the available annotation does not pin down. A statistically informative variable can remain biologically undetermined, and the two questions are reported separately rather than merged into one verdict.

---

## What is unchanged

Every threshold, margin, comparator and fold of the frozen H2a protocol. The destroyed-signal diagnostic, which does not depend on the identity of any probe beyond the target gene, and whose finding — that an outcome-blind ordering keeps a variable whose information has been destroyed, on every endpoint, while a selector reading the endpoint was never once fooled in ten folds — is unaffected on the nine endpoints whose target probe is human.
