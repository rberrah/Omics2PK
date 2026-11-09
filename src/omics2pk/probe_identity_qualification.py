"""What is established about a probe's molecular identity, and how strongly.

Why the evidence is graded rather than binary
---------------------------------------------
An earlier version of this qualification returned `identified_human` for two
situations that are not equally supported, and an external review was right to
reject the conflation:

* two independent annotations pointing at the same human gene;
* **one** annotation pointing at a human gene, the other silent.

The second is a single source, not a concordance. Counting both as identified
turned 22,348 *annotation-qualified candidates* into a claim about 22,348
*molecular identities*, which the evidence does not support. The verdicts below
are therefore ordered by the strength of what is actually shown.

A second correction of the same kind: **two different accessions are not a
conflict if they designate transcripts of the same gene.** The old rule compared
accession strings and called any difference a contradiction, which manufactured
conflicts out of an archive citing one transcript and a platform citing another.
Accessions are resolved to genes before anything is compared.

Two further corrections, from a review of the graded version
-----------------------------------------------------------
* **A verdict must not attribute a measurement to the wrong gene.** The graded
  version returned `concordant_human` with `symbol_matches_gene = False`, and
  the pool then grouped that feature under the symbol the analysis uses. Where
  that symbol is the *approved symbol or a listed synonym of a different human
  gene*, the two claims are mutually exclusive and the feature is refused
  (`symbol_names_another_gene`). Where the symbol is merely unknown to the
  current gene reference -- a clone or locus placeholder such as an `FLJ`,
  `MGC`, `DKFZP` or `LOC` name, or a symbol a spreadsheet mangled into a date
  -- it carries no competing claim, the accessions decide, and the feature
  stays admissible with the gene's current symbol recorded beside it. Measured
  on this array: 632 contradictions against 5,268 stale aliases, so refusing
  both would have discarded five thousand real measurements for a
  naming-vintage reason.
* **The platform record can contradict itself.** Its accession and its
  `EntrezGeneID` sometimes name different genes. The old order of tests let the
  archive agree with one of the two and called the result a concordance between
  sources. Internal coherence is now checked first, and a self-contradictory
  platform record cannot confirm anything (`platform_self_contradictory`): 610
  features on this array, 607 of which the previous rule admitted.

What this module cannot establish, whatever the verdict
-------------------------------------------------------
Neither annotation source distributes the oligonucleotide sequence, so
hybridisation specificity is never checked. A probe designed from a non-human
coding sequence may hybridise the human orthologue with unknown efficiency, and
a probe with impeccable annotation may still cross-hybridise. Every verdict here
is about *annotation*, and the word "identified" always means "identified by
annotation".

And the corollary that matters for the science: a probe whose identity is
undetermined is **not** a probe whose measurement is wrong. It is a measurement
whose biological meaning the available annotation does not pin down. A
statistically informative variable can remain biologically undetermined, and the
two questions are reported separately rather than merged into one verdict.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Mapping, Sequence

CONTROL_PATTERN = re.compile(
    r"^(NegCont|PosCont|Corner|DarkCorner|SpikeIn|GE_Bright|E1A_r60|ETG)",
    re.IGNORECASE)

# Ordered from strongest to weakest. A pool that admits a level admits every
# stronger one.
VERDICTS = (
    "concordant_human",      # two sources, same human gene
    "single_source_human",   # one source names a human gene, the other is silent
    "symbol_names_another_gene",    # the accessions resolve to one gene while
                                    # the analysis's symbol is the approved
                                    # symbol or a listed synonym of another
    "platform_self_contradictory",  # the platform's accession and its gene
                                    # assignment name different genes, so the
                                    # platform cannot confirm anything
    "conflicting",           # the sources name different genes
    "undetermined",          # no source yields a gene
    "technical_control",     # the array's own control feature
)

ADMISSIBLE_BY_DEFAULT = ("concordant_human", "single_source_human")


class ProbeIdentityError(ValueError):
    pass


@dataclass(frozen=True)
class ProbeIdentity:
    """What is known about one feature, with the chain that established it."""

    feature_id: str
    symbol: str
    verdict: str
    reason: str
    archive_accession: str
    archive_genes: tuple[str, ...]
    platform_accession: str
    platform_accession_genes: tuple[str, ...]
    platform_entrez: str
    human_gene: str | None
    human_symbol: str | None
    symbol_matches_gene: bool | None

    @property
    def admissible(self) -> bool:
        return self.verdict in ADMISSIBLE_BY_DEFAULT

    @property
    def symbol_is_a_stale_alias(self) -> bool:
        """Admissible, but the analysis names it something the gene no longer is.

        The measurement stands; what does not stand is the assumption that the
        analysis's column name is the gene's current symbol. A mechanistic
        ordering keyed on current symbols cannot reach these features, which is
        a limit of the array's naming rather than of the measurement.
        """
        return self.admissible and self.symbol_matches_gene is False


def _base(accession: str) -> str:
    return accession.split(".")[0].strip()


def symbol_to_gene_index(human_gene_symbols: Mapping[str, tuple[str, frozenset]],
                         ) -> dict[str, tuple[str, ...]]:
    """Which human genes carry each approved symbol or listed synonym.

    Built from the same reference the verdicts use, so a symbol cleared here is
    cleared against exactly the gene set the accessions are resolved against.
    Synonyms are included: a symbol listed as another gene's former name is a
    competing claim, not a free pass.
    """
    index: dict[str, list[str]] = {}
    for gene, (approved, synonyms) in human_gene_symbols.items():
        for name in (approved, *synonyms):
            if not name:
                continue
            carriers = index.setdefault(name, [])
            if gene not in carriers:
                carriers.append(gene)
    return {name: tuple(carriers) for name, carriers in index.items()}


def qualify_probe(feature_id: str, *, symbol: str, archive_accession: str,
                  platform_accession: str, platform_entrez: str,
                  accession_to_human_genes: Mapping[str, Sequence[str]],
                  human_gene_symbols: Mapping[str, tuple[str, frozenset]],
                  symbol_to_human_genes: Mapping[str, Sequence[str]] | None = None,
                  ) -> ProbeIdentity:
    """Grade one probe. Every argument is annotation; none is measurement.

    `symbol_to_human_genes` maps a current approved symbol or a listed synonym
    to the human genes carrying it. It is what separates a symbol that names
    another gene from a symbol the reference no longer knows at all. Without it
    the stricter branch applies and every mismatch is refused, because a rule
    that cannot tell the two apart must not assume the harmless case.
    """
    symbol = (symbol or "").strip()
    archive_accession = _base(archive_accession or "")
    platform_accession = _base(platform_accession or "")
    platform_entrez = (platform_entrez or "").strip()

    def blank(verdict, reason, **extra):
        fields = {
            "feature_id": feature_id, "symbol": symbol, "verdict": verdict,
            "reason": reason, "archive_accession": archive_accession,
            "archive_genes": (), "platform_accession": platform_accession,
            "platform_accession_genes": (), "platform_entrez": platform_entrez,
            "human_gene": None, "human_symbol": None,
            "symbol_matches_gene": None,
        }
        fields.update(extra)
        return ProbeIdentity(**fields)

    if symbol and CONTROL_PATTERN.match(symbol):
        return blank("technical_control",
                     "the array's own control feature; it carries no transcript")

    archive_genes = tuple(accession_to_human_genes.get(archive_accession, ()))
    platform_genes = tuple(accession_to_human_genes.get(platform_accession, ()))

    # Internal coherence of the platform record comes FIRST. Its accession and
    # its gene assignment sometimes name different genes, and a record that
    # contradicts itself cannot confirm what the archive says: letting the
    # archive agree with one of its two fields would report a concordance
    # between sources that does not exist.
    if platform_genes and platform_entrez and platform_entrez not in platform_genes:
        return blank(
            "platform_self_contradictory",
            (f"the platform's accession {platform_accession} resolves to human "
             f"gene(s) {sorted(platform_genes)} while the same record assigns "
             f"gene {platform_entrez}; a record that contradicts itself cannot "
             "confirm the archive"),
            archive_genes=archive_genes, platform_accession_genes=platform_genes)

    # The platform's own gene assignment counts as evidence only when the gene
    # it names is a human gene; otherwise it is evidence of a conflict.
    platform_gene_is_human = platform_entrez in human_gene_symbols
    platform_claims = set(platform_genes)
    if platform_gene_is_human:
        platform_claims.add(platform_entrez)
    archive_claims = set(archive_genes)

    def describe(gene):
        human_symbol = human_gene_symbols[gene][0]
        matches = (symbol == human_symbol
                   or symbol in human_gene_symbols[gene][1]) if symbol else None
        return human_symbol, matches

    def claims_another_gene(gene):
        """Which other human genes, if any, does the analysis's symbol name?"""
        if symbol_to_human_genes is None:
            return ()
        return tuple(sorted(set(symbol_to_human_genes.get(symbol, ())) - {gene}))

    def refuse_wrong_gene(gene, human_symbol, others):
        detail = (f"the symbol is carried by human gene(s) {list(others)}"
                  if others else
                  "no symbol-to-gene reference was supplied, so the symbol "
                  "cannot be cleared")
        return blank(
            "symbol_names_another_gene",
            (f"the accessions resolve to human gene {gene} ({human_symbol}) but "
             f"{detail}; grouping this measurement under {symbol!r} would "
             "attribute it to the wrong gene"),
            archive_genes=archive_genes, platform_accession_genes=platform_genes)

    shared = sorted(archive_claims & platform_claims)
    if shared:
        gene = shared[0]
        human_symbol, matches = describe(gene)
        if matches is False:
            others = claims_another_gene(gene)
            if others or symbol_to_human_genes is None:
                return refuse_wrong_gene(gene, human_symbol, others)
        return ProbeIdentity(
            feature_id=feature_id, symbol=symbol, verdict="concordant_human",
            reason=(f"both annotations resolve to human gene {gene} "
                    f"({human_symbol})"),
            archive_accession=archive_accession, archive_genes=archive_genes,
            platform_accession=platform_accession,
            platform_accession_genes=platform_genes,
            platform_entrez=platform_entrez, human_gene=gene,
            human_symbol=human_symbol, symbol_matches_gene=matches)

    if archive_claims and platform_claims:
        return blank(
            "conflicting",
            (f"the archive resolves to human gene(s) {sorted(archive_claims)} "
             f"while the platform points at {sorted(platform_claims)}; the "
             "identity is unestablished"),
            archive_genes=archive_genes, platform_accession_genes=platform_genes)

    if platform_entrez and not platform_gene_is_human and archive_claims:
        return blank(
            "conflicting",
            (f"the archive resolves to human gene(s) {sorted(archive_claims)} "
             f"while the platform assigns gene {platform_entrez}, which is not "
             "a human gene"),
            archive_genes=archive_genes, platform_accession_genes=platform_genes)

    only = archive_claims or platform_claims
    if only:
        gene = sorted(only)[0]
        human_symbol, matches = describe(gene)
        if matches is False:
            others = claims_another_gene(gene)
            if others or symbol_to_human_genes is None:
                return refuse_wrong_gene(gene, human_symbol, others)
        which = "the archive" if archive_claims else "the platform"
        return ProbeIdentity(
            feature_id=feature_id, symbol=symbol,
            verdict="single_source_human",
            reason=(f"only {which} yields an identity, human gene {gene} "
                    f"({human_symbol}); there is no second annotation to agree "
                    "or disagree with it"),
            archive_accession=archive_accession, archive_genes=archive_genes,
            platform_accession=platform_accession,
            platform_accession_genes=platform_genes,
            platform_entrez=platform_entrez, human_gene=gene,
            human_symbol=human_symbol, symbol_matches_gene=matches)

    if platform_entrez and not platform_gene_is_human:
        return blank(
            "conflicting",
            (f"the platform assigns gene {platform_entrez}, which is not a "
             "human gene, and nothing else resolves this feature"))

    return blank("undetermined",
                 "no annotation resolves this feature to a human gene")


def build_probe_pool(identities: Sequence[ProbeIdentity], *,
                     admit: Sequence[str] = ADMISSIBLE_BY_DEFAULT,
                     restrict_to_genes: Sequence[str] | None = None,
                     ) -> dict:
    """Group admissible probes by the symbol the analysis uses.

    The pool is built at the **probe** level and only then grouped, so a symbol
    with some admissible and some inadmissible probes keeps the admissible ones
    instead of being dropped whole. Excluding such a gene would discard a real
    measurement because a *different* probe of the same gene is badly annotated.
    """
    admit = tuple(admit)
    for verdict in admit:
        if verdict not in VERDICTS:
            raise ProbeIdentityError(f"unknown verdict: {verdict}")
    if "technical_control" in admit:
        raise ProbeIdentityError(
            "a technical control carries no transcript and is never admissible")

    allowed_genes = set(restrict_to_genes) if restrict_to_genes is not None else None
    by_symbol: dict[str, list[str]] = {}
    dropped = {verdict: 0 for verdict in VERDICTS}
    partial = []
    for identity in identities:
        if identity.verdict not in admit:
            dropped[identity.verdict] += 1
            continue
        if not identity.symbol:
            dropped["undetermined"] += 1
            continue
        if allowed_genes is not None and identity.symbol not in allowed_genes:
            continue
        by_symbol.setdefault(identity.symbol, []).append(identity.feature_id)

    all_by_symbol: dict[str, int] = {}
    for identity in identities:
        if identity.symbol:
            all_by_symbol[identity.symbol] = all_by_symbol.get(identity.symbol, 0) + 1
    for symbol, members in by_symbol.items():
        if len(members) < all_by_symbol.get(symbol, 0):
            partial.append(symbol)

    return {
        "admitted_verdicts": list(admit),
        "restricted_to_gene_set": allowed_genes is not None,
        "symbols": len(by_symbol),
        "probes": sum(len(v) for v in by_symbol.values()),
        "by_symbol": by_symbol,
        "symbols_keeping_only_some_probes": sorted(partial),
        "probes_dropped_by_verdict": dropped,
        "note": ("Built probe by probe, then grouped. A symbol with both "
                 "admissible and inadmissible probes keeps the admissible "
                 "ones; dropping the symbol whole would discard a real "
                 "measurement because another probe of the same gene is badly "
                 "annotated."),
    }


def summarise(identities: Sequence[ProbeIdentity]) -> dict:
    counts = {verdict: 0 for verdict in VERDICTS}
    symbol_mismatch = 0
    for identity in identities:
        counts[identity.verdict] = counts.get(identity.verdict, 0) + 1
        if identity.symbol_matches_gene is False:
            symbol_mismatch += 1
    return {
        "probes": len(identities),
        "counts_by_verdict": counts,
        "admissible_by_default": sum(
            counts[v] for v in ADMISSIBLE_BY_DEFAULT),
        "symbol_disagrees_with_the_resolved_gene": symbol_mismatch,
        "admissible_under_a_stale_alias": sum(
            1 for identity in identities if identity.symbol_is_a_stale_alias),
        "what_a_stale_alias_means": (
            "The analysis's column name is not the gene's current symbol, and "
            "no other human gene carries that name. The measurement stands and "
            "the gene is recorded beside it, but a mechanistic ordering keyed "
            "on current symbols cannot reach the feature under that name. That "
            "is a limit of the array's naming, not a defect of the "
            "measurement, and it is not repaired by renaming anything here: "
            "the analysis universe is keyed on the array's own symbols."),
        "what_identified_does_not_mean": (
            "No verdict checks hybridisation specificity, because neither "
            "source distributes the oligonucleotide sequence. 'Identified' "
            "means identified by annotation. An undetermined probe is not a "
            "wrong measurement; it is a measurement whose biological meaning "
            "the available annotation does not pin down."),
    }
