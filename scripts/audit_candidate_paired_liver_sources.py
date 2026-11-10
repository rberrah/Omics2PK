"""Qualify candidate paired liver sources against a rule declared before use.

Seven candidates were put forward on 2026-11-10 — three by the commissioning
note and four found while checking them. This script applies one rule to all
seven and records why each one passes or fails, so that "we looked and there is
nothing" is a verdict with evidence behind it rather than an impression.

Two kinds of fact are handled differently, on purpose:

* **Access and licence** are re-queried from Europe PMC at every run, because
  they can change, and the run **fails** if what comes back disagrees with what
  was recorded. An access claim that is not re-checked is an access claim that
  has quietly expired.
* **Whether per-donor values exist** cannot be read from an index: it needs the
  supplement opened and looked at. That determination is therefore recorded
  here as evidence — which file was retrieved, its SHA-256, and what was in it
  — and the script applies the rule to it rather than re-deriving it.

The script fits nothing, predicts nothing and downloads no supplement. Run it
with ``python -I``.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

EUROPE_PMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

#: Donor floor, declared. Leave-one-donor-out with a one-column baseline and an
#: inner leave-one-out for the penalty leaves n-2 donors in the inner training
#: set; below ten donors that is under eight, and neither the penalty choice nor
#: the interval is interpretable. A source under the floor can still be a
#: published comparator.
DONOR_FLOOR = 10

#: The outcome kinds the project's readiness gate accepts, plus the Axis A
#: escape hatch for a declared non-pharmacological target.
ADMISSIBLE_OUTCOMES = (
    "drug_specific_clint",
    "drug_specific_metabolic_activity",
    "clinical_pk",
    "declared_non_pharmacological_target",
)

REASONS = (
    "ABSTAIN_ACCESS",
    "ABSTAIN_NO_PER_DONOR_VALUES",
    "ABSTAIN_TOO_FEW_DONORS",
    "ABSTAIN_NO_ADMISSIBLE_ENDPOINT",
)

#: Every candidate, with the evidence gathered on 2026-11-10. ``expected``
#: holds the Europe PMC fields as read that day; the run fails if they moved.
CANDIDATES = {
    "drozdzik_2018_liver_intestine_9": {
        "label": "Drozdzik 2018, paired liver and intestine, 9 donors",
        "citation": (
            "Drozdzik M, Busch D, Lapczuk J, Muller J, Ostrowski M, Kurzawski "
            "M, Oswald S. Clin Pharmacol Ther 2018;104(3):515-524."
        ),
        "query": "EXT_ID:29205295",
        "expected": {"isOpenAccess": "N", "inPMC": "N", "hasSuppl": "N"},
        "n_donors": 9,
        "layers": ("mRNA (qRT-PCR)", "protein (LC-MS/MS)"),
        "same_donor_layers": True,
        "outcome_kind": None,
        "per_donor_layer_and_outcome": None,
        "evidence": (
            "Europe PMC records isOpenAccess N, inPMC N and hasSuppl N, so the "
            "supplement is not in any open index. The publisher's pages answer "
            "HTTP 403 to every automated client (Cloudflare interstitial "
            "'Just a moment...'), both for the article "
            "https://ascpt.onlinelibrary.wiley.com/doi/10.1002/cpt.967 and for "
            "the downloadSupplement endpoint. A web summary asserted that "
            "Table S1 and S2 hold per-donor protein and mRNA for the nine "
            "donors; that assertion could NOT be checked against the file and "
            "is therefore not treated as a fact."
        ),
        "reason": "ABSTAIN_ACCESS",
        "usable_as": "citation_only",
        "why_it_matters": (
            "This is the most conceptually valuable candidate of the seven: "
            "mRNA AND protein for CYPs and UGTs in the SAME nine donors would "
            "make the project's RNA-versus-protein comparison within-cohort "
            "instead of cross-cohort. It fails on access alone."
        ),
    },
    "wenzel_2024_jejunum_liver_11": {
        "label": "Wenzel 2024, intra-subject jejunum and liver, 11 donors",
        "citation": (
            "Wenzel C, Lapczuk-Romanska J, Malinowski D, Ostrowski M, "
            "Drozdzik M, Oswald S. Clin Pharmacol Ther 2024;115(2):221-230."
        ),
        "query": "EXT_ID:37739780",
        "expected": {"isOpenAccess": "N", "hasSuppl": "N"},
        "n_donors": 11,
        "layers": ("mRNA (qRT-PCR)", "protein (targeted proteomics)"),
        "same_donor_layers": True,
        "outcome_kind": None,
        "per_donor_layer_and_outcome": None,
        "evidence": (
            "Europe PMC records isOpenAccess N and hasSuppl N; the licence is "
            "cc by-nc-nd, which would forbid redistributing a derivative even "
            "if the file were reachable. Same publisher, same Cloudflare 403."
        ),
        "reason": "ABSTAIN_ACCESS",
        "usable_as": "citation_only",
        "why_it_matters": (
            "The same two layers in the same subjects as Drozdzik 2018, with "
            "eleven donors instead of nine, and the same access barrier."
        ),
    },
    "almajdoub_2026_kcat_matched_4": {
        "label": "Al-Majdoub 2026, kcat in matched liver and intestine, 4 donors",
        "citation": (
            "Al-Majdoub ZM, et al. Drug Metab Dispos 2026;54(4). "
            "PMID 41855594, PMCID PMC13197946."
        ),
        "query": "PMCID:PMC13197946",
        "expected": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "n_donors": 4,
        "layers": (
            "protein abundance (POOLED commercial lots)",
            "kcat / specificity constant (per donor)",
        ),
        "same_donor_layers": False,
        "outcome_kind": "drug_specific_metabolic_activity",
        "per_donor_layer_and_outcome": False,
        "evidence": (
            "Open access, CC BY, supplement retrieved and opened: mmc1.docx "
            "(934,810 bytes, SHA-256 prefix 80b321c6f197214c9c741464). "
            "Supplemental Tables 8 and 9 DO give intrinsic clearance and Vmax "
            "'from matched individual tissue donors', so the OUTCOME is per "
            "donor. But every abundance table -- Supplemental Tables 11, 12, "
            "13 and 14, for liver, intestine and kidney -- reports 'pooled "
            "donor' microsomes identified by LOT NUMBER. The "
            "abundance-to-activity pairing per donor therefore does not exist "
            "in this source at all, independent of its size. FOUR paired "
            "donors, and the article's own main comparison excludes one for "
            "low intestinal activity, leaving THREE. The commissioning note "
            "described this as 'only three donors for some individual "
            "comparisons'; the correct reading is that three to four is the "
            "whole cohort, not a subset."
        ),
        "secondary_disqualification": (
            "ABSTAIN_TOO_FEW_DONORS: four donors, three after the article's "
            "own exclusion, against a declared floor of ten. Either "
            "disqualification alone is sufficient."
        ),
        "reason": "ABSTAIN_NO_PER_DONOR_VALUES",
        "usable_as": "published_comparator",
        "why_it_matters": (
            "Scientifically the closest paper to this project's own question — "
            "it measures activity PER UNIT of enzyme in matched tissues and "
            "concludes that scaling activity through abundance needs care. It "
            "is a citation for the argument and cannot be a dataset: no "
            "contrast, interval or cross-validation is estimable at n = 3."
        ),
    },
    "zhang_2015_mppgl_128": {
        "label": "Zhang 2015, MPPGL and ten CYP activities, 128 livers",
        "citation": (
            "Zhang H, Gao N, Tian X, et al. Sci Rep 2015;5:17671. "
            "PMCID PMC4669488. CC BY."
        ),
        "query": "EXT_ID:26635233",
        "expected": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "n_donors": 128,
        "layers": ("MPPGL", "ten CYP activities in microsomes AND in tissue"),
        "same_donor_layers": True,
        "outcome_kind": "drug_specific_metabolic_activity",
        "per_donor_layer_and_outcome": False,
        "evidence": (
            "Supplement retrieved from the Springer CDN "
            "(41598_2015_BFsrep17671_MOESM1_ESM.pdf, 1,004,095 bytes, SHA-256 "
            "23c1ceba3e8e142df50e449cb02b42f506efc2737945b0228553c211efb5bf7b) "
            "and read in full: 41 lines of text. It holds ONE table of "
            "aggregates — ranges, means, SD and 95 % prediction intervals for "
            "five scaling parameters on the n = 78 prediction subset — and "
            "three figure captions. NO per-donor value for MPPGL or for any "
            "CYP activity. Recovering the individual points would mean "
            "digitising scatter figures, which this project forbids."
        ),
        "reason": "ABSTAIN_NO_PER_DONOR_VALUES",
        "usable_as": "published_comparator",
        "why_it_matters": (
            "At 128 livers it is by far the largest candidate, and its "
            "published direction CONCORDS with this project's own Axis B "
            "finding: MPPGL correlates strongly with CYP activity measured in "
            "liver TISSUE (r >= 0.6, P < 0.001 for most CYPs) and not at all "
            "with activity measured in MICROSOMES (no significant "
            "correlation). Concordance in direction between two cohorts is "
            "worth recording; it is not a replication, because no per-donor "
            "value is available to recompute anything."
        ),
    },
    "drozdzik_2023_hcv_58": {
        "label": "Drozdzik 2023, hepatitis C livers, 58 samples",
        "citation": (
            "Drozdzik M, Lapczuk-Romanska J, Wenzel C, et al. Int J Mol Sci "
            "2023;24(5):4543. PMCID PMC10002520. CC BY."
        ),
        "query": "DOI:\"10.3390/ijms24054543\"",
        "expected": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "n_donors": 58,
        "layers": ("mRNA (qRT-PCR)", "protein (LC-MS/MS)"),
        "same_donor_layers": True,
        "outcome_kind": None,
        "per_donor_layer_and_outcome": False,
        "evidence": (
            "Supplement retrieved and opened: the Europe PMC archive carries "
            "ijms-24-04543-s001.zip (133,511 bytes, SHA-256 prefix "
            "537c94b7d60348d9ca689dce), which contains exactly one file, "
            "ijms-2146685-supplementary.pdf (141,723 bytes). Supplemental "
            "Table S1 reports, per enzyme and per Child-Pugh group, the "
            "median, min, max, mean, SD and CV % — GROUP-LEVEL aggregates "
            "only. No per-donor row for either layer."
        ),
        "reason": "ABSTAIN_NO_PER_DONOR_VALUES",
        "usable_as": "published_comparator",
        "why_it_matters": (
            "CC BY, 58 same-donor livers with both layers measured: on paper "
            "the candidate that would have closed the project's biggest gap. "
            "It fails because the individual values are not published, which "
            "is a sharing decision rather than a measurement limit."
        ),
    },
    "drozdzik_2021_pathology_97": {
        "label": "Drozdzik 2021, liver pathologies, about 97 livers",
        "citation": (
            "Drozdzik M, Lapczuk-Romanska J, Wenzel C, et al. Pharmaceutics "
            "2021;13(9):1334. PMCID PMC8471929. CC BY."
        ),
        "query": "EXT_ID:34575411",
        "expected": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "n_donors": 97,
        "layers": ("mRNA (qRT-PCR)", "protein (LC-MS/MS)"),
        "same_donor_layers": True,
        "outcome_kind": None,
        "per_donor_layer_and_outcome": False,
        "evidence": (
            "Read from the Europe PMC full text. The paper never states a "
            "total in one sentence; the group sizes sum to 97 (controls 20, "
            "HCV 21, ALD 20, AIH 20, PBC 10, PSC 6). Its **Data Availability "
            "Statement is 'Not applicable'**, and the supplementary tables it "
            "names are subject characteristics and BETWEEN-GROUP differences "
            "in abundance and mRNA — group-level by description. Both layers "
            "come from separate aliquots of the same harvest, which the "
            "Methods imply rather than state."
        ),
        "reason": "ABSTAIN_NO_PER_DONOR_VALUES",
        "usable_as": "published_comparator",
        "why_it_matters": (
            "It publishes per-enzyme mRNA-to-protein correlations in liver at "
            "about 97 samples, which is comparator content for the project's "
            "layer question — the same role Grangeon plays for "
            "protein-to-activity."
        ),
    },
    "song_2018_hlm_16": {
        "label": "Song 2018, individual liver microsomes, 16 donors",
        "citation": (
            "Song W, Yu L, Peng Z. Data Brief 2018;18:831-834. "
            "PMCID PMC5996730. CC BY."
        ),
        "query": "PMCID:PMC5996730",
        "expected": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "n_donors": 16,
        "layers": ("microsomal proteomics", "CYP450, MEH and MGST activities"),
        "same_donor_layers": True,
        "outcome_kind": "drug_specific_metabolic_activity",
        "per_donor_layer_and_outcome": False,
        "evidence": (
            "Read from the Europe PMC full text. Table 1 is headed 'Individual "
            "human liver microsomes (n = 16)' but reports only the maximum, "
            "the minimum, the max/min ratio and group means by age and sex — "
            "no per-donor row. The supplementary files are Mascot-processed "
            "MS identifications in .doc form (S1-S4, one per organ), which are "
            "protein identifications rather than per-donor quantification."
        ),
        "reason": "ABSTAIN_NO_PER_DONOR_VALUES",
        "usable_as": "published_comparator",
        "why_it_matters": (
            "Sixteen individual donors with both proteomics and activity would "
            "have been a second Wegler. The activities are summarised across "
            "donors, so no donor-level pairing survives."
        ),
    },
}


class AuditError(RuntimeError):
    pass


def europe_pmc(query, timeout=60):
    url = (
        EUROPE_PMC
        + "?"
        + urllib.parse.urlencode(
            {"query": query, "resultType": "core", "format": "json"}
        )
    )
    request = urllib.request.Request(
        url, headers={"User-Agent": "Omics2PK-research/1.0"}
    )
    with urllib.request.urlopen(request, timeout=timeout) as handle:
        payload = json.load(handle)
    results = payload.get("resultList", {}).get("result", [])
    if not results:
        raise AuditError(f"Europe PMC returned no record for {query!r}")
    return results[0]


def qualify(entry):
    """Apply the declared rule. Order matters: access is checked first,
    because a source nobody can open cannot be assessed on anything else."""
    if entry["reason"] == "ABSTAIN_ACCESS":
        verdict = "ABSTAIN"
        reason = "ABSTAIN_ACCESS"
    elif entry["per_donor_layer_and_outcome"] is not True:
        verdict = "ABSTAIN"
        reason = "ABSTAIN_NO_PER_DONOR_VALUES"
    elif entry["n_donors"] < DONOR_FLOOR:
        verdict = "ABSTAIN"
        reason = "ABSTAIN_TOO_FEW_DONORS"
    elif entry["outcome_kind"] not in ADMISSIBLE_OUTCOMES:
        verdict = "ABSTAIN"
        reason = "ABSTAIN_NO_ADMISSIBLE_ENDPOINT"
    else:
        verdict = "ACQUIRED"
        reason = None
    if reason is not None and reason != entry["reason"]:
        raise AuditError(
            f"{entry['label']}: the rule returns {reason} but the recorded "
            f"reason is {entry['reason']}. One of the two is wrong and the "
            "disagreement must be resolved rather than reported."
        )
    return verdict, reason


def audit(check_network=True):
    report = {
        "audit": "candidate paired liver sources",
        "date": "2026-11-10",
        "fits_nothing": True,
        "rule": {
            "a_source_is_usable_only_if_all_four_hold": [
                "retrievable without subscription, authorisation request or "
                "collaboration",
                "at least one measured layer AND the outcome published or "
                "deposited per individual donor; a per-donor outcome paired "
                "with a POOLED layer does not satisfy this",
                f"at least {DONOR_FLOOR} donors carrying the required layers",
                "an admissible outcome kind, or a declared non-pharmacological "
                "target for Axis A",
            ],
            "donor_floor": DONOR_FLOOR,
            "why_the_floor": (
                "Leave-one-donor-out with an inner leave-one-out for the "
                "penalty leaves n-2 donors in the inner training set; under "
                "ten donors that is under eight, and neither the penalty "
                "choice nor the interval is interpretable."
            ),
            "admissible_outcomes": list(ADMISSIBLE_OUTCOMES),
            "reason_codes": list(REASONS),
            "a_source_failing_only_on_per_donor_values": (
                "may still be a published_comparator, never a dataset"
            ),
        },
        "how_the_two_kinds_of_fact_are_treated": (
            "Access and licence are re-queried from Europe PMC at every run "
            "and a disagreement with the recorded value fails the audit. "
            "Whether per-donor values exist cannot be read from an index, so "
            "it is recorded as evidence naming the file, its checksum and its "
            "contents."
        ),
        "candidates": {},
    }
    counts = {}
    for key, entry in CANDIDATES.items():
        verdict, reason = qualify(entry)
        row = {
            "label": entry["label"],
            "citation": entry["citation"],
            "n_donors": entry["n_donors"],
            "layers": list(entry["layers"]),
            "same_donor_layers": entry["same_donor_layers"],
            "outcome_kind": entry["outcome_kind"],
            "per_donor_layer_and_outcome": entry["per_donor_layer_and_outcome"],
            "verdict": verdict,
            "reason": reason,
            "usable_as": entry["usable_as"],
            "evidence": entry["evidence"],
            "why_it_matters": entry["why_it_matters"],
            "secondary_disqualification": entry.get("secondary_disqualification"),
        }
        if check_network:
            record = europe_pmc(entry["query"])
            observed = {
                field: record.get(field) for field in entry["expected"]
            }
            if observed != entry["expected"]:
                raise AuditError(
                    f"{entry['label']}: Europe PMC now reports {observed}, "
                    f"recorded {entry['expected']}. The access verdict rests "
                    "on these fields, so the audit stops rather than "
                    "reporting a stale one."
                )
            row["europe_pmc_reverified"] = observed
            row["licence_as_indexed"] = record.get("license")
            row["title_as_indexed"] = record.get("title")
        counts[reason or verdict] = counts.get(reason or verdict, 0) + 1
        report["candidates"][key] = row

    report["summary"] = {
        "n_candidates": len(CANDIDATES),
        "n_usable_as_datasets": sum(
            1 for row in report["candidates"].values()
            if row["verdict"] == "ACQUIRED"
        ),
        "by_reason": counts,
        "finding": (
            "Of seven candidates, NONE qualifies as a dataset. Two fail on "
            "access and FIVE on the absence of a per-donor layer paired with a "
            "per-donor outcome. SIX of the seven measured both a molecular "
            "layer and either the second layer or an activity in the SAME "
            "donors, so the measurements exist; what is missing is the "
            "individual-level sharing. That is a different obstacle from 'the "
            "data has never been collected', and it is the one the funding "
            "case should name. The seventh, Al-Majdoub, is the opposite case "
            "and worth separating: it PUBLISHES per-donor activity but "
            "measured abundance on pooled commercial lots, so no donor-level "
            "pairing was ever created."
        ),
        "what_this_does_not_claim": [
            "Not that these studies are flawed; publishing group-level "
            "summaries is a normal and legitimate choice.",
            "Not that no further source exists; this is seven candidates "
            "examined, not an exhaustive census.",
            "Not that a paywalled supplement lacks per-donor values: for "
            "Drozdzik 2018 and Wenzel 2024 the contents could not be checked "
            "at all, and the verdict is about access and nothing else.",
        ],
    }
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=None)
    parser.add_argument(
        "--offline",
        action="store_true",
        help="apply the rule without re-querying Europe PMC",
    )
    args = parser.parse_args(argv)
    try:
        sys.stdout.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass
    try:
        report = audit(check_network=not args.offline)
    except AuditError as exc:
        parser.exit(2, f"INVALID: {exc}\n")
    text = json.dumps(report, indent=2, ensure_ascii=False, sort_keys=True)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
        print(f"written: {args.out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
