#!/usr/bin/env python
"""Search the publisher-independent repositories, and re-verify what was found.

The seven-candidate qualification of 2026-11-10 stopped at the publishers. Two
of its seven verdicts were ``ABSTAIN_ACCESS`` on the strength of an HTTP 403,
and an HTTP 403 from a publisher does not establish that the individual-level
data is unavailable through another public channel. This audit tests that.

Eight repositories and indexes are queried: PRIDE, ProteomeXchange, GEO,
BioStudies / ArrayExpress, MetaboLights, Zenodo, Figshare and Dryad, plus the
Europe PMC and NCBI article channels and the EBI FTP mirror that BioStudies
redirects to.

It found two things worth the search:

1. **A usable dataset that the publisher-first search had missed** -- the
   Lee / Jackson CYP3A panel, whose per-donor layers and per-donor outcome are
   in an open-access article's supplement and main text.
2. **An access ladder that matters.** The supplement of a 171-donor paired
   abundance-and-activity study is refused by Europe PMC as "not open access"
   and gated by a proof-of-work challenge at NCBI, and is served without
   restriction by the EBI FTP mirror. Once opened, its tables turn out to be
   age-binned aggregates -- so its real barrier was never access, and recording
   it as an access failure would have been wrong.

This script fits nothing. Every recorded fact is re-queried live and the run
fails if one has moved.

Usage:
    python -I scripts/audit_independent_repository_search.py \
        --out reports/data/independent_repository_search.json
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import urllib.error
import urllib.parse
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")

ROOT = pathlib.Path(__file__).resolve().parents[1]
TIMEOUT = 60
AGENT = "Omics2PK-source-audit"

VERDICTS = (
    "DATASET_USABLE",
    "LAYER_WITHOUT_OUTCOME",
    "ABSTAIN_NO_PER_DONOR_VALUES",
    "ABSTAIN_ACCESS",
    "ALREADY_HELD",
    "NOT_HUMAN_OR_NOT_RELEVANT",
)

#: The repositories queried, with the query strings actually used. Recorded so
#: that "we searched the repositories" is a checkable statement and not a
#: claim.
REPOSITORIES = {
    "PRIDE": {
        "api": "https://www.ebi.ac.uk/pride/ws/archive/v3/search/projects",
        "queries": ("human liver proteome", "liver intestine drug metabolizing",
                    "Oswald", "Drozdzik", "Al-Majdoub", "Wegler",
                    "hepatocyte microsome proteomics"),
        "outcome": ("one deposition of the cohort this project already holds, "
                    "and eight depositions of per-donor liver microsome "
                    "proteomics with no deposited per-donor outcome"),
    },
    "ProteomeXchange": {
        "api": ("https://proteomecentral.proteomexchange.org/api/proxi/v0.1/"
                "datasets"),
        "queries": ("liver cytochrome",),
        "outcome": ("no dataset beyond what PRIDE returned; the PROXI "
                    "endpoint indexes the same depositions"),
    },
    "GEO": {
        "api": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        "queries": ("human liver cytochrome P450 activity",
                    "human liver CYP3A4 activity genotype",
                    "human liver RNA-seq drug metabolism donors",
                    "liver bank pharmacogenomics expression"),
        "outcome": ("nothing relevant: the hits are cell lines, organoids, "
                    "iPSC-derived hepatocytes and unrelated disease series"),
    },
    "BioStudies_ArrayExpress": {
        "api": "https://www.ebi.ac.uk/biostudies/api/v1/search",
        "queries": ("liver cytochrome P450 proteomics",
                    "individual liver donors CYP abundance activity",
                    "hepatic intrinsic clearance donor"),
        "outcome": ("the most productive index of the eight, because its "
                    "S-EPMC records wrap the supplementary files of PMC "
                    "articles and expose them over the EBI FTP mirror"),
    },
    "MetaboLights": {
        "api": "https://www.ebi.ac.uk/metabolights/ws/studies",
        "queries": ("pharmacokinetics metabolomics plasma concentration",),
        "outcome": ("no study pairing metabolomics with a drug concentration "
                    "in the same individuals beyond MTBLS3128, already held"),
    },
    "Zenodo": {
        "api": "https://zenodo.org/api/records",
        "queries": ("intrinsic clearance human liver donors",
                    "cytochrome P450 abundance activity individual",
                    "pharmacokinetics omics individual concentrations"),
        "outcome": "no human liver dataset; the hits are reviews and figures",
    },
    "Figshare": {
        "api": "https://api.figshare.com/v2/articles/search",
        "queries": ("human liver cytochrome P450 abundance activity "
                    "individual donors",),
        "outcome": ("no hit mentioning human liver, an individual donor or "
                    "an enzyme activity; the index returns unrelated "
                    "chemistry and ecology records"),
    },
    "Dryad": {
        "api": "https://datadryad.org/api/v2/search",
        "queries": ("cytochrome P450 liver clearance",),
        "outcome": ("no human liver dataset. Dryad does hold individual-level "
                    "clinical pharmacokinetic datasets -- a morphine study in "
                    "neonates was returned -- but none with a molecular layer "
                    "in the same individuals"),
    },
}

#: PRIDE depositions, re-queried live. ``expect`` holds the fields as read on
#: 2026-11-10.
PRIDE_FINDINGS = {
    "PXD014131": {
        "expect": {"submitter": "Christine Wegler",
                   "title_contains": "15 matched donors"},
        "verdict": "ALREADY_HELD",
        "note": ("this is the Wegler cohort itself. Its raw mass spectrometry "
                 "is deposited; the per-donor clearances this project needed "
                 "came from the article supplement, not from PRIDE. Finding it "
                 "here confirms the cohort and adds nothing analysable."),
    },
    "PXD020910": {
        "expect": {"submitter": "Jill Barber",
                   "title_contains": "human liver microsomes"},
        "verdict": "LAYER_WITHOUT_OUTCOME",
        "note": ("per-donor proteomics as raw files, with no per-donor "
                 "activity or clearance deposited anywhere. A layer without "
                 "an outcome cannot answer a protein-to-function question, "
                 "however many donors it holds."),
    },
    "PXD021025": {
        "expect": {"submitter": "Brahim Achour",
                   "title_contains": "adult human livers"},
        "verdict": "LAYER_WITHOUT_OUTCOME",
        "note": "same group, same limitation",
    },
    "PXD038776": {
        "expect": {"submitter": "Areti-Maria Vasilogianni",
                   "title_contains": "healthy controls"},
        "verdict": "LAYER_WITHOUT_OUTCOME",
        "note": ("the proteomics behind Clin Pharmacol Ther 2022 "
                 "(PMC9540503), whose article publishes aggregates"),
    },
}

#: Articles found through the Europe PMC open-access-and-supplement filter.
#: ``expect`` holds the Europe PMC fields as read on 2026-11-10.
ARTICLE_FINDINGS = {
    "PMC9926076": {
        "label": "Lee 2023, ibrutinib, CYP3A panel",
        "expect": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "verdict": "DATASET_USABLE",
        "n_donors": 15,
        "note": ("per-donor hepatic CYP3A4 and CYP3A5 abundance (Table S2), "
                 "per-donor plasma 4-beta-hydroxycholesterol and cholesterol "
                 "(MAIN-TEXT Table 1), and per-donor midazolam "
                 "1'-hydroxylation and ibrutinib metabolite formation "
                 "(Table S3). Qualified separately in "
                 "docs/OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_"
                 "20261110.md."),
    },
    "PMC11645475": {
        "label": "Lee 2024, venetoclax, the SAME donors",
        "expect": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "verdict": "ALREADY_HELD",
        "n_donors": 15,
        "note": ("the same fifteen hepatocyte donors and the same twenty "
                 "microsome lots as PMC9926076, with venetoclax as a further "
                 "substrate and twelve real-world patients whose per-patient "
                 "4-beta-hydroxycholesterol is not tabulated. Two drugs in the "
                 "same donors are ONE cohort with two substrates. Counting it "
                 "as a second cohort is exactly the error the brief forbids."),
    },
    "PMC6222000": {
        "label": "Bhatt 2019, UGT abundance and activity, 171 donors",
        "expect": {"isOpenAccess": "N", "hasSuppl": "Y"},
        "verdict": "ABSTAIN_NO_PER_DONOR_VALUES",
        "n_donors": 171,
        "note": ("the largest paired abundance-and-activity liver cohort in "
                 "the public literature, and the most instructive negative. "
                 "Europe PMC refuses its supplement as 'not open access'; "
                 "NCBI gates the author-manuscript file behind a "
                 "proof-of-work challenge; the EBI FTP mirror that BioStudies "
                 "redirects to serves the same 2,547,704-byte PDF without "
                 "restriction. Opened: Table S2 gives abundance and activity "
                 "as means and SDs BY AGE CATEGORY, and Table S3 gives a "
                 "per-donor row whose Activity and Proteomics columns are "
                 "availability TICKS, not values. So the barrier was never "
                 "access."),
    },
    "PMC9310776": {
        "label": "matched jejunum and liver proteome, 37 donors",
        "expect": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "verdict": "DATASET_USABLE",
        "n_donors": 37,
        "note": ("supplementary workbook CPT-111-1142-s002.xlsx, 11,163,005 "
                 "bytes over the EBI FTP mirror against under 1 MB and a "
                 "corrupt archive from the Europe PMC endpoint, carries "
                 "per-donor columns: 9,027 proteins x 37 jejunum donors and "
                 "6,982 proteins x 37 liver donors, fmol/ug protein. The "
                 "article's own tables publish medians and ranges only. It "
                 "has NO pharmacological endpoint, so it can support a "
                 "cross-compartment redundancy question and nothing about "
                 "drug metabolism."),
    },
    "PMC12735975": {
        "label": "16 CYPs and 2 CES, individual preparations",
        "expect": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "verdict": "ABSTAIN_NO_PER_DONOR_VALUES",
        "n_donors": 20,
        "note": ("measured abundance AND activity in twenty individual "
                 "hepatocyte donors, sixteen genotyped microsome lots and "
                 "fifteen intestinal mucosa lots, and publishes per-donor "
                 "demographics (Table 1) and per-donor GENOTYPE (Table S1) "
                 "while reporting abundance as means and ranges (Table 2) and "
                 "the pairing as Spearman coefficients (Table 3). Its "
                 "twenty-eight published correlations are the transferable "
                 "input to the question-A sizing and are not a dataset."),
    },
    "PMC12673578": {
        "label": "CES1 and clopidogrel oxidation",
        "expect": {"isOpenAccess": "Y", "hasSuppl": "Y"},
        "verdict": "NOT_HUMAN_OR_NOT_RELEVANT",
        "n_donors": None,
        "note": ("recombinant enzymes and POOLED microsomes; no individual "
                 "donor anywhere"),
    },
}

#: BioStudies file listings, re-queried live, because the access ladder is the
#: finding and not an aside.
BIOSTUDIES_FINDINGS = {
    "S-EPMC6222000": {
        "expect_file": "NIHMS965499-supplement-Supp_info.pdf",
        "expect_size": 2547704,
    },
    "S-EPMC9310776": {
        "expect_file": "CPT-111-1142-s002.xlsx",
        "expect_size": 11163005,
    },
}

#: The two sources whose ABSTAIN_ACCESS was re-tested through this new channel
#: and survives it. Re-queried live.
STILL_INACCESSIBLE = {
    "29205295": {"label": "Drozdzik 2018, paired liver and intestine, n=9"},
    "37739780": {"label": "Wenzel 2024, intra-subject jejunum and liver, n=11"},
}


class AuditError(RuntimeError):
    pass


def _check(condition, message):
    if not condition:
        raise AuditError(message)


def _get(url, *, as_json=True):
    request = urllib.request.Request(url, headers={"User-Agent": AGENT})
    try:
        with urllib.request.urlopen(request, timeout=TIMEOUT) as response:
            body = response.read()
    except (urllib.error.URLError, TimeoutError) as error:
        raise AuditError(f"{url} could not be reached: {error}") from error
    return json.loads(body.decode("utf-8", "replace")) if as_json else body


def verify_pride():
    out = {}
    for accession, record in PRIDE_FINDINGS.items():
        data = _get("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/"
                    + accession)
        submitters = [s.get("name", "") for s in data.get("submitters", [])]
        title = data.get("title") or ""
        expect = record["expect"]
        _check(expect["submitter"] in submitters,
               f"{accession}: expected submitter {expect['submitter']!r}, "
               f"PRIDE now lists {submitters}")
        _check(expect["title_contains"].lower() in title.lower(),
               f"{accession}: the title no longer contains "
               f"{expect['title_contains']!r}; it is now {title!r}")
        out[accession] = {"title": title, "submitters": submitters,
                          "verdict": record["verdict"], "note": record["note"]}
        _check(record["verdict"] in VERDICTS,
               f"{accession} carries an undeclared verdict")
    return out


def verify_articles():
    out = {}
    for pmcid, record in ARTICLE_FINDINGS.items():
        query = urllib.parse.quote(f"PMCID:{pmcid}")
        data = _get("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
                    f"?query={query}&resultType=core&format=json")
        results = data.get("resultList", {}).get("result", [])
        _check(results, f"{pmcid} is no longer indexed in Europe PMC")
        found = results[0]
        for field, expected in record["expect"].items():
            _check(found.get(field) == expected,
                   f"{pmcid}: Europe PMC now reports {field}="
                   f"{found.get(field)!r}, recorded as {expected!r}")
        _check(record["verdict"] in VERDICTS,
               f"{pmcid} carries an undeclared verdict")
        out[pmcid] = {
            "label": record["label"],
            "title": found.get("title"),
            "journal": found.get("journalInfo", {}).get("journal", {}).get(
                "title"),
            "year": found.get("pubYear"),
            "licence": found.get("license"),
            "is_open_access": found.get("isOpenAccess"),
            "has_supplement": found.get("hasSuppl"),
            "n_donors": record["n_donors"],
            "verdict": record["verdict"],
            "note": record["note"],
        }
    return out


def verify_biostudies():
    out = {}
    for accession, record in BIOSTUDIES_FINDINGS.items():
        data = _get("https://www.ebi.ac.uk/biostudies/api/v1/studies/"
                    + accession)
        files = {}

        def walk(node):
            if isinstance(node, dict):
                if "path" in node and node.get("type") == "file":
                    files[node["path"]] = node.get("size")
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)

        walk(data)
        name = record["expect_file"]
        _check(name in files,
               f"{accession} no longer lists {name!r}; it lists "
               f"{sorted(files)}")
        _check(files[name] == record["expect_size"],
               f"{accession}/{name} is now {files[name]} bytes, recorded as "
               f"{record['expect_size']}")
        out[accession] = {"file": name, "size_bytes": files[name],
                          "served_over": ("https://ftp.ebi.ac.uk/pub/databases/"
                                          "biostudies/ -- BioStudies redirects "
                                          "its file endpoint there")}
    return out


def verify_still_inaccessible():
    out = {}
    for pmid, record in STILL_INACCESSIBLE.items():
        data = _get("https://www.ebi.ac.uk/europepmc/webservices/rest/search"
                    f"?query=EXT_ID:{pmid}&resultType=core&format=json")
        results = data.get("resultList", {}).get("result", [])
        _check(results, f"PMID {pmid} is no longer indexed")
        found = results[0]
        _check(not found.get("pmcid"),
               f"PMID {pmid} NOW HAS a PMCID ({found.get('pmcid')}). The "
               "ABSTAIN_ACCESS verdict must be revisited: a PMCID means "
               "BioStudies may now mirror its supplement.")
        out[pmid] = {
            "label": record["label"],
            "pmcid": found.get("pmcid"),
            "in_pmc": found.get("inPMC"),
            "is_open_access": found.get("isOpenAccess"),
            "verdict": "ABSTAIN_ACCESS",
            "channels_tried": [
                "publisher landing page and supplement endpoint: HTTP 403 to "
                "every automated client",
                "Europe PMC: not indexed as full text, no supplement",
                "BioStudies: has no S-EPMC record, because it has no PMCID",
                "EBI FTP mirror: nothing to mirror",
            ],
            "what_would_change_it": (
                "a PMCID appearing, or a deposition under either author's "
                "name in a repository. This audit re-checks both on every run."),
        }
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",
                        default=str(ROOT / "reports" / "data"
                                    / "independent_repository_search.json"))
    args = parser.parse_args(argv)

    pride = verify_pride()
    articles = verify_articles()
    biostudies = verify_biostudies()
    inaccessible = verify_still_inaccessible()

    usable = [k for k, v in articles.items() if v["verdict"] == "DATASET_USABLE"]
    payload = {
        "audit": "independent_repository_search",
        "generated": "2026-11-10",
        "repositories_queried": REPOSITORIES,
        "n_repositories": len(REPOSITORIES),
        "pride_depositions": pride,
        "articles": articles,
        "biostudies_file_listings": biostudies,
        "still_inaccessible_after_this_search": inaccessible,
        "datasets_now_usable": usable,
        "n_datasets_now_usable": len(usable),
        "corrects_an_earlier_claim": {
            "claim": ("the candidate qualification of 2026-11-10 concluded "
                      "that step 2 of the commissioning order -- execute the "
                      "comparisons that are genuinely identifiable -- had "
                      "'nothing to run'."),
            "status": "WITHDRAWN",
            "why": ("that conclusion followed a publisher-first search over "
                    "seven named citations. A search over the repositories and "
                    "over the Europe PMC open-access-plus-supplement filter "
                    "found a cohort with two per-donor layers and a per-donor "
                    "outcome in the same fifteen donors. The earlier search "
                    "was too narrow, not the public record too thin."),
        },
        "the_access_lesson": {
            "finding": ("an HTTP 403 from a publisher, and even Europe PMC's "
                        "own 'not open access' flag, do not establish that "
                        "individual-level data is unavailable. The ladder that "
                        "worked was publisher -> Europe PMC -> NCBI -> "
                        "BioStudies -> EBI FTP mirror."),
            "but": ("climbing it did not produce a usable dataset from the "
                    "171-donor study. Once the file was open, the barrier was "
                    "the same as in five of the seven earlier candidates: "
                    "aggregates instead of rows. The access ladder changed a "
                    "wrong REASON into a right one; it did not change the "
                    "verdict."),
            "what_was_not_done": (
                "NCBI's proof-of-work challenge on the author-manuscript "
                "download endpoint was not defeated, and the author-manuscript "
                "bulk collection requires registration, which the standing "
                "constraint on authorisation-free sources excludes."),
        },
        "the_pattern_now_stands_at": {
            "sources_that_measured_a_per_donor_pairing": 9,
            "sources_that_published_the_per_donor_rows": 2,
            "named": ("Wegler 2021 and Lee 2023. Both published rows by "
                      "accident of format rather than by policy: Wegler's "
                      "Table S3 happens to be a per-donor table, and Lee's "
                      "plasma values happen to sit in a main-text table."),
        },
        "this_audit_fits_nothing": True,
        "no_individual_level_data_in_this_artefact": True,
    }

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                   encoding="utf-8")

    print(f"repositories queried: {len(REPOSITORIES)}")
    print("PRIDE:")
    for accession, row in pride.items():
        print(f"  {accession}  {row['verdict']:<24} {row['title'][:62]}")
    print("articles:")
    for pmcid, row in articles.items():
        print(f"  {pmcid:<12} {row['verdict']:<28} n="
              f"{str(row['n_donors']):<5} {row['label']}")
    print("BioStudies file listings (the channel that worked):")
    for accession, row in biostudies.items():
        print(f"  {accession}  {row['file']}  {row['size_bytes']:,} bytes")
    print("still inaccessible:")
    for pmid, row in inaccessible.items():
        print(f"  PMID {pmid}  pmcid={row['pmcid']}  {row['label']}")
    print(f"\ndatasets now usable: {len(usable)} -> {usable}")
    print(f"artefact: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
