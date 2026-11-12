"""Tests for the Lee / Jackson panel qualification.

The audit's whole value is that it refuses to proceed on a mis-parsed table, so
these tests mostly check that the refusals fire.
"""
from __future__ import annotations

import importlib.util
import json
import pathlib
import re
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_lee_jackson_cyp3a_panel.py"
ARTEFACT = ROOT / "reports" / "data" / "lee_jackson_cyp3a_panel.json"
DOC = (ROOT / "docs"
       / "OMICS2PK_LEE_JACKSON_CYP3A_PANEL_QUALIFICATION_20261110.md")


@pytest.fixture(scope="module")
def mod():
    spec = importlib.util.spec_from_file_location("lj_audit", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["lj_audit"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def payload():
    if not ARTEFACT.is_file():
        pytest.skip("the audit has not been run in this checkout")
    return json.loads(ARTEFACT.read_text(encoding="utf-8"))


def _rows(n=15, **overrides):
    out = {}
    for index in range(n):
        donor = f"D{index:02d}"
        out[donor] = {
            "age_years": 20.0 + index, "sex": "M", "race": "C",
            "bmi": 25.0,
            "cyp3a4_abundance": 10.0 + index,
            "cyp3a5_abundance": 1.0 + index * 0.1,
            "midazolam_1_hydroxylation": 5.0 + index * 2,
            "ibrutinib_m37": 4.0 + index * 1.5,
            "ibrutinib_m34": 0.1 + index * 0.01,
            "ibrutinib_m25": 0.2 + index * 0.03,
            "plasma_4b_hc_ng_per_ml": 10.0 + index,
            "plasma_cholesterol_mg_per_ml": 1.0,
            "plasma_4b_hc_cholesterol_ratio": 0.1 + index * 0.02,
        }
    for donor, patch in overrides.items():
        out[donor].update(patch)
    return out


# --------------------------------------------------------------------------
# 1. the statistics the audit carries itself
# --------------------------------------------------------------------------
class TestSpearman:
    def test_a_perfect_monotone_pair_is_one(self, mod):
        assert mod.spearman([1, 2, 3, 4], [10, 20, 30, 40]) == pytest.approx(1.0)

    def test_a_reversed_pair_is_minus_one(self, mod):
        assert mod.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == pytest.approx(-1.0)

    def test_it_is_invariant_to_a_monotone_transform(self, mod):
        left = [1.0, 2.0, 7.0, 9.0, 11.0]
        right = [3.0, 1.0, 8.0, 4.0, 20.0]
        import math
        assert mod.spearman(left, right) == pytest.approx(
            mod.spearman([math.log(v) for v in left], right))

    def test_ties_get_the_average_rank(self, mod):
        assert mod.ranks([5.0, 5.0, 9.0]) == [1.5, 1.5, 3.0]

    def test_a_constant_side_is_refused(self, mod):
        with pytest.raises(mod.AuditError, match="variation on both sides"):
            mod.spearman([1, 1, 1, 1], [1, 2, 3, 4])

    def test_too_few_pairs_is_refused(self, mod):
        with pytest.raises(mod.AuditError, match="at least three"):
            mod.spearman([1, 2], [3, 4])


# --------------------------------------------------------------------------
# 2. the join refuses what it should refuse
# --------------------------------------------------------------------------
class TestJoin:
    def test_the_same_fifteen_donors_join(self, mod):
        rows = _rows()
        abundance = {d: {k: v for k, v in r.items()
                         if k in {"age_years", "sex", "race", "bmi",
                                  "cyp3a4_abundance", "cyp3a5_abundance"}}
                     for d, r in rows.items()}
        activity = {d: {k: v for k, v in r.items()
                        if k.startswith(("midazolam", "ibrutinib"))}
                    for d, r in rows.items()}
        plasma = {d: {k: v for k, v in r.items() if k.startswith("plasma")}
                  for d, r in rows.items()}
        joined = mod.join_on_the_same_donors(abundance, activity, plasma)
        assert len(joined) == 15

    def test_a_different_donor_set_is_refused(self, mod):
        rows = _rows()
        abundance = {d: {"age_years": 1.0, "sex": "M", "race": "C",
                         "bmi": 2.0, "cyp3a4_abundance": 3.0,
                         "cyp3a5_abundance": 4.0} for d in rows}
        activity = dict(abundance)
        activity.pop("D00")
        with pytest.raises(mod.AuditError, match="same donors"):
            mod.join_on_the_same_donors(abundance, activity, abundance)

    def test_a_wrong_donor_count_is_refused(self, mod):
        small = {f"D{i}": {"age_years": 1.0, "sex": "M", "race": "C",
                           "bmi": 2.0, "cyp3a4_abundance": 3.0,
                           "cyp3a5_abundance": 4.0} for i in range(9)}
        with pytest.raises(mod.AuditError, match="hepatocyte donors"):
            mod.join_on_the_same_donors(small, small, small)

    def test_a_different_row_order_is_refused_not_silently_reordered(self, mod):
        base = {f"D{i:02d}": {"age_years": 1.0, "sex": "M", "race": "C",
                              "bmi": 2.0, "cyp3a4_abundance": 3.0,
                              "cyp3a5_abundance": 4.0} for i in range(15)}
        shuffled = {key: base[key] for key in reversed(list(base))}
        with pytest.raises(mod.AuditError, match="DIFFERENT order"):
            mod.join_on_the_same_donors(base, shuffled, base)

    def test_a_missing_cell_is_refused(self, mod):
        base = {f"D{i:02d}": {"age_years": 1.0, "sex": "M", "race": "C",
                              "bmi": 2.0, "cyp3a4_abundance": 3.0,
                              "cyp3a5_abundance": 4.0} for i in range(15)}
        base["D03"]["cyp3a4_abundance"] = None
        with pytest.raises(mod.AuditError, match="no value for"):
            mod.join_on_the_same_donors(base, base, base)


# --------------------------------------------------------------------------
# 3. the published anchors are the gate
# --------------------------------------------------------------------------
class TestAnchors:
    def test_three_anchors_are_declared(self, mod):
        assert len(mod.PUBLISHED_ANCHORS) == 3

    def test_two_anchors_exclude_the_infant_donor_the_article_excludes(self, mod):
        excluding = [a for a in mod.PUBLISHED_ANCHORS
                     if mod.INFANT_DONOR in a["exclude"]]
        assert len(excluding) == 2

    def test_the_tolerance_leaves_room_for_rounding_and_no_more(self, mod):
        assert 0.005 < mod.ANCHOR_TOLERANCE <= 0.02

    def test_a_table_that_misses_an_anchor_is_refused(self, mod):
        rows = _rows()
        # Break the relation the anchors measure: shuffle one column by one
        # position, which leaves every marginal distribution intact.
        donors = list(rows)
        values = [rows[d]["plasma_4b_hc_cholesterol_ratio"] for d in donors]
        for donor, value in zip(donors, values[1:] + values[:1]):
            rows[donor]["plasma_4b_hc_cholesterol_ratio"] = value
        with pytest.raises(mod.AuditError, match="does not reproduce"):
            mod.verify_published_anchors(rows)

    def test_the_real_run_reproduced_all_three(self, payload):
        for row in payload["published_anchors_verified"]:
            assert row["within_tolerance"] is True
            assert row["absolute_gap"] <= 0.02


# --------------------------------------------------------------------------
# 4. the per-donor table never lands in the repository
# --------------------------------------------------------------------------
class TestPerDonorTable:
    def test_writing_inside_the_repository_is_refused(self, mod, tmp_path):
        with pytest.raises(mod.AuditError, match="must not be written inside"):
            mod.write_per_donor_table(_rows(), str(ROOT / "reports"))

    def test_writing_outside_the_repository_succeeds(self, mod, tmp_path):
        target = mod.write_per_donor_table(_rows(), str(tmp_path))
        assert target.is_file()
        assert target.read_text(encoding="utf-8").count("\n") == 16

    def test_the_run_wrote_it_outside(self, payload):
        where = pathlib.Path(payload[
            "per_donor_table_is_outside_this_repository"])
        assert ROOT not in where.parents and ROOT != where.parent


# --------------------------------------------------------------------------
# 5. the artefact and the document
# --------------------------------------------------------------------------
class TestArtefactAndDocument:
    def test_the_audit_fits_nothing(self, payload):
        assert payload["this_audit_fits_nothing"] is True
        assert payload[
            "no_donor_identifiers_values_or_predictions_exported"] is True

    def test_it_holds_no_long_numeric_list(self, payload):
        def arrays(obj):
            found = []
            if isinstance(obj, dict):
                for value in obj.values():
                    found += arrays(value)
            elif isinstance(obj, list):
                numeric = [x for x in obj if isinstance(x, (int, float))]
                if len(numeric) >= 10:
                    found.append(len(numeric))
                for value in obj:
                    if not isinstance(value, (int, float)):
                        found += arrays(value)
            return found

        assert not arrays(payload)

    def test_the_donor_floor_is_met(self, payload):
        assert payload["n_donors"] == 15
        assert payload["meets_donor_floor"] is True

    def test_the_circularity_check_is_recorded_as_answered(self, payload):
        check = payload["circularity_check"]
        assert check["answer"] == "no"
        assert "DIFFERENT denominators" in check["cost_of_that_difference"]

    def test_the_venetoclax_paper_is_recorded_as_the_same_cohort(self, payload):
        text = payload["independence"]["from_the_venetoclax_paper"]
        assert "NOT independent" in text
        assert "never two cohorts" in text

    def test_the_document_carries_every_recorded_correlation(self, payload):
        text = DOC.read_text(encoding="utf-8")
        for name, value in payload[
                "rank_correlations_all_donors"].items():
            printed = f"{value:+.4f}".replace("-", "−")
            assert printed in text, (
                f"{name} is {printed} in the artefact and absent from the "
                "qualification document")

    def test_the_document_names_what_the_source_cannot_do(self, payload):
        text = re.sub(r"[\s*_`>|]+", " ", DOC.read_text(encoding="utf-8"))
        assert "No intrinsic clearance" in text
        assert "No microsomal outcome" in text
        assert "Procurement plasma" in text

    def test_the_document_states_the_anchor_gaps(self, payload):
        text = DOC.read_text(encoding="utf-8")
        for row in payload["published_anchors_verified"]:
            assert f"{row['absolute_gap']:.4f}" in text
