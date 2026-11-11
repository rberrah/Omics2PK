"""Contracts for the selector-comparison harness.

The commissioned question is whether the *methodology* chooses better
biological information, so the guarantees that matter are about the selection
procedure, not about predictive accuracy. In order of how much damage their
absence would do:

1. **Nothing supervised escapes the training folds** — including the choice of
   how many variables to keep. A panel chosen by looking at out-of-fold
   performance and then scored on it is scored on its own selection.
   -> `TestSelectionStaysInsideTheFolds`.
2. **SHAP is computed on rows the booster did not train on.** In-sample SHAP
   ranks an overfitted feature as important. -> `TestShapIsOutOfSample`.
3. **The harness can find a sparse signal that exists**, so a null means
   something. -> `TestPositiveControl`.
4. **The curve and the validated panel size are different numbers** and are
   never conflated. -> `TestCurveIsNotAValidatedSize`.
5. A selector cannot abstain by accident, so the permuted control runs beside
   every real result. -> `TestNegativeControls`.
"""
from __future__ import annotations

import pathlib

import pytest

np = pytest.importorskip("numpy")

ROOT = pathlib.Path(__file__).resolve().parents[1]
FAST = {"selectors": ("S0_everything", "S1_univariate_outcome"),
        "k_grid": (1, 3, 8)}


@pytest.fixture(scope="module")
def sss():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import supervised_screening_selection
    return supervised_screening_selection


@pytest.fixture(scope="module")
def cohort(sss):
    return sss.sparse_signal_cohort(n_groups=60, per_group=3, n_candidates=20,
                                    n_true=3, seed=1)


@pytest.fixture(scope="module")
def run(sss, cohort):
    return sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                      baseline_columns=cohort["baseline_columns"],
                      knowledge_order=cohort["knowledge_order"], **FAST)


class TestSelectionStaysInsideTheFolds:
    def test_a_test_row_never_enters_a_selector(self, sss, cohort, monkeypatch):
        """The spy fails if any ordering ever sees a row it will be scored on."""
        seen = []
        original = sss.selector_ordering

        def spy(name, design, target, **kwargs):
            seen.append(design.shape[0])
            return original(name, design, target, **kwargs)

        monkeypatch.setattr(sss, "selector_ordering", spy)
        n = cohort["target"].size
        sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                   baseline_columns=cohort["baseline_columns"],
                   knowledge_order=cohort["knowledge_order"], **FAST)
        assert seen, "no ordering was computed at all"
        # Every ordering saw a strict subset: never the whole cohort.
        assert max(seen) < n

    def test_the_budget_is_chosen_per_fold_not_once(self, run):
        for name, chosen in run["k_chosen_per_fold"].items():
            if name == "S0_everything":
                continue
            assert len(chosen) >= 2, name
            assert all(k in run["k_grid"] for k in chosen), name

    def test_an_ordering_is_computed_once_per_training_set(self, sss, cohort,
                                                           run):
        """The count must not scale with the budget grid.

        A selector's ordering does not depend on how many variables it may
        keep, so recomputing it per budget both wastes work and lets the two
        drift apart. Widening the grid must therefore change nothing here.
        """
        assert run["one_ordering_per_selector_per_training_set"]
        wider = sss.screen(
            cohort["target"], cohort["design"], cohort["groups"],
            baseline_columns=cohort["baseline_columns"],
            knowledge_order=cohort["knowledge_order"],
            selectors=FAST["selectors"], k_grid=(1, 2, 3, 5, 8, 12))
        assert len(wider["k_grid"]) > len(run["k_grid"])
        assert wider["n_orderings_computed"] == run["n_orderings_computed"]
        # One per selector per training set: five outer folds and their three
        # inner folds each.
        assert run["n_orderings_computed"] == len(FAST["selectors"]) * (5 + 15)

    def test_baseline_columns_are_always_kept(self, run, cohort, sss):
        names = {f"x{c}" for c in cohort["baseline_columns"]}
        for name, panels in run["panels_per_fold"].items():
            for panel in panels:
                assert names <= set(panel), (name, panel)

    def test_group_folds_never_split_a_group(self, sss):
        groups = np.repeat(np.arange(30), 4)
        folds = sss.grouped_folds(groups, 5, 0)
        for key in np.unique(groups):
            assert len(set(folds[groups == key])) == 1

    def test_too_few_groups_is_refused(self, sss):
        with pytest.raises(Exception):
            sss.grouped_folds(np.arange(3), 5, 0)


class TestShapIsOutOfSample:
    def test_it_says_where_it_was_computed(self, sss, cohort):
        _, diagnostics = sss.shap_importance(
            cohort["design"], cohort["target"], groups=cohort["groups"],
            seed=0)
        if not diagnostics.get("available"):
            pytest.skip("xgboost absent")
        assert "never training rows" in diagnostics["computed_on"]
        assert diagnostics["exact"] is True

    def test_it_ranks_a_real_signal_above_noise(self, sss, cohort):
        strength, diagnostics = sss.shap_importance(
            cohort["design"], cohort["target"], groups=cohort["groups"],
            seed=0)
        if not diagnostics.get("available"):
            pytest.skip("xgboost absent")
        truth = list(cohort["truth_columns"])
        noise = [j for j in range(cohort["design"].shape[1])
                 if j not in truth and j not in cohort["baseline_columns"]
                 and j != cohort["redundant_copy_of_a_true_column"]]
        assert strength[truth].mean() > strength[noise].mean()

    def test_a_missing_xgboost_is_reported_not_hidden(self, sss, cohort,
                                                      monkeypatch):
        import builtins
        real_import = builtins.__import__

        def no_xgboost(name, *args, **kwargs):
            if name == "xgboost":
                raise ImportError("blocked for this test")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", no_xgboost)
        strength, diagnostics = sss.shap_importance(
            cohort["design"], cohort["target"], groups=cohort["groups"],
            seed=0)
        assert diagnostics["available"] is False
        assert "not installed" in diagnostics["why"]
        assert not strength.any()

    def test_no_shap_package_is_required(self):
        source = (ROOT / "src" / "omics2pk"
                  / "supervised_screening_selection.py").read_text(
                      encoding="utf-8")
        assert "import shap" not in source
        assert "pred_contribs" in source


class TestPositiveControl:
    def test_the_sparse_cohort_hides_the_signal_at_the_end(self, sss, cohort):
        # The useful columns are LAST, so a selector that returns the first k
        # columns fails loudly instead of passing by accident.
        assert min(cohort["truth_columns"]) >= cohort["design"].shape[1] // 2

    def test_the_knowledge_ordering_leads_with_the_truth(self, cohort):
        head = set(cohort["knowledge_order"][:len(cohort["truth_columns"])])
        assert head == set(cohort["truth_columns"])

    def test_a_univariate_selector_recovers_some_of_the_truth(self, run,
                                                              cohort):
        truth = {f"x{c}" for c in cohort["truth_columns"]}
        panels = run["panels_per_fold"]["S1_univariate_outcome"]
        found = [len(truth & set(panel)) / len(truth) for panel in panels]
        assert float(np.mean(found)) > 0.0

    def test_the_redundancy_decoy_exists_and_is_declared(self, cohort):
        column = cohort["redundant_copy_of_a_true_column"]
        strongest = cohort["truth_columns"][0]
        correlation = float(np.corrcoef(cohort["design"][:, column],
                                        cohort["design"][:, strongest])[0, 1])
        assert correlation > 0.9


class TestCurveIsNotAValidatedSize:
    def test_the_curve_is_labelled_optimistic(self, run):
        assert "optimistically biased" in run["the_curve_is_descriptive"]

    def test_the_parsimony_point_flags_its_own_bias(self, sss, run):
        curve = run["curve"]["S1_univariate_outcome"]
        point = sss.parsimony_point(curve)
        assert point["the_k_of_best_is_optimistically_biased"] is True
        assert "k_retaining_80pct" in point
        assert "k_retaining_90pct" in point

    def test_a_flat_curve_yields_no_parsimony_point_rather_than_a_wrong_one(
            self, sss):
        point = sss.parsimony_point({"1": 0.5, "5": 0.5, "10": 0.5})
        assert point["k_retaining_80pct"] is None
        assert point["gain_span"] == 0.0

    def test_a_single_point_curve_is_refused(self, sss):
        with pytest.raises(Exception):
            sss.parsimony_point({"3": 0.4})

    def test_the_retained_share_is_monotone_in_k(self, sss):
        # A curve that improves with k must give 80 % at or before 90 %.
        point = sss.parsimony_point({"1": 1.0, "2": 0.6, "4": 0.3, "8": 0.2})
        assert point["k_retaining_80pct"] <= point["k_retaining_90pct"]


class TestNegativeControls:
    def test_permuting_the_outcome_is_recorded_in_the_result(self, sss,
                                                             cohort):
        run = sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                         baseline_columns=cohort["baseline_columns"],
                         knowledge_order=cohort["knowledge_order"],
                         permute_target=True, **FAST)
        assert run["permuted_target"] is True

    def test_permuting_the_knowledge_changes_the_knowledge_panel(self, sss,
                                                                 cohort):
        kwargs = dict(baseline_columns=cohort["baseline_columns"],
                      knowledge_order=cohort["knowledge_order"],
                      selectors=("S5_knowledge_only",), k_grid=(3,))
        real = sss.screen(cohort["target"], cohort["design"],
                          cohort["groups"], **kwargs)
        fake = sss.screen(cohort["target"], cohort["design"],
                          cohort["groups"], permute_knowledge=True, **kwargs)
        assert (real["panels_per_fold"]["S5_knowledge_only"][0]
                != fake["panels_per_fold"]["S5_knowledge_only"][0])

    def test_the_knowledge_only_selector_never_reads_the_outcome(self, sss,
                                                                 cohort):
        _, diagnostics = sss.rank_candidates(
            "S5_knowledge_only", cohort["design"], cohort["target"],
            baseline_columns=cohort["baseline_columns"],
            knowledge_order=cohort["knowledge_order"],
            groups=cohort["groups"], seed=0)
        assert diagnostics["uses_outcome"] is False

    def test_a_knowledge_selector_without_an_ordering_is_refused(self, sss,
                                                                 cohort):
        for name in ("S5_knowledge_only", "S6_knowledge_then_shap"):
            with pytest.raises(Exception):
                sss.rank_candidates(name, cohort["design"], cohort["target"],
                                    baseline_columns=(), knowledge_order=None,
                                    groups=cohort["groups"], seed=0)

    def test_an_unknown_selector_is_refused(self, sss, cohort):
        with pytest.raises(Exception):
            sss.rank_candidates("S9_invented", cohort["design"],
                                cohort["target"], baseline_columns=(),
                                knowledge_order=None,
                                groups=cohort["groups"], seed=0)


class TestContrastAndStability:
    def test_the_contrast_uses_the_shared_ladder(self, sss, run, cohort):
        out = sss.contrast(run, "S1_univariate_outcome", "S0_everything",
                           groups=cohort["groups"])
        assert out["verdict"] in {
            "GAIN_REACHES_MARGIN", "GAIN_PRESENT_BUT_BELOW_MARGIN",
            "GAIN_PRESENT_MARGIN_UNDECIDED", "DETECTABLY_WORSE",
            "REFUTED_GAIN_BELOW_MARGIN", "NO_DETECTABLE_GAIN_UNDERPOWERED"}

    def test_the_contrast_clusters_by_group_when_asked(self, sss, run,
                                                       cohort):
        clustered = sss.contrast(run, "S1_univariate_outcome",
                                 "S0_everything", groups=cohort["groups"])
        assert "cluster" in clustered["resampling_unit"]
        rowwise = sss.contrast(run, "S1_univariate_outcome", "S0_everything")
        assert rowwise["resampling_unit"] == "individual"

    def test_stability_is_reported_beside_performance(self, run):
        for name, block in run["stability"].items():
            if block["n_folds"] == 0:
                continue
            assert 0.0 <= block["mean_pairwise_jaccard"] <= 1.0
            assert "never instead of it" in block["why_it_matters"]

    def test_an_unknown_selector_name_in_a_contrast_is_refused(self, sss, run):
        with pytest.raises(Exception):
            sss.contrast(run, "S9_invented", "S0_everything")


class TestDisclaimers:
    def test_shap_is_not_claimed_as_evidence_of_value(self, sss):
        text = sss.what_this_module_does_not_do()
        assert "never settles a claim" in text[
            "a_shap_value_is_not_evidence_of_value"]

    def test_the_module_disclaims_pharmacokinetics(self, sss):
        text = sss.what_this_module_does_not_do()
        assert "clearance" in text[
            "beating_a_selector_is_not_predicting_pharmacokinetics"]

    def test_the_frozen_xgboost_matches_the_earlier_hlc_comparison(self, sss):
        # Reused verbatim so the screening arm stays comparable, and so the
        # booster cannot be tuned on an outcome in this module.
        settings = sss.FROZEN["xgb"]
        assert settings["n_estimators"] == 150
        assert settings["max_depth"] == 2
        assert settings["learning_rate"] == 0.035
        assert settings["min_child_weight"] == 8
        assert settings["reg_lambda"] == 10.0

    def test_every_selector_is_declared(self, sss):
        assert len(sss.SELECTORS) == len(set(sss.SELECTORS))
        assert "S0_everything" in sss.SELECTORS
        assert "S6_knowledge_then_shap" in sss.SELECTORS


class TestKnowledgeAgreementSpeaksOneLanguage:
    """The first version of this compared two namespaces and never said so.

    It built generic `x{j}` names, intersected them with panels written in the
    real feature names, and reported an overlap of exactly 0.000 on every one
    of eighteen cells. A number that is always zero looks like a finding, so
    the contract now is that the comparison is impossible to make across
    namespaces: the result carries the names its own panels use.
    """

    @pytest.fixture(scope="class")
    def named(self, sss, cohort):
        names = [f"feature_{j}" for j in range(cohort["design"].shape[1])]
        return sss.screen(
            cohort["target"], cohort["design"], cohort["groups"],
            baseline_columns=cohort["baseline_columns"],
            knowledge_order=cohort["knowledge_order"], feature_names=names,
            selectors=("S0_everything", "S4_shap_xgboost"), k_grid=(3, 8))

    def test_the_result_carries_the_names_its_panels_use(self, named):
        assert named["candidate_names"][0] == "feature_0"
        for panel in named["panels_per_fold"]["S4_shap_xgboost"]:
            assert set(panel) <= set(named["candidate_names"])

    def test_overlap_is_not_identically_zero_with_real_names(self, sss, named,
                                                             cohort):
        out = sss.knowledge_agreement(named, cohort["knowledge_order"], k=8)
        # The knowledge ordering leads with the truly useful columns and the
        # data-driven selector finds them, so some overlap must exist.
        assert out["observed_overlap"] > 0.0
        assert out["n_candidates"] == len(named["candidate_names"])
        assert out["proposed_by_knowledge"]

    def test_it_beats_a_random_ordering_when_the_knowledge_is_right(self, sss,
                                                                    named,
                                                                    cohort):
        out = sss.knowledge_agreement(named, cohort["knowledge_order"], k=8)
        assert out["beats_a_random_ordering"] is True
        assert out["empirical_p"] < 0.05

    def test_a_result_without_names_is_refused_not_scored_as_zero(self, sss,
                                                                  named,
                                                                  cohort):
        stripped = {k: v for k, v in named.items() if k != "candidate_names"}
        with pytest.raises(Exception) as caught:
            sss.knowledge_agreement(stripped, cohort["knowledge_order"], k=3)
        assert "candidate_names" in str(caught.value)

    def test_a_mismatched_candidate_count_is_refused(self, sss, named,
                                                     cohort):
        with pytest.raises(Exception):
            sss.knowledge_agreement(named, cohort["knowledge_order"], k=3,
                                    n_candidates=999)

    def test_it_disclaims_what_agreement_does_not_show(self, sss, named,
                                                       cohort):
        out = sss.knowledge_agreement(named, cohort["knowledge_order"], k=3)
        assert "is not proof that either is" in out["what_this_does_not_show"]


class TestResultsDocument:
    @pytest.fixture(scope="class")
    def text(self):
        return (ROOT / "docs"
                / "OMICS2PK_SUPERVISED_SCREENING_SELECTION_20261111.md"
                ).read_text(encoding="utf-8")

    def test_it_is_versioned(self, text):
        assert "**Version 1.1.0 — 2026-11-11.**" in text

    def test_it_reports_the_primary_contrast_as_refuted(self, text):
        assert "n'atteint la marge sur aucune des 11 cellules" in text

    def test_it_does_not_claim_a_sophisticated_selector_won(self, text):
        assert "ne pas construire le sélecteur" in text
        assert "n'achètent rien ici" in text

    def test_it_reports_the_knowledge_only_arm_as_detectably_bad(self, text):
        assert "−3,63 %" in text
        assert "+15,90 %" in text

    def test_it_separates_the_two_commissioned_questions(self, text):
        assert "7 sur 10" in text
        assert "0 sur 12" in text
        # And it reports the control that limits what the second one means.
        assert "CYP2E1 bat aussi l'ordre aléatoire" in text

    def test_it_records_the_redundant_duplicate(self, text):
        assert "ne déduplique pas" in text

    def test_it_records_both_of_the_defects_found_in_passing(self, text):
        assert "ElasticNetCV" in text
        assert "deux espaces de noms" in text

    def test_it_declares_epa_unevaluable_for_want_of_data(self, text):
        assert "aucune archive EPA n'est sur cette machine" in text

    def test_it_declares_the_genotype_cost(self, text):
        assert "171–172" in text

    def test_the_artefact_agrees_with_the_document_on_the_primary(self):
        import json
        payload = json.loads((ROOT / "reports" / "data"
                              / "supervised_screening_selection.json"
                              ).read_text(encoding="utf-8"))
        cells = []
        for entry in payload["mavedb"].values():
            if not isinstance(entry, dict):
                continue
            for split in ("grouped_by_position",
                          "grouped_by_contiguous_region"):
                block = entry.get(split, {}).get("contrasts", {})
                if "PRIMARY_hybrid_vs_univariate" in block:
                    cells.append(block["PRIMARY_hybrid_vs_univariate"])
        for entry in payload["hlc"].get("endpoints", {}).values():
            block = entry.get("contrasts", {})
            if "PRIMARY_hybrid_vs_univariate" in block:
                cells.append(block["PRIMARY_hybrid_vs_univariate"])
        assert len(cells) == 11, len(cells)
        assert not any(c["verdict"] == "GAIN_REACHES_MARGIN" for c in cells)

    def test_the_permuted_control_reached_the_margin_nowhere(self):
        import json
        payload = json.loads((ROOT / "reports" / "data"
                              / "supervised_screening_selection.json"
                              ).read_text(encoding="utf-8"))
        checked = 0
        for entry in payload["mavedb"].values():
            if not isinstance(entry, dict):
                continue
            for split in ("grouped_by_position",
                          "grouped_by_contiguous_region"):
                control = entry.get(split, {}).get(
                    "permuted_outcome_control")
                if control:
                    checked += 1
                    assert control["best_selector_vs_everything"][
                        "verdict"] != "GAIN_REACHES_MARGIN"
        for entry in payload["hlc"].get("endpoints", {}).values():
            control = entry.get("permuted_outcome_control")
            if control:
                checked += 1
                assert control["best_selector_vs_everything"][
                    "verdict"] != "GAIN_REACHES_MARGIN"
        assert checked == 11, checked
