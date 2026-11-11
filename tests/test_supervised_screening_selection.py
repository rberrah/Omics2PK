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


class TestAbstentionIsAvailable:
    """`k = 0` means keep the baseline covariates and nothing else.

    A method that cannot choose it is not automatic, it is obedient. On HLC
    CYP2E1 every selector was DETECTABLY_WORSE than taking everything, so
    "select nothing" was the right answer and no arm could express it.
    """

    def test_zero_is_in_the_frozen_grid(self, sss):
        assert 0 in sss.FROZEN["k_grid"]

    def test_an_empty_panel_predicts_the_training_mean(self, sss, cohort):
        run = sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                         baseline_columns=(),
                         knowledge_order=cohort["knowledge_order"],
                         selectors=("S1_univariate_outcome",), k_grid=(0, 3))
        assert run["abstention_was_available"] is True
        # With no baseline and no selected column there is nothing to fit, so
        # the loss must be that of an intercept-only model: worse than k=3 on
        # a cohort that has signal, and finite rather than an error.
        curve = run["curve"]["S1_univariate_outcome"]
        assert curve["0"] > curve["3"]
        assert np.isfinite(curve["0"])

    def test_with_a_baseline_k_zero_keeps_exactly_the_baseline(self, sss,
                                                               cohort):
        run = sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                         baseline_columns=cohort["baseline_columns"],
                         knowledge_order=cohort["knowledge_order"],
                         selectors=("S1_univariate_outcome",), k_grid=(0,))
        for panel in run["panels_per_fold"]["S1_univariate_outcome"]:
            assert set(panel) == {f"x{c}"
                                  for c in cohort["baseline_columns"]}

    def test_abstention_is_counted_per_selector(self, sss, cohort):
        run = sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                         baseline_columns=cohort["baseline_columns"],
                         knowledge_order=cohort["knowledge_order"],
                         selectors=("S1_univariate_outcome",), k_grid=(0, 3))
        counts = run["n_folds_that_abstained"]
        assert "S1_univariate_outcome" in counts
        assert 0 <= counts["S1_univariate_outcome"] <= 5

    def test_a_cohort_with_real_signal_does_not_abstain_everywhere(self, sss,
                                                                   cohort):
        # The guard that matters in the other direction: abstention must not
        # become the answer to everything once it is available.
        run = sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                         baseline_columns=cohort["baseline_columns"],
                         knowledge_order=cohort["knowledge_order"],
                         selectors=("S1_univariate_outcome",),
                         k_grid=(0, 3, 8))
        assert run["n_folds_that_abstained"]["S1_univariate_outcome"] < 5

    def test_a_negative_budget_is_refused(self, sss, cohort):
        with pytest.raises(Exception):
            sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                       selectors=("S0_everything",), k_grid=(-1, 3))


class TestDeduplication:
    def test_exactly_one_of_a_near_identical_pair_survives(self, sss, cohort):
        decoy = cohort["redundant_copy_of_a_true_column"]
        twin = cohort["truth_columns"][0]
        kept, dropped = sss.drop_redundant(
            np.arange(cohort["design"].shape[1]), cohort["design"],
            keep_columns=cohort["baseline_columns"])
        # WHICH of the two survives is not a contract: they carry the same
        # information, so keeping either is correct. That exactly one survives
        # IS the contract, because a budget of k must buy k pieces of
        # information.
        assert len({decoy, twin} & set(kept.tolist())) == 1
        assert any(abs(r) >= 0.9 for _, _, r in dropped)

    def test_the_baseline_is_never_dropped(self, sss, cohort):
        kept, _ = sss.drop_redundant(
            np.arange(cohort["design"].shape[1]), cohort["design"],
            keep_columns=cohort["baseline_columns"])
        assert set(cohort["baseline_columns"]) <= set(kept.tolist())

    def test_it_drops_nothing_when_nothing_is_redundant(self, sss):
        rng = np.random.default_rng(0)
        design = rng.normal(size=(200, 6))
        kept, dropped = sss.drop_redundant(np.arange(6), design)
        assert kept.size == 6
        assert dropped == []

    def test_a_zero_variance_column_correlates_with_nothing(self, sss):
        design = np.column_stack([np.ones(100), np.arange(100.0)])
        kept, dropped = sss.drop_redundant(np.arange(2), design)
        assert kept.size == 2
        assert dropped == []

    def test_it_respects_the_ordering_it_is_given(self, sss, cohort):
        decoy = cohort["redundant_copy_of_a_true_column"]
        twin = cohort["truth_columns"][0]
        # Put the true column first and the decoy must go, not the reverse.
        order = [twin] + [j for j in range(cohort["design"].shape[1])
                          if j != twin]
        kept, _ = sss.drop_redundant(order, cohort["design"])
        assert twin in kept.tolist()
        assert decoy not in kept.tolist()

    def test_the_deduplicated_selector_reports_what_it_dropped(self, sss,
                                                               cohort):
        order, diagnostics = sss.selector_ordering(
            "S7_shap_deduplicated", cohort["design"], cohort["target"],
            baseline_columns=cohort["baseline_columns"],
            knowledge_order=cohort["knowledge_order"],
            groups=cohort["groups"], seed=0, max_k=8)
        if not diagnostics.get("available", True):
            pytest.skip("xgboost absent")
        assert diagnostics["deduplicated"] is True
        assert diagnostics["n_dropped_as_redundant"] >= 1
        assert diagnostics["redundancy_threshold"] == 0.90

    def test_the_deduplicated_selector_is_declared(self, sss):
        assert "S7_shap_deduplicated" in sss.SELECTORS


class TestMatchedRandomPools:
    """The control that separates knowledge from dimension reduction.

    Comparing 31 mechanistic probes against 30 218 candidates changes two
    things at once: the biology used to choose, and the dimension searched. A
    matched random pool holds the dimension fixed and removes the biology.
    """

    @pytest.fixture(scope="class")
    def array_like(self):
        rng = np.random.default_rng(3)
        n, p = 200, 3000
        level = rng.normal(8, 2, p)
        spread = np.exp(rng.normal(-1, 1, p))
        design = level + rng.normal(size=(n, p)) * spread
        pool = np.asarray(sorted(rng.choice(p, 20, replace=False)))
        return design, pool

    def test_it_draws_the_requested_size_without_overlap(self, sss,
                                                          array_like):
        design, pool = array_like
        drawn, report = sss.matched_random_pool(design, pool, seed=0)
        assert drawn.size == pool.size == report["drawn"]
        assert not (set(drawn.tolist()) & set(pool.tolist()))
        assert len(set(drawn.tolist())) == drawn.size

    def test_it_matches_level_and_spread_closely(self, sss, array_like):
        design, pool = array_like
        _, report = sss.matched_random_pool(design, pool, seed=0)
        assert abs(report["matched_mean_level"]
                   - report["target_mean_level"]) < 0.5
        assert abs(report["matched_mean_log_variance"]
                   - report["target_mean_log_variance"]) < 0.2

    def test_an_unmatched_draw_would_be_detectably_different(self, sss,
                                                             array_like):
        """Shows the matching is doing work, not decorating the output."""
        design, pool = array_like
        _, report = sss.matched_random_pool(design, pool, seed=0)
        rng = np.random.default_rng(0)
        naive = rng.choice(design.shape[1], size=pool.size, replace=False)
        naive_level = float(design.mean(axis=0)[naive].mean())
        matched_error = abs(report["matched_mean_level"]
                            - report["target_mean_level"])
        naive_error = abs(naive_level - report["target_mean_level"])
        # Matching must be at least as close as an unmatched draw, and on a
        # design with a wide level range it is usually much closer.
        assert matched_error <= naive_error + 1e-9

    def test_it_excludes_the_baseline(self, sss, array_like):
        design, pool = array_like
        drawn, _ = sss.matched_random_pool(design, pool, exclude=(0, 1),
                                           seed=0)
        assert 0 not in drawn.tolist()
        assert 1 not in drawn.tolist()

    def test_different_seeds_give_different_pools(self, sss, array_like):
        design, pool = array_like
        first, _ = sss.matched_random_pool(design, pool, seed=0)
        second, _ = sss.matched_random_pool(design, pool, seed=1)
        assert first.tolist() != second.tolist()

    def test_redundancy_is_reported_and_declared_unmatched(self, sss,
                                                            array_like):
        design, pool = array_like
        _, report = sss.matched_random_pool(design, pool, seed=0)
        assert "target_within_pool_redundancy" in report
        assert "matched_within_pool_redundancy" in report
        assert "REPORTED not matched" in report["matched_on"]
        assert "co-regulated" in report["why_redundancy_is_only_reported"]

    def test_the_runner_declares_the_number_of_pools(self):
        script = (ROOT / "scripts"
                  / "run_supervised_screening_selection.py").read_text(
                      encoding="utf-8")
        assert "N_MATCHED_POOLS = 40" in script
        assert "matched_pool_null" in script


class TestNonInferiority:
    """A non-significant difference is not equivalence."""

    @pytest.fixture(scope="class")
    def run(self, sss, cohort):
        return sss.screen(cohort["target"], cohort["design"],
                          cohort["groups"],
                          baseline_columns=cohort["baseline_columns"],
                          knowledge_order=cohort["knowledge_order"],
                          selectors=("S0_everything",
                                     "S1_univariate_outcome"),
                          k_grid=(0, 3, 8))

    def test_an_identical_arm_is_non_inferior_to_itself(self, sss, run,
                                                        cohort):
        out = sss.non_inferiority(run, "S1_univariate_outcome",
                                  "S1_univariate_outcome",
                                  groups=cohort["groups"])
        assert out["verdict"] == "NON_INFERIOR"
        assert out["relative_gap"] == pytest.approx(0.0)

    def test_a_much_worse_arm_is_inferior(self, sss, run, cohort):
        # S0 keeps every column including 17 noise ones, so on this cohort it
        # is substantially worse than selecting.
        out = sss.non_inferiority(run, "S0_everything",
                                  "S1_univariate_outcome",
                                  groups=cohort["groups"])
        assert out["verdict"] == "INFERIOR"

    def test_the_verdict_set_is_closed(self, sss, run, cohort):
        for candidate in ("S0_everything", "S1_univariate_outcome"):
            out = sss.non_inferiority(run, candidate, "S1_univariate_outcome",
                                       groups=cohort["groups"])
            assert out["verdict"] in {"NON_INFERIOR", "INFERIOR", "UNDECIDED"}

    def test_undecided_is_never_called_equivalence(self, sss, run, cohort):
        out = sss.non_inferiority(run, "S1_univariate_outcome",
                                  "S0_everything", groups=cohort["groups"])
        assert out["an_undecided_result_is_not_equivalence"] is True
        assert "not the same claim" in out["how_to_read"]

    def test_it_is_a_separate_test_from_the_superiority_ladder(self, sss):
        source = (ROOT / "src" / "omics2pk"
                  / "supervised_screening_selection.py").read_text(
                      encoding="utf-8")
        # The non-inferiority margin is applied to MINUS the threshold, which
        # the six-rung ladder never does.
        assert "low > -threshold" in source
        assert "high < -threshold" in source

    def test_a_tighter_margin_is_harder_to_satisfy(self, sss, run, cohort):
        wide = sss.non_inferiority(run, "S0_everything",
                                   "S1_univariate_outcome",
                                   groups=cohort["groups"], margin=0.5)
        tight = sss.non_inferiority(run, "S0_everything",
                                    "S1_univariate_outcome",
                                    groups=cohort["groups"], margin=0.001)
        assert wide["non_inferiority_margin"] > tight[
            "non_inferiority_margin"]
        assert not (wide["verdict"] == "INFERIOR"
                    and tight["verdict"] == "NON_INFERIOR")

    def test_an_unknown_selector_is_refused(self, sss, run):
        with pytest.raises(Exception):
            sss.non_inferiority(run, "S9_invented", "S0_everything")


class TestRedundancyMatchedPools:
    """The secondary control: is the advantage in WHICH genes, or in HOW they covary?

    Rejection sampling cannot build this and that was measured, not assumed: a
    pool with internal mean |r| of 0.88 is co-regulated, and 24 level-matched
    draws from an array whose typical probe pair is uncorrelated reached 0.05.
    Calling that a redundancy-matched control would have been false. The pool
    is grown to the target instead.
    """

    @pytest.fixture(scope="class")
    def modular_array(self):
        """An array with MANY co-expression modules, as a real one has."""
        rng = np.random.default_rng(0)
        n, p = 200, 2500
        level = rng.normal(8, 2, p)
        spread = np.exp(rng.normal(-1, 1, p))
        design = level + rng.normal(size=(n, p)) * spread
        for _ in range(30):
            members = rng.choice(p, 50, replace=False)
            design[:, members] += rng.normal(size=(n, 1)) * 1.5
        pool = np.asarray(sorted(rng.choice(p, 20, replace=False)))
        design[:, pool] += rng.normal(size=(n, 1)) * 1.5
        return design, pool

    @pytest.fixture(scope="class")
    def single_module_array(self):
        """An array where the pool is the ONLY correlated module."""
        rng = np.random.default_rng(1)
        n, p = 200, 2000
        design = rng.normal(8, 1, size=(n, p))
        pool = np.asarray(sorted(rng.choice(p, 20, replace=False)))
        design[:, pool] += rng.normal(size=(n, 1)) * 2.0
        return design, pool

    def test_it_matches_redundancy_when_modules_exist(self, sss,
                                                       modular_array):
        design, pool = modular_array
        _, report = sss.matched_random_pool(design, pool, seed=2,
                                            match_redundancy=True)
        target = report["target_within_pool_redundancy"]
        achieved = report["matched_within_pool_redundancy"]
        assert target > 0.3, "the fixture must plant a co-regulated pool"
        assert abs(achieved - target) < 0.1 * target

    def test_level_matching_alone_does_not_match_redundancy(self, sss,
                                                             modular_array):
        design, pool = modular_array
        _, report = sss.matched_random_pool(design, pool, seed=2,
                                            match_redundancy=False)
        # This is the whole reason the secondary control needs its own
        # sampler: the primary control's draws are NOT redundancy-matched.
        assert (report["matched_within_pool_redundancy"]
                < 0.5 * report["target_within_pool_redundancy"])
        assert "REPORTED not matched" in report["matched_on"]

    def test_it_reports_a_shortfall_instead_of_claiming_a_match(
            self, sss, single_module_array):
        """When no comparable module exists the control is not achievable.

        That is a fact about the array, not a bug, and the achieved number
        must show it rather than the report asserting success.
        """
        design, pool = single_module_array
        _, report = sss.matched_random_pool(design, pool, seed=2,
                                            match_redundancy=True)
        assert report["redundancy_matched"] is True
        assert (report["matched_within_pool_redundancy"]
                < 0.5 * report["target_within_pool_redundancy"])

    def test_it_declares_that_level_and_spread_drift(self, sss,
                                                     modular_array):
        design, pool = modular_array
        _, report = sss.matched_random_pool(design, pool, seed=2,
                                            match_redundancy=True)
        assert "level and spread then drift" in report["matched_on"]

    def test_the_grown_pool_has_the_requested_size_and_no_duplicates(
            self, sss, modular_array):
        design, pool = modular_array
        drawn, _ = sss.matched_random_pool(design, pool, seed=2,
                                           match_redundancy=True)
        assert drawn.size == pool.size
        assert len(set(drawn.tolist())) == drawn.size

    def test_the_grown_pool_respects_the_exclusion(self, sss, modular_array):
        design, pool = modular_array
        drawn, _ = sss.matched_random_pool(design, pool, exclude=(0, 1),
                                           seed=2, match_redundancy=True)
        assert not ({0, 1} & set(drawn.tolist()))


class TestInductivePoolMatching:
    """A pool matched on every donor is transductive, not inductive."""

    def test_the_matched_selector_draws_from_training_rows_only(self, sss,
                                                                cohort):
        seen = []
        original = sss.matched_random_pool

        def spy(design, pool, **kwargs):
            seen.append(design.shape[0])
            return original(design, pool, **kwargs)

        import unittest.mock as mock
        pool = np.asarray(cohort["truth_columns"], dtype=int)
        with mock.patch.object(sss, "matched_random_pool", spy):
            sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                       baseline_columns=cohort["baseline_columns"],
                       knowledge_order=cohort["knowledge_order"],
                       selectors=("S9_matched_pool_then_univariate",),
                       k_grid=(0, 2), pool_size=pool.size,
                       match_target_pool=pool)
        assert seen, "the matched selector never drew a pool"
        # Never the whole cohort: every draw saw a strict training subset.
        assert max(seen) < cohort["target"].size

    def test_the_result_records_that_matching_was_inductive(self, sss,
                                                             cohort):
        pool = np.asarray(cohort["truth_columns"], dtype=int)
        run = sss.screen(cohort["target"], cohort["design"], cohort["groups"],
                         baseline_columns=cohort["baseline_columns"],
                         knowledge_order=cohort["knowledge_order"],
                         selectors=("S9_matched_pool_then_univariate",),
                         k_grid=(0, 2), pool_size=pool.size,
                         match_target_pool=pool)
        assert run["pool_matching_was_inductive"] is True
        assert run["pool_redundancy_matched"] is False

    def test_the_matched_selector_without_a_target_pool_is_refused(self, sss,
                                                                    cohort):
        with pytest.raises(Exception):
            sss.rank_candidates("S9_matched_pool_then_univariate",
                                cohort["design"], cohort["target"],
                                baseline_columns=(), knowledge_order=None,
                                groups=cohort["groups"], seed=0)

    def test_the_matched_selector_is_declared(self, sss):
        assert "S9_matched_pool_then_univariate" in sss.SELECTORS

    def test_the_confirmatory_count_can_clear_bonferroni(self, sss):
        """40 pools cannot; 1000 can, and the arithmetic is the reason."""
        exploratory = sss.FROZEN["n_matched_pools"]
        confirmatory = sss.FROZEN["n_matched_pools_confirmatory"]
        bonferroni = 0.05 / 3  # three prespecified endpoints
        assert 1.0 / (1 + exploratory) > bonferroni
        assert 1.0 / (1 + confirmatory) < bonferroni


class TestFullArrayArtefact:
    """Guards on the high-dimension artefact, including the project's own rules.

    The pre-registration declared one refutation condition that fired: the
    variance filter costs a margin-reaching gain on CYP2D6, so the arms that
    need it are not interpretable there. That has to stay visible in the
    artefact, not only in prose.
    """

    @pytest.fixture(scope="class")
    def payload(self):
        import json
        path = (ROOT / "reports" / "data"
                / "supervised_screening_fullarray.json")
        if not path.is_file():
            pytest.skip("the high-dimension artefact has not been generated")
        return json.loads(path.read_text(encoding="utf-8"))

    def test_no_numeric_list_of_ten_or_more_entries(self, payload):
        """The project's rule against per-row data in a committed artefact.

        It exists to stop per-individual material entering the repository.
        k_grid is not that, but a guard that holds only when convenient is not
        a guard, so the declared grids are written as text.
        """
        offenders = []

        def walk(node, path):
            if isinstance(node, list):
                if node and all(isinstance(v, (int, float)) for v in node) \
                        and len(node) >= 10:
                    offenders.append((path, len(node)))
                for index, value in enumerate(node):
                    walk(value, f"{path}[{index}]")
            elif isinstance(node, dict):
                for key, value in node.items():
                    walk(value, f"{path}.{key}")

        walk(payload, "")
        assert not offenders, offenders

    def test_it_carries_no_donor_level_material(self, payload):
        import json
        raw = json.dumps(payload)
        assert payload["no_individual_level_data_in_this_artefact"] is True
        # 30 596 probe names per endpoint made a 6 MB artefact for a list
        # nothing reads; the count is kept instead.
        assert len(raw) < 2_000_000, len(raw)

    def test_it_names_the_protocol_it_followed(self, payload):
        block = payload["hlc_full_array"]
        assert block["protocol"] == \
            "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
        assert block["prefilter_top_variance"] == 2000
        assert block["mechanistic_panel_size"] == 31

    def test_the_completeness_filter_shrank_the_mechanistic_pool(self,
                                                                 payload):
        # 6 of 31 mechanistic probes are not finite across donors. A knowledge
        # pool that shrinks in silence would invalidate the arm that uses it.
        for entry in payload["hlc_full_array"]["endpoints"].values():
            assert entry["n_mechanistic_probes_dropped_as_incomplete"] >= 1
            assert entry["pool_size"] < block_size(payload)

    def test_the_prefilter_cost_is_reported_for_every_endpoint(self, payload):
        for name, entry in payload["hlc_full_array"]["endpoints"].items():
            cost = entry["cost_of_the_variance_prefilter"]
            assert "S1_univariate_outcome" in cost, name
            assert "verdict" in cost["S1_univariate_outcome"], name

    def test_the_refutation_condition_that_fired_is_visible(self, payload):
        """CYP2D6: the whole array beats the filtered array by ~14 %."""
        cost = payload["hlc_full_array"]["endpoints"]["CYP2D6_activity"][
            "cost_of_the_variance_prefilter"]["S1_univariate_outcome"]
        assert cost["verdict"] == "GAIN_REACHES_MARGIN"
        assert cost["relative_gain"] > 0.10

    def test_the_matched_pool_control_ran_on_every_endpoint(self, payload):
        for name, entry in payload["hlc_full_array"]["endpoints"].items():
            pools = entry["matched_random_pools"]
            assert pools["n_pools"] >= 30, name
            assert "empirical_p" in pools, name

    def test_the_empirical_p_floor_is_respected(self, payload):
        """No p may be below 1/(1+n_pools); the floor is arithmetic."""
        for name, entry in payload["hlc_full_array"]["endpoints"].items():
            pools = entry["matched_random_pools"]
            floor = 1.0 / (1 + pools["n_pools"])
            assert pools["empirical_p"] >= floor - 1e-12, name

    def test_the_permuted_control_reached_the_margin_nowhere(self, payload):
        for name, entry in payload["hlc_full_array"]["endpoints"].items():
            control = entry["permuted_outcome_control"][
                "best_selector_vs_everything"]
            assert control["verdict"] != "GAIN_REACHES_MARGIN", name

    def test_abstention_was_available_and_used(self, payload):
        negative = payload["hlc_full_array"]["endpoints"]["CYP2E1_activity"]
        block = negative["variance_prefiltered"]
        assert block["abstention_was_available"] is True
        # The tree-based arms abstain fully on the endpoint where selecting
        # hurt; the univariate ones do not, which is reported, not hidden.
        assert block["n_folds_that_abstained"]["S4_shap_xgboost"] >= 3


def block_size(payload):
    return payload["hlc_full_array"]["mechanistic_panel_size"]


class TestVarianceProfileDiagnostic:
    """The secondary diagnostic added to protocol 1b, and its own limits.

    It exists because the variance filter cost +13.98 % on CYP2D6, the very
    endpoint where the mechanistic pool beat matched chance. The first reading
    -- that the pool wins by bypassing a filter the other arms suffer -- is
    false for the matched comparison, because those pools are matched ON log
    variance. These tests pin that down so the claim cannot drift back.
    """

    @pytest.fixture(scope="class")
    def low_variance_pool(self):
        rng = np.random.default_rng(0)
        n, p = 200, 4000
        spread = np.exp(rng.normal(-1, 1.2, p))
        design = rng.normal(8, 1, size=(n, p)) * spread
        order = np.argsort(design.var(axis=0))
        pool = np.asarray(sorted(rng.choice(order[:1000], 20, replace=False)))
        return design, pool

    def test_matched_pools_cannot_differ_in_variance_by_construction(
            self, sss, low_variance_pool):
        """The measurement behind the correction, not an assertion about it."""
        design, pool = low_variance_pool
        matched = [sss.matched_random_pool(design, pool, exclude=(0, 1),
                                           seed=100 + i)[0]
                   for i in range(12)]
        profile = sss.variance_profile(design, pool, prefilter_n=2000,
                                       matched_pools=matched, exclude=(0, 1))
        assert abs(profile["delta_median_variance_percentile"]) < 2.0
        assert profile["pool_is_lower_variance_than_its_matches"] is False

    def test_it_detects_a_pool_the_filter_would_discard(self, sss,
                                                         low_variance_pool):
        design, pool = low_variance_pool
        profile = sss.variance_profile(design, pool, prefilter_n=2000,
                                       exclude=(0, 1))
        assert profile["pool"]["share_excluded_by_the_filter"] > 0.5
        assert profile["pool"]["median_variance_percentile"] < 50.0

    def test_a_high_variance_pool_survives_the_filter(self, sss,
                                                       low_variance_pool):
        design, _ = low_variance_pool
        order = np.argsort(design.var(axis=0))[::-1]
        pool = np.asarray(sorted(order[:20]))
        profile = sss.variance_profile(design, pool, prefilter_n=2000,
                                       exclude=(0, 1))
        assert profile["pool"]["share_excluded_by_the_filter"] == 0.0
        assert profile["pool"]["median_variance_percentile"] > 90.0

    def test_it_reads_no_outcome(self, sss, low_variance_pool):
        design, pool = low_variance_pool
        profile = sss.variance_profile(design, pool, prefilter_n=2000,
                                       exclude=(0, 1))
        assert profile["this_reads_no_outcome"] is True
        # The signature takes no target, which is the real guarantee.
        import inspect
        parameters = inspect.signature(sss.variance_profile).parameters
        assert "target" not in parameters
        assert "activity" not in parameters

    def test_it_states_what_it_cannot_show(self, sss, low_variance_pool):
        design, pool = low_variance_pool
        profile = sss.variance_profile(design, pool, prefilter_n=2000,
                                       exclude=(0, 1))
        assert "cannot explain the mechanistic pool's advantage over THEM" \
            in profile["what_it_cannot_show"]

    def test_the_predictive_split_partitions_the_pool(self, sss,
                                                       low_variance_pool):
        design, pool = low_variance_pool
        kept, dropped = sss.split_pool_by_filter(design, pool,
                                                 prefilter_n=2000,
                                                 exclude=(0, 1))
        assert kept.size + dropped.size == pool.size
        assert not (set(kept.tolist()) & set(dropped.tolist()))
        assert set(kept.tolist()) | set(dropped.tolist()) == set(pool.tolist())

    def test_the_split_reads_no_outcome_either(self, sss):
        import inspect
        parameters = inspect.signature(sss.split_pool_by_filter).parameters
        assert "target" not in parameters

    def test_a_wide_filter_keeps_everything(self, sss, low_variance_pool):
        design, pool = low_variance_pool
        kept, dropped = sss.split_pool_by_filter(
            design, pool, prefilter_n=design.shape[1], exclude=(0, 1))
        assert dropped.size == 0
        assert kept.size == pool.size

    def test_the_amendment_is_recorded_as_prior_to_execution(self):
        """The chronology table, not a sentence in the banner.

        A banner claiming "nothing here postdates a result" stops being true
        the moment results are added, and it was quietly false for one
        revision. A per-commit table cannot go stale that way: each row says
        what was declared and when, and the rows above the results sections
        are the ones that have to predate them.
        """
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "Chronologie" in text
        for commit in ("b1b4f23", "6e3549f", "fcd2a2a", "83aa8c2"):
            assert commit in text, commit
        # The two facts that make the chronology meaningful.
        assert "arr\u00eat\u00e9 et sa sortie partielle supprim\u00e9e sans \u00eatre lue" in text
        assert "**avant** l'ex\u00e9cution de 1b" in text
        assert "le protocole 1 n'est pas relanc\u00e9 pour lui" in text
        # And it records that the first reading of the problem was wrong.
        assert "Cette lecture est fausse pour la comparaison principale" in text

    def test_the_banner_no_longer_claims_what_became_false(self):
        """It said no result of any protocol had been examined. Two have."""
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "Aucun r\u00e9sultat d'aucun protocole" not in text
        assert "PROTOCOLES 1 (\u00a75) ET 1b (\u00a76) EX\u00c9CUT\u00c9S" in text
        # Whatever has not run must still say so. Protocol 2 reached its
        # declared NOT_EVALUABLE verdict without any download, and protocol 3
        # is still wholly prior.
        assert "Protocole 3 d\u00e9clar\u00e9, aucune donn\u00e9e t\u00e9l\u00e9charg\u00e9e" in text


class TestCyp2d6NetworkAblation:
    """The exploratory ablation that bounded the CYP2D6 result.

    Protocol 1b found the mechanistic pool beating 0 of 1000 matched pools,
    and the panel selected inside the folds contained POR, the obligate
    electron donor of microsomal CYPs. Biologically coherent -- and the
    ablation showed it adds nothing. These guards keep that bound visible.
    """

    @pytest.fixture(scope="class")
    def payload(self):
        import json
        path = ROOT / "reports" / "data" / "cyp2d6_network_ablation.json"
        if not path.is_file():
            pytest.skip("the ablation artefact has not been generated")
        return json.loads(path.read_text(encoding="utf-8"))

    def test_it_is_labelled_exploratory_in_the_artefact(self, payload):
        assert payload["status"] == "EXPLORATORY"
        assert payload["not_part_of_any_frozen_protocol"] is True
        assert "same cohort and the same donors" in payload["why_exploratory"]

    def test_the_single_transcript_arm_exists_and_is_one_probe(self, payload):
        arm = payload["arms"]["CYP2D6_alone"]
        assert arm["n_probes"] == 1
        assert arm["genes"] == ["CYP2D6"]

    def test_the_pool_does_not_beat_the_single_transcript(self, payload):
        out = payload["contrasts"]["full_pool_vs_CYP2D6_alone"]
        assert out["verdict"] == "REFUTED_GAIN_BELOW_MARGIN"
        assert out["relative_gain"] < 0.0

    def test_the_single_transcript_is_non_inferior_to_the_pool(self, payload):
        out = payload["non_inferiority"][
            "CYP2D6_alone_non_inferior_to_full_pool"]
        assert out["verdict"] == "NON_INFERIOR"

    def test_por_adds_nothing(self, payload):
        out = payload["contrasts"]["CYP2D6_plus_POR_vs_CYP2D6_alone"]
        assert out["verdict"] == "REFUTED_GAIN_BELOW_MARGIN"

    def test_the_network_without_the_target_carries_nothing_detectable(
            self, payload):
        out = payload["contrasts"]["pool_without_CYP2D6_vs_baseline_only"]
        assert out["verdict"] != "GAIN_REACHES_MARGIN"
        assert out["interval"][0] < 0.0 < out["interval"][1]

    def test_the_transcript_still_beats_age_and_sex_by_a_margin(self, payload):
        out = payload["contrasts"]["CYP2D6_alone_vs_baseline_only"]
        assert out["verdict"] == "GAIN_REACHES_MARGIN"

    def test_it_records_why_the_ablation_followed(self, payload):
        prior = payload["from_protocol_1b"]
        assert prior["folds_that_chose_k_equal_1"] == 2
        assert prior["folds_containing_POR"] == 2
        assert prior["folds_containing_CYP2D6"] == 5
        assert "not evidence of added information" in \
            prior["why_this_ablation_follows"]

    def test_the_document_retracts_the_word_confirmatory(self):
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "\u00ab Confirmatoire \u00bb est retir\u00e9" in text
        assert "AJOUT POST\u00c9RIEUR" in text
        # And the declared protocol text is NOT rewritten, only annotated.
        assert "## 1bis. Protocole 1b \u2014 la r\u00e9plication confirmatoire" in text
        assert "exige de **nouveaux donneurs**" in text

    def test_the_document_states_the_result_is_one_transcript(self):
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "est, en entier, le transcrit CYP2D6" in text
        # And keeps the one thing knowledge demonstrably did.
        assert "transcrit unique conserver quand un filtre" in text

    def test_the_per_fold_correction_is_recorded(self):
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "Compt\u00e9 en PLIS, pas en sondes" in text
        assert "CYP2D6 seul**" in text


class TestComtVerdict:
    """Protocol 2 came back NOT_EVALUABLE, by its own declared condition.

    The temptation it was written to resist: take protein abundance as the
    endpoint and RNA as the predictor. That is predicting one omic from
    another, which is a different question and says nothing about choosing a
    measurement in order to predict a FUNCTION.
    """

    @pytest.fixture(scope="class")
    def text(self):
        return (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")

    def test_the_verdict_is_recorded(self, text):
        assert "PROTOCOLE 2 : `NOT_EVALUABLE`" in text
        assert "## 7. Protocole 2 \u2014 `NOT_EVALUABLE`" in text

    def test_it_is_the_declared_condition_not_a_revision(self, text):
        assert "verdict d\u00e9clar\u00e9 en \u00a72.2, pas une r\u00e9vision" in text
        # The condition itself must still be readable in §2.2.
        assert "Un jeu d'activit\u00e9 COMT doit exister dans le m\u00eame d\u00e9p\u00f4t" in text

    def test_the_four_score_sets_are_named(self, text):
        for urn in ("00000661-a-1", "00000661-b-1", "00000661-c-1",
                    "00000661-d-1"):
            assert urn in text, urn
        assert "Les quatre sont des mesures d'abondance" in text

    def test_it_refuses_both_substitutions(self, text):
        # Protein abundance as endpoint, and polysome occupancy as a proxy.
        assert "remplacer la question par" in text
        assert "n'est pas promue en proxy d'activit\u00e9" in text

    def test_the_inventory_is_reported_with_numbers(self, text):
        assert "| g\u00e8nes avec **les deux** | **4** |" in text
        assert "ARAF et BRAF, et eux seuls" in text

    def test_the_weaker_question_is_named_as_weaker(self, text):
        assert "m\u00eame couche omique dans trois conditions" in text
        assert "version plus faible de la question" in text

    def test_no_further_data_was_downloaded(self, text):
        assert "Aucune donn\u00e9e ARAF ou BRAF suppl\u00e9mentaire n'est "\
            "t\u00e9l\u00e9charg\u00e9e" in text

    def test_the_decision_is_left_open_rather_than_taken(self, text):
        assert "le choix n'est pas le mien" in text
        assert "Deux voies possibles" in text
