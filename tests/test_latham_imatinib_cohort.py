"""Non-regression for the extraction bug that nearly produced a false claim.

The published table marks a missing value with `U+2212 MINUS SIGN`. A
completeness check built on a hand-written token list read three absent
concentrations as present, and the conclusion drawn was that the article
contradicted itself. It did not. Donor GKP has no CYP3A4, CYP3A5 or total
CYP3A concentration, which is exactly what the article's text says.

In order of how much damage their absence would do:

1. **A dash is missing and a zero is not.** -> `TestDashIsNotZero`.
2. **GKP is excluded by the data, never by a choice.** -> `TestGkpRegression`.
3. **n is per block, not per table.** -> `TestNIsPerBlock`.
4. **Everything learned is learned on training rows only**, because at
   thirteen donors a singleton category is absent from training in the fold
   that holds its carrier out. -> `TestEncodingLearnsInFoldOnly`.
"""
from __future__ import annotations

import pathlib

import pytest

np = pytest.importorskip("numpy")

ROOT = pathlib.Path(__file__).resolve().parents[1]
TABLE = pathlib.Path(r"C:\Users\abdel\Documents\Omics"
                     r"\Latham2025_imatinib_PHH"
                     r"\latham2025_table4_phh_per_donor.csv")
MINUS = "\u2212"


@pytest.fixture(scope="module")
def latham():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import latham_imatinib_cohort
    return latham_imatinib_cohort


@pytest.fixture(scope="module")
def loaded(latham):
    if not TABLE.is_file():
        pytest.skip("the per-donor table is not present on this machine")
    return latham.load_table4(TABLE)


class TestDashIsNotZero:
    @pytest.mark.parametrize("dash", [
        "-",            # hyphen-minus
        "\u2212",       # MINUS SIGN, the one that caused the bug
        "\u2013",       # en dash
        "\u2014",       # em dash
        "\u2012",       # figure dash
        "\u2011",       # non-breaking hyphen
        "\u2010",       # hyphen
    ])
    def test_every_dash_is_missing(self, latham, dash):
        assert latham.is_missing(dash) is True
        assert latham.parse_number(dash) is None

    def test_zero_is_a_measurement_not_a_missing_value(self, latham):
        assert latham.is_missing("0") is False
        assert latham.parse_number("0") == 0.0
        assert latham.parse_number("0.0") == 0.0

    def test_a_missing_value_never_becomes_zero(self, latham):
        # The failure mode this guards: a dash coerced to 0.0 would read as a
        # measured concentration of zero, which is a different claim.
        for dash in ("-", MINUS, "\u2013"):
            assert latham.parse_number(dash) is not 0.0  # noqa: F632
            assert latham.parse_number(dash) is None

    def test_an_unrecognised_token_raises_rather_than_vanishing(self, latham):
        for token in ("~1.3", "<0.5", "about 2", "1.3 or 1.4"):
            with pytest.raises(Exception):
                latham.parse_number(token)

    def test_declared_words_are_missing(self, latham):
        for word in ("", "NA", "n/a", "ND", "NR", "."):
            assert latham.is_missing(word) is True

    def test_the_bug_character_is_named_in_the_source(self, latham):
        source = (ROOT / "src" / "omics2pk"
                  / "latham_imatinib_cohort.py").read_text(encoding="utf-8")
        assert "MINUS_SIGN" in source
        assert "U+2212" in source
        # And the fix is by category, not by a longer hand-written list.
        assert "unicodedata.category" in source


class TestGkpRegression:
    def test_gkp_has_no_cyp3a_concentrations(self, loaded):
        parsed = loaded["parsed"]
        i = parsed["donor"].index("GKP")
        assert parsed["protein_cyp3a4"][i] is None
        assert parsed["protein_cyp3a5"][i] is None
        assert parsed["protein_cyp3a_total"][i] is None

    def test_gkp_does_have_its_cyp2c8_concentration(self, loaded):
        parsed = loaded["parsed"]
        i = parsed["donor"].index("GKP")
        assert parsed["protein_cyp2c8"][i] == pytest.approx(1.3)

    def test_gkp_has_activities_and_a_clearance(self, loaded):
        """So it stays available for the arms that do not need the panel."""
        parsed = loaded["parsed"]
        i = parsed["donor"].index("GKP")
        for role in ("activity_cyp2c8", "activity_cyp3a_midazolam",
                     "clearance", "metabolite", "age", "bmi"):
            assert parsed[role][i] is not None, role

    def test_gkp_is_the_only_donor_missing_anything(self, loaded):
        missing = {d for block in loaded["completeness"].values()
                   for d in block["donors_missing"]}
        assert missing == {"GKP"}

    def test_the_primary_cohort_is_thirteen_and_excludes_gkp(self, latham,
                                                             loaded):
        rows, donors = latham.primary_cohort(loaded)
        assert len(rows) == 13
        assert "GKP" not in donors
        assert len(set(donors)) == 13

    def test_gkp_is_excluded_by_the_data_not_by_a_list(self, latham):
        """No donor identifier is hard-coded as an exclusion anywhere."""
        source = (ROOT / "src" / "omics2pk"
                  / "latham_imatinib_cohort.py").read_text(encoding="utf-8")
        body = source.split("def what_this_module_does_not_do")[0]
        code = "\n".join(line for line in body.splitlines()
                         if not line.strip().startswith("#"))
        # GKP may be named in prose, never in a filter.
        assert 'primary_cohort' in code
        assert '"GKP"' not in code
        assert "'GKP'" not in code


class TestNIsPerBlock:
    def test_every_block_reports_its_own_n(self, loaded):
        expected = {"clinical": 14, "activities": 14, "proteomics": 13,
                    "clearance": 14, "metabolite": 14}
        actual = {k: v["n_complete"]
                  for k, v in loaded["completeness"].items()}
        assert actual == expected

    def test_the_article_text_is_reconciled_not_contradicted(self, loaded):
        # 14 for activity, 13 for protein quantitation: exactly the text.
        assert loaded["completeness"]["activities"]["n_complete"] == 14
        assert loaded["completeness"]["proteomics"]["n_complete"] == 13
        assert "is correct" in loaded["why_n_is_per_block"]

    def test_the_table_has_fourteen_rows(self, loaded):
        assert loaded["n_rows"] == 14

    def test_every_numeric_column_converted_or_is_declared_missing(self,
                                                                   loaded):
        """The check the first extraction lacked: verify after conversion."""
        parsed = loaded["parsed"]
        numeric = ("age", "bmi", "activity_cyp2c8",
                   "activity_cyp3a_midazolam",
                   "activity_cyp3a_testosterone", "protein_cyp2c8",
                   "protein_cyp3a4", "protein_cyp3a5",
                   "protein_cyp3a_total", "clearance", "metabolite")
        for role in numeric:
            for value in parsed[role]:
                assert value is None or isinstance(value, float), role

    def test_no_role_resolved_to_the_same_column_as_another(self, loaded):
        indices = list(loaded["columns"].values())
        assert len(set(indices)) == len(indices)

    def test_cyp3a_total_did_not_capture_cyp3a4(self, loaded):
        columns = loaded["columns"]
        assert columns["protein_cyp3a_total"] != columns["protein_cyp3a4"]
        assert columns["protein_cyp3a_total"] != columns["protein_cyp3a5"]


class TestEncodingLearnsInFoldOnly:
    def test_standardisation_uses_training_rows_only(self, latham, loaded):
        rows, _ = latham.primary_cohort(loaded)
        full = latham.encode_blocks(loaded, rows, train_rows=rows)
        partial = latham.encode_blocks(loaded, rows, train_rows=rows[:-1])
        # Dropping one donor from the training set must change the centring,
        # or the standardisation is not being learned in-fold at all.
        assert not np.allclose(full["clinical"], partial["clinical"])

    def test_the_cyp3a5_ordering_is_fixed_a_priori(self, latham):
        assert latham.CYP3A5_ORDINAL == {"PM": 0.0, "IM": 1.0, "NM": 2.0}

    def test_a_singleton_category_gives_a_constant_column_not_an_error(
            self, latham, loaded):
        """The fold that holds the only *1/*15 carrier out must still run."""
        parsed = loaded["parsed"]
        rows, _ = latham.primary_cohort(loaded)
        carriers = [i for i in rows
                    if parsed["genotype_cyp3a4"][i] != "*1/*1"]
        if not carriers:
            pytest.skip("no non-reference CYP3A4 carrier in the cohort")
        train = [i for i in rows if i not in carriers]
        encoded = latham.encode_blocks(loaded, rows, train_rows=train)
        column = encoded["pgx"][:, 1]
        # Constant across the training rows, which is what ridge needs to see.
        train_positions = [rows.index(i) for i in train]
        assert len(set(column[train_positions].tolist())) == 1

    def test_an_unmapped_metaboliser_status_raises(self, latham, loaded):
        import copy
        broken = copy.deepcopy(loaded)
        broken["parsed"]["metabolizer_cyp3a5"][0] = "UNKNOWN"
        rows, _ = latham.primary_cohort(broken)
        with pytest.raises(Exception) as caught:
            latham.encode_blocks(broken, rows, train_rows=rows)
        assert "CYP3A5_ORDINAL" in str(caught.value)

    def test_it_declares_what_is_learned_and_what_is_fixed(self, latham,
                                                            loaded):
        rows, _ = latham.primary_cohort(loaded)
        encoded = latham.encode_blocks(loaded, rows, train_rows=rows)
        assert encoded["learned_in_fold"] == ("centring and scaling only",)
        assert any("PM<IM<NM" in s for s in encoded["fixed_a_priori"])

    def test_training_rows_must_be_a_subset(self, latham, loaded):
        rows, _ = latham.primary_cohort(loaded)
        with pytest.raises(Exception):
            latham.encode_blocks(loaded, rows, train_rows=[999])


class TestDisclaimers:
    def test_it_states_that_a_dash_is_not_a_zero(self, latham):
        text = latham.what_this_module_does_not_do()
        assert "U+2212" in text["a_dash_is_not_a_zero"]
        assert "contradicted itself" in text["a_dash_is_not_a_zero"]

    def test_it_states_the_endpoint_is_not_a_clearance(self, latham):
        text = latham.what_this_module_does_not_do()
        assert "no renal excretion" in text[
            "an_in_vitro_clearance_is_not_a_clearance"]

    def test_it_states_a_singleton_is_not_estimable(self, latham):
        text = latham.what_this_module_does_not_do()
        assert "does not make the category" in text[
            "a_singleton_category_is_not_an_effect"]


class TestAbstentionControl:
    """`k = 0` must be chosen when it IS the right answer.

    My §5.1 said that no fold choosing k = 0 would be an implementation
    defect. That was wrong: it may simply mean the model keeps the extra
    variables, which is an answer and not a fault. Abstention is therefore
    verified on a synthetic cohort where abstaining is correct -- the extra
    block is pure noise -- and never by observing whether it happens on real
    data.
    """

    @pytest.fixture(scope="class")
    def screening(self):
        import sys
        sys.path.insert(0, str(ROOT / "src"))
        from omics2pk import supervised_screening_selection
        return supervised_screening_selection

    def test_abstention_is_chosen_when_the_extra_block_is_noise(self,
                                                                screening):
        """Thirteen rows, a real baseline, and candidates that carry nothing."""
        rng = np.random.default_rng(0)
        n = 13
        baseline = rng.normal(size=(n, 2))
        noise = rng.normal(size=(n, 3))
        target = 1.4 * baseline[:, 0] - 0.9 * baseline[:, 1] \
            + rng.normal(0, 0.25, n)
        design = np.column_stack([baseline, noise])
        out = screening.pool_arm_loss(
            target, design, np.arange(n), baseline_columns=(0, 1),
            pool_for_fold=lambda rows, fold: np.array([2, 3, 4]),
            k_grid=(0, 1, 2, 3), n_outer=n, n_inner=4,
            min_scored=13, min_train=6)
        # Leave-one-donor-out at this size, with nothing to gain: most folds
        # must keep the baseline alone.
        assert out["n_folds_that_abstained"] >= n // 2, \
            out["k_chosen_per_fold"]

    def test_abstention_is_not_chosen_when_the_block_carries_signal(
            self, screening):
        """The other direction: abstaining must not become the answer to all."""
        rng = np.random.default_rng(1)
        n = 13
        baseline = rng.normal(size=(n, 2))
        useful = rng.normal(size=(n, 3))
        target = 0.3 * baseline[:, 0] + 2.5 * useful[:, 0] \
            - 2.0 * useful[:, 1] + rng.normal(0, 0.2, n)
        design = np.column_stack([baseline, useful])
        out = screening.pool_arm_loss(
            target, design, np.arange(n), baseline_columns=(0, 1),
            pool_for_fold=lambda rows, fold: np.array([2, 3, 4]),
            k_grid=(0, 1, 2, 3), n_outer=n, n_inner=4,
            min_scored=13, min_train=6)
        assert out["n_folds_that_abstained"] < n, out["k_chosen_per_fold"]

    def test_the_protocol_records_that_the_old_rule_was_wrong(self):
        text = (ROOT / "docs"
                / "OMICS2PK_LATHAM_IMATINIB_PREREGISTRATION_20261112.md"
                ).read_text(encoding="utf-8")
        assert "n'est pas un défaut" in text
        assert "contrôle synthétique" in text
        # And the old clause is struck rather than deleted.
        assert "défaut d'implémentation.~~" in text

    def test_the_protocol_keeps_the_wrong_first_interpretation_visible(self):
        text = (ROOT / "docs"
                / "OMICS2PK_LATHAM_IMATINIB_PREREGISTRATION_20261112.md"
                ).read_text(encoding="utf-8")
        assert "FAUX, voir §7.4" in text
        assert "effacer une interprétation fausse" in text
        assert "U+2212 MINUS SIGN" in text
