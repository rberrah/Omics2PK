"""Contract tests for the Mueller 2026 potato-meal source qualification.

Offline: they exercise the pinned provenance constants and the pure logic,
never the network and never the workbook. The published identifiers are
pinned here so a change to them is a deliberate, reviewable edit.

The rule these tests exist to defend: this cohort administers no drug, so it
can never support H1 or any pharmacokinetic claim, however favourable its
biomarker behaviour turns out to be.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "qualify_muller2026_potato_meal_source.py"
)


def _load():
    spec = importlib.util.spec_from_file_location("qualify_muller2026", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def mod():
    return _load()


def test_provenance_is_pinned(mod):
    assert mod.EXPECTED_SHA256 == (
        "68d53071db79f961a63e433137a8a3e7dbd40daf8c8ad0f9f7f1b268bb3a8e07"
    )
    assert mod.EXPECTED_BYTES == 20965
    assert "JulianPeterMueller" in mod.DATA_URL
    assert mod.DATA_URL.startswith("https://")


def test_checksum_guard_rejects_a_wrong_file(mod, tmp_path):
    decoy = tmp_path / "not_the_workbook.xlsx"
    decoy.write_bytes(b"not the published workbook")
    with pytest.raises(mod.QualificationError, match="checksum mismatch"):
        mod.verify_source(decoy)


def test_sampling_structure_is_the_published_one(mod):
    assert mod.TIMEPOINTS_H == (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 29.0)
    assert len(mod.TIMEPOINTS_H) == 9
    assert mod.N_PARTICIPANTS == 20


def test_three_analyte_blocks_at_distinct_offsets(mod):
    assert set(mod.PLASMA_BLOCKS) == {
        "solanidine_ng_ml",
        "oh4_solanidine_ng_ml",
        "ssda_counts",
    }
    offsets = sorted(mod.PLASMA_BLOCKS.values())
    assert offsets == [3, 26, 49]
    # Blocks must not overlap: 20 participants plus two header rows apart.
    for lower, upper in zip(offsets, offsets[1:]):
        assert upper - lower >= mod.N_PARTICIPANTS


def test_ssda_is_declared_semi_quantitative(mod):
    """SSDA is in counts, the same weakness already logged for ARIANA."""
    assert "ssda_counts" in mod.SEMI_QUANTITATIVE_ANALYTES
    assert "solanidine_ng_ml" not in mod.SEMI_QUANTITATIVE_ANALYTES
    assert "oh4_solanidine_ng_ml" not in mod.SEMI_QUANTITATIVE_ANALYTES


def test_activity_score_convention_conflict_is_recorded(mod):
    """*41 is 0.25 here and the Rudesheim grid predates that CPIC change."""
    assert 0.25 not in mod.RUDESHEIM_AS_GRID
    assert mod.RUDESHEIM_AS_GRID == (0.0, 0.5, 1.25, 1.5, 2.0, 3.0)


def test_source_loq_rule_is_pinned(mod):
    """The source imputes below-LOQ values; the rule must stay quotable."""
    rule = mod.SOURCE_LOQ_RULE.lower()
    assert "limit of quantification" in rule
    assert "0.01 ng/ml" in rule and "0.02 ng/ml" in rule


class _FakeFont:
    def __init__(self, bold):
        self.bold = bold


class _FakeCell:
    def __init__(self, value, bold=False):
        self.value = value
        self.font = _FakeFont(bold)


class _FakeSheet:
    """Minimal worksheet stub: one analyte block of 20 x 9 cells."""

    def __init__(self, bold_map=None, missing=(), first_row=3):
        self.bold_map = bold_map or {}
        self.missing = set(missing)
        self.first_row = first_row

    def cell(self, row, column):
        if column == 1:
            return _FakeCell(row - self.first_row + 1)
        index = column - 2
        if (row, index) in self.missing:
            return _FakeCell(None)
        return _FakeCell(1.0, bold=self.bold_map.get((row, index), False))


def test_bold_cells_are_counted_as_loq_imputed(mod, monkeypatch):
    monkeypatch.setattr(mod, "PLASMA_BLOCKS", {"only": 3})
    bold = {(3, 0): True, (3, 1): True, (4, 0): True}
    result = mod.qualify_plasma_blocks(_FakeSheet(bold_map=bold))["only"]
    assert result["loq_imputed_cells"] == 3
    assert result["cells_present"] == 180
    assert result["cells_missing"] == 0
    assert result["complete"] is True
    assert result["participant_ids_are_1_to_n"] is True
    assert result["loq_imputed_percent"] == pytest.approx(100.0 * 3 / 180, abs=0.05)


def test_a_participant_entirely_below_loq_is_flagged(mod, monkeypatch):
    """The real cohort's only poor metaboliser is fully imputed; that must surface."""
    monkeypatch.setattr(mod, "PLASMA_BLOCKS", {"only": 3})
    bold = {(3, index): True for index in range(len(mod.TIMEPOINTS_H))}
    result = mod.qualify_plasma_blocks(_FakeSheet(bold_map=bold))["only"]
    assert result["participants_fully_below_loq"] == [1]
    assert result["most_imputed_participants"][0]["n_imputed"] == 9


def test_missing_cells_make_the_block_incomplete(mod, monkeypatch):
    monkeypatch.setattr(mod, "PLASMA_BLOCKS", {"only": 3})
    result = mod.qualify_plasma_blocks(_FakeSheet(missing={(5, 2)}))["only"]
    assert result["cells_missing"] == 1
    assert result["cells_present"] == 179
    assert result["complete"] is False


def test_plain_value_read_would_miss_the_imputation(mod, monkeypatch):
    """Imputed cells hold ordinary numbers; only the font distinguishes them.

    This is the reason the qualification uses openpyxl rather than pandas.
    """
    monkeypatch.setattr(mod, "PLASMA_BLOCKS", {"only": 3})
    sheet = _FakeSheet(bold_map={(3, 0): True})
    imputed_cell = sheet.cell(row=3, column=2)
    ordinary_cell = sheet.cell(row=4, column=2)
    assert imputed_cell.value == ordinary_cell.value
    assert imputed_cell.font.bold is not ordinary_cell.font.bold
    assert mod.qualify_plasma_blocks(sheet)["only"]["loq_imputed_cells"] == 1
