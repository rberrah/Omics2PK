"""Contracts for the 3D structural annotation layer.

The guarantee that matters most is the **numbering guard**. A crystal
structure's residue numbers need not match the reference sequence, and an
offset would produce plausible wrong distances rather than an error — the worst
possible failure mode, because nothing would look broken. The guard asserts a
fact that is true only under correct numbering: UniProt calls Cys435 CYP2C9's
axial heme ligand, so its sulphur must lie within a coordinate bond of the
iron. `TestNumberingGuard` checks that the guard fires when numbering is wrong.

These tests run without the structures present: the ones that need a real PDB
file skip if it is absent, because the structures live outside the repository.
"""
from __future__ import annotations

import pathlib

import pytest

np = pytest.importorskip("numpy")

ROOT = pathlib.Path(__file__).resolve().parents[1]
PDB_DIR = pathlib.Path(
    "C:/Users/abdel/Documents/Omics/PDB_public_structures")


@pytest.fixture(scope="module")
def sa():
    import sys
    sys.path.insert(0, str(ROOT / "src"))
    from omics2pk import structural_annotation
    return structural_annotation


def _atom(record, serial, atom, resname, chain, resseq, x, y, z, element):
    """Build a fixed-column PDB line, so the parser is tested on real format."""
    return (f"{record:<6}{serial:>5} {atom:^4}{'':1}{resname:>3} {chain}"
            f"{resseq:>4}{'':1}   {x:>8.3f}{y:>8.3f}{z:>8.3f}"
            f"{1.0:>6.2f}{0.0:>6.2f}          {element:>2}  ")


@pytest.fixture
def toy():
    """Two residues and a metal, at distances chosen so the geometry is known."""
    lines = [
        _atom("ATOM", 1, "CA", "CYS", "A", 435, 0.0, 0.0, 0.0, "C"),
        _atom("ATOM", 2, "SG", "CYS", "A", 435, 2.0, 0.0, 0.0, "S"),
        _atom("ATOM", 3, "CA", "THR", "A", 301, 10.0, 0.0, 0.0, "C"),
        _atom("HETATM", 4, "FE", "HEM", "A", 500, 4.0, 0.0, 0.0, "FE"),
        _atom("HETATM", 5, "C1", "HEM", "A", 500, 5.0, 0.0, 0.0, "C"),
        _atom("HETATM", 6, "C1", "FLP", "A", 501, 12.0, 0.0, 0.0, "C"),
        _atom("HETATM", 7, "O", "HOH", "A", 600, 1.0, 1.0, 1.0, "O"),
        _atom("ATOM", 8, "H", "CYS", "A", 435, 2.5, 0.0, 0.0, "H"),
    ]
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# 1. the parser
# --------------------------------------------------------------------------
class TestParser:
    def test_it_reads_residues_and_ligands_separately(self, sa, toy):
        residues, ligands = sa.parse_pdb(toy)
        assert set(residues) == {435, 301}
        assert set(ligands) == {"HEM", "FLP"}

    def test_water_is_excluded(self, sa, toy):
        _, ligands = sa.parse_pdb(toy)
        assert "HOH" not in ligands

    def test_hydrogen_is_excluded(self, sa, toy):
        residues, _ = sa.parse_pdb(toy)
        assert "H" not in residues[435]["atoms"], (
            "hydrogens are absent from most structures, so including them "
            "would make a distance depend on the depositor")

    def test_residue_names_are_kept(self, sa, toy):
        residues, _ = sa.parse_pdb(toy)
        assert residues[435]["name"] == "CYS"
        assert residues[301]["name"] == "THR"

    def test_an_empty_file_is_refused_not_returned_empty(self, sa):
        with pytest.raises(sa.StructuralAnnotationError,
                           match="no polymer residue"):
            sa.parse_pdb("HEADER    nothing here\n")

    def test_insertion_codes_are_skipped(self, sa):
        line = _atom("ATOM", 1, "CA", "ALA", "A", 50, 0.0, 0.0, 0.0, "C")
        inserted = line[:26] + "B" + line[27:]
        good = _atom("ATOM", 2, "CA", "GLY", "A", 51, 1.0, 0.0, 0.0, "C")
        residues, _ = sa.parse_pdb(inserted + "\n" + good + "\n")
        assert set(residues) == {51}, (
            "an insertion code is not a reference-sequence position")

    def test_only_the_first_chain_is_used(self, sa):
        a = _atom("ATOM", 1, "CA", "ALA", "A", 10, 0.0, 0.0, 0.0, "C")
        b = _atom("ATOM", 2, "CA", "GLY", "B", 11, 1.0, 0.0, 0.0, "C")
        residues, _ = sa.parse_pdb(a + "\n" + b + "\n")
        assert set(residues) == {10}


# --------------------------------------------------------------------------
# 2. distances
# --------------------------------------------------------------------------
class TestDistances:
    def test_distances_are_minimum_over_heavy_atoms(self, sa, toy):
        residues, ligands = sa.parse_pdb(toy)
        centre = {"cofactor_ligands": ("HEM",), "substrate_ligands": ("FLP",),
                  "metal_atom": "FE", "pdb": "TOY"}
        out = sa.residue_distances(residues, ligands, centre)
        # SG at x=2 is 2.0 from FE at x=4; CA at x=0 is 4.0. Minimum is 2.0.
        assert out[435]["to_cofactor"] == pytest.approx(2.0)
        assert out[435]["to_metal"] == pytest.approx(2.0)
        # THR CA at x=10 is 2.0 from FLP at x=12.
        assert out[301]["to_substrate"] == pytest.approx(2.0)

    def test_a_missing_cofactor_is_refused(self, sa, toy):
        residues, ligands = sa.parse_pdb(toy)
        centre = {"cofactor_ligands": ("ZZZ",), "substrate_ligands": ("FLP",),
                  "metal_atom": "FE", "pdb": "TOY"}
        with pytest.raises(sa.StructuralAnnotationError,
                           match="no cofactor atoms"):
            sa.residue_distances(residues, ligands, centre)

    def test_a_missing_substrate_is_refused(self, sa, toy):
        residues, ligands = sa.parse_pdb(toy)
        centre = {"cofactor_ligands": ("HEM",), "substrate_ligands": ("ZZZ",),
                  "metal_atom": "FE", "pdb": "TOY"}
        with pytest.raises(sa.StructuralAnnotationError,
                           match="no substrate atoms"):
            sa.residue_distances(residues, ligands, centre)


# --------------------------------------------------------------------------
# 3. the numbering guard — the one that prevents a silent wrong answer
# --------------------------------------------------------------------------
class TestNumberingGuard:
    @staticmethod
    def _centre(**overrides):
        base = {"cofactor_ligands": ("HEM",), "substrate_ligands": ("FLP",),
                "metal_atom": "FE", "axial_residue": 435, "axial_atom": "SG",
                "pdb": "TOY", "uniprot": "P11712"}
        base.update(overrides)
        return base

    def test_correct_geometry_passes_and_reports_the_bond(self, sa, toy):
        residues, ligands = sa.parse_pdb(toy)
        report = sa.verify_catalytic_geometry(residues, ligands,
                                              self._centre())
        assert report["verified"] is True
        assert report["axial_bond_angstrom"] == pytest.approx(2.0)

    def test_a_renumbered_structure_is_caught(self, sa, toy):
        """Shift the axial residue's number: the guard must fire."""
        shifted = toy.replace(" 435 ", " 436 ")
        residues, ligands = sa.parse_pdb(shifted)
        with pytest.raises(sa.StructuralAnnotationError,
                           match="no residue 435"):
            sa.verify_catalytic_geometry(residues, ligands, self._centre())

    def test_a_wrong_residue_identity_is_caught(self, sa, toy):
        wrong = toy.replace("CYS", "SER")
        residues, ligands = sa.parse_pdb(wrong)
        with pytest.raises(sa.StructuralAnnotationError, match="not CYS"):
            sa.verify_catalytic_geometry(residues, ligands, self._centre())

    def test_a_distance_too_long_for_a_bond_is_caught(self, sa, toy):
        """Move the iron far away: the axial claim becomes false."""
        far = toy.replace(f"{4.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
                          f"{1.0:>6.2f}{0.0:>6.2f}          FE",
                          f"{40.0:>8.3f}{0.0:>8.3f}{0.0:>8.3f}"
                          f"{1.0:>6.2f}{0.0:>6.2f}          FE")
        residues, ligands = sa.parse_pdb(far)
        with pytest.raises(sa.StructuralAnnotationError,
                           match="coordinate bond allows"):
            sa.verify_catalytic_geometry(residues, ligands, self._centre())

    def test_a_missing_metal_is_caught(self, sa, toy):
        residues, ligands = sa.parse_pdb(toy)
        with pytest.raises(sa.StructuralAnnotationError, match="no FE atom"):
            sa.verify_catalytic_geometry(residues, ligands,
                                         self._centre(metal_atom="FE",
                                                      cofactor_ligands=("FLP",)))

    def test_the_motif_variant_of_the_guard_fires_too(self, sa, toy):
        centre = self._centre(axial_residue=None, axial_atom=None,
                              motif_span=(300, 302))
        residues, ligands = sa.parse_pdb(toy)
        # THR301 is 6 A from FE, beyond the 5 A contact shell.
        with pytest.raises(sa.StructuralAnnotationError,
                           match="beyond contact range"):
            sa.verify_catalytic_geometry(residues, ligands, centre)


# --------------------------------------------------------------------------
# 4. the feature block
# --------------------------------------------------------------------------
class TestFeatureBlock:
    def test_the_block_has_the_declared_width(self, sa):
        row = sa.structural_features(10, {10: {"to_cofactor": 3.0,
                                               "to_substrate": 9.0,
                                               "to_metal": 3.0,
                                               "residue": "A"}})
        assert len(row) == len(sa.STRUCTURAL_NAMES)

    def test_an_unresolved_position_is_flagged_not_faked(self, sa):
        row = sa.structural_features(999, {})
        assert row[sa.STRUCTURAL_NAMES.index("unresolved_in_structure")] == 1.0

    def test_contact_and_second_shell_are_distinguished(self, sa):
        close = sa.structural_features(1, {1: {"to_cofactor": 3.0,
                                               "to_substrate": 20.0,
                                               "to_metal": 3.0,
                                               "residue": "A"}})
        middle = sa.structural_features(1, {1: {"to_cofactor": 7.0,
                                                "to_substrate": 20.0,
                                                "to_metal": 7.0,
                                                "residue": "A"}})
        assert close[sa.STRUCTURAL_NAMES.index("in_cofactor_contact")] == 1.0
        assert middle[sa.STRUCTURAL_NAMES.index("in_cofactor_contact")] == 0.0
        assert middle[sa.STRUCTURAL_NAMES.index("in_second_shell")] == 1.0

    def test_build_structural_block_keeps_the_unresolved_flag(self, sa):
        distances = {1: {"to_cofactor": 3.0, "to_substrate": 4.0,
                         "to_metal": 3.0, "residue": "A"},
                     2: {"to_cofactor": 9.0, "to_substrate": 11.0,
                         "to_metal": 9.0, "residue": "G"}}
        block = sa.build_structural_block([1, 2, 999], distances)
        flag = sa.STRUCTURAL_NAMES.index("unresolved_in_structure")
        assert block[2, flag] == 1.0
        assert np.isfinite(block).all(), (
            "an unresolved position must get a finite filled value, so a "
            "model can use the flag rather than crash on a NaN")

    def test_the_block_is_finite_for_every_resolved_position(self, sa):
        distances = {p: {"to_cofactor": float(p), "to_substrate": float(p) + 1,
                         "to_metal": float(p), "residue": "A"}
                     for p in range(1, 40)}
        block = sa.build_structural_block(list(range(1, 40)), distances)
        assert np.isfinite(block).all()


# --------------------------------------------------------------------------
# 5. coverage — the number that decided the earlier diagnosis
# --------------------------------------------------------------------------
class TestCoverage:
    @staticmethod
    def _distances(n=100):
        return {p: {"to_cofactor": float(abs(p - 50)),
                    "to_substrate": float(abs(p - 50)) + 1.0,
                    "to_metal": float(abs(p - 50)), "residue": "A"}
                for p in range(1, n + 1)}

    def test_coverage_counts_the_discordant_variants_reached(self, sa):
        distances = self._distances()
        positions = list(range(1, 101))
        discordant = [abs(p - 50) <= 4 for p in positions]
        out = sa.coverage(positions, discordant, distances, cutoff=8.0)
        assert out["share_of_discordant_variants_reached"] == 1.0

    def test_a_ruler_that_misses_everything_reports_zero(self, sa):
        distances = self._distances()
        positions = list(range(1, 101))
        discordant = [p <= 3 for p in positions]
        out = sa.coverage(positions, discordant, distances, cutoff=8.0)
        assert out["share_of_discordant_variants_reached"] == 0.0

    def test_enrichment_is_a_ratio_of_rates(self, sa):
        distances = self._distances()
        positions = list(range(1, 101))
        # Discordance must occur on BOTH sides of the cutoff, or the ratio has
        # a zero denominator and is unbounded rather than large.
        discordant = [abs(p - 50) <= 4 or p == 1 for p in positions]
        out = sa.coverage(positions, discordant, distances, cutoff=8.0)
        assert out["enrichment_ratio"] > 1.0
        assert out["discordance_rate_within"] > out["discordance_rate_beyond"]

    def test_an_unbounded_ratio_says_why_instead_of_returning_a_bare_nan(
            self, sa):
        distances = self._distances()
        positions = list(range(1, 101))
        discordant = [abs(p - 50) <= 4 for p in positions]
        out = sa.coverage(positions, discordant, distances, cutoff=8.0)
        assert np.isnan(out["enrichment_ratio"])
        assert "unbounded rather than large" in             out["enrichment_unavailable_because"]

    def test_unresolved_positions_are_counted_and_not_called_near(self, sa):
        distances = self._distances(50)
        positions = list(range(1, 101))
        discordant = [False] * 100
        out = sa.coverage(positions, discordant, distances, cutoff=8.0)
        assert out["n_unresolved"] == 50

    def test_the_interval_is_returned_or_declared_unavailable(self, sa):
        distances = self._distances()
        positions = list(range(1, 101))
        discordant = [abs(p - 50) <= 4 for p in positions]
        low, high = sa.enrichment_interval(positions, discordant, distances,
                                           n_boot=400)
        assert np.isnan(low) or low > 0.0


# --------------------------------------------------------------------------
# 6. the declared centres, and what the module refuses to claim
# --------------------------------------------------------------------------
class TestDeclarations:
    def test_every_gene_declares_a_centre(self, sa):
        assert set(sa.CATALYTIC_CENTRES) == {"CYP2C9", "NUDT15", "BRAF", "ARAF"}

    def test_each_centre_names_its_structure_and_a_substrate(self, sa):
        for gene, centre in sa.CATALYTIC_CENTRES.items():
            assert centre["pdb"], gene
            assert centre["cofactor_ligands"], gene
            assert centre["uniprot"], gene
            # A substrate is defined by a ligand OR by a protein chain, and
            # exactly one of the two. Neither would make `to_substrate` a
            # silent infinity and turn the ruler cofactor-only without saying
            # so; both would make it ambiguous which one was measured.
            ligand = bool(centre["substrate_ligands"])
            chain = bool(centre.get("substrate_chain"))
            assert ligand != chain, gene

    def test_braf_is_the_pre_registered_test_and_names_its_chain(self, sa):
        centre = sa.CATALYTIC_CENTRES["BRAF"]
        assert centre["pdb"] == "6PP9"
        assert centre["uniprot"] == "P15056"
        # MEK1 carries its own AMP-PNP and Mg, so the chain is load-bearing.
        assert centre["chain"] == "A"
        assert centre["substrate_chain"] == "B"
        assert centre["substrate_is_a_protein"] == "MEK1"

    def test_brafs_guard_is_identities_plus_a_nucleotide_contact(self, sa):
        centre = sa.CATALYTIC_CENTRES["BRAF"]
        assert centre["axial_residue"] is None
        identities = centre["residue_identities"]
        # The DFG motif and the catalytic lysine, at the numbers UniProt uses.
        assert identities[483] == "LYS"
        assert (identities[594], identities[595], identities[596]) == \
            ("ASP", "PHE", "GLY")
        assert centre["nucleotide_contact"]["residue"] == 483
        assert centre["nucleotide_contact"]["max_angstrom"] <= 4.0

    def test_the_salt_bridge_was_deliberately_not_the_guard(self, sa):
        # The obvious kinase guard, Lys483 to Glu501, is BROKEN at 12.69 A in
        # 6PP9's inactive aC-out conformation. Declaring it would have failed
        # on correctly numbered coordinates.
        source = (ROOT / "src" / "omics2pk"
                  / "structural_annotation.py").read_text(encoding="utf-8")
        assert "why_not_the_salt_bridge" in source

    def test_cyp2c9_uses_the_axial_cysteine_as_its_guard(self, sa):
        assert sa.CATALYTIC_CENTRES["CYP2C9"]["axial_residue"] == 435

    def test_nudt15_uses_the_nudix_box_as_its_guard(self, sa):
        centre = sa.CATALYTIC_CENTRES["NUDT15"]
        assert centre["axial_residue"] is None
        assert centre["motif_span"] == (48, 69)

    def test_it_refuses_to_read_a_distance_as_causal(self, sa):
        out = sa.what_this_module_does_not_do()
        assert "is a hypothesis" in out["a_distance_is_not_a_cause"]
        assert "never" in out["a_distance_is_not_a_cause"]

    def test_it_declares_the_single_conformation_limit(self, sa):
        out = sa.what_this_module_does_not_do()
        assert "one conformation" in out["one_structure_is_one_conformation"]
        assert "membrane anchor truncated" in \
            out["one_structure_is_one_conformation"]

    def test_it_declares_the_bound_ligand_limit(self, sa):
        out = sa.what_this_module_does_not_do()
        assert "one pocket" in out["a_bound_ligand_is_not_every_substrate"]

    def test_it_says_numbering_is_verified_not_trusted(self, sa):
        out = sa.what_this_module_does_not_do()
        assert "plausible wrong distances" in \
            out["numbering_is_verified_not_trusted"]


# --------------------------------------------------------------------------
# 7. against the real structures, when they are present
# --------------------------------------------------------------------------
@pytest.mark.skipif(not PDB_DIR.is_dir(),
                    reason="the structures live outside the repository")
class TestAgainstTheRealStructures:
    @pytest.mark.parametrize("gene", ["CYP2C9", "NUDT15"])
    def test_the_numbering_guard_passes_on_the_declared_structure(self, sa,
                                                                   gene):
        centre = sa.CATALYTIC_CENTRES[gene]
        path = PDB_DIR / f"{centre['pdb']}.pdb"
        if not path.is_file():
            pytest.skip(f"{path.name} not downloaded")
        residues, ligands = sa.parse_pdb(
            path.read_text(encoding="utf-8", errors="replace"))
        report = sa.verify_catalytic_geometry(residues, ligands, centre)
        assert report["verified"] is True

    def test_the_cyp2c9_axial_bond_is_a_coordinate_bond(self, sa):
        centre = sa.CATALYTIC_CENTRES["CYP2C9"]
        path = PDB_DIR / f"{centre['pdb']}.pdb"
        if not path.is_file():
            pytest.skip("1R9O not downloaded")
        residues, ligands = sa.parse_pdb(
            path.read_text(encoding="utf-8", errors="replace"))
        report = sa.verify_catalytic_geometry(residues, ligands, centre)
        assert 2.0 <= report["axial_bond_angstrom"] <= 2.8, (
            "Cys435 SG to heme Fe should be a coordinate bond near 2.3 A; "
            f"got {report['axial_bond_angstrom']}")

    def test_the_residues_the_sequence_ruler_missed_are_close_in_3d(self, sa):
        """CYP2C9 Thr301 is 134 residues from the only annotation and in
        contact with the heme. That gap is the whole reason this layer exists."""
        centre = sa.CATALYTIC_CENTRES["CYP2C9"]
        path = PDB_DIR / f"{centre['pdb']}.pdb"
        if not path.is_file():
            pytest.skip("1R9O not downloaded")
        residues, ligands = sa.parse_pdb(
            path.read_text(encoding="utf-8", errors="replace"))
        distances = sa.residue_distances(residues, ligands, centre)
        assert distances[301]["to_cofactor"] < 5.0
        assert abs(301 - centre["axial_residue"]) > 100


class TestArafCentreAndMmcifReader:
    """ARAF, and the mmCIF reader that had to exist for it.

    ARAF has seven RCSB entries and exactly one carries an ATP analogue --
    9AXM -- which RCSB serves only as mmCIF. The declared selection rule said
    "legacy PDB format", which was a limitation of this module's reader rather
    than a scientific criterion, so the clause was lifted and recorded as a
    deviation. The rule picks 9AXM either way, so the lift could not steer the
    choice.
    """

    def test_araf_is_declared_with_its_chain_and_protein_substrate(self, sa):
        centre = sa.CATALYTIC_CENTRES["ARAF"]
        assert centre["pdb"] == "9AXM"
        assert centre["uniprot"] == "P10398"
        # MEK1 carries its own AMP-PNP and Mg, as in 6PP9.
        assert centre["chain"] == "B"
        assert centre["substrate_chain"] == "A"
        assert centre["substrate_is_a_protein"] == "MEK1"

    def test_the_araf_guard_is_identities_plus_a_nucleotide_contact(self, sa):
        centre = sa.CATALYTIC_CENTRES["ARAF"]
        assert centre["axial_residue"] is None
        identities = centre["residue_identities"]
        assert identities[336] == "LYS"      # UniProt catalytic lysine
        assert identities[429] == "ASP"      # UniProt proton acceptor
        assert (identities[447], identities[448], identities[449]) == \
            ("ASP", "PHE", "GLY")            # DFG, found by scanning
        assert centre["nucleotide_contact"]["residue"] == 336
        assert centre["nucleotide_contact"]["max_angstrom"] <= 4.0

    def test_the_mmcif_reader_uses_author_numbering(self, sa):
        source = (ROOT / "src" / "omics2pk"
                  / "structural_annotation.py").read_text(encoding="utf-8")
        # label_seq_id restarts at 1 per entity and would shift every residue.
        assert "auth_seq_id" in source
        assert "label_seq_id` restarts at 1" in source

    def test_the_dispatcher_detects_the_format(self, sa):
        assert hasattr(sa, "parse_structure")
        assert hasattr(sa, "chain_points")
        source = (ROOT / "src" / "omics2pk"
                  / "structural_annotation.py").read_text(encoding="utf-8")
        assert '"_atom_site." in text' in source

    def test_an_mmcif_without_an_atom_site_loop_is_refused(self, sa):
        with pytest.raises(Exception):
            sa.parse_mmcif("data_XXXX\n_entry.id XXXX\n")

    def test_both_readers_agree_where_both_formats_exist(self, sa):
        """The validation that keeps the reader from rewriting BRAF silently."""
        pdb = pathlib.Path(
            r"C:\Users\abdel\Documents\Omics\PDB_public_structures\6PP9.pdb")
        cif = pathlib.Path(r"C:\Users\abdel\Documents\Omics"
                           r"\PDB_public_structures_araf\6PP9_validation.cif")
        if not (pdb.is_file() and cif.is_file()):
            pytest.skip("both formats of 6PP9 are not present locally")
        centre = sa.CATALYTIC_CENTRES["BRAF"]
        from_pdb = sa.parse_pdb(pdb.read_text(encoding="utf-8",
                                              errors="replace"), chain="A")
        from_cif = sa.parse_mmcif(cif.read_text(encoding="utf-8",
                                                errors="replace"), chain="A")
        assert set(from_pdb[0]) == set(from_cif[0])
        assert {k: len(v) for k, v in from_pdb[1].items()} == \
            {k: len(v) for k, v in from_cif[1].items()}
        guard_pdb = sa.verify_catalytic_geometry(*from_pdb, centre)
        guard_cif = sa.verify_catalytic_geometry(*from_cif, centre)
        assert guard_pdb["nucleotide_contact_angstrom"] == \
            guard_cif["nucleotide_contact_angstrom"]

    def test_the_deviations_are_declared_in_the_protocol(self):
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        # No .lower(): the search string carries "PDB" in capitals, so
        # folding the document would guarantee a miss.
        assert "clause \u00ab format PDB h\u00e9ritable \u00bb est lev\u00e9e" in text
        assert "ne peut donc pas avoir \u00e9t\u00e9 dirig\u00e9e par le r\u00e9sultat" in text
        assert "n'a pas d'abondance au d\u00e9part" in text
        assert "2,979" in text
        assert "rompu \u00e0 10,49" in text


class TestArafResults:
    """ARAF's outcome, including the prediction that failed and cost me a theory.

    The primary contrast held -- the structural ablation is refuted under the
    region split, as predicted -- but three of five predictions failed, and
    one of them was the only prediction that could fail alone. So the reason
    for predicting the primary was wrong, and the record has to say so.
    """

    @pytest.fixture(scope="class")
    def araf(self):
        import json
        path = ROOT / "reports" / "data" / "variant_mechanism_benchmark.json"
        if not path.is_file():
            pytest.skip("the benchmark artefact has not been generated")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if "ARAF" not in payload.get("genes", {}):
            pytest.skip("ARAF has not been run")
        return payload["genes"]["ARAF"]

    @pytest.fixture(scope="class")
    def hsp90(self):
        import json
        path = ROOT / "reports" / "data" / "araf_hsp90_proxy.json"
        if not path.is_file():
            pytest.skip("the HSP90 artefact has not been generated")
        return json.loads(path.read_text(encoding="utf-8"))

    def test_the_structure_actually_covers_the_assayed_positions(self, araf):
        """The N-terminal construct had zero overlap; this must not recur."""
        coverage = araf["structural_layer"]["coverage"]
        assert coverage["n_within"] > 0
        assert coverage["n_unresolved"] < araf["n_shared_parsed"] / 2

    def test_the_primary_contrast_is_refuted(self, araf):
        out = araf["runs"]["grouped_by_contiguous_region"][
            "with_structural_layer"]["ablation_of_the_structural_block"][
                "M4_fixed_union"]
        assert out["verdict"] == "REFUTED_GAIN_BELOW_MARGIN"

    def test_the_structural_block_helps_under_the_easier_split(self, araf):
        """Fourth protein with the same split-dependence."""
        out = araf["runs"]["grouped_by_position"]["with_structural_layer"][
            "ablation_of_the_structural_block"]["M4_fixed_union"]
        assert out["interval"][0] > 0.0

    def test_araf_is_loss_of_function_unlike_braf(self, hsp90):
        polarity = hsp90["assay_polarity"]
        assert polarity["missense_above_synonymous_in_activity"] is False
        assert polarity["activity"]["missense"] < \
            polarity["activity"]["synonymous"]
        assert polarity["activity"]["nonsense"] < \
            polarity["activity"]["missense"]

    def test_the_abundance_correlation_is_positive_not_negative(self, araf):
        # BRAF was -0.290. The prediction that ARAF would match it failed.
        assert araf["spearman_abundance_activity"] > 0.0

    def test_the_hsp90_prediction_failed_in_the_opposite_direction(self,
                                                                   hsp90):
        out = hsp90["contrasts"]["hsp90i_better_than_dmso"]
        assert out["verdict"] == "DETECTABLY_WORSE"
        assert out["relative_gain"] < 0.0
        # And the inner folds rejected the predictor outright.
        assert hsp90["k_chosen_per_fold"]["abundance_hsp90i"] == [0] * 5

    def test_the_two_abundances_are_not_measuring_different_things(self,
                                                                   hsp90):
        # rho 0.714 between them, so the failure is not a measurement problem.
        assert hsp90["spearman"]["dmso_vs_hsp90i"] > 0.5

    def test_the_document_records_that_the_reason_was_wrong(self):
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "ma raison de l'avoir prédit est fausse" in text
        assert "n'est pas une confirmation" in text
        assert "L'explication HSP90 de l'échec BRAF est réfutée" in text

    def test_the_document_revises_the_scope_claim(self):
        text = (ROOT / "docs"
                / "OMICS2PK_HIGH_DIMENSION_PREREGISTRATION_20261111.md"
                ).read_text(encoding="utf-8")
        assert "enzyme métabolique contre kinase de signalisation" in text
        # And states why ARAF is the evidence for it.
        assert "ARAF le montre : perte de" in text
