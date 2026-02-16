"""Tests for Faulon displacement (ported from displacement.test.ts)"""
import pytest
import json
from app.core.molecule import MoleculeGraph
from app.core.displacement import attempt_displacement
from app.core.random import SeededRandom


class TestBasicValidation:
    """Test basic displacement validation"""

    def test_null_for_less_than_4_atoms(self):
        """Molecules with < 4 atoms return None"""
        methane = MoleculeGraph.create_linear_alkane(1)  # CH4, 1 atom
        ethane = MoleculeGraph.create_linear_alkane(2)  # C2H6, 2 atoms
        propane = MoleculeGraph.create_linear_alkane(3)  # C3H8, 3 atoms
        rng = SeededRandom(42)

        assert attempt_displacement(methane, rng) is None
        assert attempt_displacement(ethane, rng) is None
        assert attempt_displacement(propane, rng) is None

    def test_preserves_atom_count(self):
        """100 displacements on hexane all have 6 atoms"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(12345)

        for _ in range(100):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                assert result.get_atom_count() == hexane.get_atom_count()

    def test_preserves_atom_types(self):
        """All atoms remain carbon after displacement"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(54321)

        for i in range(100):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                for j in range(result.get_atom_count()):
                    assert result.get_atom_element(j) == 'C'

    def test_does_not_mutate_original(self):
        """Original graph unchanged after displacement"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        # Capture original bond orders
        original_bonds = {}
        for i in range(6):
            for j in range(i + 1, 6):
                original_bonds[(i, j)] = hexane.get_bond_order(i, j)

        rng = SeededRandom(99999)
        attempt_displacement(hexane, rng)

        # Verify original unchanged
        for i in range(6):
            for j in range(i + 1, 6):
                assert hexane.get_bond_order(i, j) == original_bonds[(i, j)]


class TestBondOrderConstraints:
    """Test bond order constraints"""

    def test_bond_orders_in_range(self):
        """200 displacements, all bond orders in [0, 3]"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(777)

        for _ in range(200):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                n = result.get_atom_count()
                for i in range(n):
                    for j in range(n):
                        order = result.get_bond_order(i, j)
                        assert 0 <= order <= 3

    def test_handles_existing_double_bonds(self):
        """Molecule with C=C-C-C handles displacement"""
        # Create via SMILES: C=CCC
        from rdkit import Chem
        mol = Chem.MolFromSmiles('C=CCC')
        graph = MoleculeGraph(mol)

        rng = SeededRandom(333)

        for _ in range(100):
            result = attempt_displacement(graph, rng)
            if result is not None:
                # All bonds should still be in valid range
                n = result.get_atom_count()
                for i in range(n):
                    for j in range(n):
                        order = result.get_bond_order(i, j)
                        assert 0 <= order <= 3


class TestValencePreservation:
    """Test that Faulon displacement preserves valence (equations 7-9)"""

    def test_always_valid_valences(self):
        """200 displacements, result always has valid valences (Faulon preserves valence)"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(444)

        for _ in range(200):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                assert result.has_valid_valences()

    def test_produces_some_displacements(self):
        """At least some displacements succeed on hexane"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(555)
        success_count = 0

        for _ in range(100):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                success_count += 1

        assert success_count > 0

    def test_may_produce_disconnected(self):
        """Displacement may produce disconnected graphs (SA engine handles reconnection)"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(888)

        connected_count = 0
        disconnected_count = 0
        for _ in range(200):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                if result.is_connected():
                    connected_count += 1
                else:
                    disconnected_count += 1

        # Should have both connected and disconnected results
        assert connected_count > 0
        assert disconnected_count > 0


class TestReproducibility:
    """Test reproducibility with seeded RNG"""

    def test_reproducible_same_seed(self):
        """Seed 12345 produces identical results twice"""
        hexane = MoleculeGraph.create_linear_alkane(6)

        # Run 1
        rng1 = SeededRandom(12345)
        results1 = []
        for _ in range(10):
            result = attempt_displacement(hexane, rng1)
            if result is not None:
                results1.append(result.to_smiles())
            else:
                results1.append(None)

        # Run 2
        rng2 = SeededRandom(12345)
        results2 = []
        for _ in range(10):
            result = attempt_displacement(hexane, rng2)
            if result is not None:
                results2.append(result.to_smiles())
            else:
                results2.append(None)

        assert results1 == results2

    def test_different_seeds_differ(self):
        """Seeds 111 vs 222 produce different results"""
        hexane = MoleculeGraph.create_linear_alkane(6)

        rng1 = SeededRandom(111)
        results1 = []
        for _ in range(10):
            result = attempt_displacement(hexane, rng1)
            if result is not None:
                results1.append(result.to_smiles())
            else:
                results1.append(None)

        rng2 = SeededRandom(222)
        results2 = []
        for _ in range(10):
            result = attempt_displacement(hexane, rng2)
            if result is not None:
                results2.append(result.to_smiles())
            else:
                results2.append(None)

        # Should differ at some point
        assert results1 != results2


class TestStress500Displacements:
    """Stress test: 500 consecutive displacements"""

    def test_stress_500_displacements(self):
        """500 iterations on hexane, all preserve valence, some valid"""
        hexane = MoleculeGraph.create_linear_alkane(6)
        rng = SeededRandom(42424242)

        valid_count = 0
        for _ in range(500):
            result = attempt_displacement(hexane, rng)
            if result is not None:
                # Faulon preserves valence by construction
                assert result.has_valid_valences()
                valid_count += 1

        # Should have at least some successful displacements
        assert valid_count > 0


class TestCyclohexaneDisplacement:
    """Test displacement on cyclic structure"""

    def test_cyclohexane_displacement(self):
        """Cyclic structure displacement works and preserves valence"""
        cyclohexane = MoleculeGraph.create_cyclohexane()
        rng = SeededRandom(777)

        success_count = 0
        for _ in range(100):
            result = attempt_displacement(cyclohexane, rng)
            if result is not None:
                assert result.has_valid_valences()
                success_count += 1

        # Should have at least some successful displacements
        assert success_count > 0


class TestBondConservationEquations:
    """Test that displacement equations produce valid results"""

    def test_bond_conservation_equations(self):
        """4-atom molecule, equations preserve atom count and valence"""
        butane = MoleculeGraph.create_linear_alkane(4)
        rng = SeededRandom(999)

        for _ in range(50):
            result = attempt_displacement(butane, rng)
            if result is not None:
                # Faulon preserves atom count and valence
                assert result.get_atom_count() == 4
                assert result.has_valid_valences()


class TestUnsaturatedDisplacement:
    """Test displacement on unsaturated molecules (benzene, naphthalene)."""

    def test_benzene_displacement_no_crash(self):
        """Displacement on benzene (C6H6) doesn't crash."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles('c1ccccc1')
        graph = MoleculeGraph(mol)
        rng = SeededRandom(42)

        success_count = 0
        for _ in range(100):
            result = attempt_displacement(graph, rng)
            if result is not None:
                assert result.get_atom_count() == 6
                success_count += 1

        assert success_count > 0

    def test_naphthalene_displacement_no_crash(self):
        """Displacement on naphthalene (C10H8) doesn't crash."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles('c1ccc2ccccc2c1')
        graph = MoleculeGraph(mol)
        rng = SeededRandom(42)

        success_count = 0
        for _ in range(100):
            result = attempt_displacement(graph, rng)
            if result is not None:
                assert result.get_atom_count() == 10
                success_count += 1

        assert success_count > 0

    def test_unsaturated_bond_orders_valid(self):
        """Bond orders stay in [0, 3] for unsaturated molecules."""
        from rdkit import Chem
        mol = Chem.MolFromSmiles('c1ccccc1')
        graph = MoleculeGraph(mol)
        rng = SeededRandom(777)

        for _ in range(100):
            result = attempt_displacement(graph, rng)
            if result is not None:
                n = result.get_atom_count()
                for i in range(n):
                    for j in range(n):
                        order = result.get_bond_order(i, j)
                        assert 0 <= order <= 3
