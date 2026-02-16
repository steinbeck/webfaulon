"""Tests for initial structure generation (ported from initialStructure.test.ts)"""
import pytest
from app.core.initial_structure import generate_initial_structure


class TestSaturatedAlkanes:
    """Test saturated alkane generation"""

    def test_hexane_structure(self):
        """C6H14 -> 6 atoms, all C, connected, valid valences, 14 implicit H"""
        graph = generate_initial_structure('C6H14')

        # Check atom count
        assert graph.get_atom_count() == 6

        # Check all atoms are carbon
        for i in range(6):
            assert graph.get_atom_element(i) == 'C'

        # Check connectivity
        assert graph.is_connected()

        # Check valid valences
        assert graph.has_valid_valences()

        # Check total implicit H matches formula
        total_h = sum(graph.get_implicit_h(i) for i in range(6))
        assert total_h == 14

    def test_butane_structure(self):
        """C4H10 -> 4 atoms, 10 implicit H"""
        graph = generate_initial_structure('C4H10')

        assert graph.get_atom_count() == 4
        assert graph.is_connected()
        assert graph.has_valid_valences()

        total_h = sum(graph.get_implicit_h(i) for i in range(4))
        assert total_h == 10

    def test_methane_structure(self):
        """CH4 -> 1 atom, 4 implicit H"""
        graph = generate_initial_structure('CH4')

        assert graph.get_atom_count() == 1
        assert graph.is_connected()
        assert graph.has_valid_valences()
        assert graph.get_implicit_h(0) == 4


class TestUnsaturation:
    """Test unsaturated molecule generation"""

    def test_unsaturation_hdi_1(self):
        """C6H12 -> 6 atoms, 12 implicit H, at least one double bond"""
        graph = generate_initial_structure('C6H12')

        assert graph.get_atom_count() == 6
        assert graph.is_connected()
        assert graph.has_valid_valences()

        total_h = sum(graph.get_implicit_h(i) for i in range(6))
        assert total_h == 12

        # Should have at least one double bond (bond order = 2)
        has_double_bond = False
        for i in range(6):
            for j in range(i + 1, 6):
                if graph.get_bond_order(i, j) == 2:
                    has_double_bond = True
                    break
            if has_double_bond:
                break
        assert has_double_bond

    def test_unsaturation_hdi_4(self):
        """C6H6 -> 6 atoms, 6 implicit H"""
        graph = generate_initial_structure('C6H6')

        assert graph.get_atom_count() == 6
        assert graph.is_connected()
        assert graph.has_valid_valences()

        total_h = sum(graph.get_implicit_h(i) for i in range(6))
        assert total_h == 6


class TestHeteroatoms:
    """Test heteroatom handling"""

    def test_heteroatom_oxygen(self):
        """C2H6O -> 3 atoms (2C + 1O), 6 implicit H"""
        graph = generate_initial_structure('C2H6O')

        assert graph.get_atom_count() == 3  # 2 carbons + 1 oxygen
        assert graph.is_connected()
        assert graph.has_valid_valences()

        # Count element types
        c_count = sum(1 for i in range(3) if graph.get_atom_element(i) == 'C')
        o_count = sum(1 for i in range(3) if graph.get_atom_element(i) == 'O')
        assert c_count == 2
        assert o_count == 1

        total_h = sum(graph.get_implicit_h(i) for i in range(3))
        assert total_h == 6

    def test_heteroatom_nitrogen(self):
        """C2H7N -> 3 atoms, 7 implicit H"""
        graph = generate_initial_structure('C2H7N')

        assert graph.get_atom_count() == 3  # 2 C + 1 N
        assert graph.is_connected()
        assert graph.has_valid_valences()

        total_h = sum(graph.get_implicit_h(i) for i in range(3))
        assert total_h == 7


class TestSimpleUnsaturatedMolecules:
    """Test simple unsaturated molecules"""

    def test_ethene(self):
        """C2H4 -> 2 atoms, double bond, 4 implicit H"""
        graph = generate_initial_structure('C2H4')

        assert graph.get_atom_count() == 2
        assert graph.is_connected()
        assert graph.has_valid_valences()

        total_h = sum(graph.get_implicit_h(i) for i in range(2))
        assert total_h == 4

        # Should have a double bond
        assert graph.get_bond_order(0, 1) == 2

    def test_ethyne(self):
        """C2H2 -> 2 atoms, triple bond, 2 implicit H"""
        graph = generate_initial_structure('C2H2')

        assert graph.get_atom_count() == 2
        assert graph.is_connected()
        assert graph.has_valid_valences()

        total_h = sum(graph.get_implicit_h(i) for i in range(2))
        assert total_h == 2

        # Should have a triple bond
        assert graph.get_bond_order(0, 1) == 3


class TestErrorHandling:
    """Test error handling"""

    def test_empty_formula_raises(self):
        """Empty string raises error"""
        with pytest.raises(Exception):
            generate_initial_structure('')

    def test_impossible_formula_raises(self):
        """C6H99 raises error (impossible hydrogen count)"""
        with pytest.raises(Exception):
            generate_initial_structure('C6H99')


class TestDeterminism:
    """Test deterministic output"""

    def test_deterministic(self):
        """Same formula produces same structure"""
        graph1 = generate_initial_structure('C6H14')
        graph2 = generate_initial_structure('C6H14')

        # Should produce identical SMILES
        assert graph1.to_smiles() == graph2.to_smiles()


class TestEdgeCases:
    """Test edge cases"""

    def test_h2_no_heavy_atoms(self):
        """H2 -> 0 heavy atoms"""
        graph = generate_initial_structure('H2')
        assert graph.get_atom_count() == 0  # No heavy atoms, H2 is just implicit
