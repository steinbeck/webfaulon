"""Tests for Wiener Index computation (ported from wiener.test.ts)"""
import pytest
from app.core.molecule import MoleculeGraph
from app.core.wiener import compute_wiener_index


class TestKnownValuesLinearAlkanes:
    """Test known Wiener Index values for linear alkanes"""

    def test_methane_wiener_0(self):
        """Single atom has Wiener Index 0"""
        graph = MoleculeGraph.create_linear_alkane(1)
        assert compute_wiener_index(graph) == 0

    def test_ethane_wiener_1(self):
        """Ethane (C2) has Wiener Index 1"""
        graph = MoleculeGraph.create_linear_alkane(2)
        assert compute_wiener_index(graph) == 1

    def test_propane_wiener_4(self):
        """Propane (C3) has Wiener Index 4"""
        graph = MoleculeGraph.create_linear_alkane(3)
        # Formula: n(n^2-1)/6 = 3*8/6 = 4
        assert compute_wiener_index(graph) == 4

    def test_butane_wiener_10(self):
        """n-Butane (C4) has Wiener Index 10"""
        graph = MoleculeGraph.create_linear_alkane(4)
        # Formula: 4*15/6 = 10
        assert compute_wiener_index(graph) == 10

    def test_pentane_wiener_20(self):
        """n-Pentane (C5) has Wiener Index 20"""
        graph = MoleculeGraph.create_linear_alkane(5)
        # Formula: 5*24/6 = 20
        assert compute_wiener_index(graph) == 20

    def test_hexane_wiener_35(self):
        """n-Hexane (C6) has Wiener Index 35"""
        graph = MoleculeGraph.create_linear_alkane(6)
        # Formula: 6*35/6 = 35
        assert compute_wiener_index(graph) == 35


class TestKnownValuesCyclicStructures:
    """Test known Wiener Index values for cyclic structures"""

    def test_cyclohexane_wiener_27(self):
        """Cyclohexane has Wiener Index 27"""
        graph = MoleculeGraph.create_cyclohexane()
        # Pairs at distance 1: 6 (adjacent)
        # Pairs at distance 2: 6 (separated by 1)
        # Pairs at distance 3: 3 (opposite vertices)
        # Total: 6*1 + 6*2 + 3*3 = 6 + 12 + 9 = 27
        assert compute_wiener_index(graph) == 27


class TestKnownValuesBranchedStructures:
    """Test known Wiener Index values for branched structures"""

    def test_isobutane_wiener_9(self):
        """Isobutane has Wiener Index 9"""
        graph = MoleculeGraph.create_branched('isobutane')
        # Central C bonded to 3 others
        # 3 pairs at distance 1 (center to each terminal)
        # 3 pairs at distance 2 (terminal to terminal through center)
        # Total: 3*1 + 3*2 = 3 + 6 = 9
        assert compute_wiener_index(graph) == 9

    def test_neopentane_wiener_16(self):
        """Neopentane has Wiener Index 16"""
        graph = MoleculeGraph.create_branched('neopentane')
        # Central C bonded to 4 others
        # 4 pairs at distance 1 (center to each terminal)
        # 6 pairs at distance 2 (4 choose 2 = 6 terminal pairs through center)
        # Total: 4*1 + 6*2 = 4 + 12 = 16
        assert compute_wiener_index(graph) == 16


class TestLinearAlkaneFormula:
    """Test linear alkane formula: Wiener Index = n(n^2-1)/6"""

    def test_linear_formula_n7(self):
        """n=7 should give (7*48)/6 = 56"""
        graph = MoleculeGraph.create_linear_alkane(7)
        expected = (7 * 48) // 6
        assert compute_wiener_index(graph) == expected

    def test_linear_formula_n8(self):
        """n=8 should give (8*63)/6 = 84"""
        graph = MoleculeGraph.create_linear_alkane(8)
        expected = (8 * 63) // 6
        assert compute_wiener_index(graph) == expected

    def test_linear_formula_n10(self):
        """n=10 should give (10*99)/6 = 165"""
        graph = MoleculeGraph.create_linear_alkane(10)
        expected = (10 * 99) // 6
        assert compute_wiener_index(graph) == expected


class TestEdgeCases:
    """Test edge cases and error handling"""

    def test_disconnected_raises(self):
        """Disconnected graph should raise error"""
        # Will test this once we can create disconnected graphs
        pass
