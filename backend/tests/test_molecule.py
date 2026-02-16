"""Tests for MoleculeGraph RDKit wrapper (ported from MolGraph.test.ts)"""
import pytest
from app.core.molecule import MoleculeGraph


class TestMoleculeGraphConstruction:
    """Test basic construction and factory methods"""

    def test_create_linear_alkane_methane(self):
        """Single carbon atom with 4 implicit hydrogens"""
        graph = MoleculeGraph.create_linear_alkane(1)
        assert graph.get_atom_count() == 1
        assert graph.get_atom_element(0) == 'C'
        assert graph.get_implicit_h(0) == 4  # No bonds, so 4 H

    def test_create_linear_alkane_ethane(self):
        """Two carbons with single bond, 3 implicit H each"""
        graph = MoleculeGraph.create_linear_alkane(2)
        assert graph.get_atom_count() == 2
        assert graph.get_bond_order(0, 1) == 1
        assert graph.get_implicit_h(0) == 3  # 1 bond, so 3 H
        assert graph.get_implicit_h(1) == 3

    def test_create_linear_alkane_pentane(self):
        """Pentane with correct implicit H distribution"""
        graph = MoleculeGraph.create_linear_alkane(5)
        assert graph.get_atom_count() == 5
        # Terminal carbons
        assert graph.get_implicit_h(0) == 3
        assert graph.get_implicit_h(4) == 3
        # Internal carbons
        assert graph.get_implicit_h(1) == 2
        assert graph.get_implicit_h(2) == 2
        assert graph.get_implicit_h(3) == 2

    def test_create_linear_alkane_decane(self):
        """Decane is connected"""
        graph = MoleculeGraph.create_linear_alkane(10)
        assert graph.get_atom_count() == 10
        assert graph.is_connected()


class TestCyclohexaneFactory:
    """Test cyclohexane factory method"""

    def test_create_cyclohexane(self):
        """6 carbons in a ring"""
        graph = MoleculeGraph.create_cyclohexane()
        assert graph.get_atom_count() == 6

    def test_cyclohexane_bond_orders(self):
        """Each carbon bonded to 2 neighbors with 2 implicit H"""
        graph = MoleculeGraph.create_cyclohexane()
        for i in range(6):
            assert graph.get_bond_order_sum(i) == 2
            assert graph.get_implicit_h(i) == 2  # 4 - 2 = 2

    def test_cyclohexane_connected(self):
        """Cyclohexane is connected"""
        graph = MoleculeGraph.create_cyclohexane()
        assert graph.is_connected()


class TestBranchedFactories:
    """Test branched structure factory methods"""

    def test_create_isobutane(self):
        """Isobutane: central carbon bonded to 3 others"""
        graph = MoleculeGraph.create_branched('isobutane')
        assert graph.get_atom_count() == 4
        assert graph.get_bond_order_sum(0) == 3  # Central carbon
        assert graph.get_implicit_h(0) == 1  # 4 - 3 = 1
        # Terminal carbons
        assert graph.get_bond_order_sum(1) == 1
        assert graph.get_implicit_h(1) == 3

    def test_create_neopentane(self):
        """Neopentane: central carbon bonded to 4 others"""
        graph = MoleculeGraph.create_branched('neopentane')
        assert graph.get_atom_count() == 5
        assert graph.get_bond_order_sum(0) == 4  # Central carbon
        assert graph.get_implicit_h(0) == 0  # 4 - 4 = 0
        # Terminal carbons
        for i in range(1, 5):
            assert graph.get_bond_order_sum(i) == 1
            assert graph.get_implicit_h(i) == 3


class TestConnectivity:
    """Test connectivity validation"""

    def test_connectivity_connected(self):
        """Linear pentane is connected"""
        graph = MoleculeGraph.create_linear_alkane(5)
        assert graph.is_connected()

    def test_connectivity_disconnected(self):
        """Two separate fragments are not connected"""
        # Create two separate ethane molecules (not implemented via factory)
        # Will test this once we have manual construction
        pass

    def test_connectivity_single_atom(self):
        """Single atom is connected"""
        graph = MoleculeGraph.create_linear_alkane(1)
        assert graph.is_connected()


class TestValidValences:
    """Test valence validation"""

    def test_valid_valences_alkane(self):
        """Linear pentane has valid valences"""
        graph = MoleculeGraph.create_linear_alkane(5)
        assert graph.has_valid_valences()


class TestSetBond:
    """Test bond mutation"""

    def test_set_bond_updates_order(self):
        """Setting bond from 1 to 2 changes bond order"""
        graph = MoleculeGraph.create_linear_alkane(3)
        assert graph.get_bond_order(0, 1) == 1

        # Change to double bond
        graph.set_bond(0, 1, 2)
        assert graph.get_bond_order(0, 1) == 2

    def test_set_bond_symmetric(self):
        """Setting bond(0,1,2) means both get_bond_order(0,1) and (1,0) return 2"""
        graph = MoleculeGraph.create_linear_alkane(3)
        graph.set_bond(0, 1, 2)
        assert graph.get_bond_order(0, 1) == 2
        assert graph.get_bond_order(1, 0) == 2  # Symmetric


class TestClone:
    """Test cloning functionality"""

    def test_clone_independent(self):
        """Modifying clone doesn't affect original"""
        original = MoleculeGraph.create_linear_alkane(3)
        clone = original.clone()

        # Modify clone
        clone.set_bond(0, 1, 2)

        # Original should be unchanged
        assert original.get_bond_order(0, 1) == 1
        assert clone.get_bond_order(0, 1) == 2


class TestToSmiles:
    """Test SMILES generation"""

    def test_to_smiles(self):
        """Linear hexane produces valid canonical SMILES"""
        graph = MoleculeGraph.create_linear_alkane(6)
        smiles = graph.to_smiles()
        assert smiles is not None
        assert len(smiles) > 0
        # Should be CCCCCC (canonical for hexane)
        assert 'C' in smiles


class TestGetImplicitH:
    """Test implicit hydrogen counts"""

    def test_get_implicit_h(self):
        """Verify implicit hydrogen counts match TypeScript values"""
        graph = MoleculeGraph.create_linear_alkane(5)
        # Terminal carbons: 3 H
        assert graph.get_implicit_h(0) == 3
        assert graph.get_implicit_h(4) == 3
        # Internal carbons: 2 H
        assert graph.get_implicit_h(1) == 2
        assert graph.get_implicit_h(2) == 2
        assert graph.get_implicit_h(3) == 2
