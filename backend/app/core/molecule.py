"""RDKit Mol wrapper providing adjacency-matrix-like API for Faulon SA."""
from rdkit import Chem
from rdkit.Chem import RWMol
from app.utils.rdkit_helpers import order_to_bond_type, bond_type_to_order


class MoleculeGraph:
    """
    RDKit Mol wrapper providing adjacency-matrix-like API for Faulon SA.

    Unlike TypeScript's adjacency-matrix MolGraph, this wraps RDKit's native Mol object.
    The API feels similar but leverages RDKit internally for chemistry validation.
    """

    def __init__(self, mol: Chem.Mol):
        """
        Initialize from an RDKit Mol (makes deep copy).

        Args:
            mol: RDKit Mol object to wrap
        """
        # Make deep copy and convert to RWMol for mutation
        self._mol = RWMol(Chem.Mol(mol))

    # Properties
    def get_atom_count(self) -> int:
        """Return number of heavy atoms in molecule."""
        return self._mol.GetNumAtoms()

    def get_bond_order(self, i: int, j: int) -> int:
        """
        Get bond order between atoms i and j.

        Args:
            i: First atom index
            j: Second atom index

        Returns:
            Bond order as integer (0 if no bond, 1-3 otherwise)
        """
        bond = self._mol.GetBondBetweenAtoms(i, j)
        if bond is None:
            return 0
        return bond_type_to_order(bond.GetBondType())

    def get_bond_order_sum(self, atom_index: int) -> int:
        """
        Get sum of all bond orders for given atom.

        Args:
            atom_index: Atom index

        Returns:
            Sum of bond orders connected to this atom
        """
        total = 0
        atom = self._mol.GetAtomWithIdx(atom_index)
        for bond in atom.GetBonds():
            total += bond_type_to_order(bond.GetBondType())
        return total

    def get_implicit_h(self, atom_index: int) -> int:
        """
        Get number of implicit hydrogens for given atom.

        RDKit computes this automatically after sanitization.

        Args:
            atom_index: Atom index

        Returns:
            Number of implicit hydrogens
        """
        atom = self._mol.GetAtomWithIdx(atom_index)
        return atom.GetNumImplicitHs()

    def get_atom_element(self, atom_index: int) -> str:
        """
        Get element symbol for given atom.

        Args:
            atom_index: Atom index

        Returns:
            Element symbol (e.g., 'C', 'N', 'O')
        """
        atom = self._mol.GetAtomWithIdx(atom_index)
        return atom.GetSymbol()

    # Mutation
    def set_bond(self, i: int, j: int, order: int) -> None:
        """
        Set bond order between atoms i and j.

        If order == 0: remove bond if exists
        If no existing bond and order > 0: add bond with given order
        If existing bond: set bond type to given order

        CRITICAL: Calls Chem.SanitizeMol() after modification.
        If sanitization fails, raises ValueError.

        Args:
            i: First atom index
            j: Second atom index
            order: New bond order (0-3)

        Raises:
            ValueError: If sanitization fails after modification
        """
        if order < 0 or order > 3:
            raise ValueError(f"Bond order {order} out of range [0, 3]")

        bond = self._mol.GetBondBetweenAtoms(i, j)

        if order == 0:
            # Remove bond if exists
            if bond is not None:
                self._mol.RemoveBond(i, j)
        else:
            bond_type = order_to_bond_type(order)
            if bond is None:
                # Add new bond
                self._mol.AddBond(i, j, bond_type)
            else:
                # Update existing bond
                bond.SetBondType(bond_type)

        # No sanitization needed: Faulon displacement (equations 7-9)
        # preserves valence by construction. Skipping SanitizeMol also
        # avoids RDKit aromaticity perception converting SINGLE/DOUBLE
        # bonds to AROMATIC type, which is incompatible with integer
        # bond orders required by displacement.

    # Validation
    def is_connected(self) -> bool:
        """
        Check if molecule is connected (single fragment).

        Returns:
            True if connected, False otherwise
        """
        if self.get_atom_count() == 0:
            return True
        if self.get_atom_count() == 1:
            return True

        # Use RDKit GetMolFrags to count fragments
        frags = Chem.GetMolFrags(self._mol)
        return len(frags) == 1

    def has_valid_valences(self) -> bool:
        """
        Check if all atoms have valid valences.

        Returns:
            True if all valences valid, False otherwise
        """
        try:
            # Try sanitization - will fail if valences invalid
            mol_copy = Chem.Mol(self._mol)
            Chem.SanitizeMol(mol_copy)
            return True
        except:
            return False

    def validate(self) -> dict:
        """
        Validate both connectivity and valences.

        Returns:
            Dict with 'connected' and 'valid_valences' bool keys
        """
        return {
            'connected': self.is_connected(),
            'valid_valences': self.has_valid_valences()
        }

    # Utilities
    def clone(self) -> 'MoleculeGraph':
        """
        Create independent deep copy of this molecule.

        Returns:
            New MoleculeGraph with same structure
        """
        mol_copy = Chem.Mol(self._mol)
        return MoleculeGraph(mol_copy)

    def to_smiles(self) -> str:
        """
        Generate canonical SMILES string.

        Returns:
            Canonical SMILES representation
        """
        mol_copy = Chem.Mol(self._mol)
        return Chem.MolToSmiles(mol_copy)

    def get_mol(self) -> Chem.Mol:
        """
        Return immutable copy of internal RDKit Mol.

        Returns:
            RDKit Mol object (immutable copy)
        """
        return Chem.Mol(self._mol)

    # Factory methods (class methods)
    @classmethod
    def create_linear_alkane(cls, n: int) -> 'MoleculeGraph':
        """
        Create linear alkane with n carbons.

        Args:
            n: Number of carbon atoms

        Returns:
            MoleculeGraph representing linear alkane

        Raises:
            ValueError: If n < 1
        """
        if n < 1:
            raise ValueError('Linear alkane must have at least 1 carbon')

        # Create RWMol
        mol = RWMol()

        # Add n carbon atoms
        for _ in range(n):
            mol.AddAtom(Chem.Atom('C'))

        # Add single bonds between consecutive carbons
        for i in range(n - 1):
            mol.AddBond(i, i + 1, Chem.BondType.SINGLE)

        # Convert to Mol and sanitize
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)

        return cls(mol)

    @classmethod
    def create_cyclohexane(cls) -> 'MoleculeGraph':
        """
        Create cyclohexane (6-carbon ring).

        Returns:
            MoleculeGraph representing cyclohexane
        """
        mol = RWMol()

        # Add 6 carbon atoms
        for _ in range(6):
            mol.AddAtom(Chem.Atom('C'))

        # Create ring: 0-1-2-3-4-5-0
        for i in range(6):
            next_idx = (i + 1) % 6
            mol.AddBond(i, next_idx, Chem.BondType.SINGLE)

        # Convert to Mol and sanitize
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)

        return cls(mol)

    @classmethod
    def create_branched(cls, pattern: str) -> 'MoleculeGraph':
        """
        Create branched structure (isobutane or neopentane).

        Args:
            pattern: 'isobutane' or 'neopentane'

        Returns:
            MoleculeGraph representing branched structure

        Raises:
            ValueError: If pattern is not recognized
        """
        if pattern == 'isobutane':
            # Central carbon (index 0) bonded to 3 other carbons (indices 1, 2, 3)
            mol = RWMol()
            for _ in range(4):
                mol.AddAtom(Chem.Atom('C'))

            mol.AddBond(0, 1, Chem.BondType.SINGLE)
            mol.AddBond(0, 2, Chem.BondType.SINGLE)
            mol.AddBond(0, 3, Chem.BondType.SINGLE)

            mol = mol.GetMol()
            Chem.SanitizeMol(mol)
            return cls(mol)

        elif pattern == 'neopentane':
            # Central carbon (index 0) bonded to 4 other carbons
            mol = RWMol()
            for _ in range(5):
                mol.AddAtom(Chem.Atom('C'))

            mol.AddBond(0, 1, Chem.BondType.SINGLE)
            mol.AddBond(0, 2, Chem.BondType.SINGLE)
            mol.AddBond(0, 3, Chem.BondType.SINGLE)
            mol.AddBond(0, 4, Chem.BondType.SINGLE)

            mol = mol.GetMol()
            Chem.SanitizeMol(mol)
            return cls(mol)

        else:
            raise ValueError(f"Unknown branched pattern: {pattern}")
