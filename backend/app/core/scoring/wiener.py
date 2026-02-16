"""Wiener Index scoring component."""
from rdkit import Chem
from app.core.molecule import MoleculeGraph


class WienerIndexComponent:
    """
    Wiener Index scoring component (sum of pairwise distances).

    The Wiener Index is the sum of all pairwise distances between heavy atoms
    in the molecular graph. It's a topological descriptor used in QSAR studies.

    This component wraps the same algorithm as the standalone compute_wiener_index()
    function, refactored into the ScoringComponent protocol interface.
    """
    name = "wiener_index"

    def compute(self, mol_graph: MoleculeGraph) -> float:
        """
        Compute Wiener Index using RDKit GetDistanceMatrix.

        Args:
            mol_graph: MoleculeGraph to compute Wiener Index for

        Returns:
            Wiener Index (sum of all pairwise distances)

        Raises:
            ValueError: If graph is disconnected
        """
        mol = mol_graph.get_mol()
        n = mol.GetNumAtoms()

        # Edge case: single atom has Wiener index 0
        if n <= 1:
            return 0

        # Use RDKit's built-in distance matrix computation
        dist_matrix = Chem.GetDistanceMatrix(mol)

        # Check for disconnected graph (distance matrix will have large/inf values)
        for i in range(n):
            for j in range(n):
                if dist_matrix[i][j] > 1000:  # Effectively infinite
                    raise ValueError('Graph is disconnected')

        # Sum upper triangle (avoid double-counting)
        wiener = sum(
            dist_matrix[i][j]
            for i in range(n)
            for j in range(i + 1, n)
        )

        return wiener
