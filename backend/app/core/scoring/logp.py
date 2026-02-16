"""LogP (partition coefficient) scoring component."""
from app.core.molecule import MoleculeGraph
from rdkit import Chem
from rdkit.Chem import Descriptors


class LogPComponent:
    """
    Wildman-Crippen LogP scoring component (partition coefficient).

    LogP varies across constitutional isomers because it depends on atom
    connectivity and local environment, not just molecular formula.
    This makes it suitable for SA optimization, unlike molecular weight
    which is constant for all isomers of a given formula.
    """
    name = "logp"

    def compute(self, mol_graph: MoleculeGraph) -> float:
        """
        Compute LogP using RDKit Descriptors.MolLogP (Wildman-Crippen method).

        Args:
            mol_graph: MoleculeGraph to compute LogP for

        Returns:
            LogP value (dimensionless partition coefficient)
        """
        mol = mol_graph.get_mol()
        # Ensure implicit valences are calculated (required by MolLogP)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        return Descriptors.MolLogP(mol)
