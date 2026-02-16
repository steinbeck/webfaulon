"""Protocol interface for scoring components."""
from typing import Protocol, runtime_checkable
from app.core.molecule import MoleculeGraph


@runtime_checkable
class ScoringComponent(Protocol):
    """
    Structural interface for scoring components.

    All scoring components must provide:
    - name: str attribute identifying the component
    - compute(mol_graph) -> float method that scores a molecule
    """
    name: str

    def compute(self, mol_graph: MoleculeGraph) -> float:
        """
        Compute score for the given molecular graph.

        Args:
            mol_graph: MoleculeGraph to score

        Returns:
            Float score value
        """
        ...
