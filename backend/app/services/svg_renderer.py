"""SVG rendering service for molecular structures.

Uses RDKit MolDraw2DSVG to generate 2D depictions from SMILES strings.
"""

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def generate_molecule_svg(smiles: str, width: int = 300, height: int = 300) -> str:
    """Generate 2D SVG depiction of molecule from SMILES.

    Args:
        smiles: SMILES notation string
        width: SVG width in pixels (default 300)
        height: SVG height in pixels (default 300)

    Returns:
        SVG string (inline, embeddable in HTML or SSE events)

    Raises:
        ValueError: If SMILES is invalid or empty
    """
    if not smiles:
        raise ValueError("Empty SMILES string")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    return drawer.GetDrawingText()
