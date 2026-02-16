"""
Generate deterministic initial molecular structure from a molecular formula.

Algorithm:
1. Parse formula to get heavy atom counts (C, N, O, etc.)
2. Create RDKit RWMol with heavy atoms in deterministic order
3. Connect all heavy atoms in a linear chain with single bonds
4. Compute HDI from formula
5. Add unsaturation (double/triple bonds) to satisfy HDI
6. Validate that implicit H matches formula H count

This produces a valid, connected starting structure for SA optimization.
The specific arrangement doesn't matter - SA will rearrange via displacement.
"""
from rdkit import Chem
from rdkit.Chem import RWMol
from app.core.molecule import MoleculeGraph
from app.core.formula_parser import parse_formula, compute_hdi

# Standard valences for common elements
STANDARD_VALENCE = {
    'H': 1,
    'C': 4,
    'N': 3,
    'O': 2,
    'S': 2,
    'P': 3,
    'F': 1,
    'Cl': 1,
    'Br': 1,
    'I': 1,
}


def generate_initial_structure(formula: str) -> MoleculeGraph:
    """
    Generate a deterministic initial molecular structure from a molecular formula.

    Args:
        formula: Molecular formula string (e.g., "C6H14")

    Returns:
        Connected MoleculeGraph with valid valences

    Raises:
        ValueError: If formula is invalid or impossible to satisfy
    """
    # Parse formula
    formula_map = parse_formula(formula)
    target_h = formula_map.get('H', 0)
    hdi = compute_hdi(formula_map)

    # Extract heavy atoms (all except H)
    heavy_atoms = []

    # Add atoms in deterministic order: C, N, O, S, P, then halogens
    element_order = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']

    for element in element_order:
        count = formula_map.get(element, 0)
        for _ in range(count):
            heavy_atoms.append(element)

    # Handle edge case: no heavy atoms (e.g., H2)
    if len(heavy_atoms) == 0:
        # Return empty graph (all hydrogens are implicit)
        mol = Chem.MolFromSmiles('')  # Empty molecule
        return MoleculeGraph(mol)

    # Handle edge case: single atom
    if len(heavy_atoms) == 1:
        mol = RWMol()
        mol.AddAtom(Chem.Atom(heavy_atoms[0]))
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)
        graph = MoleculeGraph(mol)

        # Verify hydrogen count
        actual_h = graph.get_implicit_h(0)
        if actual_h != target_h:
            raise ValueError(
                f"Cannot generate structure for {formula}: expected {target_h} H, got {actual_h}"
            )

        return graph

    # Create RWMol and add heavy atoms
    n = len(heavy_atoms)
    mol = RWMol()
    for element in heavy_atoms:
        mol.AddAtom(Chem.Atom(element))

    # Step 1: Connect atoms in a linear chain with single bonds
    for i in range(n - 1):
        mol.AddBond(i, i + 1, Chem.BondType.SINGLE)

    # Step 2: Add unsaturation to satisfy HDI
    # Each double bond adds 1 HDI, each triple bond adds 2 HDI
    # Strategy: upgrade bonds iteratively until HDI is satisfied
    remaining_hdi = hdi

    # Keep upgrading bonds until we satisfy HDI or can't upgrade any more
    while remaining_hdi > 0:
        upgraded = False

        # Try to upgrade each bond in the chain
        for i in range(n - 1):
            if remaining_hdi == 0:
                break

            j = i + 1
            bond = mol.GetBondBetweenAtoms(i, j)
            current_bond_order = _bond_type_to_order(bond.GetBondType())

            # Can't upgrade beyond triple bond
            if current_bond_order >= 3:
                continue

            # Calculate bond order sums for both atoms
            def get_bond_order_sum(mol, atom_index):
                total = 0
                atom = mol.GetAtomWithIdx(atom_index)
                for b in atom.GetBonds():
                    total += _bond_type_to_order(b.GetBondType())
                return total

            atom1_element = mol.GetAtomWithIdx(i).GetSymbol()
            atom2_element = mol.GetAtomWithIdx(j).GetSymbol()

            valence1 = STANDARD_VALENCE.get(atom1_element)
            valence2 = STANDARD_VALENCE.get(atom2_element)

            if valence1 is None or valence2 is None:
                raise ValueError(f"Unknown element valence")

            # Try upgrading by 1 (e.g., 1->2 or 2->3)
            new_bond_order = current_bond_order + 1

            # Calculate what the bond sums would be after upgrade
            current_bond_sum1 = get_bond_order_sum(mol, i)
            current_bond_sum2 = get_bond_order_sum(mol, j)

            new_bond_sum1 = current_bond_sum1 + 1  # Adding 1 to the bond
            new_bond_sum2 = current_bond_sum2 + 1

            implicit_h1 = valence1 - new_bond_sum1
            implicit_h2 = valence2 - new_bond_sum2

            # Check if both atoms can handle the upgrade
            if implicit_h1 >= 0 and implicit_h2 >= 0:
                bond.SetBondType(_order_to_bond_type(new_bond_order))
                remaining_hdi -= 1
                upgraded = True

        # If we couldn't upgrade any bond, we can't satisfy HDI
        if not upgraded:
            break

    # If we couldn't satisfy HDI, the formula might be impossible
    if remaining_hdi > 0:
        raise ValueError(
            f"Cannot generate structure for {formula}: unable to satisfy HDI={hdi}"
        )

    # Convert to Mol and sanitize
    mol = mol.GetMol()
    Chem.SanitizeMol(mol)

    # Create the graph
    graph = MoleculeGraph(mol)

    # Step 3: Verify hydrogen count
    actual_h = sum(graph.get_implicit_h(i) for i in range(n))

    if actual_h != target_h:
        raise ValueError(
            f"Cannot generate structure for {formula}: expected {target_h} H, got {actual_h}"
        )

    # Step 4: Validate structure
    if not graph.is_connected():
        raise ValueError(f"Generated structure is not connected")

    if not graph.has_valid_valences():
        raise ValueError(f"Generated structure has invalid valences")

    return graph


def _bond_type_to_order(bond_type: Chem.BondType) -> int:
    """Convert RDKit BondType to integer order."""
    if bond_type == Chem.BondType.SINGLE:
        return 1
    elif bond_type == Chem.BondType.DOUBLE:
        return 2
    elif bond_type == Chem.BondType.TRIPLE:
        return 3
    else:
        return 0


def _order_to_bond_type(order: int) -> Chem.BondType:
    """Convert integer order to RDKit BondType."""
    if order == 1:
        return Chem.BondType.SINGLE
    elif order == 2:
        return Chem.BondType.DOUBLE
    elif order == 3:
        return Chem.BondType.TRIPLE
    else:
        return Chem.BondType.UNSPECIFIED
