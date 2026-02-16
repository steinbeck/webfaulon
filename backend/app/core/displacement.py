"""
Faulon displacement operation (equations 7-11 from 1996 paper).

The core SA mutation: select 4 atoms, redistribute bond orders
while preserving valences and connectivity.

Reference: J. Chem. Inf. Comput. Sci. 1996, 36, 731-740, page 733
"""
from app.core.molecule import MoleculeGraph
from app.core.random import SeededRandom

# Maximum bond order (triple bond)
MAX_BOND_ORDER = 3

# Minimum number of atoms required for displacement
MIN_ATOMS_FOR_DISPLACEMENT = 4

# Maximum attempts to find a non-trivial displacement before giving up.
# Sparse graphs (long chains) often produce no-op displacements where
# b == a because the randomly selected atoms share no bonds.
MAX_DISPLACEMENT_ATTEMPTS = 10


def compute_displacement_bonds(
    a11: int, a12: int, a21: int, a22: int, rng: SeededRandom
) -> dict | None:
    """
    Compute new bond orders via Faulon equations 7-11.

    Args:
        a11: Current bond order between x1 and y1
        a12: Current bond order between y1 and y2
        a21: Current bond order between x1 and x2
        a22: Current bond order between x2 and y2
        rng: Random number generator for choosing b11

    Returns:
        Dict with b11, b12, b21, b22 keys, or None if no valid displacement
    """
    # Equation 10: b11 >= MAX(0, a11-a22, a11-a12, a11+a12-3, a11+a21-3)
    b11_min = max(
        0,
        a11 - a22,
        a11 - a12,
        a11 + a12 - MAX_BOND_ORDER,
        a11 + a21 - MAX_BOND_ORDER
    )

    # Equation 11: b11 <= MIN(3, a11+a12, a11+a21, a11-a22+3)
    b11_max = min(
        MAX_BOND_ORDER,
        a11 + a12,
        a11 + a21,
        a11 - a22 + MAX_BOND_ORDER
    )

    # If no valid range, displacement is impossible for this atom selection
    if b11_min > b11_max:
        return None

    # Choose b11 randomly from valid range
    b11 = rng.next_int(b11_min, b11_max)

    # Compute remaining new bond orders using equations 7-9
    # Equation 7: b12 = a11 + a12 - b11
    # Equation 8: b21 = a11 + a21 - b11
    # Equation 9: b22 = a22 - a11 + b11
    b12 = a11 + a12 - b11
    b21 = a11 + a21 - b11
    b22 = a22 - a11 + b11

    # Verify all bond orders are in valid range [0, MAX_BOND_ORDER]
    # (Should be guaranteed by constraint equations, but verify as safety check)
    if (
        b11 < 0 or b11 > MAX_BOND_ORDER or
        b12 < 0 or b12 > MAX_BOND_ORDER or
        b21 < 0 or b21 > MAX_BOND_ORDER or
        b22 < 0 or b22 > MAX_BOND_ORDER
    ):
        raise RuntimeError(
            f"Bond order out of range after displacement: b11={b11}, b12={b12}, "
            f"b21={b21}, b22={b22}. This indicates a bug in the displacement equations."
        )

    return {"b11": b11, "b12": b12, "b21": b21, "b22": b22}


def attempt_displacement(
    mol_graph: MoleculeGraph,
    rng: SeededRandom
) -> MoleculeGraph | None:
    """
    Attempt a Faulon displacement on the given graph.

    Retries up to MAX_DISPLACEMENT_ATTEMPTS times to find a non-trivial
    displacement (one that actually changes bond orders). Sparse graphs
    like long chains often produce no-op displacements where the randomly
    selected 4 atoms share no bonds.

    Faulon equations 7-9 preserve valence by construction, so no
    sanitization is needed after bond changes. Connectivity is NOT
    checked here -- the SA engine handles disconnected structures
    by persisting them and continuing to apply operators until
    reconnection occurs.

    Args:
        mol_graph: The molecular graph to displace
        rng: Seeded random number generator for reproducibility

    Returns:
        A new graph with bonds redistributed, or None if move is invalid
    """
    atom_count = mol_graph.get_atom_count()

    # Need at least 4 atoms for displacement
    if atom_count < MIN_ATOMS_FOR_DISPLACEMENT:
        return None

    for _ in range(MAX_DISPLACEMENT_ATTEMPTS):
        # Select 4 distinct atoms (Faulon paper p.733, Table 2, step 2.1)
        selected = rng.select_n_distinct(4, atom_count)
        x1, y1, x2, y2 = selected[0], selected[1], selected[2], selected[3]

        # Read current bond orders (a notation from paper)
        # Paper p.733: "Let a11, a12, a21, and a22 be the orders of
        # the bonds [x1,y1], [x1,y2], [x2,y1], and [x2,y2]"
        a11 = mol_graph.get_bond_order(x1, y1)
        a12 = mol_graph.get_bond_order(y1, y2)
        a21 = mol_graph.get_bond_order(x1, x2)
        a22 = mol_graph.get_bond_order(x2, y2)

        # Compute new bond orders using Faulon equations 7-11
        new_bonds = compute_displacement_bonds(a11, a12, a21, a22, rng)

        # If no valid displacement exists for this atom selection, retry
        if new_bonds is None:
            continue

        # Skip no-op displacements (b == a, nothing changes)
        if (new_bonds["b11"] == a11 and new_bonds["b12"] == a12
                and new_bonds["b21"] == a21 and new_bonds["b22"] == a22):
            continue

        # Clone the graph (don't mutate original)
        new_graph = mol_graph.clone()

        # Apply new bond orders (no sanitization -- Faulon preserves valence)
        new_graph.set_bond(x1, y1, new_bonds["b11"])
        new_graph.set_bond(y1, y2, new_bonds["b12"])
        new_graph.set_bond(x1, x2, new_bonds["b21"])
        new_graph.set_bond(x2, y2, new_bonds["b22"])

        return new_graph

    # All attempts produced no-ops or invalid ranges
    return None
