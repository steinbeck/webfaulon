"""RDKit helper utilities for bond type conversion and sanitization."""
from rdkit import Chem
from rdkit.Chem import BondType


def order_to_bond_type(order: int) -> BondType:
    """
    Convert Faulon integer bond order (0-3) to RDKit BondType enum.

    Args:
        order: Bond order as integer (0=no bond, 1=single, 2=double, 3=triple)

    Returns:
        RDKit BondType enum value

    Raises:
        ValueError: If order is not in range [0, 3]
    """
    # Explicit dict mapping - NOT direct cast
    order_map = {
        0: BondType.UNSPECIFIED,  # No bond
        1: BondType.SINGLE,
        2: BondType.DOUBLE,
        3: BondType.TRIPLE,
    }

    if order not in order_map:
        raise ValueError(f"Invalid bond order: {order}. Must be in range [0, 3]")

    return order_map[order]


def bond_type_to_order(bond_type: BondType) -> int:
    """
    Convert RDKit BondType enum to Faulon integer bond order.

    Args:
        bond_type: RDKit BondType enum value

    Returns:
        Bond order as integer (0-3)

    Raises:
        ValueError: If bond_type cannot be mapped to integer order
    """
    type_map = {
        BondType.UNSPECIFIED: 0,
        BondType.SINGLE: 1,
        BondType.DOUBLE: 2,
        BondType.TRIPLE: 3,
    }

    if bond_type not in type_map:
        raise ValueError(f"Cannot map bond type {bond_type} to integer order")

    return type_map[bond_type]
