"""
Molecular formula parser and HDI (Hydrogen Deficiency Index) computation.

This module provides utilities for parsing molecular formulas and computing
the degree of unsaturation (HDI), which indicates the number of rings and/or
multiple bonds in a molecule.
"""

import re
from typing import Dict

# Type alias for formula maps
FormulaMap = Dict[str, int]

# Known elements supported by the parser
KNOWN_ELEMENTS = {'H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'}

# Halogen elements (count as hydrogen equivalents in HDI calculation)
HALOGENS = {'F', 'Cl', 'Br', 'I'}


def parse_formula(formula: str) -> FormulaMap:
    """
    Parse a molecular formula string into a map of element symbols to counts.

    Examples:
        - "C6H14" -> {"C": 6, "H": 14}
        - "CH4" -> {"C": 1, "H": 4}
        - "C2H6O" -> {"C": 2, "H": 6, "O": 1}

    Args:
        formula: Molecular formula string (e.g., "C6H14")

    Returns:
        Dictionary mapping element symbols to counts

    Raises:
        ValueError: If formula is empty or contains unknown elements
    """
    if not formula or formula.strip() == '':
        raise ValueError('Empty formula')

    result: FormulaMap = {}

    # Regex to match element symbol (capital + optional lowercase) followed by optional digits
    # Example: C6H14 -> [("C", "6"), ("H", "14")]
    # Example: CH4 -> [("C", ""), ("H", "4")]
    pattern = r'([A-Z][a-z]?)(\d*)'
    matches = re.findall(pattern, formula)

    for element, count_str in matches:
        # Skip empty matches (happens at end of string)
        if not element:
            continue

        # Validate element is known
        if element not in KNOWN_ELEMENTS:
            raise ValueError(f'Unknown element: {element}')

        # Parse count (default to 1 if not specified)
        count = int(count_str) if count_str else 1

        # Add to result (accumulate if element appears multiple times)
        result[element] = result.get(element, 0) + count

    return result


def compute_hdi(formula: FormulaMap) -> int:
    """
    Compute the Hydrogen Deficiency Index (HDI) for a molecular formula.

    HDI = (2C + 2 + N - H - Halogens) / 2

    HDI indicates the number of rings and/or multiple bonds:
    - HDI = 0: saturated (no rings or double bonds)
    - HDI = 1: one ring or one double bond
    - HDI = 2: two rings, one triple bond, or one ring + one double bond
    - etc.

    Examples:
        - C6H14 (hexane): HDI = 0
        - C6H12 (cyclohexane or 1-hexene): HDI = 1
        - C6H6 (benzene): HDI = 4

    Args:
        formula: Parsed formula map from parse_formula()

    Returns:
        Hydrogen Deficiency Index (non-negative integer)

    Note:
        Oxygen and sulfur don't affect HDI (they're divalent like C in the formula).
        Halogens count as hydrogen equivalents.
    """
    C = formula.get('C', 0)
    H = formula.get('H', 0)
    N = formula.get('N', 0)

    # Count halogens
    halogens = sum(formula.get(element, 0) for element in HALOGENS)

    # HDI = (2C + 2 + N - H - Halogens) / 2
    # Oxygen and sulfur don't affect HDI (they're divalent)
    hdi = (2 * C + 2 + N - H - halogens) / 2

    # HDI can't be negative (clamp to 0)
    return max(0, int(hdi))
