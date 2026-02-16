"""Tests for formula parser - ported from src/core/__tests__/formulaParser.test.ts"""

import pytest
import sys
from pathlib import Path

# Add backend directory to path to import app modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.core.formula_parser import parse_formula, compute_hdi


class TestParseFormula:
    """Test parse_formula function"""

    def test_parse_simple_alkanes(self):
        """Parse simple alkane formulas"""
        assert parse_formula('C6H14') == {'C': 6, 'H': 14}
        assert parse_formula('C8H10') == {'C': 8, 'H': 10}
        assert parse_formula('CH4') == {'C': 1, 'H': 4}
        assert parse_formula('C4H10') == {'C': 4, 'H': 10}

    def test_parse_implicit_count_1(self):
        """Handle implicit count of 1"""
        assert parse_formula('CH4') == {'C': 1, 'H': 4}
        assert parse_formula('CO2') == {'C': 1, 'O': 2}

    def test_parse_heteroatoms(self):
        """Parse formulas with heteroatoms"""
        assert parse_formula('C2H6O') == {'C': 2, 'H': 6, 'O': 1}
        assert parse_formula('C6H5Cl') == {'C': 6, 'H': 5, 'Cl': 1}
        assert parse_formula('C3H7N') == {'C': 3, 'H': 7, 'N': 1}
        assert parse_formula('C2H6S') == {'C': 2, 'H': 6, 'S': 1}

    def test_parse_multiple_heteroatoms(self):
        """Parse complex formulas with multiple heteroatoms"""
        assert parse_formula('C2H5NO2') == {'C': 2, 'H': 5, 'N': 1, 'O': 2}
        assert parse_formula('C3H9N3') == {'C': 3, 'H': 9, 'N': 3}

    def test_parse_empty_raises(self):
        """Empty string raises ValueError"""
        with pytest.raises(ValueError, match="Empty formula"):
            parse_formula('')

    def test_parse_unknown_element_raises(self):
        """Unknown element raises ValueError with element name"""
        with pytest.raises(ValueError, match="Unknown element: X"):
            parse_formula('XYZ')
        with pytest.raises(ValueError, match="Unknown element: Q"):
            parse_formula('C6Q14')

    def test_parse_double_digit_counts(self):
        """Handle formulas with double-digit counts"""
        assert parse_formula('C10H22') == {'C': 10, 'H': 22}
        assert parse_formula('C12H26') == {'C': 12, 'H': 26}

    def test_parse_order_independence(self):
        """Different orderings produce same result"""
        result1 = parse_formula('C6H14')
        result2 = parse_formula('H14C6')
        assert result1 == {'C': 6, 'H': 14}
        assert result2 == {'C': 6, 'H': 14}


class TestComputeHDI:
    """Test compute_hdi function"""

    def test_hdi_saturated_alkanes(self):
        """HDI for saturated alkanes is 0"""
        assert compute_hdi({'C': 6, 'H': 14}) == 0  # hexane
        assert compute_hdi({'C': 4, 'H': 10}) == 0  # butane
        assert compute_hdi({'C': 1, 'H': 4}) == 0   # methane

    def test_hdi_unsaturated(self):
        """HDI for unsaturated hydrocarbons"""
        assert compute_hdi({'C': 6, 'H': 12}) == 1  # one ring or double bond
        assert compute_hdi({'C': 6, 'H': 6}) == 4   # benzene
        assert compute_hdi({'C': 2, 'H': 4}) == 1   # ethene
        assert compute_hdi({'C': 2, 'H': 2}) == 2   # ethyne

    def test_hdi_with_oxygen(self):
        """HDI for molecules with oxygen"""
        assert compute_hdi({'C': 2, 'H': 6, 'O': 1}) == 0  # ethanol
        assert compute_hdi({'C': 2, 'H': 4, 'O': 1}) == 1  # acetaldehyde

    def test_hdi_with_nitrogen(self):
        """HDI for molecules with nitrogen"""
        # HDI = (2C + 2 + N - H - Halogens) / 2
        assert compute_hdi({'C': 2, 'H': 7, 'N': 1}) == 0  # ethylamine (C2H5NH2)
        assert compute_hdi({'C': 2, 'H': 5, 'N': 1}) == 1  # acetonitrile (CH3CN)

    def test_hdi_with_halogens(self):
        """HDI for molecules with halogens"""
        assert compute_hdi({'C': 6, 'H': 5, 'Cl': 1}) == 4  # chlorobenzene
        assert compute_hdi({'C': 2, 'H': 5, 'Br': 1}) == 0  # bromoethane

    def test_hdi_multiple_heteroatoms(self):
        """HDI with multiple heteroatoms"""
        assert compute_hdi({'C': 2, 'H': 5, 'N': 1, 'O': 2}) == 1  # glycine

    def test_hdi_edge_cases(self):
        """Edge cases for HDI"""
        assert compute_hdi({'C': 1, 'H': 4}) == 0  # methane
        assert compute_hdi({'H': 2}) == 0  # H2
