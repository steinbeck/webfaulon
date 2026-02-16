"""Tests for SVG renderer service.

Following TDD RED-GREEN-REFACTOR pattern.
These tests are written FIRST and will initially fail.
"""

import pytest
import time

from app.services.svg_renderer import generate_molecule_svg


# Valid SMILES tests
def test_generates_svg_from_valid_smiles():
    """Should generate SVG string from valid SMILES."""
    svg = generate_molecule_svg("CCCCCC")

    assert isinstance(svg, str)
    assert "<svg" in svg
    assert "</svg>" in svg


def test_svg_contains_drawing_elements():
    """SVG should contain drawing elements (not empty canvas)."""
    svg = generate_molecule_svg("CCCCCC")

    # RDKit SVG should contain path or ellipse elements
    assert "path" in svg or "ellipse" in svg


def test_custom_dimensions():
    """Should generate SVG with custom dimensions."""
    svg = generate_molecule_svg("CCCCCC", width=400, height=400)

    # Check for dimension attributes in SVG
    assert "400" in svg


def test_default_dimensions():
    """Should use 300x300 when no dimensions provided."""
    svg = generate_molecule_svg("CCCCCC")

    # Check for default 300x300 dimensions
    assert "300" in svg


def test_aromatic_smiles():
    """Should handle aromatic SMILES (benzene)."""
    svg = generate_molecule_svg("c1ccccc1")

    assert isinstance(svg, str)
    assert "<svg" in svg
    assert "</svg>" in svg


def test_branched_molecule():
    """Should handle branched molecules."""
    svg = generate_molecule_svg("CC(C)CC")

    assert isinstance(svg, str)
    assert "<svg" in svg
    assert "</svg>" in svg


# Error handling tests
def test_invalid_smiles_raises_valueerror():
    """Should raise ValueError for invalid SMILES."""
    with pytest.raises(ValueError, match="Invalid SMILES"):
        generate_molecule_svg("INVALID")


def test_empty_smiles_raises_valueerror():
    """Should raise ValueError for empty SMILES."""
    with pytest.raises(ValueError, match="Empty SMILES"):
        generate_molecule_svg("")


# Performance test
def test_svg_generation_under_50ms():
    """SVG generation should be fast enough for SSE streaming."""
    start = time.time()
    svg = generate_molecule_svg("CCCCCC")
    elapsed = time.time() - start

    assert elapsed < 0.05  # 50ms
    assert "<svg" in svg  # Sanity check
