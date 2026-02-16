"""
Integration tests for multi-component SA engine.

Tests that SAEngine and SAParams work correctly with the scoring
component framework (protocol, registry, weighted sums).
"""

import pytest
from pydantic import ValidationError
from rdkit import Chem
from rdkit.Chem import Descriptors

from app.core.sa_engine import SAEngine
from app.models.sa_params import SAParams


class TestSAParamsComponentWeights:
    """Test SAParams component_weights field validation."""

    def test_default_weights(self):
        """SAParams with no component_weights defaults to wiener_index=1.0."""
        params = SAParams(formula="C6H14")
        assert params.component_weights == {"wiener_index": 1.0}

    def test_custom_weights(self):
        """SAParams accepts custom component_weights dict."""
        params = SAParams(
            formula="C6H14",
            component_weights={"wiener_index": 0.7, "logp": 0.3}
        )
        assert params.component_weights == {"wiener_index": 0.7, "logp": 0.3}

    def test_negative_weight_rejected(self):
        """SAParams rejects negative weights."""
        with pytest.raises(ValidationError) as exc_info:
            SAParams(
                formula="C6H14",
                component_weights={"wiener_index": -1.0}
            )
        assert "negative weight" in str(exc_info.value).lower()

    def test_all_zero_weights_rejected(self):
        """SAParams rejects all-zero weights."""
        with pytest.raises(ValidationError) as exc_info:
            SAParams(
                formula="C6H14",
                component_weights={"wiener_index": 0.0}
            )
        assert "non-zero" in str(exc_info.value).lower()

    def test_empty_weights_rejected(self):
        """SAParams rejects empty component_weights dict."""
        with pytest.raises(ValidationError) as exc_info:
            SAParams(
                formula="C6H14",
                component_weights={}
            )
        assert "at least one" in str(exc_info.value).lower()


class TestSAEngineWeightedEnergy:
    """Test SAEngine energy computation with weighted components."""

    def test_wiener_only_weight_1(self):
        """Default weights produce initial_energy == 35 (linear hexane Wiener)."""
        params = SAParams(
            formula="C6H14",
            component_weights={"wiener_index": 1.0},
            seed=42
        )
        engine = SAEngine(params)
        engine.init()
        assert engine._current_energy == 35

    def test_wiener_weight_2_doubles_energy(self):
        """Wiener weight of 2.0 doubles the energy."""
        params = SAParams(
            formula="C6H14",
            component_weights={"wiener_index": 2.0},
            seed=42
        )
        engine = SAEngine(params)
        engine.init()
        assert engine._current_energy == 70  # 35 * 2.0

    def test_two_components_combined(self):
        """Two components with weight 1.0 each sum correctly."""
        # First compute expected LogP for linear hexane
        # Formula C6H14 starts as linear hexane
        smiles = "CCCCCC"
        mol = Chem.MolFromSmiles(smiles)
        expected_logp = Descriptors.MolLogP(mol)

        params = SAParams(
            formula="C6H14",
            component_weights={"wiener_index": 1.0, "logp": 1.0},
            seed=42
        )
        engine = SAEngine(params)
        engine.init()

        # Energy should be 35 (Wiener) + LogP
        expected_energy = 35.0 + expected_logp
        assert engine._current_energy == pytest.approx(expected_energy, abs=0.01)

    def test_zero_weight_component_excluded(self):
        """Components with zero weight don't contribute to energy."""
        params = SAParams(
            formula="C6H14",
            component_weights={"wiener_index": 1.0, "logp": 0.0},
            seed=42
        )
        engine = SAEngine(params)
        engine.init()
        # Only Wiener contributes
        assert engine._current_energy == 35

    def test_default_weights_match_original(self):
        """Default params produce same results as hardcoded engine."""
        # This test verifies backward compatibility by checking that
        # the initial energy with default weights matches the expected
        # Wiener Index for linear hexane
        params = SAParams(
            formula="C6H14",
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=100,
            num_cycles=2,
            optimization_mode="MINIMIZE",
            seed=42
        )
        engine = SAEngine(params)
        result = engine.run()

        # Linear hexane has Wiener Index of 35
        assert result.initial_energy == 35
        # Verify accounting sums to total steps
        total_steps = 100 * 2
        assert (result.accepted_moves + result.rejected_moves +
                result.invalid_moves) == total_steps


class TestSAEngineValidation:
    """Test SAEngine validation of component_weights."""

    def test_unknown_component_raises(self):
        """Unknown component name raises ValueError."""
        params = SAParams(
            formula="C6H14",
            component_weights={"nonexistent": 1.0}
        )
        with pytest.raises(ValueError) as exc_info:
            SAEngine(params)

        error_msg = str(exc_info.value).lower()
        assert "unknown" in error_msg
        assert "component" in error_msg

    def test_unknown_component_lists_available(self):
        """ValueError for unknown component lists available components."""
        params = SAParams(
            formula="C6H14",
            component_weights={"nonexistent": 1.0}
        )
        with pytest.raises(ValueError) as exc_info:
            SAEngine(params)

        error_msg = str(exc_info.value).lower()
        # Should mention both available components
        assert "wiener_index" in error_msg
        assert "logp" in error_msg


class TestSAEngineBackwardCompatibility:
    """Test that existing code patterns still work."""

    def test_existing_tests_pattern_still_works(self):
        """SAParams without component_weights works (backward compat)."""
        params = SAParams(
            formula="C6H14",
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=100,
            num_cycles=2,
            optimization_mode="MINIMIZE",
            seed=42
        )
        engine = SAEngine(params)
        result = engine.run()

        # Should succeed and produce valid result
        assert result.best_smiles is not None
        total_steps = 100 * 2
        assert (result.accepted_moves + result.rejected_moves +
                result.invalid_moves) == total_steps

    def test_run_completes_with_multi_component(self):
        """Full SA run with multiple components completes successfully."""
        params = SAParams(
            formula="C6H14",
            component_weights={"wiener_index": 0.5, "logp": 0.5},
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=100,
            num_cycles=2,
            optimization_mode="MINIMIZE",
            seed=42
        )
        engine = SAEngine(params)
        result = engine.run()

        # Verify result is valid
        assert result.best_smiles is not None
        assert Chem.MolFromSmiles(result.best_smiles) is not None

        # Verify accounting
        total_steps = 100 * 2
        assert (result.accepted_moves + result.rejected_moves +
                result.invalid_moves) == total_steps

        # Verify acceptance ratio is in valid range
        acceptance_ratio = result.accepted_moves / total_steps
        assert 0 <= acceptance_ratio <= 1
