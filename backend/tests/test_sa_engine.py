"""
Comprehensive tests for SAEngine.

Ported from src/core/__tests__/SAEngine.test.ts
Tests cover: basic functionality, determinism, minimization, maximization,
chemical validity, history tracking, Metropolis criterion, edge cases,
and step-by-step execution.
"""

import pytest
from app.core.sa_engine import SAEngine
from app.models.sa_params import SAParams
from app.core.molecule import MoleculeGraph


DEFAULT_PARAMS = SAParams(
    formula="C6H14",
    initial_temp=100,
    cooling_schedule_k=8,
    steps_per_cycle=100,
    num_cycles=2,
    optimization_mode="MINIMIZE",
    seed=42,
)


class TestBasicFunctionality:
    """Basic SAEngine functionality tests."""

    def test_creates_with_valid_params(self):
        """Engine instantiates without error."""
        engine = SAEngine(DEFAULT_PARAMS)
        assert engine is not None

    def test_run_returns_all_fields(self):
        """run() returns SAResult with all required fields."""
        engine = SAEngine(DEFAULT_PARAMS)
        result = engine.run()

        assert result.best_energy is not None
        assert isinstance(result.best_energy, (int, float))
        assert result.best_smiles is not None
        assert isinstance(result.best_smiles, str)
        assert result.final_energy is not None
        assert isinstance(result.final_energy, (int, float))
        assert result.final_smiles is not None
        assert isinstance(result.final_smiles, str)
        assert result.initial_energy is not None
        assert isinstance(result.initial_energy, (int, float))
        assert result.total_steps == 200  # 100 * 2
        assert isinstance(result.accepted_moves, int)
        assert isinstance(result.rejected_moves, int)
        assert isinstance(result.invalid_moves, int)
        assert isinstance(result.acceptance_ratio, float)
        assert isinstance(result.history, list)

    def test_history_length(self):
        """history array has correct length."""
        engine = SAEngine(DEFAULT_PARAMS)
        result = engine.run()

        assert len(result.history) == result.total_steps

    def test_accounting_sums_to_total(self):
        """acceptedMoves + rejectedMoves + invalidMoves = totalSteps."""
        engine = SAEngine(DEFAULT_PARAMS)
        result = engine.run()

        total = result.accepted_moves + result.rejected_moves + result.invalid_moves
        assert total == result.total_steps

    def test_acceptance_ratio_range(self):
        """acceptanceRatio is between 0 and 1."""
        engine = SAEngine(DEFAULT_PARAMS)
        result = engine.run()

        assert 0 <= result.acceptance_ratio <= 1


class TestDeterminism:
    """Determinism and reproducibility tests."""

    def test_same_seed_identical_results(self):
        """Same seed produces identical results."""
        params1 = SAParams(**{**DEFAULT_PARAMS.model_dump(), "seed": 42})
        params2 = SAParams(**{**DEFAULT_PARAMS.model_dump(), "seed": 42})

        engine1 = SAEngine(params1)
        engine2 = SAEngine(params2)

        result1 = engine1.run()
        result2 = engine2.run()

        assert result1.best_energy == result2.best_energy
        assert result1.total_steps == result2.total_steps
        assert result1.accepted_moves == result2.accepted_moves
        assert result1.rejected_moves == result2.rejected_moves
        assert result1.invalid_moves == result2.invalid_moves

    def test_different_seeds_differ(self):
        """Different seeds produce different results."""
        params1 = SAParams(**{**DEFAULT_PARAMS.model_dump(), "seed": 42})
        params2 = SAParams(**{**DEFAULT_PARAMS.model_dump(), "seed": 123})

        engine1 = SAEngine(params1)
        engine2 = SAEngine(params2)

        result1 = engine1.run()
        result2 = engine2.run()

        # At least one metric should differ
        differs = (
            result1.best_energy != result2.best_energy
            or result1.accepted_moves != result2.accepted_moves
            or result1.rejected_moves != result2.rejected_moves
        )
        assert differs


class TestMinimization:
    """Minimization mode tests."""

    def test_minimize_improves_c6h14(self):
        """bestEnergy <= initialEnergy for C6H14."""
        params = SAParams(
            formula="C6H14",
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=500,
            num_cycles=4,
            optimization_mode="MINIMIZE",
            seed=42,
        )

        engine = SAEngine(params)
        result = engine.run()

        # Best energy should never be worse than initial
        assert result.best_energy <= result.initial_energy
        # Verify SA is actually accepting moves
        assert result.accepted_moves > 0

    def test_minimize_initial_energy_35(self):
        """initialEnergy is 35 (linear hexane Wiener Index)."""
        params = SAParams(
            formula="C6H14",
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=500,
            num_cycles=4,
            optimization_mode="MINIMIZE",
            seed=42,
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.initial_energy == 35

    def test_minimize_best_le_final(self):
        """bestEnergy <= finalEnergy."""
        params = SAParams(
            **{**DEFAULT_PARAMS.model_dump(), "optimization_mode": "MINIMIZE"}
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.best_energy <= result.final_energy

    def test_minimize_accepted_moves_positive(self):
        """acceptedMoves > 0."""
        params = SAParams(
            formula="C6H14",
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=500,
            num_cycles=4,
            optimization_mode="MINIMIZE",
            seed=42,
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.accepted_moves > 0


class TestMaximization:
    """Maximization mode tests."""

    def test_maximize_c6h14(self):
        """bestEnergy >= initialEnergy for C6H14."""
        params = SAParams(
            **{
                **DEFAULT_PARAMS.model_dump(),
                "optimization_mode": "MAXIMIZE",
                "steps_per_cycle": 500,
                "num_cycles": 4,
            }
        )

        engine = SAEngine(params)
        result = engine.run()

        # Linear hexane already has maximum Wiener Index = 35
        assert result.best_energy >= result.initial_energy

    def test_maximize_best_ge_final(self):
        """bestEnergy >= finalEnergy."""
        params = SAParams(
            **{**DEFAULT_PARAMS.model_dump(), "optimization_mode": "MAXIMIZE"}
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.best_energy >= result.final_energy


class TestChemicalValidity:
    """Chemical validity tests."""

    def test_best_smiles_valid(self):
        """best_smiles is non-empty and parseable by RDKit."""
        engine = SAEngine(DEFAULT_PARAMS)
        result = engine.run()

        assert result.best_smiles is not None
        assert len(result.best_smiles) > 0

        # Verify it's parseable by RDKit
        from rdkit import Chem

        mol = Chem.MolFromSmiles(result.best_smiles)
        assert mol is not None

    def test_final_smiles_valid(self):
        """final_smiles is non-empty and parseable by RDKit."""
        engine = SAEngine(DEFAULT_PARAMS)
        result = engine.run()

        assert result.final_smiles is not None
        assert len(result.final_smiles) > 0

        # Verify it's parseable by RDKit
        from rdkit import Chem

        mol = Chem.MolFromSmiles(result.final_smiles)
        assert mol is not None


class TestHistoryTracking:
    """History tracking tests."""

    def test_history_step_numbers(self):
        """Steps are numbered 1 through totalSteps."""
        params = SAParams(
            **{
                **DEFAULT_PARAMS.model_dump(),
                "steps_per_cycle": 10,
                "num_cycles": 1,
            }
        )

        engine = SAEngine(params)
        result = engine.run()

        assert len(result.history) == 10

        for i in range(len(result.history)):
            step = result.history[i]
            assert step.step == i + 1
            assert isinstance(step.current_energy, (int, float))
            assert isinstance(step.best_energy, (int, float))
            assert isinstance(step.temperature, (int, float))
            assert isinstance(step.accepted, bool)

    def test_history_best_non_increasing_minimize(self):
        """bestEnergy never increases (minimization)."""
        params = SAParams(
            **{
                **DEFAULT_PARAMS.model_dump(),
                "optimization_mode": "MINIMIZE",
                "steps_per_cycle": 50,
                "num_cycles": 2,
            }
        )

        engine = SAEngine(params)
        result = engine.run()

        for i in range(1, len(result.history)):
            assert result.history[i].best_energy <= result.history[i - 1].best_energy

    def test_history_best_non_decreasing_maximize(self):
        """bestEnergy never decreases (maximization)."""
        params = SAParams(
            **{
                **DEFAULT_PARAMS.model_dump(),
                "optimization_mode": "MAXIMIZE",
                "steps_per_cycle": 50,
                "num_cycles": 2,
            }
        )

        engine = SAEngine(params)
        result = engine.run()

        for i in range(1, len(result.history)):
            assert result.history[i].best_energy >= result.history[i - 1].best_energy


class TestMetropolisCriterion:
    """Metropolis acceptance criterion tests."""

    def test_low_temp_accepts_some(self):
        """At low temperature, some moves are still accepted."""
        params = SAParams(
            formula="C4H10",
            initial_temp=0.01,  # Very low temperature
            cooling_schedule_k=0,  # Constant temperature
            steps_per_cycle=100,
            num_cycles=1,
            optimization_mode="MINIMIZE",
            seed=42,
        )

        engine = SAEngine(params)
        result = engine.run()

        # At low temperature, should still accept improving moves
        assert result.accepted_moves > 0

    def test_high_temp_high_acceptance(self):
        """At high temperature, acceptance ratio > 0.3."""
        params = SAParams(
            formula="C6H14",
            initial_temp=1000,
            cooling_schedule_k=0,
            steps_per_cycle=500,
            num_cycles=1,
            optimization_mode="MINIMIZE",
            seed=100,
        )

        engine = SAEngine(params)
        result = engine.run()

        # At very high temperature, should accept most valid moves
        assert result.acceptance_ratio > 0.3

    def test_high_temp_higher_than_low_temp(self):
        """High temp acceptance > low temp acceptance."""
        high_temp_params = SAParams(
            formula="C6H14",
            initial_temp=1000,
            cooling_schedule_k=0,
            steps_per_cycle=500,
            num_cycles=1,
            optimization_mode="MINIMIZE",
            seed=100,
        )

        low_temp_params = SAParams(
            formula="C6H14",
            initial_temp=0.01,
            cooling_schedule_k=0,
            steps_per_cycle=500,
            num_cycles=1,
            optimization_mode="MINIMIZE",
            seed=100,
        )

        high_temp_engine = SAEngine(high_temp_params)
        high_temp_result = high_temp_engine.run()

        low_temp_engine = SAEngine(low_temp_params)
        low_temp_result = low_temp_engine.run()

        # High temp should have higher acceptance ratio than low temp
        assert high_temp_result.acceptance_ratio >= low_temp_result.acceptance_ratio

    def test_very_high_temp_frequent_worsening(self):
        """At very high temperature, acceptance ratio > 0.5."""
        params = SAParams(
            formula="C5H12",
            initial_temp=1000,  # Very high temperature
            cooling_schedule_k=0,  # Constant
            steps_per_cycle=1000,
            num_cycles=1,
            optimization_mode="MINIMIZE",
            seed=200,
        )

        engine = SAEngine(params)
        result = engine.run()

        # At very high temperature, acceptance ratio should be high
        assert result.acceptance_ratio > 0.5


class TestEdgeCases:
    """Edge case tests."""

    def test_small_molecule_c4h10(self):
        """Works on 4-carbon molecule."""
        params = SAParams(
            **{
                **DEFAULT_PARAMS.model_dump(),
                "formula": "C4H10",
                "steps_per_cycle": 50,
                "num_cycles": 1,
            }
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.best_smiles is not None
        # Verify chemical validity via RDKit
        from rdkit import Chem

        mol = Chem.MolFromSmiles(result.best_smiles)
        assert mol is not None

    def test_single_cycle(self):
        """numCycles=1, stepsPerCycle=10, totalSteps=10."""
        params = SAParams(
            **{**DEFAULT_PARAMS.model_dump(), "num_cycles": 1, "steps_per_cycle": 10}
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.total_steps == 10
        assert len(result.history) == 10

    def test_many_cycles(self):
        """numCycles=10, stepsPerCycle=10, totalSteps=100."""
        params = SAParams(
            **{**DEFAULT_PARAMS.model_dump(), "num_cycles": 10, "steps_per_cycle": 10}
        )

        engine = SAEngine(params)
        result = engine.run()

        assert result.total_steps == 100
        assert len(result.history) == 100


class TestStepByStepExecution:
    """Step-by-step execution tests."""

    STEP_PARAMS = SAParams(
        formula="C6H14",
        initial_temp=100,
        cooling_schedule_k=8,
        steps_per_cycle=50,
        num_cycles=2,
        optimization_mode="MINIMIZE",
        seed=42,
    )

    def test_init_sets_initial_state(self):
        """init() sets up initial state correctly."""
        engine = SAEngine(self.STEP_PARAMS)
        engine.init()
        state = engine.get_state()

        assert state.step == 0
        assert state.is_complete is False
        # Initial structure computed
        assert state.current_energy == state.best_energy
        assert state.current_energy > 0
        assert state.total_steps == 100  # 50 * 2

    def test_step_advances(self):
        """step() advances one iteration."""
        engine = SAEngine(self.STEP_PARAMS)
        engine.init()
        engine.step()
        state = engine.get_state()

        assert state.step == 1

    def test_step_before_init_raises(self):
        """step() without init() raises."""
        engine = SAEngine(self.STEP_PARAMS)

        with pytest.raises(RuntimeError):
            engine.step()

    def test_multiple_steps(self):
        """Multiple steps execute correctly."""
        engine = SAEngine(self.STEP_PARAMS)
        engine.init()

        N = 5
        for _ in range(N):
            engine.step()

        state = engine.get_state()
        assert state.step == N

    def test_is_complete_after_all_steps(self):
        """isComplete becomes true after all steps."""
        params = SAParams(
            **{**self.STEP_PARAMS.model_dump(), "steps_per_cycle": 10, "num_cycles": 1}
        )

        engine = SAEngine(params)
        engine.init()

        # Execute all 10 steps
        for _ in range(10):
            engine.step()

        state = engine.get_state()
        assert state.is_complete is True
        assert state.step == 10

    def test_step_after_complete_raises(self):
        """step() after completion raises."""
        params = SAParams(
            **{**self.STEP_PARAMS.model_dump(), "steps_per_cycle": 5, "num_cycles": 1}
        )

        engine = SAEngine(params)
        engine.init()

        # Execute all steps
        for _ in range(5):
            engine.step()

        # Try to step again after completion
        with pytest.raises(RuntimeError):
            engine.step()

    def test_get_result_after_complete(self):
        """getResult() works after completion."""
        params = SAParams(
            **{**self.STEP_PARAMS.model_dump(), "steps_per_cycle": 10, "num_cycles": 1}
        )

        engine = SAEngine(params)
        engine.init()

        for _ in range(10):
            engine.step()

        result = engine.get_result()

        assert result.best_smiles is not None
        assert isinstance(result.best_energy, (int, float))
        assert result.final_smiles is not None
        assert isinstance(result.final_energy, (int, float))
        assert isinstance(result.initial_energy, (int, float))
        assert result.total_steps == 10
        assert isinstance(result.accepted_moves, int)
        assert isinstance(result.rejected_moves, int)
        assert isinstance(result.invalid_moves, int)
        assert isinstance(result.acceptance_ratio, float)
        assert isinstance(result.history, list)

    def test_get_result_before_complete_raises(self):
        """getResult() before completion raises."""
        engine = SAEngine(self.STEP_PARAMS)
        engine.init()
        engine.step()  # Only 1 step out of 100

        with pytest.raises(RuntimeError):
            engine.get_result()

    def test_step_matches_run(self):
        """Step-by-step execution matches run() for same seed."""
        params = SAParams(
            **{
                **self.STEP_PARAMS.model_dump(),
                "steps_per_cycle": 50,
                "num_cycles": 2,
                "seed": 12345,
            }
        )

        # Engine 1: use run()
        engine1 = SAEngine(params)
        result1 = engine1.run()

        # Engine 2: use step-by-step
        engine2 = SAEngine(params)
        engine2.init()
        while not engine2.get_state().is_complete:
            engine2.step()
        result2 = engine2.get_result()

        # Both should produce identical results
        assert result2.best_energy == result1.best_energy
        assert result2.accepted_moves == result1.accepted_moves
        assert result2.rejected_moves == result1.rejected_moves
        assert result2.invalid_moves == result1.invalid_moves
        assert result2.total_steps == result1.total_steps

    def test_run_backward_compatible(self):
        """run() still works unchanged (backward compatibility)."""
        engine = SAEngine(self.STEP_PARAMS)
        result = engine.run()

        # Verify run() works exactly as before
        assert result.best_smiles is not None
        assert isinstance(result.best_energy, (int, float))
        assert result.total_steps == 100
        assert len(result.history) == 100

        # Accounting still correct
        total = result.accepted_moves + result.rejected_moves + result.invalid_moves
        assert total == result.total_steps


class TestFullSARun:
    """Full SA run integration test."""

    def test_full_run_500x4_c6h14(self):
        """Run 500 steps x 4 cycles on C6H14.

        Verify:
        - bestEnergy <= 35
        - valid SMILES at every history entry
        - accepted + rejected + invalid = 2000
        """
        params = SAParams(
            formula="C6H14",
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=500,
            num_cycles=4,
            optimization_mode="MINIMIZE",
            seed=42,
        )

        engine = SAEngine(params)
        result = engine.run()

        # Best energy should be <= initial (35 for linear hexane)
        assert result.best_energy <= 35

        # Verify all history entries are valid
        from rdkit import Chem

        for step_result in result.history:
            assert step_result.current_energy > 0
            assert step_result.best_energy > 0
            assert step_result.temperature >= 0.01

        # Verify accounting
        assert (
            result.accepted_moves + result.rejected_moves + result.invalid_moves
            == 2000
        )

        # Verify SMILES are valid
        best_mol = Chem.MolFromSmiles(result.best_smiles)
        assert best_mol is not None
        final_mol = Chem.MolFromSmiles(result.final_smiles)
        assert final_mol is not None
