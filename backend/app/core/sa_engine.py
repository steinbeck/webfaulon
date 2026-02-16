"""
Simulated Annealing Engine for molecular graph optimization.

Implements the SA algorithm from Faulon's 1996 paper (Table 2).
Optimizes Wiener Index by iteratively applying displacement operations
with Metropolis acceptance criterion.
"""

import math
from typing import Optional

from app.models.sa_params import SAParams, SAResult, SAStepResult, SAEngineState
from app.core.molecule import MoleculeGraph
from app.core.initial_structure import generate_initial_structure
from app.core.wiener import compute_wiener_index
from app.core.displacement import attempt_displacement
from app.core.random import SeededRandom
from app.core.cooling import compute_temperature


class SAEngine:
    """Simulated Annealing Engine for molecular graph optimization.

    Usage:
        engine = SAEngine(SAParams(
            formula='C6H14',
            initial_temp=100,
            cooling_schedule_k=8,
            steps_per_cycle=500,
            num_cycles=4,
            optimization_mode='MINIMIZE',
            seed=42
        ))
        result = engine.run()
    """

    def __init__(self, params: SAParams):
        """Initialize SAEngine with parameters.

        Args:
            params: Simulated annealing configuration parameters
        """
        self._params = params
        self._rng = SeededRandom(params.seed)

        # These are set in init()
        self._current_graph: Optional[MoleculeGraph] = None
        self._current_energy: float = 0
        self._best_graph: Optional[MoleculeGraph] = None
        self._best_energy: float = 0
        self._initial_energy: float = 0
        self._accepted_moves: int = 0
        self._rejected_moves: int = 0
        self._invalid_moves: int = 0
        self._history: list[SAStepResult] = []
        self._total_steps: int = 0
        self._global_step: int = 0
        self._current_temperature: float = 0
        self._initialized: bool = False
        self._completed: bool = False

    def init(self) -> None:
        """Initialize the SA engine for step-by-step execution.

        Sets up initial structure, computes initial energy, and prepares state.
        Must be called before step().
        """
        # Generate initial structure
        self._current_graph = generate_initial_structure(self._params.formula)
        self._current_energy = compute_wiener_index(self._current_graph)
        self._initial_energy = self._current_energy
        self._best_graph = self._current_graph.clone()
        self._best_energy = self._current_energy

        # Reset state
        self._total_steps = self._params.steps_per_cycle * self._params.num_cycles
        self._global_step = 0
        self._accepted_moves = 0
        self._rejected_moves = 0
        self._invalid_moves = 0
        self._history = []
        self._current_temperature = self._params.initial_temp
        self._initialized = True
        self._completed = False

    def step(self) -> None:
        """Execute a single SA iteration.

        Advances the algorithm by one step. Can be called repeatedly with
        arbitrary delays between calls (enabling pause/resume).

        Raises:
            RuntimeError: If init() hasn't been called or execution is complete
        """
        if not self._initialized:
            raise RuntimeError("SAEngine.step() called before init()")

        if self._completed:
            raise RuntimeError("SAEngine.step() called after completion")

        # Increment step counter
        self._global_step += 1

        # Compute current temperature
        self._current_temperature = compute_temperature(
            self._global_step - 1,  # 0-indexed for temperature calculation
            self._total_steps,
            self._params.initial_temp,
            self._params.cooling_schedule_k,
        )

        # Execute one SA iteration
        self._iterate(self._current_temperature, self._global_step)

        # Check if complete
        if self._global_step == self._total_steps:
            self._completed = True

    def get_state(self) -> SAEngineState:
        """Get current execution state.

        Returns a snapshot of the current SA execution state.

        Returns:
            Current state including step number, energy values, and completion status

        Raises:
            RuntimeError: If init() hasn't been called
        """
        if not self._initialized:
            raise RuntimeError("SAEngine.get_state() called before init()")

        cycle = (
            0
            if self._global_step == 0
            else (self._global_step - 1) // self._params.steps_per_cycle + 1
        )

        return SAEngineState(
            step=self._global_step,
            total_steps=self._total_steps,
            cycle=cycle,
            current_energy=self._current_energy,
            best_energy=self._best_energy,
            best_smiles=self._best_graph.to_smiles(),
            temperature=self._current_temperature,
            accepted_moves=self._accepted_moves,
            rejected_moves=self._rejected_moves,
            invalid_moves=self._invalid_moves,
            is_complete=self._completed,
        )

    def get_result(self) -> SAResult:
        """Get final optimization result.

        Returns the complete SAResult with best/final graphs and statistics.
        Can only be called after all steps are complete.

        Returns:
            Final optimization result

        Raises:
            RuntimeError: If execution is not complete
        """
        if not self._completed:
            raise RuntimeError("SAEngine.get_result() called before completion")

        return SAResult(
            best_energy=self._best_energy,
            best_smiles=self._best_graph.to_smiles(),
            final_energy=self._current_energy,
            final_smiles=self._current_graph.to_smiles(),
            initial_energy=self._initial_energy,
            total_steps=self._total_steps,
            accepted_moves=self._accepted_moves,
            rejected_moves=self._rejected_moves,
            invalid_moves=self._invalid_moves,
            acceptance_ratio=self._accepted_moves / self._total_steps,
            history=self._history,
        )

    def run(self) -> SAResult:
        """Run the full SA optimization to completion.

        Convenience method that delegates to init(), step(), and get_result().
        Maintains backward compatibility with existing code.

        Returns:
            Complete optimization result
        """
        self.init()
        while not self._completed:
            self.step()
        return self.get_result()

    def _iterate(self, temperature: float, step_number: int) -> None:
        """Execute a single SA iteration.

        Args:
            temperature: Current temperature (kT)
            step_number: Current step number (1-indexed for display)
        """
        # Attempt displacement
        proposed_graph = attempt_displacement(self._current_graph, self._rng)

        # If displacement returned None, it's invalid
        if proposed_graph is None:
            self._invalid_moves += 1
            self._record_step(step_number, temperature, False)
            return

        # Compute energy of proposed graph
        proposed_energy = compute_wiener_index(proposed_graph)

        # Compute energy delta based on optimization mode
        # For MINIMIZE: positive delta = worsening move
        # For MAXIMIZE: positive delta = worsening move
        if self._params.optimization_mode == "MINIMIZE":
            delta_e = proposed_energy - self._current_energy
        else:
            delta_e = self._current_energy - proposed_energy

        # Metropolis acceptance criterion
        accepted = self._metropolis_accept(delta_e, temperature)

        if accepted:
            self._accepted_moves += 1
            self._current_graph = proposed_graph
            self._current_energy = proposed_energy

            # Update best if this is better
            if self._is_better(self._current_energy, self._best_energy):
                self._best_graph = self._current_graph.clone()
                self._best_energy = self._current_energy
        else:
            self._rejected_moves += 1

        self._record_step(step_number, temperature, accepted)

    def _is_better(self, proposed_energy: float, current_best: float) -> bool:
        """Check if proposed energy is better than current best.

        Args:
            proposed_energy: Energy to test
            current_best: Current best energy

        Returns:
            True if proposed is better than current best
        """
        if self._params.optimization_mode == "MINIMIZE":
            return proposed_energy < current_best
        return proposed_energy > current_best

    def _metropolis_accept(self, delta_e: float, temperature: float) -> bool:
        """Metropolis acceptance criterion.

        Always accepts improving moves (delta_e <= 0).
        Probabilistically accepts worsening moves: P = exp(-delta_e / temperature).

        Args:
            delta_e: Energy change (positive = worsening for both modes)
            temperature: Current temperature (kT)

        Returns:
            True if move should be accepted
        """
        # Always accept improving moves
        if delta_e <= 0:
            return True

        # Probabilistically accept worsening moves
        probability = math.exp(-delta_e / temperature)
        return self._rng.next() < probability

    def _record_step(
        self, step_number: int, temperature: float, accepted: bool
    ) -> None:
        """Record step result in history array.

        Args:
            step_number: Step number (1-indexed)
            temperature: Current temperature
            accepted: Whether move was accepted
        """
        self._history.append(
            SAStepResult(
                step=step_number,
                current_energy=self._current_energy,
                best_energy=self._best_energy,
                temperature=temperature,
                accepted=accepted,
            )
        )
