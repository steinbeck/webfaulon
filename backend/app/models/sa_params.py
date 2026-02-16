"""
Pydantic models for Simulated Annealing parameters and results.

These models define the API contract for SA engine configuration and output.
"""

from pydantic import BaseModel, Field
from typing import Literal


class SAParams(BaseModel):
    """Simulated annealing configuration parameters."""

    formula: str = Field(
        ..., pattern=r"^([A-Z][a-z]?\d*)+$", examples=["C6H14"]
    )
    initial_temp: float = Field(
        100.0, gt=0, description="Initial temperature kT_0"
    )
    cooling_schedule_k: float = Field(
        8.0, ge=0, description="Cooling rate parameter k"
    )
    steps_per_cycle: int = Field(
        500, gt=0, description="Steps per temperature cycle"
    )
    num_cycles: int = Field(4, gt=0, description="Number of cooling cycles")
    optimization_mode: Literal["MINIMIZE", "MAXIMIZE"] = "MINIMIZE"
    seed: int = Field(42, description="Random seed for reproducibility")


class SAStepResult(BaseModel):
    """Result of a single SA iteration."""

    step: int
    current_energy: float
    best_energy: float
    temperature: float
    accepted: bool


class SAEngineState(BaseModel):
    """Snapshot of SA engine state for progress reporting."""

    step: int
    total_steps: int
    cycle: int
    current_energy: float
    best_energy: float
    best_smiles: str  # SMILES not MolBlock (memory per research pitfall #3)
    temperature: float
    accepted_moves: int
    rejected_moves: int
    invalid_moves: int
    is_complete: bool


class SAResult(BaseModel):
    """Complete optimization result."""

    best_energy: float
    best_smiles: str
    final_energy: float
    final_smiles: str
    initial_energy: float
    total_steps: int
    accepted_moves: int
    rejected_moves: int
    invalid_moves: int
    acceptance_ratio: float
    history: list[SAStepResult]
