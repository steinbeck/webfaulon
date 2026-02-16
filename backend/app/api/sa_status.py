"""SA status endpoint - query session state."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from app.dependencies import session_manager
from app.services.svg_renderer import generate_molecule_svg

router = APIRouter(prefix="/api/sa", tags=["SA Status"])


class StatusResponse(BaseModel):
    """Response model for status endpoint."""

    session_id: str
    session_state: str
    step: int
    total_steps: int
    cycle: int
    current_energy: float
    best_energy: float
    best_smiles: str
    best_svg: str
    temperature: float
    accepted_moves: int
    rejected_moves: int
    invalid_moves: int
    is_complete: bool


@router.get("/{session_id}/status", response_model=StatusResponse)
async def get_session_status(session_id: str) -> StatusResponse:
    """Get current state of SA session.

    Args:
        session_id: Session ID to query

    Returns:
        Complete session state snapshot with SVG

    Raises:
        HTTPException: 404 if session not found
    """
    session = session_manager.get(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")

    # Get engine state
    state = session.engine.get_state()

    # Generate SVG for best molecule
    svg = generate_molecule_svg(state.best_smiles)

    return StatusResponse(
        session_id=session_id,
        session_state=session.state,
        step=state.step,
        total_steps=state.total_steps,
        cycle=state.cycle,
        current_energy=state.current_energy,
        best_energy=state.best_energy,
        best_smiles=state.best_smiles,
        best_svg=svg,
        temperature=state.temperature,
        accepted_moves=state.accepted_moves,
        rejected_moves=state.rejected_moves,
        invalid_moves=state.invalid_moves,
        is_complete=state.is_complete,
    )
