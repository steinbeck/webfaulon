"""SA control endpoints - start, pause, reset sessions."""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel

from app.dependencies import session_manager

router = APIRouter(prefix="/api/sa", tags=["SA Control"])


class ControlResponse(BaseModel):
    """Response model for control endpoints."""

    session_id: str
    status: str


@router.post("/{session_id}/start", response_model=ControlResponse)
async def start_session(session_id: str) -> ControlResponse:
    """Start SA execution for a session.

    Args:
        session_id: Session ID to start

    Returns:
        Session ID and status

    Raises:
        HTTPException: 404 if session not found, 400 if invalid state
    """
    session = session_manager.get(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")

    if session.state == "running":
        raise HTTPException(status_code=400, detail="Session already running")

    if session.state == "complete":
        raise HTTPException(status_code=400, detail="Session already complete")

    session.state = "running"
    return ControlResponse(session_id=session_id, status="running")


@router.post("/{session_id}/pause", response_model=ControlResponse)
async def pause_session(session_id: str) -> ControlResponse:
    """Pause SA execution for a session.

    Args:
        session_id: Session ID to pause

    Returns:
        Session ID and status

    Raises:
        HTTPException: 404 if session not found, 400 if not running
    """
    session = session_manager.get(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")

    if session.state != "running":
        raise HTTPException(status_code=400, detail="Session not running")

    session.state = "paused"
    return ControlResponse(session_id=session_id, status="paused")


@router.post("/{session_id}/reset", response_model=ControlResponse)
async def reset_session(session_id: str) -> ControlResponse:
    """Reset SA session to idle state.

    Args:
        session_id: Session ID to reset

    Returns:
        Session ID and status

    Raises:
        HTTPException: 404 if session not found
    """
    session = session_manager.get(session_id)
    if session is None:
        raise HTTPException(status_code=404, detail="Session not found or expired")

    # Re-initialize the engine
    session.engine.init()
    session.state = "idle"
    return ControlResponse(session_id=session_id, status="idle")
