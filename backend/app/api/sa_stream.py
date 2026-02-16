"""SA streaming endpoint - SSE progress updates."""

from sse_starlette.sse import EventSourceResponse
from fastapi import APIRouter, Request, HTTPException
import asyncio
import json

from app.dependencies import session_manager
from app.services.svg_renderer import generate_molecule_svg

router = APIRouter(prefix="/api/sa", tags=["SA Streaming"])


@router.get("/{session_id}/stream")
async def stream_progress(session_id: str, request: Request):
    """Stream SA optimization progress via SSE.

    The SA engine runs steps inside the async generator. Each step yields
    a progress event. The generator respects pause/running state and
    checks for client disconnection.

    Events emitted:
    - "progress": step-by-step state updates (id = step number)
    - "complete": final result with best molecule SVG
    - "error": optimization error
    - "waiting": heartbeat when paused or idle (keeps connection alive)

    Headers:
    - X-Accel-Buffering: no (disables nginx buffering for real-time delivery)
    """
    session = session_manager.get(session_id)
    if not session:
        raise HTTPException(404, "Session not found or expired")

    async def event_generator():
        try:
            idle_count = 0
            while True:
                # Check client disconnect
                if await request.is_disconnected():
                    break

                # Session might have been deleted
                current_session = session_manager.get(session_id)
                if current_session is None:
                    break

                # Only step if running
                if current_session.state == "running":
                    idle_count = 0

                    # Check if engine is complete
                    if current_session.engine._completed:
                        # Send completion event
                        state = current_session.engine.get_state()
                        svg = generate_molecule_svg(state.best_smiles)
                        current_session.state = "complete"

                        yield {
                            "event": "complete",
                            "data": json.dumps({
                                "step": state.step,
                                "total_steps": state.total_steps,
                                "best_energy": state.best_energy,
                                "best_smiles": state.best_smiles,
                                "best_svg": svg,
                                "accepted_moves": state.accepted_moves,
                                "rejected_moves": state.rejected_moves,
                                "invalid_moves": state.invalid_moves,
                                "disconnected_moves": state.disconnected_moves,
                            })
                        }
                        break

                    # Execute one SA step
                    current_session.engine.step()
                    state = current_session.engine.get_state()

                    # Generate SVG for best molecule
                    svg = generate_molecule_svg(state.best_smiles)

                    yield {
                        "id": str(state.step),
                        "event": "progress",
                        "data": json.dumps({
                            "step": state.step,
                            "total_steps": state.total_steps,
                            "cycle": state.cycle,
                            "temperature": state.temperature,
                            "current_energy": state.current_energy,
                            "best_energy": state.best_energy,
                            "best_smiles": state.best_smiles,
                            "best_svg": svg,
                            "accepted_moves": state.accepted_moves,
                            "rejected_moves": state.rejected_moves,
                            "invalid_moves": state.invalid_moves,
                                "disconnected_moves": state.disconnected_moves,
                            "is_complete": state.is_complete,
                        })
                    }

                    # Check if just completed
                    if state.is_complete:
                        current_session.state = "complete"
                        yield {
                            "event": "complete",
                            "data": json.dumps({
                                "step": state.step,
                                "total_steps": state.total_steps,
                                "best_energy": state.best_energy,
                                "best_smiles": state.best_smiles,
                                "best_svg": svg,
                                "accepted_moves": state.accepted_moves,
                                "rejected_moves": state.rejected_moves,
                                "invalid_moves": state.invalid_moves,
                                "disconnected_moves": state.disconnected_moves,
                            })
                        }
                        break

                    # Yield control to event loop (~100 events/sec throttling)
                    await asyncio.sleep(0.01)

                elif current_session.state == "complete":
                    # Already complete, send final state and close
                    state = current_session.engine.get_state()
                    svg = generate_molecule_svg(state.best_smiles)
                    yield {
                        "event": "complete",
                        "data": json.dumps({
                            "step": state.step,
                            "total_steps": state.total_steps,
                            "best_energy": state.best_energy,
                            "best_smiles": state.best_smiles,
                            "best_svg": svg,
                        })
                    }
                    break

                else:
                    # Paused or idle -- send heartbeat to keep connection alive
                    idle_count += 1
                    if idle_count % 10 == 0:  # Every ~1 second
                        yield {
                            "event": "waiting",
                            "data": json.dumps({
                                "session_state": current_session.state,
                                "message": f"Session is {current_session.state}"
                            })
                        }
                    await asyncio.sleep(0.1)  # Check state 10x/sec when idle

        except Exception as e:
            yield {
                "event": "error",
                "data": json.dumps({"error": str(e)})
            }

    return EventSourceResponse(
        event_generator(),
        headers={"X-Accel-Buffering": "no"}
    )
