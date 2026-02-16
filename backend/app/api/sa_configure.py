"""SA configuration endpoint - creates new SA sessions."""

from fastapi import APIRouter
from pydantic import BaseModel

from app.models.sa_params import SAParams
from app.dependencies import session_manager

router = APIRouter(prefix="/api/sa", tags=["SA Configuration"])


class ConfigureResponse(BaseModel):
    """Response model for configure endpoint."""

    session_id: str
    status: str


@router.post("/configure", response_model=ConfigureResponse)
async def configure_session(params: SAParams) -> ConfigureResponse:
    """Create a new SA session with provided parameters.

    Args:
        params: SA configuration parameters

    Returns:
        Session ID and status
    """
    session_id = session_manager.create(params)
    return ConfigureResponse(session_id=session_id, status="configured")
