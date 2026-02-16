"""Error response models for structured API error handling."""

from pydantic import BaseModel


class ErrorResponse(BaseModel):
    """Standard error response structure for all API errors."""

    error: str
    detail: str | None = None
    status_code: int
