"""Shared dependencies for API endpoints."""

from app.services.session_manager import SessionManager
from app.config import Settings

settings = Settings()
session_manager = SessionManager(ttl_seconds=settings.session_ttl_seconds)
