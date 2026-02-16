"""Pytest configuration and fixtures for WebFaulon backend tests."""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture
def client():
    """Create a TestClient for the FastAPI application."""
    return TestClient(app)
