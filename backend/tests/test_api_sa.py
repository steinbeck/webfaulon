"""Tests for SA API endpoints."""

import pytest
from fastapi.testclient import TestClient

from app.main import app
from app.dependencies import session_manager


@pytest.fixture(autouse=True)
def clear_sessions():
    """Clear all sessions between tests."""
    session_manager._sessions.clear()
    yield
    session_manager._sessions.clear()


@pytest.fixture
def client():
    """Test client fixture."""
    return TestClient(app)


# ===== Configure endpoint tests =====


def test_configure_returns_session_id(client):
    """POST /api/sa/configure with valid params returns session_id and status."""
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})

    assert response.status_code == 200
    data = response.json()
    assert "session_id" in data
    assert "status" in data
    assert data["status"] == "configured"
    # Verify UUID format (36 chars with hyphens)
    assert len(data["session_id"]) == 36
    assert data["session_id"].count("-") == 4


def test_configure_with_custom_params(client):
    """POST with all params returns session_id."""
    response = client.post(
        "/api/sa/configure",
        json={
            "formula": "C6H14",
            "initial_temp": 150.0,
            "cooling_schedule_k": 10.0,
            "steps_per_cycle": 1000,
            "num_cycles": 5,
            "optimization_mode": "MAXIMIZE",
            "seed": 123,
        },
    )

    assert response.status_code == 200
    data = response.json()
    assert "session_id" in data
    assert data["status"] == "configured"


def test_configure_invalid_formula(client):
    """POST with empty formula returns 422 validation error."""
    response = client.post("/api/sa/configure", json={"formula": ""})

    assert response.status_code == 422


def test_configure_invalid_params(client):
    """POST with invalid params returns 422."""
    response = client.post(
        "/api/sa/configure", json={"formula": "C6H14", "initial_temp": -1}
    )

    assert response.status_code == 422


# ===== Control endpoint tests =====


def test_start_session(client):
    """POST /api/sa/{id}/start returns 200 with status running."""
    # Create session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]

    # Start session
    response = client.post(f"/api/sa/{session_id}/start")

    assert response.status_code == 200
    data = response.json()
    assert data["session_id"] == session_id
    assert data["status"] == "running"


def test_start_already_running(client):
    """Starting a running session returns 400."""
    # Create and start session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]
    client.post(f"/api/sa/{session_id}/start")

    # Try to start again
    response = client.post(f"/api/sa/{session_id}/start")

    assert response.status_code == 400


def test_start_nonexistent_session(client):
    """Starting nonexistent session returns 404."""
    response = client.post("/api/sa/fake-id-12345/start")

    assert response.status_code == 404


def test_pause_running_session(client):
    """POST /api/sa/{id}/pause returns 200 with status paused."""
    # Create and start session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]
    client.post(f"/api/sa/{session_id}/start")

    # Pause session
    response = client.post(f"/api/sa/{session_id}/pause")

    assert response.status_code == 200
    data = response.json()
    assert data["session_id"] == session_id
    assert data["status"] == "paused"


def test_pause_not_running(client):
    """Pausing an idle session returns 400."""
    # Create session (idle state)
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]

    # Try to pause idle session
    response = client.post(f"/api/sa/{session_id}/pause")

    assert response.status_code == 400


def test_reset_session(client):
    """POST /api/sa/{id}/reset returns 200 with status idle."""
    # Create and start session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]
    client.post(f"/api/sa/{session_id}/start")

    # Reset session
    response = client.post(f"/api/sa/{session_id}/reset")

    assert response.status_code == 200
    data = response.json()
    assert data["session_id"] == session_id
    assert data["status"] == "idle"


def test_reset_nonexistent(client):
    """Resetting nonexistent session returns 404."""
    response = client.post("/api/sa/fake-id-12345/reset")

    assert response.status_code == 404


# ===== Status endpoint tests =====


def test_status_returns_state(client):
    """GET /api/sa/{id}/status returns current state with all fields."""
    # Create session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]

    # Get status
    response = client.get(f"/api/sa/{session_id}/status")

    assert response.status_code == 200
    data = response.json()
    # Verify all required fields
    assert data["session_id"] == session_id
    assert data["session_state"] == "idle"
    assert "step" in data
    assert "total_steps" in data
    assert "cycle" in data
    assert "current_energy" in data
    assert "best_energy" in data
    assert "best_smiles" in data
    assert "best_svg" in data
    assert "temperature" in data
    assert "accepted_moves" in data
    assert "rejected_moves" in data
    assert "invalid_moves" in data
    assert "is_complete" in data


def test_status_includes_svg(client):
    """Status response contains best_svg field with SVG content."""
    # Create session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]

    # Get status
    response = client.get(f"/api/sa/{session_id}/status")

    assert response.status_code == 200
    data = response.json()
    assert "best_svg" in data
    assert "<svg" in data["best_svg"]


def test_status_nonexistent(client):
    """Getting status for nonexistent session returns 404."""
    response = client.get("/api/sa/fake-id-12345/status")

    assert response.status_code == 404


def test_status_after_start(client):
    """Status shows session_state as running after start."""
    # Create and start session
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    session_id = response.json()["session_id"]
    client.post(f"/api/sa/{session_id}/start")

    # Get status
    response = client.get(f"/api/sa/{session_id}/status")

    assert response.status_code == 200
    data = response.json()
    assert data["session_state"] == "running"


# ===== Integration tests =====


def test_health_check_still_works(client):
    """GET /health returns 200 (regression check)."""
    response = client.get("/health")

    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "healthy"
    assert data["version"] == "2.0.0"


def test_full_configure_start_pause_reset_flow(client):
    """Configure -> start -> pause -> reset flow works sequentially."""
    # Configure
    response = client.post("/api/sa/configure", json={"formula": "C6H14"})
    assert response.status_code == 200
    session_id = response.json()["session_id"]
    assert response.json()["status"] == "configured"

    # Status shows idle
    response = client.get(f"/api/sa/{session_id}/status")
    assert response.status_code == 200
    assert response.json()["session_state"] == "idle"

    # Start
    response = client.post(f"/api/sa/{session_id}/start")
    assert response.status_code == 200
    assert response.json()["status"] == "running"

    # Status shows running
    response = client.get(f"/api/sa/{session_id}/status")
    assert response.status_code == 200
    assert response.json()["session_state"] == "running"

    # Pause
    response = client.post(f"/api/sa/{session_id}/pause")
    assert response.status_code == 200
    assert response.json()["status"] == "paused"

    # Status shows paused
    response = client.get(f"/api/sa/{session_id}/status")
    assert response.status_code == 200
    assert response.json()["session_state"] == "paused"

    # Reset
    response = client.post(f"/api/sa/{session_id}/reset")
    assert response.status_code == 200
    assert response.json()["status"] == "idle"

    # Status shows idle again
    response = client.get(f"/api/sa/{session_id}/status")
    assert response.status_code == 200
    assert response.json()["session_state"] == "idle"
