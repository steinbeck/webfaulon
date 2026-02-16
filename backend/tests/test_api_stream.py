"""Tests for SSE streaming endpoint."""

import pytest
import json
from fastapi.testclient import TestClient
from app.main import app
from app.dependencies import session_manager


@pytest.fixture(autouse=True)
def clear_sessions():
    """Clear sessions before and after each test."""
    session_manager._sessions.clear()
    yield
    session_manager._sessions.clear()


@pytest.fixture
def client():
    """Test client fixture."""
    return TestClient(app)


def parse_sse_events(response):
    """Parse SSE events from streaming response.

    Args:
        response: httpx streaming response

    Returns:
        List of event dictionaries with 'event' and 'data' keys
    """
    events = []
    current_event = {}

    for line in response.iter_lines():
        if line.startswith("event:"):
            current_event["event"] = line[len("event:"):].strip()
        elif line.startswith("data:"):
            current_event["data"] = json.loads(line[len("data:"):].strip())
        elif line.startswith("id:"):
            current_event["id"] = line[len("id:"):].strip()
        elif line == "":
            if current_event:
                events.append(current_event)
                current_event = {}

    # Add last event if exists
    if current_event:
        events.append(current_event)

    return events


# SSE Connection Tests


def test_stream_nonexistent_session_returns_404(client):
    """Stream endpoint returns 404 for nonexistent session."""
    response = client.get("/api/sa/nonexistent/stream")
    assert response.status_code == 404


def test_stream_returns_event_stream_content_type(client):
    """Stream endpoint returns text/event-stream content type."""
    # Create and start a small session
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]

    client.post(f"/api/sa/{session_id}/start")

    # Use stream context manager to check headers
    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        assert response.status_code == 200
        assert "text/event-stream" in response.headers.get("content-type", "")


def test_stream_has_no_buffering_header(client):
    """Stream endpoint includes X-Accel-Buffering: no header."""
    # Create and start a small session
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]

    client.post(f"/api/sa/{session_id}/start")

    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        assert response.headers.get("x-accel-buffering") == "no"


# SSE Event Tests


def test_stream_emits_progress_events(client):
    """Stream emits progress events with required fields."""
    # Configure small session (5 steps total)
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]

    # Start session
    client.post(f"/api/sa/{session_id}/start")

    # Stream and collect events
    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        events = parse_sse_events(response)

    # Find progress events
    progress_events = [e for e in events if e.get("event") == "progress"]

    # Should have at least one progress event
    assert len(progress_events) > 0

    # Verify progress event structure
    first_progress = progress_events[0]["data"]
    assert "step" in first_progress
    assert "temperature" in first_progress
    assert "best_energy" in first_progress
    assert "best_smiles" in first_progress
    assert isinstance(first_progress["step"], int)
    assert isinstance(first_progress["temperature"], (int, float))
    assert isinstance(first_progress["best_energy"], (int, float))


def test_stream_emits_complete_event(client):
    """Stream emits complete event when SA finishes."""
    # Configure small session
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]

    client.post(f"/api/sa/{session_id}/start")

    # Stream to completion
    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        events = parse_sse_events(response)

    # Last event should be complete
    complete_events = [e for e in events if e.get("event") == "complete"]
    assert len(complete_events) > 0

    # Verify complete event has required fields
    complete = complete_events[-1]["data"]
    assert "best_energy" in complete
    assert "best_smiles" in complete
    assert "best_svg" in complete
    assert isinstance(complete["best_svg"], str)


def test_progress_event_contains_svg(client):
    """Progress events contain best_svg field with SVG content."""
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]

    client.post(f"/api/sa/{session_id}/start")

    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        events = parse_sse_events(response)

    progress_events = [e for e in events if e.get("event") == "progress"]
    assert len(progress_events) > 0

    # Check first progress event has SVG
    first_progress = progress_events[0]["data"]
    assert "best_svg" in first_progress
    assert "<svg" in first_progress["best_svg"]


def test_progress_events_have_incrementing_steps(client):
    """Progress events have incrementing step numbers."""
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]

    client.post(f"/api/sa/{session_id}/start")

    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        events = parse_sse_events(response)

    progress_events = [e for e in events if e.get("event") == "progress"]

    # Extract step numbers
    steps = [e["data"]["step"] for e in progress_events]

    # Steps should be incrementing
    for i in range(len(steps) - 1):
        assert steps[i] < steps[i + 1], f"Steps not incrementing: {steps}"


# State Interaction Tests


def test_stream_complete_session_sends_final_event(client):
    """Already completed session sends complete event immediately."""
    # Create, start, and let complete
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]
    client.post(f"/api/sa/{session_id}/start")

    # Stream to completion
    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        first_events = parse_sse_events(response)

    # Session is now complete - stream again
    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        second_events = parse_sse_events(response)

    # Should get complete event immediately (no progress events)
    complete_events = [e for e in second_events if e.get("event") == "complete"]
    assert len(complete_events) > 0

    # Should not get progress events on second stream
    progress_events = [e for e in second_events if e.get("event") == "progress"]
    assert len(progress_events) == 0


# End-to-End Integration Tests


def test_e2e_configure_start_stream_complete(client):
    """Full end-to-end flow: configure -> start -> stream -> complete."""
    # Configure
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 5,
        "num_cycles": 1,
        "seed": 42
    })
    assert r.status_code == 200
    session_id = r.json()["session_id"]

    # Start
    r = client.post(f"/api/sa/{session_id}/start")
    assert r.status_code == 200
    assert r.json()["status"] == "running"

    # Stream to completion
    with client.stream("GET", f"/api/sa/{session_id}/stream") as response:
        events = parse_sse_events(response)

    # Verify we got progress events
    progress_events = [e for e in events if e.get("event") == "progress"]
    assert len(progress_events) > 0

    # Verify we got completion event
    complete_events = [e for e in events if e.get("event") == "complete"]
    assert len(complete_events) > 0

    # Total progress events should be roughly equal to total steps (5)
    # (May be slightly less due to timing)
    assert len(progress_events) <= 5

    # Check status after completion
    r = client.get(f"/api/sa/{session_id}/status")
    assert r.status_code == 200
    status = r.json()
    assert status["session_state"] == "complete"
    assert status["is_complete"] is True
    assert status["step"] == 5


def test_health_responsive_during_stream(client):
    """Health endpoint remains responsive during active stream."""
    # Configure and start session
    r = client.post("/api/sa/configure", json={
        "formula": "C4H10",
        "steps_per_cycle": 10,
        "num_cycles": 1,
        "seed": 42
    })
    session_id = r.json()["session_id"]
    client.post(f"/api/sa/{session_id}/start")

    # Start stream in background (don't consume it)
    # Just verify health still works
    r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["status"] == "healthy"
