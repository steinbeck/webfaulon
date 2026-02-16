"""Tests for SessionManager service.

Following TDD RED-GREEN-REFACTOR pattern.
These tests are written FIRST and will initially fail.
"""

import pytest
from datetime import datetime, timedelta
from uuid import UUID

from app.services.session_manager import SessionManager, SASession
from app.models.sa_params import SAParams


@pytest.fixture
def sample_params():
    """Sample SAParams for testing."""
    return SAParams(
        formula="C6H14",
        initial_temp=100.0,
        cooling_schedule_k=8.0,
        steps_per_cycle=500,
        num_cycles=4,
        optimization_mode="MINIMIZE",
        seed=42
    )


@pytest.fixture
def session_manager():
    """SessionManager instance with short TTL for testing."""
    return SessionManager(ttl_seconds=3600)


# Creation tests
def test_create_returns_uuid_string(session_manager, sample_params):
    """create() should return a valid UUID string."""
    session_id = session_manager.create(sample_params)

    # Should be a string
    assert isinstance(session_id, str)

    # Should be a valid UUID
    uuid_obj = UUID(session_id)
    assert str(uuid_obj) == session_id


def test_create_stores_session(session_manager, sample_params):
    """After create(), get() should return the session."""
    session_id = session_manager.create(sample_params)
    session = session_manager.get(session_id)

    assert session is not None
    assert session.session_id == session_id


def test_create_initializes_engine(session_manager, sample_params):
    """Session engine should be initialized (get_state() works)."""
    session_id = session_manager.create(sample_params)
    session = session_manager.get(session_id)

    # get_state() should not raise an error
    state = session.engine.get_state()
    assert state is not None
    assert state.step == 0
    assert state.is_complete is False


def test_create_sets_idle_state(session_manager, sample_params):
    """Session state should be 'idle' after creation."""
    session_id = session_manager.create(sample_params)
    session = session_manager.get(session_id)

    assert session.state == "idle"


def test_create_multiple_unique_ids(session_manager, sample_params):
    """Creating multiple sessions should return different IDs."""
    id1 = session_manager.create(sample_params)
    id2 = session_manager.create(sample_params)
    id3 = session_manager.create(sample_params)

    assert id1 != id2
    assert id2 != id3
    assert id1 != id3


# Retrieval tests
def test_get_existing_session(session_manager, sample_params):
    """get() should return SASession with correct session_id."""
    session_id = session_manager.create(sample_params)
    session = session_manager.get(session_id)

    assert isinstance(session, SASession)
    assert session.session_id == session_id
    assert session.engine is not None
    assert isinstance(session.created_at, datetime)
    assert isinstance(session.last_accessed, datetime)


def test_get_nonexistent_returns_none(session_manager):
    """get() should return None for unknown session ID."""
    result = session_manager.get("nonexistent-id")
    assert result is None


def test_get_updates_last_accessed(session_manager, sample_params):
    """Accessing a session should update last_accessed timestamp."""
    session_id = session_manager.create(sample_params)
    session1 = session_manager.get(session_id)
    initial_accessed = session1.last_accessed

    # Small delay
    import time
    time.sleep(0.01)

    session2 = session_manager.get(session_id)
    assert session2.last_accessed > initial_accessed


# State management tests
def test_session_state_transitions(session_manager, sample_params):
    """Session state should be mutable."""
    session_id = session_manager.create(sample_params)
    session = session_manager.get(session_id)

    assert session.state == "idle"

    session.state = "running"
    assert session.state == "running"

    session.state = "paused"
    assert session.state == "paused"

    session.state = "idle"
    assert session.state == "idle"


def test_delete_session(session_manager, sample_params):
    """delete() should remove session."""
    session_id = session_manager.create(sample_params)

    # Session exists
    assert session_manager.get(session_id) is not None

    # Delete it
    result = session_manager.delete(session_id)
    assert result is True

    # Should be gone
    assert session_manager.get(session_id) is None

    # Deleting again should return False
    result2 = session_manager.delete(session_id)
    assert result2 is False


# TTL cleanup tests
def test_cleanup_removes_expired(session_manager, sample_params):
    """cleanup_expired() should remove sessions past TTL."""
    session_id = session_manager.create(sample_params)
    session = session_manager.get(session_id)

    # Manually set last_accessed to past
    session.last_accessed = datetime.now() - timedelta(seconds=7200)  # 2 hours ago

    # Run cleanup
    removed_count = session_manager.cleanup_expired()

    assert removed_count == 1
    assert session_manager.get(session_id) is None


def test_cleanup_keeps_active(session_manager, sample_params):
    """cleanup_expired() should NOT remove recently accessed sessions."""
    session_id = session_manager.create(sample_params)

    # Run cleanup immediately
    removed_count = session_manager.cleanup_expired()

    assert removed_count == 0
    assert session_manager.get(session_id) is not None


def test_active_session_count(session_manager, sample_params):
    """active_session_count property should return correct count."""
    assert session_manager.active_session_count == 0

    id1 = session_manager.create(sample_params)
    assert session_manager.active_session_count == 1

    id2 = session_manager.create(sample_params)
    assert session_manager.active_session_count == 2

    session_manager.delete(id1)
    assert session_manager.active_session_count == 1

    session_manager.delete(id2)
    assert session_manager.active_session_count == 0
