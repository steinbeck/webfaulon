"""Session management service for SA execution.

Manages session lifecycle (create, get, delete, expire) with TTL-based cleanup.
Each session wraps an SAEngine instance with state tracking.
"""

import uuid
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Optional

from app.core.sa_engine import SAEngine
from app.models.sa_params import SAParams


@dataclass
class SASession:
    """Session wrapping an SAEngine instance.

    Attributes:
        session_id: Unique session identifier (UUID string)
        engine: SAEngine instance for this session
        created_at: Session creation timestamp
        last_accessed: Last access timestamp (updated on get())
        state: Session state ("idle", "running", "paused", "complete")
    """

    session_id: str
    engine: SAEngine
    created_at: datetime
    last_accessed: datetime
    state: str  # "idle", "running", "paused", "complete"


class SessionManager:
    """Manages SA execution sessions with TTL-based expiration.

    Usage:
        manager = SessionManager(ttl_seconds=3600)
        session_id = manager.create(params)
        session = manager.get(session_id)
        session.engine.step()
        manager.cleanup_expired()
    """

    def __init__(self, ttl_seconds: int = 3600):
        """Initialize SessionManager.

        Args:
            ttl_seconds: Session TTL in seconds (default 1 hour)
        """
        self._ttl = timedelta(seconds=ttl_seconds)
        self._sessions: dict[str, SASession] = {}

    def create(self, params: SAParams) -> str:
        """Create new session from SAParams.

        Args:
            params: SA configuration parameters

        Returns:
            Unique session ID (UUID string)
        """
        # Generate unique ID
        session_id = str(uuid.uuid4())

        # Create and initialize engine
        engine = SAEngine(params)
        engine.init()

        # Create session
        now = datetime.now()
        session = SASession(
            session_id=session_id,
            engine=engine,
            created_at=now,
            last_accessed=now,
            state="idle"
        )

        # Store session
        self._sessions[session_id] = session

        return session_id

    def get(self, session_id: str) -> Optional[SASession]:
        """Retrieve session by ID.

        Updates last_accessed timestamp on successful retrieval.

        Args:
            session_id: Session ID to retrieve

        Returns:
            SASession if exists, None otherwise
        """
        session = self._sessions.get(session_id)

        if session is not None:
            # Update last_accessed timestamp
            session.last_accessed = datetime.now()

        return session

    def delete(self, session_id: str) -> bool:
        """Delete session by ID.

        Args:
            session_id: Session ID to delete

        Returns:
            True if session existed and was deleted, False otherwise
        """
        if session_id in self._sessions:
            del self._sessions[session_id]
            return True
        return False

    def cleanup_expired(self) -> int:
        """Remove sessions that exceed TTL.

        Returns:
            Number of sessions removed
        """
        now = datetime.now()
        expired_ids = [
            session_id
            for session_id, session in self._sessions.items()
            if now - session.last_accessed > self._ttl
        ]

        for session_id in expired_ids:
            del self._sessions[session_id]

        return len(expired_ids)

    @property
    def active_session_count(self) -> int:
        """Get count of active sessions.

        Returns:
            Number of sessions currently stored
        """
        return len(self._sessions)
