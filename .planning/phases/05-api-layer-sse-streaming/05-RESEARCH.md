# Phase 5: API Layer & SSE Streaming - Research

**Researched:** 2026-02-16
**Domain:** FastAPI REST + SSE hybrid API, real-time streaming
**Confidence:** HIGH

## Summary

Phase 5 implements a REST + SSE hybrid API for configuring and streaming real-time SA optimization progress. The standard stack combines FastAPI's native async capabilities with sse-starlette for production-ready Server-Sent Events, FastAPI BackgroundTasks for lightweight background execution, and in-memory dictionary-based session management with TTL cleanup.

FastAPI's async foundation makes it exceptionally well-suited for SSE streaming. The sse-starlette library (production-ready, W3C SSE spec compliant) provides EventSourceResponse for streaming, automatic client disconnect detection, and multi-threaded safety. For this classroom/demo use case, FastAPI's built-in BackgroundTasks is appropriate—Celery/Redis would be over-engineering for ephemeral, single-server optimization sessions.

The architecture separates concerns cleanly: POST endpoints configure sessions and return IDs, SSE endpoints stream real-time progress via async generators, control endpoints (start/pause/reset) manipulate session state, and status endpoints enable client reconnection. RDKit's MolDraw2DSVG generates 2D molecule SVGs directly in Python, eliminating need for external rendering services.

**Primary recommendation:** Use sse-starlette's EventSourceResponse with async generators yielding SAEngineState snapshots; run SA in same event loop (not separate process) using SAEngine's step-by-step API; store session state in dict with datetime-based TTL tracking; leverage request.is_disconnected() for cleanup; set X-Accel-Buffering: no header for nginx compatibility.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| sse-starlette | 2.x+ | Production SSE streaming | W3C spec compliant, automatic disconnect detection, thread-safe, FastAPI native |
| FastAPI | 0.129+ | REST API + async runtime | Already in use, native async/await, auto OpenAPI docs, excellent SSE support |
| asyncio | stdlib | Event loop management | Python standard library, FastAPI foundation, no external deps |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| RDKit | (existing) | SVG molecule rendering via MolDraw2DSVG | Already in use, generates inline SVG strings |
| Pydantic v2 | 2.12+ | Request/response validation | Already in use, SAParams/SAResult models exist |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| sse-starlette | fastapi-sse (PyPI) | Less mature, fewer features, no multi-threading support documented |
| BackgroundTasks | Celery + Redis | Over-engineering for classroom use; adds infrastructure complexity, deployment burden |
| In-memory dict | Redis session store | Overkill for ephemeral sessions; no persistence needed, single-server deployment |
| asyncio tasks | Threading/multiprocessing | Complicates state sharing; async already sufficient for I/O-bound SSE streaming |

**Installation:**
```bash
poetry add sse-starlette
```

## Architecture Patterns

### Recommended Project Structure
```
backend/app/
├── api/                  # API endpoints
│   ├── __init__.py
│   ├── sa_configure.py   # POST /api/sa/configure → session_id
│   ├── sa_stream.py      # GET /api/sa/{session_id}/stream → SSE
│   ├── sa_control.py     # POST /api/sa/{session_id}/start|pause|reset
│   └── sa_status.py      # GET /api/sa/{session_id}/status → state snapshot
├── services/             # Business logic
│   ├── __init__.py
│   ├── session_manager.py  # Session CRUD + TTL cleanup
│   └── sa_runner.py      # SA execution orchestration
├── models/               # (existing) Pydantic models
│   └── sa_params.py      # SAParams, SAResult, SAEngineState
└── main.py               # FastAPI app with router registration
```

### Pattern 1: SSE Streaming with EventSourceResponse
**What:** Async generator yielding SSE events, wrapped in EventSourceResponse
**When to use:** Streaming SA progress, real-time updates from background computation
**Example:**
```python
# Source: https://github.com/sysid/sse-starlette
from sse_starlette.sse import EventSourceResponse
from fastapi import Request

async def sa_progress_stream(request: Request, session_id: str):
    """Stream SA optimization progress via SSE."""
    async def event_generator():
        session = session_manager.get(session_id)
        while not session.engine.is_complete():
            # Check for client disconnect
            if await request.is_disconnected():
                break

            # Execute one SA step
            session.engine.step()
            state = session.engine.get_state()

            # Yield SSE event with id, event type, and data
            yield {
                "id": str(state.step),
                "event": "progress",
                "data": state.model_dump_json()
            }

            await asyncio.sleep(0.01)  # Throttle to ~100 events/sec

        # Send completion event
        result = session.engine.get_result()
        yield {
            "event": "complete",
            "data": result.model_dump_json()
        }

    return EventSourceResponse(event_generator())
```

### Pattern 2: Session State Management with TTL
**What:** In-memory dict storing session state with timestamp-based expiration
**When to use:** Ephemeral session management, classroom/demo deployments
**Example:**
```python
from datetime import datetime, timedelta
from typing import Dict, Optional
from dataclasses import dataclass

@dataclass
class SASession:
    session_id: str
    engine: SAEngine
    created_at: datetime
    last_accessed: datetime
    state: str  # "idle", "running", "paused", "complete"

class SessionManager:
    def __init__(self, ttl_seconds: int = 3600):
        self._sessions: Dict[str, SASession] = {}
        self._ttl = timedelta(seconds=ttl_seconds)

    def create(self, params: SAParams) -> str:
        """Create new SA session, return session ID."""
        session_id = str(uuid.uuid4())
        engine = SAEngine(params)
        engine.init()  # Initialize SA state

        self._sessions[session_id] = SASession(
            session_id=session_id,
            engine=engine,
            created_at=datetime.now(),
            last_accessed=datetime.now(),
            state="idle"
        )
        return session_id

    def get(self, session_id: str) -> Optional[SASession]:
        """Get session, update last_accessed, return None if expired."""
        session = self._sessions.get(session_id)
        if session:
            session.last_accessed = datetime.now()
            return session
        return None

    def cleanup_expired(self):
        """Remove sessions exceeding TTL (call periodically)."""
        now = datetime.now()
        expired = [
            sid for sid, sess in self._sessions.items()
            if now - sess.last_accessed > self._ttl
        ]
        for sid in expired:
            del self._sessions[sid]
```

### Pattern 3: Control Endpoints for State Manipulation
**What:** POST endpoints that modify session state (start/pause/reset)
**When to use:** User-controlled SA execution flow
**Example:**
```python
@router.post("/{session_id}/start")
async def start_optimization(session_id: str):
    """Start SA optimization (transitions from idle/paused to running)."""
    session = session_manager.get(session_id)
    if not session:
        raise HTTPException(404, "Session not found")

    if session.state == "running":
        raise HTTPException(400, "Already running")

    session.state = "running"
    return {"status": "started"}

@router.post("/{session_id}/pause")
async def pause_optimization(session_id: str):
    """Pause SA optimization (halts iteration)."""
    session = session_manager.get(session_id)
    if not session:
        raise HTTPException(404, "Session not found")

    session.state = "paused"
    return {"status": "paused"}

@router.post("/{session_id}/reset")
async def reset_optimization(session_id: str):
    """Reset SA to initial state."""
    session = session_manager.get(session_id)
    if not session:
        raise HTTPException(404, "Session not found")

    session.engine.init()  # Re-initialize SA
    session.state = "idle"
    return {"status": "reset"}
```

### Pattern 4: Reconnection Support with Status Endpoint
**What:** GET endpoint returning current SA state for client reconnection
**When to use:** Client network interruption, page refresh, manual reconnect
**Example:**
```python
@router.get("/{session_id}/status")
async def get_status(session_id: str):
    """Get current SA state (for reconnection)."""
    session = session_manager.get(session_id)
    if not session:
        raise HTTPException(404, "Session not found or expired")

    state = session.engine.get_state()
    return {
        "session_id": session_id,
        "state": session.state,
        "current_step": state.step,
        "total_steps": state.total_steps,
        "best_energy": state.best_energy,
        "best_smiles": state.best_smiles,
        "is_complete": state.is_complete
    }
```

### Pattern 5: SVG Generation for Molecule Depiction
**What:** RDKit MolDraw2DSVG generates inline SVG from SMILES
**When to use:** Streaming molecule depictions to frontend
**Example:**
```python
# Source: https://www.rdkit.org/docs/GettingStartedInPython.html
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def generate_molecule_svg(smiles: str, width: int = 300, height: int = 300) -> str:
    """Generate 2D SVG depiction of molecule from SMILES.

    Args:
        smiles: SMILES notation string
        width: SVG width in pixels
        height: SVG height in pixels

    Returns:
        SVG string (inline, can be embedded in HTML or sent via SSE)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()

    return drawer.GetDrawingText()
```

### Anti-Patterns to Avoid
- **Separate processes for SA execution:** SA is CPU-bound but short-lived (seconds); async event loop handles streaming overhead, no need for multiprocessing complexity
- **Storing MoleculeGraph objects in SSE events:** Convert to SMILES first (memory efficiency, JSON serialization)
- **Missing X-Accel-Buffering header:** Causes nginx buffering, breaks real-time streaming
- **Infinite SSE loops without disconnect check:** Leaks connections when clients navigate away
- **Global session dict without locks:** Race conditions on concurrent access (use single-threaded async event loop, no need for threading.Lock)

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SSE event formatting | Custom SSE message serializer | sse-starlette EventSourceResponse | Handles id/event/data fields, retry, comments, multi-line data, disconnect detection automatically |
| Session expiration | Manual datetime checks everywhere | Centralized SessionManager.cleanup_expired() | Consolidates TTL logic, prevents scattered expiration checks |
| Client disconnect detection | Polling request state | await request.is_disconnected() | Native Starlette support, async-friendly, prevents busy-waiting |
| SVG molecule rendering | Canvas/D3.js depiction | RDKit MolDraw2DSVG | Handles stereochemistry, aromaticity, wedge bonds—deceptively complex domain |
| Event stream reconnection | Custom Last-Event-ID parsing | SSE protocol's built-in id field | Browsers automatically send Last-Event-ID header on reconnect |

**Key insight:** SSE protocol has subtle edge cases (multi-line data escaping, comment handling, retry timing, browser reconnection). sse-starlette abstracts W3C spec complexity; RDKit abstracts chemoinformatics visualization; native FastAPI/Starlette features handle disconnect detection.

## Common Pitfalls

### Pitfall 1: Nginx Buffering Breaks Real-Time Streaming
**What goes wrong:** SSE events don't appear in frontend until ~16KB accumulates, defeating real-time purpose
**Why it happens:** nginx buffers HTTP responses by default for efficiency; SSE needs immediate delivery
**How to avoid:** Set `X-Accel-Buffering: no` header in EventSourceResponse or configure nginx with `proxy_buffering off`
**Warning signs:** Events arrive in large batches after long delays; localhost works but production doesn't

### Pitfall 2: Missing Client Disconnect Detection Causes Memory Leaks
**What goes wrong:** Background generators keep running after client closes tab/navigates away, accumulating memory
**Why it happens:** Async generators don't automatically stop when client disconnects
**How to avoid:** Check `await request.is_disconnected()` in event loop; clean up session state on disconnect
**Warning signs:** Server memory grows over time; active SSE connections increase without bound

### Pitfall 3: BackgroundTasks Unsuitable for Long-Running SA Optimization
**What goes wrong:** BackgroundTasks tie up event loop; server becomes unresponsive during SA execution
**Why it happens:** BackgroundTasks run in same event loop AFTER response is sent, not truly concurrent
**How to avoid:** Use pattern from Pattern 1—run SA steps INSIDE SSE generator with await asyncio.sleep() for yielding control
**Warning signs:** Health endpoint times out during SA execution; concurrent requests block

### Pitfall 4: Storing Engine State in SSE Events Instead of Session
**What goes wrong:** No way to reconnect—state only exists in SSE stream, lost on disconnect
**Why it happens:** Treating SSE as the state store rather than a projection of session state
**How to avoid:** Store SAEngine in SessionManager; SSE generator calls engine.get_state() each iteration
**Warning signs:** Reconnection requires starting over; pause/resume impossible

### Pitfall 5: Missing id Field in SSE Events Breaks Reconnection
**What goes wrong:** Browser reconnects but server can't determine where client left off
**Why it happens:** Omitting id field in yielded dict; Last-Event-ID header never set by browser
**How to avoid:** Always include `"id": str(state.step)` in yielded SSE events
**Warning signs:** EventSource shows reconnecting but stream starts from beginning

### Pitfall 6: Race Condition on Session State Access
**What goes wrong:** Control endpoints modify session.state while SSE generator reads it, causing inconsistent behavior
**Why it happens:** Multiple coroutines (SSE stream, control endpoints) access same session object concurrently
**How to avoid:** Use asyncio.Lock per session OR design state machine so transitions are atomic (e.g., only SSE generator modifies engine state)
**Warning signs:** Pause doesn't stop iteration; reset doesn't clear history

### Pitfall 7: Forgetting to Initialize SAEngine Before Streaming
**What goes wrong:** engine.step() raises RuntimeError because init() wasn't called
**Why it happens:** Session created but engine.init() never invoked
**How to avoid:** Call engine.init() in SessionManager.create() immediately after instantiation
**Warning signs:** First SSE connection fails with "step() called before init()"

## Code Examples

Verified patterns from official sources:

### SSE Event Format with All Fields
```python
# Source: https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events
# SSE events support id, event, data, and retry fields

async def event_generator():
    yield {
        "id": "1",              # Browser stores as Last-Event-ID
        "event": "progress",    # Client listens via addEventListener("progress", ...)
        "data": json.dumps({    # Payload (can be multi-line, auto-escaped)
            "step": 1,
            "energy": 42.5
        }),
        "retry": 3000          # Reconnection interval in ms (optional)
    }
```

### Periodic Session Cleanup (Background Task)
```python
# Called on FastAPI startup event
import asyncio

async def periodic_cleanup(session_manager: SessionManager, interval: int = 300):
    """Remove expired sessions every 5 minutes."""
    while True:
        await asyncio.sleep(interval)
        session_manager.cleanup_expired()

@app.on_event("startup")
async def startup_event():
    asyncio.create_task(periodic_cleanup(session_manager))
```

### Complete SSE Endpoint with Throttling and Disconnect Handling
```python
from sse_starlette.sse import EventSourceResponse
from fastapi import APIRouter, Request, HTTPException
import asyncio

router = APIRouter(prefix="/api/sa", tags=["SA Streaming"])

@router.get("/{session_id}/stream")
async def stream_progress(session_id: str, request: Request):
    """Stream SA optimization progress via SSE.

    Events:
    - progress: step-by-step state updates (id = step number)
    - complete: final result with full history
    - error: optimization error
    """
    session = session_manager.get(session_id)
    if not session:
        raise HTTPException(404, "Session not found or expired")

    async def event_generator():
        try:
            while not session.engine.is_complete():
                # Check client disconnect
                if await request.is_disconnected():
                    break

                # Only step if running (respect pause state)
                if session.state == "running":
                    session.engine.step()
                    state = session.engine.get_state()

                    # Generate SVG for best molecule
                    svg = generate_molecule_svg(state.best_smiles)

                    yield {
                        "id": str(state.step),
                        "event": "progress",
                        "data": json.dumps({
                            "step": state.step,
                            "temperature": state.temperature,
                            "best_energy": state.best_energy,
                            "best_svg": svg,
                            "is_complete": state.is_complete
                        })
                    }

                # Throttle to ~100 events/sec
                await asyncio.sleep(0.01)

            # Send completion event
            if session.engine.is_complete():
                result = session.engine.get_result()
                svg = generate_molecule_svg(result.best_smiles)

                yield {
                    "event": "complete",
                    "data": json.dumps({
                        "best_energy": result.best_energy,
                        "best_smiles": result.best_smiles,
                        "best_svg": svg,
                        "acceptance_ratio": result.acceptance_ratio
                    })
                }

        except Exception as e:
            yield {
                "event": "error",
                "data": json.dumps({"error": str(e)})
            }

    return EventSourceResponse(
        event_generator(),
        headers={"X-Accel-Buffering": "no"}  # Disable nginx buffering
    )
```

### Configure Endpoint (Session Creation)
```python
from pydantic import BaseModel
import uuid

class ConfigureResponse(BaseModel):
    session_id: str
    status: str

@router.post("/configure", response_model=ConfigureResponse)
async def configure_sa(params: SAParams):
    """Create new SA session with configuration.

    Returns session ID for subsequent operations.
    """
    session_id = session_manager.create(params)
    return ConfigureResponse(
        session_id=session_id,
        status="configured"
    )
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| WebSockets for server push | SSE for unidirectional streaming | ~2020-2021 | SSE simpler for one-way updates; automatic reconnection; HTTP/2 eliminates 6-connection limit |
| Celery for all background tasks | FastAPI BackgroundTasks for lightweight tasks | ~2021-2022 | Reduced infrastructure for simple use cases; Celery still needed for distributed/persistent queues |
| Manual SSE formatting | sse-starlette library | ~2020 | Production-ready W3C compliance; automatic disconnect detection |
| Storing full MoleculeGraph in session | Storing SMILES strings | Phase 4 decision | Memory efficiency; JSON serialization; reconstruct graph when needed |
| Synchronous RDKit drawing | Same (RDKit is sync) | N/A | RDKit not async-native; but fast enough (<10ms) to call directly in async generators |

**Deprecated/outdated:**
- **fastapi-sse-starlette (PyPI):** Early SSE library, superseded by sse-starlette with better features
- **StreamingResponse with manual SSE formatting:** sse-starlette's EventSourceResponse handles edge cases correctly
- **HTTP/1.1 6-connection limit concern:** HTTP/2 (default in modern deployments) supports 100+ simultaneous streams

## Open Questions

1. **Should we implement Last-Event-ID resume logic?**
   - What we know: SSE protocol supports Last-Event-ID header; browser sends it on reconnect
   - What's unclear: Whether classroom use case needs true resume (vs. just showing current state)
   - Recommendation: Start with status endpoint for reconnection (simpler); add Last-Event-ID replay if users request it

2. **What's the optimal throttling interval for SSE events?**
   - What we know: await asyncio.sleep(0.01) = 100 events/sec; too fast overwhelms browser rendering
   - What's unclear: Ideal balance between smooth animation and performance
   - Recommendation: Start with 0.01 (100/sec); make configurable via query param; monitor browser CPU usage

3. **Should sessions persist across server restarts?**
   - What we know: Current design uses in-memory dict (ephemeral)
   - What's unclear: Whether classroom use case tolerates losing state on deploy/restart
   - Recommendation: Start ephemeral; add pickle-based persistence if requested (YAGNI principle)

4. **How should we handle concurrent start/pause/reset requests?**
   - What we know: Multiple browser tabs could control same session; race conditions possible
   - What's unclear: Whether to use locks, reject concurrent operations, or allow "last write wins"
   - Recommendation: Use asyncio.Lock per session OR design idempotent operations (e.g., pause sets flag, SSE generator checks flag)

## Sources

### Primary (HIGH confidence)
- [sse-starlette GitHub](https://github.com/sysid/sse-starlette) - EventSourceResponse API, disconnect detection, threading support
- [FastAPI Background Tasks Official](https://fastapi.tiangolo.com/tutorial/background-tasks/) - BackgroundTasks usage, limitations, when to use Celery
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html) - MolDraw2DSVG workflow, API methods
- [MDN Server-Sent Events](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events) - SSE event format, Last-Event-ID, retry mechanism

### Secondary (MEDIUM confidence)
- [FastAPI and SSE: How to Build Streamable MCP Servers](https://www.aubergine.co/insights/a-guide-to-building-streamable-mcp-servers-with-fastapi-and-sse) - Session ID patterns, authentication considerations
- [How to use Server-Sent Events with FastAPI and React](https://www.softgrade.org/sse-with-fastapi-react-langgraph/) - Nginx X-Accel-Buffering configuration
- [FastAPI Concurrency Documentation](https://fastapi.tiangolo.com/async/) - Async/await patterns, event loop behavior
- [ExpiringDict GitHub](https://github.com/mailgun/expiringdict) - TTL dict patterns (reference for manual implementation)
- [Understanding client disconnection in FastAPI](https://fastapiexpert.com/blog/2024/06/06/understanding-client-disconnection-in-fastapi/) - request.is_disconnected() behavior

### Tertiary (LOW confidence)
- [Server-Sent Events: A Comprehensive Guide](https://medium.com/@moali314/server-sent-events-a-comprehensive-guide-e4b15d147576) - SSE overview, general concepts
- [FastAPI Best Practices 2026](https://fastlaunchapi.dev/blog/fastapi-best-practices-production-2026) - Production deployment patterns

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - sse-starlette is mature, well-documented, W3C compliant; FastAPI async native
- Architecture: HIGH - Patterns verified against official docs, GitHub examples, production deployments
- Pitfalls: MEDIUM-HIGH - nginx buffering verified in official docs; disconnect detection confirmed; BackgroundTasks limitation in official FastAPI docs; race condition pitfall inferred from async patterns

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (30 days - stable ecosystem, FastAPI/Starlette mature)
