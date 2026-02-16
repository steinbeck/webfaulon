---
phase: 05
plan: 02
subsystem: api-layer
tags: [fastapi, rest-api, session-control, endpoints, tdd]
dependency_graph:
  requires:
    - 05-01 (SessionManager and SVG Renderer)
    - 04-04 (SAEngine step-by-step API)
  provides:
    - POST /api/sa/configure endpoint for session creation
    - POST /api/sa/{id}/start|pause|reset for session control
    - GET /api/sa/{id}/status for state snapshots
  affects:
    - 05-03 (SSE streaming will use same session management)
tech_stack:
  added: []
  patterns:
    - Dependency injection via shared singleton
    - Pydantic response models for type safety
    - FastAPI router modular organization
    - Periodic background cleanup on startup
key_files:
  created:
    - backend/app/api/__init__.py
    - backend/app/api/sa_configure.py
    - backend/app/api/sa_control.py
    - backend/app/api/sa_status.py
    - backend/app/dependencies.py
    - backend/tests/test_api_sa.py
  modified:
    - backend/app/main.py
decisions:
  - Shared SessionManager singleton via dependencies.py (simpler than FastAPI app.state)
  - 5-minute periodic cleanup interval (balance between memory usage and overhead)
  - Error handling returns 404 for missing sessions, 400 for invalid state transitions
  - Status endpoint includes SVG inline (not separate endpoint, simplifies client)
metrics:
  duration: 116s
  completed_date: 2026-02-16
  tasks_completed: 2
  tests_added: 17
  tests_total: 185
  commits: 2
---

# Phase 5 Plan 02: SA REST API Endpoints Summary

**One-liner:** FastAPI endpoints for SA session configuration, control (start/pause/reset), and status queries with SVG rendering.

## What Was Built

### API Endpoints

**Configuration Endpoint** (`backend/app/api/sa_configure.py`):
- **POST /api/sa/configure** - Creates new SA session from SAParams
- Accepts SAParams request body (Pydantic auto-validates)
- Calls `session_manager.create(params)` to initialize session
- Returns JSON: `{"session_id": "<uuid>", "status": "configured"}`
- Uses `ConfigureResponse` Pydantic model for type safety

**Control Endpoints** (`backend/app/api/sa_control.py`):
- **POST /api/sa/{session_id}/start** - Transitions session to running state
  - Validates session exists (404 if not found)
  - Validates state transitions (400 if already running/complete)
  - Sets `session.state = "running"`

- **POST /api/sa/{session_id}/pause** - Pauses running session
  - Validates session is currently running (400 otherwise)
  - Sets `session.state = "paused"`

- **POST /api/sa/{session_id}/reset** - Reinitializes session to idle
  - Calls `session.engine.init()` to reset engine state
  - Sets `session.state = "idle"`

**Status Endpoint** (`backend/app/api/sa_status.py`):
- **GET /api/sa/{session_id}/status** - Returns full state snapshot for reconnection
- Gets engine state via `session.engine.get_state()`
- Generates SVG for best molecule via `generate_molecule_svg(state.best_smiles)`
- Returns comprehensive JSON with:
  - Session metadata: `session_id`, `session_state`
  - Progress: `step`, `total_steps`, `cycle`, `is_complete`
  - Energy: `current_energy`, `best_energy`
  - Molecule: `best_smiles`, `best_svg` (inline SVG)
  - Temperature: `temperature`
  - Statistics: `accepted_moves`, `rejected_moves`, `invalid_moves`
- Uses `StatusResponse` Pydantic model for type safety

### Shared Dependencies

**Dependencies Module** (`backend/app/dependencies.py`):
- Creates singleton `SessionManager` instance
- Loads `Settings` from config
- Shared by all API routers (simpler than FastAPI dependency injection)

**Main App Wiring** (`backend/app/main.py`):
- Registers all three routers: `configure_router`, `control_router`, `status_router`
- Adds startup event for periodic TTL cleanup:
  - Runs every 5 minutes (`asyncio.sleep(300)`)
  - Calls `session_manager.cleanup_expired()`
  - Prevents unbounded memory growth from abandoned sessions
- Preserves all existing middleware (CORS first) and exception handlers
- Health check endpoint unaffected

## Test Coverage

### API Test Suite (`backend/tests/test_api_sa.py`)

**17 comprehensive tests** organized by endpoint:

**Configure endpoint (4 tests):**
1. `test_configure_returns_session_id` - Valid params return UUID session_id and "configured" status
2. `test_configure_with_custom_params` - All SAParams fields accepted
3. `test_configure_invalid_formula` - Empty formula returns 422 validation error
4. `test_configure_invalid_params` - Negative temperature returns 422

**Control endpoints (7 tests):**
1. `test_start_session` - Start idle session returns "running" status
2. `test_start_already_running` - Starting running session returns 400
3. `test_start_nonexistent_session` - Nonexistent session returns 404
4. `test_pause_running_session` - Pause running session returns "paused" status
5. `test_pause_not_running` - Pausing idle session returns 400
6. `test_reset_session` - Reset reinitializes to "idle" status
7. `test_reset_nonexistent` - Nonexistent session returns 404

**Status endpoint (4 tests):**
1. `test_status_returns_state` - Returns all 14 required fields
2. `test_status_includes_svg` - `best_svg` contains `<svg>` tag
3. `test_status_nonexistent` - Nonexistent session returns 404
4. `test_status_after_start` - `session_state` reflects current state

**Integration tests (2 tests):**
1. `test_health_check_still_works` - Regression check for /health endpoint
2. `test_full_configure_start_pause_reset_flow` - End-to-end state transitions

**Test infrastructure:**
- `clear_sessions` autouse fixture clears session manager between tests (prevents state leakage)
- Uses FastAPI `TestClient` for synchronous testing
- All tests use fresh sessions per test case

### Test Results

- **New tests:** 17
- **Total tests:** 185 (168 existing + 17 new)
- **Status:** All pass, no regressions
- **Coverage:** Configure, control, status endpoints fully tested with error cases

## Verification Results

**Unit test verification:**
```
cd backend && python -m pytest tests/test_api_sa.py -v
17 passed in 0.06s
```

**Full suite regression check:**
```
cd backend && python -m pytest tests/ -v
185 passed in 1.02s
```

**Manual flow verification:**
```python
# Configure -> Status -> Start -> Pause -> Reset -> Health
Configure: 200 {'session_id': '...', 'status': 'configured'}
Status: 200, has svg: True
Start: 200 {'session_id': '...', 'status': 'running'}
Pause: 200 {'session_id': '...', 'status': 'paused'}
Reset: 200 {'session_id': '...', 'status': 'idle'}
Health: 200 {'status': 'healthy', 'version': '2.0.0'}
```

All endpoints return expected responses. Full flow works correctly.

## Deviations from Plan

None. Plan executed exactly as written.

## Integration Points

### Upstream Dependencies

**05-01 SessionManager** (Plan 01):
- `session_manager.create(params)` in configure endpoint
- `session_manager.get(session_id)` in all control/status endpoints
- TTL cleanup scheduled in main.py startup event

**05-01 SVG Renderer** (Plan 01):
- `generate_molecule_svg(state.best_smiles)` in status endpoint
- Inline SVG included in status response (not separate endpoint)

**04-04 SAEngine** (Plan 04):
- `session.engine.get_state()` for status snapshots
- `session.engine.init()` for session reset

### Downstream Impact

**05-03 SSE Streaming** (Plan 03):
- Will use same SessionManager to retrieve sessions
- Will use GET /api/sa/{id}/status pattern for reconnection
- Will generate SVGs via same renderer for real-time updates

## Key Decisions

1. **Shared singleton via dependencies.py** - Simpler than FastAPI app.state or dependency injection function. All routers import directly.

2. **5-minute periodic cleanup** - Balances memory management with overhead. Prevents unbounded growth from abandoned sessions without excessive cleanup runs.

3. **404 for missing sessions, 400 for invalid states** - Standard HTTP semantics. 404 = resource not found. 400 = bad request due to invalid state transition.

4. **Inline SVG in status response** - Simplifies client (no separate endpoint). SVGs are small (~2KB), acceptable for JSON payload.

5. **Pydantic response models** - Type safety and auto-generated OpenAPI docs. ConfigureResponse, ControlResponse, StatusResponse all explicitly typed.

6. **State validation in control endpoints** - Enforces valid transitions (can't start running session, can't pause idle session). Prevents invalid states.

## Success Criteria

All success criteria met:

- [x] POST /api/sa/configure creates session and returns session_id
- [x] POST /api/sa/{id}/start|pause|reset control session state transitions
- [x] GET /api/sa/{id}/status returns full state snapshot with SVG
- [x] All endpoints return 404 for nonexistent sessions
- [x] Control endpoints return 400 for invalid state transitions
- [x] Health check unaffected by new endpoints (regression test passes)
- [x] Periodic TTL cleanup configured on startup
- [x] 17 comprehensive API tests pass
- [x] Full test suite (185 tests) passes with no regressions
- [x] Manual flow verification successful

## Next Steps

Phase 5 Plan 03 will implement SSE streaming:
- POST /api/sa/{id}/stream endpoint for Server-Sent Events
- Real-time progress updates with SVG rendering
- Automatic execution in background task
- Client reconnection support via status endpoint
- Event stream format: `data: {"step": N, "best_svg": "...", ...}`

## Self-Check: PASSED

**Files created:**
- backend/app/api/__init__.py: FOUND
- backend/app/api/sa_configure.py: FOUND
- backend/app/api/sa_control.py: FOUND
- backend/app/api/sa_status.py: FOUND
- backend/app/dependencies.py: FOUND
- backend/tests/test_api_sa.py: FOUND

**Files modified:**
- backend/app/main.py: FOUND (routers registered, startup event added)

**Commits:**
- cb6fb56 (Task 1 - API endpoints): FOUND
- c8c555c (Task 2 - Tests and wiring): FOUND

**Imports verified:**
```python
from app.api.sa_configure import router as configure_router  # OK
from app.api.sa_control import router as control_router      # OK
from app.api.sa_status import router as status_router        # OK
```

**Routes verified:**
- Configure: ['/api/sa/configure']
- Control: ['/api/sa/{session_id}/start', '/api/sa/{session_id}/pause', '/api/sa/{session_id}/reset']
- Status: ['/api/sa/{session_id}/status']

**Tests verified:**
- 17 new tests pass
- 185 total tests pass
- Manual flow verification successful

All artifacts verified successfully.
