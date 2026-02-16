---
phase: 05
plan: 01
subsystem: api-layer-services
tags: [session-management, svg-rendering, rdkit, tdd, services]
dependency_graph:
  requires:
    - 04-04 (SAEngine step-by-step API)
    - 04-02 (MoleculeGraph with RDKit)
  provides:
    - SessionManager with TTL expiration
    - SVG renderer using RDKit MolDraw2DSVG
  affects:
    - 05-02 (will use SessionManager)
    - 05-03 (will use SVG renderer for SSE events)
tech_stack:
  added:
    - sse-starlette: SSE support for FastAPI
  patterns:
    - TDD RED-GREEN-REFACTOR workflow
    - Session lifecycle management with TTL cleanup
    - Inline SVG generation for web streaming
key_files:
  created:
    - backend/app/services/__init__.py
    - backend/app/services/session_manager.py
    - backend/app/services/svg_renderer.py
    - backend/tests/test_session_manager.py
    - backend/tests/test_svg_renderer.py
  modified:
    - backend/pyproject.toml (added sse-starlette dependency)
    - backend/app/config.py (added session_ttl_seconds setting)
decisions:
  - Use UUID strings for session IDs (standard, secure, collision-resistant)
  - Store session state as mutable string field instead of enum (flexibility for future states)
  - TTL cleanup is manual (cleanup_expired()) not automatic (avoid background threads for YAGNI)
  - SVG renderer uses default 300x300 dimensions (matches typical UI preview size)
  - Error handling uses ValueError for invalid SMILES (consistent with Python conventions)
metrics:
  duration: 192s
  completed_date: 2026-02-16
  tasks_completed: 2
  tests_added: 22
  tests_total: 168
  commits: 4 (2 RED + 2 GREEN)
---

# Phase 5 Plan 01: SessionManager and SVG Renderer Services Summary

**One-liner:** SessionManager with TTL expiration for SA engine lifecycle management, plus RDKit MolDraw2DSVG renderer for inline molecule SVGs.

## What Was Built

### SessionManager Service

A session lifecycle manager that wraps SAEngine instances with TTL-based expiration:

- **Create sessions:** `create(params: SAParams) -> str` generates UUID session IDs, initializes SAEngine (calls `engine.init()`), and stores session with "idle" state
- **Retrieve sessions:** `get(session_id: str) -> Optional[SASession]` returns session and updates `last_accessed` timestamp
- **Delete sessions:** `delete(session_id: str) -> bool` removes session by ID
- **TTL cleanup:** `cleanup_expired() -> int` removes sessions past TTL, returns count removed
- **Active count:** `active_session_count` property for monitoring

**SASession dataclass** tracks:
- `session_id`: UUID string
- `engine`: SAEngine instance
- `created_at`: creation timestamp
- `last_accessed`: last access timestamp (updated on get())
- `state`: mutable state string ("idle", "running", "paused", "complete")

**Configuration:**
- TTL defaults to 3600s (1 hour)
- Configurable via `WEBFAULON_SESSION_TTL_SECONDS` env var

### SVG Renderer Service

RDKit-based SVG generator for molecule visualization:

- **Function:** `generate_molecule_svg(smiles: str, width: int = 300, height: int = 300) -> str`
- **Returns:** Inline SVG string (embeddable in HTML or SSE events)
- **Validation:** Raises `ValueError` for invalid/empty SMILES
- **Performance:** <50ms per molecule (verified via test, suitable for SSE streaming)

Uses `rdkit.Chem.Draw.rdMolDraw2D.MolDraw2DSVG` for 2D depiction.

### Dependencies

- **sse-starlette ^2.0** added to pyproject.toml (needed by Plans 02/03 for SSE endpoints)
- Already installed (version 3.0.4 available in environment)

## Test Coverage

### SessionManager Tests (13 tests)

**Creation tests:**
- UUID string generation and validation
- Session storage and retrieval
- Engine initialization (get_state() works without error)
- Initial state is "idle"
- Multiple sessions have unique IDs

**Retrieval tests:**
- Get existing session returns SASession
- Get nonexistent returns None
- Get updates last_accessed timestamp

**State management:**
- State transitions (idle → running → paused → idle)
- Delete removes session

**TTL cleanup:**
- Expired sessions are removed
- Active sessions survive cleanup
- active_session_count tracks correctly

### SVG Renderer Tests (9 tests)

**Valid SMILES:**
- Generates SVG with `<svg>` tags
- Contains drawing elements (path/ellipse)
- Custom dimensions (400x400)
- Default dimensions (300x300)
- Aromatic SMILES (benzene)
- Branched molecules

**Error handling:**
- Invalid SMILES raises ValueError
- Empty SMILES raises ValueError

**Performance:**
- Generation completes under 50ms

### Test Suite Status

- **New tests:** 22 (13 SessionManager + 9 SVG renderer)
- **Total tests:** 168 (146 existing + 22 new)
- **Status:** All pass, no regressions

## TDD Workflow

Followed RED-GREEN-REFACTOR pattern for both tasks:

### Task 1: SessionManager
1. **RED:** Write 13 failing tests → commit
2. **GREEN:** Implement SessionManager to pass all tests → commit
3. **REFACTOR:** Not needed (implementation was clean)

### Task 2: SVG Renderer
1. **RED:** Write 9 failing tests → commit
2. **GREEN:** Implement generate_molecule_svg to pass all tests → commit
3. **REFACTOR:** Not needed (implementation was straightforward)

**Total commits:** 4 (2 RED + 2 GREEN)

## Deviations from Plan

None. Plan executed exactly as written.

## Integration Points

### Phase 04 Dependencies
- **04-04 SAEngine:** SessionManager calls `engine.init()` on creation, stores initialized engine for step-by-step execution
- **04-02 MoleculeGraph:** SVG renderer uses RDKit integration already present in molecule module

### Phase 05 Forward Dependencies
- **05-02 (Plan 02):** Will use SessionManager to create/retrieve sessions for API endpoints
- **05-03 (Plan 03):** Will use generate_molecule_svg to include molecule visuals in SSE events

## Key Decisions

1. **UUID session IDs:** Standard, secure, collision-resistant
2. **Mutable state string:** Flexibility for future states without enum changes
3. **Manual TTL cleanup:** Avoid background threads (YAGNI), let endpoints call cleanup_expired()
4. **Default 300x300 SVG:** Matches typical UI preview size
5. **ValueError for invalid SMILES:** Consistent with Python conventions

## Success Criteria

- [x] SessionManager creates sessions from SAParams
- [x] SessionManager retrieves by ID with last_accessed update
- [x] SessionManager deletes sessions
- [x] SessionManager cleans up expired sessions based on TTL
- [x] SVG renderer generates valid inline SVG from SMILES
- [x] SVG renderer handles invalid/empty SMILES with ValueError
- [x] All 22 new tests pass
- [x] Full test suite has no regressions (168 total pass)
- [x] sse-starlette dependency installed and available
- [x] Session TTL configurable via Settings

## Next Steps

Phase 5 Plan 02 will build the FastAPI endpoints:
- POST /sessions to create sessions (uses SessionManager)
- GET /sessions/{id} to retrieve session state
- POST /sessions/{id}/start to begin SA execution
- SSE streaming endpoint for real-time updates (uses SVG renderer)

## Self-Check: PASSED

**Files created:**
- backend/app/services/__init__.py: FOUND
- backend/app/services/session_manager.py: FOUND
- backend/app/services/svg_renderer.py: FOUND
- backend/tests/test_session_manager.py: FOUND
- backend/tests/test_svg_renderer.py: FOUND

**Commits:**
- e86c837 (RED - SessionManager tests): FOUND
- e7cf614 (GREEN - SessionManager implementation): FOUND
- 3ccd412 (RED - SVG renderer tests): FOUND
- feef912 (GREEN - SVG renderer implementation): FOUND

**Imports:**
- SessionManager importable: VERIFIED
- generate_molecule_svg generates valid SVG: VERIFIED (2224 chars, has <svg> tag)

**Dependencies:**
- sse-starlette installed: VERIFIED (version 3.0.4)

All artifacts verified successfully.
