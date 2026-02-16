# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time -- making the abstract algorithm from the paper tangible and intuitive.
**Current focus:** v2.0 Phase 6 -- Frontend Integration (In Progress)

## Current Position

Phase: 6 of 8 (Frontend Integration)
Plan: 1 of 2 in current phase
Status: In Progress
Last activity: 2026-02-16 -- Completed 06-01 (API Client Infrastructure)

Progress: [##################░░░░░░] 12/12 v1.0 plans complete, 4/4 v2.0 Phase 4 complete, 3/3 v2.0 Phase 5 complete, 1/2 v2.0 Phase 6 complete

## Performance Metrics

**Velocity (v1.0):**
- Total plans completed: 12
- Total execution time: ~2.5 hours

**By Phase (v1.0):**

| Phase | Plans | Status |
|-------|-------|--------|
| 1. MolGraph & SA Core | 4 | Complete |
| 2. Browser Integration | 3 | Complete |
| 3. Visualization & UX | 3 | Complete |
| 3.1. README & Deploy | 2 | Complete |

**By Phase (v2.0):**

| Phase | Plan | Duration | Tasks | Files |
|-------|------|----------|-------|-------|
| 04 | 01 | 340s (5m 40s) | 2 | 9 |
| 04 | 02 | 280s (4m 40s) | 2 | 10 |
| 04 | 03 | 347s (5m 47s) | 2 | 9 |
| 04 | 04 | 251s (4m 11s) | 2 | 3 |
| 05 | 01 | 192s (3m 12s) | 2 | 5 |
| 05 | 02 | 116s (1m 56s) | 2 | 7 |
| 05 | 03 | 461s (7m 41s) | 2 | 3 |
| 06 | 01 | 144s (2m 24s) | 2 | 5 |

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Carried forward from v1.0:

- MOL block pipeline proved more reliable than custom SMILES generation
- RDKit.js set_new_coords() required before draw_to_canvas() (zero-coord crash)
- unpkg CDN for RDKit.js in production (node_modules unavailable in deployed build)

v2.0 architectural decisions (from research):

- FastAPI + sse-starlette for backend (async, SSE native, auto OpenAPI docs)
- BackgroundTasks for SA execution (not Celery/Redis -- YAGNI)
- In-memory dict for session state (not database -- ephemeral classroom use)
- Monorepo: backend/ and frontend/ as sibling directories
- [Phase 04-01]: Poetry for dependency management with virtualenvs.in-project=true for local .venv
- [Phase 04-01]: Pydantic Settings with WEBFAULON_ env prefix for configuration
- [Phase 04-01]: CORS middleware added FIRST before other middleware
- [Phase 04-01]: StarletteHTTPException handler needed to catch default 404 responses
- [Phase 04-02]: Use pip for pytest instead of Poetry (not available in environment)
- [Phase 04-02]: Emulate JavaScript Math.imul with signed 32-bit conversion for cross-language determinism
- [Phase 04-03]: Wrap RDKit RWMol instead of pure adjacency matrix for automatic valence validation
- [Phase 04-03]: Call SanitizeMol after every bond modification for chemical correctness
- [Phase 04-03]: Use RDKit GetDistanceMatrix for Wiener Index instead of custom BFS
- [Phase 04-04]: Store SMILES strings in SAResult instead of MoleculeGraph objects (memory optimization)
- [Phase 04-04]: Implement step-by-step API (init/step/get_state/get_result) for SSE streaming support
- [Phase 04-04]: Use Pydantic v2 models for SA parameters/results for type safety and auto-validation
- [Phase 05-01]: Use UUID strings for session IDs (standard, secure, collision-resistant)
- [Phase 05-01]: Manual TTL cleanup (cleanup_expired()) not automatic to avoid background threads (YAGNI)
- [Phase 05-02]: Shared SessionManager singleton via dependencies.py (simpler than FastAPI app.state)
- [Phase 05-02]: 5-minute periodic cleanup interval (balance between memory usage and overhead)
- [Phase 05-02]: Status endpoint includes SVG inline (not separate endpoint, simplifies client)
- [Phase 05-03]: SA steps run inside SSE generator (not BackgroundTasks) for pause/resume control
- [Phase 05-03]: asyncio.sleep(0.01) yields control between steps (~100 events/sec throttling)
- [Phase 05-03]: Waiting heartbeat every ~1s during pause/idle keeps connection alive
- [Phase 05-03]: X-Accel-Buffering: no header for nginx real-time delivery
- [Phase 06-01]: Use snake_case in TypeScript types to match FastAPI JSON exactly (no camelCase conversion)
- [Phase 06-01]: SSE wrapper auto-closes on complete event to prevent EventSource leaks
- [Phase 06-01]: SVG rendering via innerHTML replaces RDKit.js WASM canvas (backend generates SVG)
- [Phase 06-01]: Vite proxy handles /api routing without URL rewriting (backend already serves /api/sa/*)

### Pending Todos

None.

### Blockers/Concerns

- ~~Faulon displacement on RDKit RWMol has no existing implementation -- validate against v1.0 test suite once ported~~ **RESOLVED: 500-displacement stress test produces zero invalid molecules, all Wiener Index values match v1.0**
- Offline usage regression: v2.0 requires backend connection (v1.0 works offline)

## Session Continuity

Last session: 2026-02-16
Stopped at: Completed 06-01-PLAN.md (API Client Infrastructure)
Resume file: None
