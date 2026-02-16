# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-15)

**Core value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time -- making the abstract algorithm from the paper tangible and intuitive.
**Current focus:** v2.0 Phase 4 -- Backend Core & RDKit Foundation

## Current Position

Phase: 4 of 8 (Backend Core & RDKit Foundation)
Plan: 1 of 4 in current phase
Status: Executing
Last activity: 2026-02-16 -- Completed 04-02 (Pure Python Utilities TDD Port)

Progress: [################░░░░░░░░] 12/12 v1.0 plans complete, 1/4 v2.0 Phase 4

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
| 04 | 02 | 280s (4m) | 2 | 10 |

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
- [Phase 04-02]: Use pip for pytest instead of Poetry (not available in environment)
- [Phase 04-02]: Emulate JavaScript Math.imul with signed 32-bit conversion for cross-language determinism

### Pending Todos

None.

### Blockers/Concerns

- Faulon displacement on RDKit RWMol has no existing implementation -- validate against v1.0 test suite once ported
- Offline usage regression: v2.0 requires backend connection (v1.0 works offline)

## Session Continuity

Last session: 2026-02-16
Stopped at: Completed 04-02-PLAN.md (Pure Python Utilities TDD Port)
Resume file: None
