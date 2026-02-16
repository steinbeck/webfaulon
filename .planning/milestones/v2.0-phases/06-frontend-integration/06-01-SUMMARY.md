---
phase: 06-frontend-integration
plan: 01
subsystem: frontend-api
tags: [typescript, sse, fetch, vite, svg]

# Dependency graph
requires:
  - phase: 05-api-layer-sse-streaming
    provides: FastAPI REST endpoints and SSE streaming endpoint with JSON payloads
provides:
  - TypeScript types mirroring backend Pydantic models (snake_case)
  - REST API client for configure/start/pause/reset/getStatus endpoints
  - SSE wrapper with typed handlers for progress/complete/waiting/error events
  - SVG-based molecule renderer replacing RDKit.js WASM canvas
  - Vite dev proxy for /api -> localhost:8000
affects: [06-02-frontend-rewiring, frontend-ui]

# Tech tracking
tech-stack:
  added: []
  patterns:
    - "snake_case TypeScript interfaces matching FastAPI JSON serialization"
    - "EventSource wrapper with auto-cleanup on complete event"
    - "innerHTML SVG rendering (no WASM, no canvas)"
    - "Vite proxy pattern for local backend development"

key-files:
  created:
    - src/api/types.ts
    - src/api/client.ts
    - src/api/sse.ts
  modified:
    - src/ui/molecule-renderer.ts
    - vite.config.ts

key-decisions:
  - "Use snake_case in TypeScript types to match FastAPI JSON exactly (no camelCase conversion)"
  - "SSE wrapper auto-closes on complete event to prevent EventSource leaks"
  - "SVG rendering via innerHTML replaces RDKit.js WASM canvas (backend generates SVG)"
  - "Vite proxy handles /api routing without URL rewriting (backend already serves /api/sa/*)"

patterns-established:
  - "API types mirror backend 1:1 with source file comments for traceability"
  - "Error handling extracts detail/error fields from FastAPI JSON responses"
  - "SSE handlers use optional callbacks (onProgress/onComplete/onWaiting/onError)"
  - "Molecule renderer guards against missing inputs (safe to call with null/undefined)"

# Metrics
duration: 144s (2m 24s)
completed: 2026-02-16
---

# Phase 06 Plan 01: API Client Infrastructure Summary

**TypeScript REST client, SSE wrapper, and SVG renderer establishing typed communication layer between frontend and Python backend**

## Performance

- **Duration:** 2 min 24 sec (144s)
- **Started:** 2026-02-16T11:01:41Z
- **Completed:** 2026-02-16T11:04:05Z
- **Tasks:** 2
- **Files modified:** 5 (3 created, 2 modified)

## Accomplishments
- Created TypeScript types for all backend models (8 interfaces with snake_case fields)
- Built REST API client covering all 5 endpoints (configure, start, pause, reset, getStatus)
- Implemented SSE wrapper with typed event handlers and automatic cleanup
- Replaced RDKit.js WASM canvas rendering with simple SVG innerHTML approach
- Added Vite dev proxy to route /api requests to localhost:8000

## Task Commits

Each task was committed atomically:

1. **Task 1: Create API types, REST client, and SSE wrapper** - `855ec5f` (feat)
   - Created src/api/types.ts with 8 interfaces mirroring backend Pydantic models
   - Created src/api/client.ts with SAAPIClient class for all 5 REST endpoints
   - Created src/api/sse.ts with SSEConnection class for typed SSE handlers
   - All types use snake_case to match FastAPI JSON serialization exactly

2. **Task 2: Rewrite molecule renderer for SVG and add Vite proxy** - `124d697` (feat)
   - Rewrote src/ui/molecule-renderer.ts to use innerHTML for backend SVG
   - Replaced renderMolecule/clearMoleculeCanvas with renderMoleculeSVG/clearMoleculeDisplay
   - Removed all RDKit.js dependencies (no WASM, no canvas operations)
   - Added Vite server proxy configuration for /api -> localhost:8000

## Files Created/Modified

### Created
- `src/api/types.ts` - TypeScript interfaces for SAConfigParams, ConfigureResponse, ControlResponse, StatusResponse, SSEProgressData, SSECompleteData, SSEWaitingData, SSEErrorData
- `src/api/client.ts` - SAAPIClient class with methods for configure, start, pause, reset, getStatus endpoints
- `src/api/sse.ts` - SSEConnection class with connect, close, isConnected methods and typed handlers

### Modified
- `src/ui/molecule-renderer.ts` - Complete rewrite from RDKit.js canvas to SVG innerHTML (88 lines removed, 42 lines added)
- `vite.config.ts` - Added server.proxy configuration for /api endpoint

## Decisions Made

1. **snake_case over camelCase for TypeScript types**
   - FastAPI serializes Pydantic models with snake_case by default
   - No runtime conversion overhead or mapping complexity
   - 1:1 correspondence with backend JSON simplifies debugging

2. **Auto-close SSE connection on complete event**
   - Prevents EventSource leaks when optimization finishes
   - No manual cleanup needed in UI code
   - Connection lifecycle matches SA session lifecycle

3. **SVG innerHTML over RDKit.js canvas**
   - Backend RDKit already generates 2D coords and SVG
   - Eliminates ~2MB WASM payload and async init complexity
   - Browser parses SVG natively with zero overhead
   - Consistent with v2.0 architecture (computation on backend)

4. **Vite proxy without URL rewriting**
   - Backend already serves at /api/sa/*, no rewrite needed
   - changeOrigin handles Host header for CORS
   - Simple one-to-one mapping for development

## Deviations from Plan

None - plan executed exactly as written.

## Issues Encountered

None - all files created and compiled successfully. Import errors in app.ts are expected (Plan 02 will update app.ts to use new renderer functions).

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for Plan 02 (Frontend Rewiring):**
- API client infrastructure complete and type-safe
- SSE wrapper ready for real-time progress updates
- SVG renderer ready to replace canvas operations
- Vite proxy configured for local backend connection
- All TypeScript types match backend contracts exactly

**Expected by Plan 02:**
- Update app.ts to use SAAPIClient instead of Web Worker
- Replace renderMolecule/clearMoleculeCanvas with renderMoleculeSVG/clearMoleculeDisplay
- Wire SSE connection for real-time progress display
- Update HTML to use <div> container instead of <canvas> for molecules

---
*Phase: 06-frontend-integration*
*Completed: 2026-02-16*

## Self-Check: PASSED

All files and commits verified:
- ✓ src/api/types.ts exists
- ✓ src/api/client.ts exists
- ✓ src/api/sse.ts exists
- ✓ Commit 855ec5f exists (Task 1)
- ✓ Commit 124d697 exists (Task 2)
