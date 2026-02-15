# Project Research Summary

**Project:** WebFaulon v2.0
**Domain:** Python Backend Migration for Molecular SA Optimization (FastAPI + RDKit + SSE)
**Researched:** 2026-02-15
**Confidence:** HIGH

## Executive Summary

WebFaulon v2.0 migrates from a browser-only architecture (Web Workers + RDKit.js WASM) to a client-server model with a FastAPI Python backend, native RDKit for cheminformatics, and Server-Sent Events for real-time progress streaming. This is a well-understood architectural transition. FastAPI 0.129.0 is the industry standard for async Python APIs in 2026, Python RDKit 2025.09.5 provides the authoritative cheminformatics toolkit, and SSE via sse-starlette 3.2.0 is the correct pattern for one-way streaming (replacing Web Worker postMessage). The frontend becomes a thin visualization layer: Alpine.js for UI state, Chart.js for charting, and SVG rendering from the backend replaces RDKit.js WASM canvas rendering. Every technology in this stack is mature, well-documented, and battle-tested.

The primary value of this migration is enabling multi-component target functions beyond Wiener Index. Python RDKit exposes 200+ molecular descriptors (logP, MW, TPSA) and opens the door to future spectroscopic scoring (NMR prediction). The pluggable component architecture -- a weighted-sum registry where each component implements `score(mol) -> float` -- fits naturally into FastAPI's dependency injection system. This architectural capability is the reason to migrate; without it, the browser-only version is sufficient.

The key risks are: (1) CPU-bound SA iterations blocking the async event loop -- mitigated by running SA in BackgroundTasks or ProcessPoolExecutor from day one, (2) RDKit sanitization failures after Faulon displacement -- mitigated by mandatory `Chem.SanitizeMol()` after every bond modification with try/catch rejection, and (3) SSE resource leaks from unclosed streams or client disconnects -- mitigated by bounded queues, adaptive reporting intervals, and disconnect detection. All three have well-documented solutions. The migration path is low-risk and the research is conclusive: proceed to implementation.

## Key Findings

### Recommended Stack

**Backend (NEW for v2.0):**
- **FastAPI 0.129.0+**: Async Python API framework with native SSE support via StreamingResponse. Requires Python 3.10+ (dropped 3.9 in 0.129.0). 5-50x faster than Flask for async workloads. Automatic OpenAPI docs at `/docs`.
- **Python RDKit 2025.09.5+**: Full cheminformatics toolkit. `Chem.GetDistanceMatrix()` for Wiener Index, `rdMolDraw2D.MolDraw2DSVG` for SVG rendering, `Chem.RWMol` for bond manipulation. Install via conda (recommended) or pip.
- **sse-starlette 3.2.0**: Production-ready SSE following W3C spec. `EventSourceResponse` for streaming progress. Released Jan 2026. Native FastAPI integration.
- **Pydantic v2**: Rust-based validation, required by FastAPI 0.128.0+ (v1 removed). Runtime type checking for request/response schemas.
- **Uvicorn 0.34.0+**: ASGI server. Single process for dev, Gunicorn + Uvicorn workers for production.
- **Poetry**: Dependency management with lockfiles. De facto Python standard in 2026.
- **pytest + pytest-asyncio**: Testing framework. `AsyncClient` for SSE endpoint tests.

**Frontend (RETAINED from v1.0, modified role):**
- **Vite 6.x + TypeScript 5.x**: Build tool and type-safe development. Unchanged.
- **Alpine.js 3.x**: UI state management. Modified to use fetch + EventSource instead of Comlink.
- **Chart.js 4.x**: SA progress charting. Now consumes SSE events instead of worker callbacks.
- **EventSource API**: Native browser SSE client. No library needed.

**Removed from frontend:** RDKit.js WASM, Comlink, Web Workers, sa-worker.ts, MolGraph.ts, SAEngine.ts.

**Critical version constraints:**
- Python >= 3.10 (FastAPI requirement)
- Pydantic v2 only (v1 removed from FastAPI 0.128.0)
- RDKit via conda for binary dependency management

### Expected Features

**Must have (table stakes for migration parity):**
- REST API for SA configuration (POST /api/sa/configure with formula, params, weights)
- SSE streaming of SA progress (maintains live chart UX from v1.x)
- Start/pause/reset via API (POST /api/sa/{session_id}/{action})
- Backend molecule validation (RDKit validates formulas and structures)
- 2D SVG molecule rendering (backend replaces RDKit.js WASM rendering)
- Session state persistence (in-memory dict with 30-min TTL)
- CORS configuration (GitHub Pages frontend calls separate backend origin)
- Health check endpoint (GET /api/health)
- Error handling with proper HTTP codes (400/404/409/500)

**Should have (migration value-adds):**
- Multi-component target function framework (pluggable components with weighted sum -- the core reason to migrate)
- Backend displacement on native RDKit RWMol (higher fidelity than adjacency matrix)
- Canonical SMILES via RDKit (authoritative, replaces custom DFS algorithm)

**Defer to v2.1+:**
- Adjustable component weights API (after framework validated)
- Component-wise score streaming in SSE events
- RDKit descriptor calculation API (expose 200+ descriptors)
- Graceful cancellation with cancel endpoint
- Concurrent session support via async task status

**Anti-features (explicitly avoid):**
- Celery/Redis task queue (overengineering for educational demo, use BackgroundTasks)
- WebSocket communication (SSE is correct for one-way streaming)
- Database persistence for sessions (in-memory dict sufficient, sessions are ephemeral)
- Authentication/authorization (students in classroom, not customers)
- GraphQL (REST with clear endpoints is simpler and more educational)
- Docker containerization at launch (deploy directly to Railway/Heroku/Fly.io)

### Architecture Approach

Three-layer backend with clear separation: API Router (FastAPI endpoints) -> Service Layer (orchestration + dependency injection) -> Core Layer (SAEngine + MoleculeService + TargetFunctionService). The frontend communicates via REST (POST to start, GET for status) + SSE (EventSource for real-time progress). SA runs in BackgroundTasks, emitting progress to asyncio.Queue, which the SSE endpoint drains as an async generator. This is the standard FastAPI pattern for long-running computation with streaming results.

**Major backend components:**
1. **FastAPI Router** -- HTTP endpoints for SA optimization (POST /optimize, GET /stream/{job_id}, GET /status/{job_id})
2. **SAService** -- Orchestrates SA execution in background tasks, manages job lifecycle and progress emission
3. **MoleculeService** -- RDKit wrapper: create molecules from formula, apply Faulon displacement, compute Wiener Index, generate SVG, convert to SMILES
4. **TargetFunctionService** -- Plugin registry for scoring components. Each component is `Callable[[Mol], float]`. Evaluation is weighted sum.
5. **SAEngine** -- Python port of TypeScript SAEngine. Stateful class with `step()` method. Uses MoleculeService for all RDKit operations.
6. **Job Store** -- In-memory Dict + asyncio.Queue per job. Tracks state, routes progress events to SSE endpoint.

**Key architectural patterns:**
- SSE for unidirectional progress streaming (EventSourceResponse + async generator)
- Dependency injection via FastAPI `Depends()` for testable services
- BackgroundTasks for non-blocking SA execution (not the event loop)
- Pydantic models for type-safe request/response + OpenAPI schema generation
- Separate POST (create job) from GET (stream progress) to prevent EventSource reconnect duplication

**Monorepo structure:** `frontend/` (Vite) and `backend/` (FastAPI) as sibling directories. Vite proxy to backend in dev mode. No shared code -- use OpenAPI schema as source of truth for types.

### Critical Pitfalls

Top 10 pitfalls identified, mapped to prevention phases:

1. **Missing RDKit sanitization after bond changes** -- RDKit does NOT auto-validate after `SetBondType()`. Must call `Chem.SanitizeMol()` after every displacement. Failure produces silently invalid molecules. **Prevention:** Try/catch sanitization, reject move on failure. Unit tests assert SMILES correctness after displacement. [Phase 1]

2. **CPU-bound SA blocking async event loop** -- SA loop (RDKit operations, Wiener computation) is CPU-bound. Running in `async def` blocks the entire server. **Prevention:** Run SA in BackgroundTasks or ProcessPoolExecutor from day one. Verify `/health` responds during optimization. [Phase 1]

3. **SSE backpressure causing memory exhaustion** -- SA generates updates faster than client consumes. Unbounded queue grows until OOM. **Prevention:** Adaptive reporting interval (max 10 events/sec), bounded asyncio.Queue (maxsize=10), `asyncio.sleep(0)` to yield control. [Phase 2]

4. **Client disconnect not canceling SA computation** -- Browser tab closes, SSE drops, but SA loop continues for minutes wasting CPU. **Prevention:** Check `request.is_disconnected()` every 100 steps. Clean up on CancelledError. [Phase 2]

5. **CORS middleware ordering breaks error responses** -- CORSMiddleware added after other middleware means exceptions bypass CORS headers. Browser shows generic CORS error instead of actual 422/500. **Prevention:** Add CORSMiddleware FIRST, before all other middleware. Use Vite proxy in dev to avoid CORS entirely. [Phase 1]

6. **RDKit Mol memory leaks in long SA runs** -- Storing Mol objects in history creates unbounded memory growth. C++ objects have complex lifetimes not immediately freed by Python GC. **Prevention:** Store SMILES strings in history, not Mol objects. Explicit `del rw_mol` after displacement. Periodic `gc.collect()`. [Phase 1]

7. **Incorrect BondType enum mapping** -- RDKit BondType enum (UNSPECIFIED=0, SINGLE=1, DOUBLE=2, TRIPLE=3, AROMATIC=12) does not match Faulon integer bond orders. **Prevention:** Explicit mapping function `faulon_order_to_rdkit_bondtype()`. Test all bond types. [Phase 1]

8. **EventSource auto-reconnect duplicating jobs** -- If SSE endpoint both starts optimization AND streams, reconnect triggers duplicate jobs. **Prevention:** Separate POST (create job, return job_id) from GET (stream progress). EventSource only on GET endpoint. [Phase 2]

9. **Frontend assumes synchronous responses** -- Web Worker pattern (`await worker.run()`) does not translate to HTTP. Single `fetch()` for entire optimization times out. **Prevention:** Build `SAClient` adapter class that mimics worker API: POST to start, EventSource for progress, Promise resolves on completion event. [Phase 3]

10. **Plugin architecture prevents pickling for ProcessPoolExecutor** -- Class-based plugins with closures/DB connections can't pickle. **Prevention:** Registry pattern with top-level functions + serializable config objects. Or use ThreadPoolExecutor (acceptable for demo). [Phase 4]

## Implications for Roadmap

### Suggested Phase Structure

#### Phase 1: Backend Core + RDKit Foundation
**Rationale:** Independent of frontend. Establishes API contract and validates RDKit integration. All downstream phases depend on this. Backend-first allows testing with curl before any frontend work.
**Delivers:** FastAPI app skeleton, CORS config, health check, Pydantic models (SAParams, SAResult, ProgressEvent), MoleculeService (RDKit wrapper), Python SAEngine (port of TypeScript), unit tests (pytest + real RDKit Mol objects).
**Features addressed:** Backend molecule validation, canonical SMILES, backend displacement on RDKit RWMol, health check, error handling, CORS.
**Pitfalls prevented:** #1 (sanitization), #2 (event loop blocking), #5 (CORS ordering), #6 (memory leaks), #7 (BondType mapping).
**Research needed:** None. RDKit API well-documented. FastAPI patterns standard.

#### Phase 2: API Layer + SSE Streaming
**Rationale:** Validates streaming independently before frontend integration. Job store + SSE is the critical integration point between computation and visualization.
**Delivers:** Job Store (in-memory + asyncio.Queue), SAService (orchestration, BackgroundTasks, progress emission), REST endpoints (POST /optimize, GET /stream/{job_id}, GET /status/{job_id}), SSE endpoint with EventSourceResponse, API tests with AsyncClient.
**Features addressed:** REST API for SA configuration, SSE streaming, start/pause/reset, session state persistence, current state query.
**Pitfalls prevented:** #3 (SSE backpressure), #4 (disconnect detection), #8 (reconnect duplication).
**Research needed:** None. sse-starlette patterns clear from official docs and 2026 examples.

#### Phase 3: Frontend Integration
**Rationale:** Lowest risk phase -- EventSource is a native browser API. Delivers end-to-end v2.0 by wiring frontend to backend. Removes all Web Worker / RDKit.js / Comlink code.
**Delivers:** API client service (fetch wrapper), SSE client service (EventSource wrapper with SAClient adapter), modified app.ts (replace worker with API), modified chart.ts (SSE events), modified molecule-renderer.ts (display SVG from backend). Removal of sa-worker.ts, MolGraph.ts, SAEngine.ts, RDKit.js, Comlink.
**Features addressed:** All table-stakes migration parity features (live chart, instant controls, molecule rendering, progress updates).
**Pitfalls prevented:** #9 (frontend sync assumption via SAClient adapter).
**Research needed:** None. EventSource is native browser API.

#### Phase 4: Multi-Component Target Function
**Rationale:** Core migration value. Not required for MVP parity but is the architectural reason for migrating. Can defer if time-constrained -- Wiener-only demonstrates SA correctly.
**Delivers:** TargetFunctionService (plugin registry), registered components (Wiener, logP, MW), updated SAParams for component weights, frontend UI for weight configuration.
**Features addressed:** Multi-component framework, adjustable weights, component-wise scoring.
**Pitfalls prevented:** #10 (plugin pickling via registry pattern).
**Research needed:** Minimal. Simple registry pattern. May want to validate weight normalization UX.

#### Phase 5: Deployment + Polish
**Rationale:** After functionality validated, deploy backend for classroom use. Error handling, loading states, and production configuration.
**Delivers:** Comprehensive error handling (invalid formula, job not found, disconnect recovery), loading states and spinners, backend deployment (Render/Railway + Docker + conda-forge base image), frontend API_URL configuration for production, updated README with v2.0 architecture.
**Features addressed:** Deployment readiness, production CORS, graceful degradation.
**Pitfalls prevented:** Security mistakes (CORS lockdown, rate limiting, formula validation).
**Research needed:** Minimal. Platform-specific deployment docs (Render/Railway) may have gotchas.

### Phase Ordering Rationale

- **Backend before frontend:** Backend can be tested independently with curl/httpx/AsyncClient. Frontend integration is the lowest-risk phase (EventSource is trivial). This ordering catches RDKit and async issues before they compound with frontend state management.
- **SSE before frontend:** Validating streaming independently reveals backpressure and disconnect issues without the noise of UI debugging. The BACKEND-INTEGRATION-SUMMARY confirms this ordering.
- **Multi-component after MVP:** Wiener-only demonstrates SA correctly. The component framework adds value but is not required for migration parity. Deferring it de-risks the critical path.
- **Deployment last:** Classroom use requires a deployed backend. But deployment complexity (conda in Docker, CORS for production domain) is isolated from functional correctness.

**Critical path:** Phase 1 -> Phase 2 -> Phase 3 -> Phase 5
**Optional branch:** Phase 3 -> Phase 4 (can run before or after Phase 5)

### Research Flags

**Phases with standard patterns (skip `/gsd:research-phase`):**
- **Phase 1:** RDKit Python API is well-documented. FastAPI project setup is standard. SAEngine port is mechanical translation.
- **Phase 2:** sse-starlette has clear examples. Job store with asyncio.Queue is a common pattern. Background tasks documented in FastAPI official docs.
- **Phase 3:** EventSource is a native browser API. Fetch wrapper is trivial. SVG display is just setting `innerHTML`.
- **Phase 5:** Railway/Render/Fly.io all have FastAPI deployment guides. Docker + conda-forge is documented.

**Phases that may benefit from validation during implementation (not full research):**
- **Phase 1:** Faulon displacement on RDKit RWMol has no existing implementation found. Validate against v1.x adjacency matrix test suite once ported.
- **Phase 4:** Multi-component scoring weight normalization UX. Simple registry pattern but the API design for weights (must sum to 1.0? auto-normalize?) needs a decision.

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | **HIGH** | FastAPI 0.129.0, RDKit 2025.09.5, sse-starlette 3.2.0 all verified via official docs, PyPI, and 2026 production guides. Version constraints confirmed. |
| Features | **HIGH** | Feature parity requirements clear from v1.x. Multi-component framework is simple plugin registry. Anti-features well-justified (YAGNI for educational demo). |
| Architecture | **HIGH** | Three-layer FastAPI pattern is industry standard. SSE + BackgroundTasks + asyncio.Queue documented in multiple 2026 sources. Monorepo structure follows FastAPI community conventions. |
| Pitfalls | **HIGH** | All 10 pitfalls sourced from official RDKit/FastAPI docs, GitHub issues, and 2024-2026 community discussions. Mitigations verified. Recovery costs assessed. |

**Overall confidence: HIGH** -- This is a well-trodden migration path with mature technologies and no novel patterns.

### Gaps to Address

Gaps are minor and can be resolved during implementation, not through additional research.

- **Faulon displacement on RDKit RWMol:** No existing Python implementation found. Port from TypeScript adjacency matrix approach and validate against v1.x test suite. If RWMol bond manipulation proves unreliable for Faulon's 4-atom displacement pattern, fall back to adjacency matrix with RDKit only for rendering/descriptors.
- **Optimal SSE event batching rate:** Research suggests 10 events/sec but actual latency depends on network conditions and SA step rate. Needs performance testing with realistic SA runs (500 steps/cycle x 4 cycles).
- **Conda in deployment:** Render/Railway support Python but conda adds complexity to Docker images. Use `continuumio/miniconda3` base image. If problematic, try `pip install rdkit` (pre-built wheels exist for Linux).
- **Offline usage regression:** v1.x works offline (static site). v2.0 requires backend connection. Document this clearly. Not a blocker but a UX change to acknowledge.

## Sources

### Primary (HIGH confidence)

**FastAPI + SSE:**
- [FastAPI Release Notes](https://fastapi.tiangolo.com/release-notes/) -- Version 0.129.0, Python 3.10+ requirement
- [FastAPI CORS](https://fastapi.tiangolo.com/tutorial/cors/) -- CORSMiddleware configuration
- [FastAPI Async Tests](https://fastapi.tiangolo.com/advanced/async-tests/) -- AsyncClient testing patterns
- [sse-starlette PyPI](https://pypi.org/project/sse-starlette/) -- Version 3.2.0, SSE implementation
- [FastAPI Concurrency](https://fastapi.tiangolo.com/async/) -- Async/await patterns

**Python RDKit:**
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html) -- Python API
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html) -- Wiener Index, molecular manipulation
- [RDKit Draw Module](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html) -- SVG rendering
- [RDKit Installation](https://www.rdkit.org/docs/Install.html) -- Conda vs pip

**Pydantic:**
- [Pydantic v2 Docs](https://docs.pydantic.dev/latest/) -- Rust-based validation, migration from v1

### Secondary (MEDIUM-HIGH confidence)

**Architecture and Patterns:**
- [FastAPI Best Practices Production 2026](https://fastlaunchapi.dev/blog/fastapi-best-practices-production-2026)
- [Modern FastAPI Architecture Patterns](https://medium.com/algomart/modern-fastapi-architecture-patterns-for-scalable-production-systems-41a87b165a8b)
- [The Complete FastAPI x pytest Guide (Feb 2026)](https://blog.greeden.me/en/2026/01/06/the-complete-fastapi-x-pytest-guide/)
- [Mastering Dependency Injection in FastAPI](https://medium.com/@azizmarzouki/mastering-dependency-injection-in-fastapi-clean-scalable-and-testable-apis-5f78099c3362)

**SSE Streaming:**
- [Implementing SSE with FastAPI](https://mahdijafaridev.medium.com/implementing-server-sent-events-sse-with-fastapi-real-time-updates-made-simple-6492f8bfc154)
- [Stop Burning CPU on Dead FastAPI Streams](https://jasoncameron.dev/posts/fastapi-cancel-on-disconnect)
- [Using Server-Sent Events -- MDN](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events)

**RDKit Pitfalls:**
- [RDKit Sanitization and File Parsing (2025)](https://greglandrum.github.io/rdkit-blog/posts/2025-06-27-sanitization-and-file-parsing.html)
- [Memory Leakage with Molecule Objects -- GitHub #3239](https://github.com/rdkit/rdkit/issues/3239)
- [Explicit Valence Error -- GitHub Discussion #8181](https://github.com/rdkit/rdkit/discussions/8181)

### Tertiary (MEDIUM confidence, validate during implementation)

- Multi-component scoring patterns inferred from scikit-learn plugin architecture
- SSE latency vs Web Worker postMessage not directly benchmarked in sources
- NMR prediction integration complexity (future, limited integration examples)
- Faulon displacement on RDKit RWMol (no existing implementation, needs validation)

---
*Research completed: 2026-02-15*
*Ready for roadmap: YES*
