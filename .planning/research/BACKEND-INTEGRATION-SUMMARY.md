# Backend Integration Research Summary

**Milestone:** v2.0 - FastAPI + Python RDKit Backend with SSE Streaming
**Researched:** 2026-02-15
**Overall Confidence:** HIGH

## Executive Summary

This research addresses the integration of a FastAPI Python backend with RDKit into the existing browser-based WebFaulon SA demo. The v1.0 architecture uses Web Workers + RDKit.js WASM; v2.0 shifts to FastAPI backend + Python RDKit + SSE streaming. This is a well-trodden architectural path with mature technologies and established patterns.

**Core transformation:**
- Web Worker computation → FastAPI backend with async Python
- RDKit.js WASM → Python RDKit (full API, better descriptors)
- Comlink IPC → REST API + Server-Sent Events
- TypeScript MolGraph → RDKit Chem.Mol objects
- Frontend rendering → Backend SVG generation

**Result:** Frontend becomes thin visualization layer, backend handles all computation/chemistry.

The research confirms this migration is **straightforward** with **low risk**:
- FastAPI 0.129.0 is industry standard for async Python APIs (native SSE support)
- Python RDKit 2025.09.5 provides complete cheminformatics toolkit
- sse-starlette 3.2.0 offers production-ready SSE streaming
- EventSource API is native browser feature (no library needed)
- Monorepo structure keeps frontend/backend in sync during development

Multi-component target function framework (new feature) fits naturally into FastAPI dependency injection as a pluggable service.

## Key Findings

### Stack Selection

**Backend:** FastAPI 0.129.0+ (Python 3.10+), RDKit 2025.09.5+, sse-starlette 3.2.0, Pydantic v2, Uvicorn 0.34.0+

**Frontend Changes:** Minimal. Replace Comlink → fetch + EventSource. Remove RDKit.js WASM, Web Worker. Retain Alpine.js, Chart.js, Vite.

**Critical Dependencies:**
- Python RDKit via conda (recommended): `conda install -c conda-forge rdkit`
- Pydantic v2 required (FastAPI 0.128.0+ dropped v1 support)
- Python 3.10+ required (FastAPI 0.129.0 dropped 3.9)

**Why FastAPI:**
- Native async/await (5-50x faster than Flask)
- Automatic OpenAPI docs (/docs endpoint)
- Built-in SSE support via StreamingResponse
- Pydantic integration for type-safe request/response
- Industry standard in 2026

**Why Python RDKit over RDKit.js:**
- Full API (WASM is subset)
- Distance matrix for Wiener Index: `Chem.GetDistanceMatrix(mol)`
- Rich descriptor library (logP, MW, TPSA) for multi-component target
- SVG rendering: `rdMolDraw2D.MolDraw2DSVG`
- Better documentation and community support

### Architecture Patterns

**Three-Layer Backend:**
```
API Router (FastAPI endpoints)
    ↓
Service Layer (orchestration + DI)
    ↓
Core Layer (SAEngine + MoleculeService + TargetFunctionService)
```

**SSE Streaming Pattern:**
1. Client POST `/api/sa/optimize` → receive job_id
2. EventSource opens `/api/sa/stream/{job_id}`
3. Backend runs SA in BackgroundTasks, emits progress to asyncio.Queue
4. SSE endpoint yields queue events: `yield {"data": json.dumps(event)}`
5. Completion event closes stream

**Dependency Injection:**
Services injected via `Depends()`. Example:
```python
@router.post("/optimize")
async def optimize(
    params: SAParams,
    sa_service: SAService = Depends(get_sa_service)
):
    job_id = sa_service.start_optimization(params)
    return {"job_id": job_id}
```

**Multi-Component Target Function:**
Plugin registry pattern. Each component = `Callable[[Mol], float]`. Evaluation = weighted sum. Client sends weights in API request.

### Integration Points

**New Components:**
- FastAPI app (main.py, CORS, startup/shutdown)
- SAService (orchestration, background tasks, progress emission)
- MoleculeService (RDKit wrapper: create, displace, evaluate, render)
- TargetFunctionService (component registry, weighted evaluation)
- Job Store (in-memory Dict + asyncio.Queue)
- SSE endpoint (EventSourceResponse streaming)
- API client (frontend: fetch wrapper)
- SSE client (frontend: EventSource wrapper)

**Modified Components:**
- app.ts: Replace worker.run() → API fetch + EventSource
- chart.ts: Accept SSE events instead of worker callbacks
- molecule-renderer.ts: Display SVG from backend (not canvas)

**Removed Components:**
- sa-worker.ts (replaced by backend)
- MolGraph.ts (replaced by Python MoleculeService)
- SAEngine.ts (replaced by Python SAEngine)
- RDKit.js WASM dependency
- Comlink dependency

### Data Flow

**Request Flow:**
```
Alpine UI → POST /api/sa/optimize → BackgroundTask.add_task(sa_service.run)
         ← {"job_id": "abc123"}

Alpine UI → EventSource /api/sa/stream/abc123
         ← event: progress, data: {"step": 10, "bestEnergy": 45, ...}
         ← event: progress, data: {"step": 20, "bestEnergy": 42, ...}
         ← event: complete, data: {"bestSmiles": "CCC(C)CC", ...}
```

**Backend Execution:**
```
BackgroundTask:
  SAEngine.init()
  for step in steps:
    SAEngine.step()
    if step % 10 == 0:
      await progress_queue.put({"type": "progress", "data": {...}})
  await progress_queue.put({"type": "complete", "data": {...}})

SSE Endpoint:
  async for event in progress_queue:
    yield {"data": json.dumps(event)}
    if event["type"] == "complete":
      break
```

## Implications for Roadmap

### Suggested Phase Structure

**Phase 1: Backend Core (4 hours)**
- FastAPI app skeleton (CORS, health check)
- Pydantic models (SAParams, SAResult, ProgressEvent)
- MoleculeService (RDKit wrapper)
- Python SAEngine (port of TypeScript version)
- Unit tests (pytest + mocked RDKit)

**Phase 2: API + SSE (3 hours)**
- Job Store (in-memory + asyncio.Queue)
- SAService (orchestration, background tasks)
- API routes (POST /optimize, GET /stream, GET /status)
- SSE endpoint with EventSourceResponse
- API tests with AsyncClient

**Phase 3: Frontend Integration (3 hours)**
- API client service (fetch wrapper)
- SSE client service (EventSource wrapper)
- Modify app.ts (replace worker with API)
- Modify chart.ts (SSE events)
- Modify molecule-renderer.ts (SVG display)
- Remove worker code, RDKit.js, Comlink

**Phase 4: Multi-Component Target (2 hours)**
- TargetFunctionService (plugin registry)
- Register Wiener, logP, MW components
- Update SAParams for component weights
- Frontend UI for weight sliders

**Phase 5: Deploy (2 hours)**
- Error handling (invalid formula, job not found, disconnect)
- Loading states (spinner)
- Backend deployment (Render/Railway + Docker + conda)
- Update frontend API_URL to production
- README with v2.0 architecture

**Total: ~14 hours (12 hours for MVP without Phase 4)**

### Phase Ordering Rationale

**Backend First (Phase 1):**
- Independent of frontend
- Can test with curl/httpx
- Establishes API contract

**API Layer Next (Phase 2):**
- Validates SSE streaming independently
- Catches backend issues before frontend integration

**Frontend Integration (Phase 3):**
- Lowest risk (EventSource is simple)
- Delivers end-to-end v2.0

**Multi-Component Additive (Phase 4):**
- Not essential for MVP
- Validates base case first (Wiener only)
- Can defer if time-constrained

### Research Flags for Phases

**Unlikely to need deeper research:**
- Phase 1: RDKit API well-documented
- Phase 2: sse-starlette patterns clear
- Phase 3: EventSource is native browser API
- Phase 4: Simple registry pattern
- Phase 5: Platform docs sufficient (Render/Railway)

**Validation during implementation:**
- RDKit Mol creation from formula (linear chain approach)
- Displacement with RWMol (GetBondBetweenAtoms, AddBond, RemoveBond)
- SanitizeMol edge cases (try/catch, return None for invalid)
- SSE reconnection behavior (test manual disconnect)

## Confidence Assessment

| Area | Confidence | Reason |
|------|------------|--------|
| **Stack** | HIGH | FastAPI, RDKit, SSE are mature (official docs, 2026 examples) |
| **Architecture** | HIGH | Standard three-layer FastAPI pattern, dependency injection best practice |
| **Features** | HIGH | Multi-component is simple plugin registry, RDKit descriptors built-in |
| **Integration** | HIGH | EventSource native, Vite proxy standard, CORS well-documented |
| **Deployment** | MEDIUM | Render/Railway support Python + conda (verified), platform-specific gotchas possible |
| **Pitfalls** | HIGH | Async/sync boundary (RDKit blocking) has known mitigations |

**Overall: HIGH** — Well-trodden path, mature technologies, no novel patterns.

## Roadmap Implications

### Critical Path (MVP)
1. Phase 1 (Backend Core) — blocks all
2. Phase 2 (API + SSE) — blocks frontend
3. Phase 3 (Frontend) — delivers v2.0
4. Phase 5 (Deploy) — classroom-ready

### Optional (Can Defer)
- Phase 4 (Multi-Component) — nice-to-have, Wiener alone demonstrates SA

### Build Order Dependencies
```
Phase 1 → Phase 2 → Phase 3 → Phase 5
                        ↓
                    Phase 4 (optional)
```

### Testing Strategy
- **Phase 1:** Unit (pytest, mock RDKit) + integration (real Mol)
- **Phase 2:** API (AsyncClient) + manual (curl + browser)
- **Phase 3:** E2E (Playwright) + manual UI
- **Phase 4:** Unit (components) + integration (SAEngine)
- **Phase 5:** Smoke (production) + monitoring

### Key Risks & Mitigations

**Risk 1: RDKit blocking event loop**
- **Mitigation:** Run SA step() loop in BackgroundTasks (not individual RDKit calls)

**Risk 2: SSE stream resource leak**
- **Mitigation:** Always emit completion/error event, client closes on completion

**Risk 3: conda deployment complexity**
- **Mitigation:** Use conda-forge Docker images (continuumio/miniconda3)

**Risk 4: CORS misconfiguration**
- **Mitigation:** Test with Vite dev server, never use `*` with credentials

## Open Questions

### Resolved (No Additional Research Needed)
- ✅ How to stream progress? SSE via EventSourceResponse
- ✅ How to structure monorepo? Separate frontend/backend dirs, Vite proxy
- ✅ How to compute Wiener Index? Chem.GetDistanceMatrix → sum
- ✅ How to render molecules? rdMolDraw2D.MolDraw2DSVG → SVG
- ✅ How to handle async RDKit? BackgroundTasks for SA loop

### Implementation Details (Validate Phase 1)
- ❓ Mol creation from formula: Linear chain or random?
  - **Assumption:** Linear (simplest), SA explores space
- ❓ SanitizeMol edge cases: How to detect invalid?
  - **Assumption:** Try/catch, return None for invalid
- ❓ SSE event format: JSON or structured?
  - **Assumption:** JSON in data field (simpler parsing)

### Future Enhancements (Post-v2.0)
- Redis job store (multi-worker)
- Celery task queue (heavy workloads)
- Database persistence (job history)
- Authentication (multi-user)
- Rate limiting (public API)

## Files Created

| File | Purpose |
|------|---------|
| `.planning/research/ARCHITECTURE.md` | System design, components, data flow, patterns (1043 lines) |
| `.planning/research/STACK.md` | Technology versions, installation, alternatives (already existed, backend-focused) |
| `.planning/research/BACKEND-INTEGRATION-SUMMARY.md` (this file) | Executive summary, roadmap implications |

**Note:** FEATURES.md and PITFALLS.md from v1.0 research are still relevant. Multi-component target function covered in ARCHITECTURE.md Pattern 2.

## Summary for Roadmap Creator

**Green light for implementation.** All critical questions answered. Follow phase ordering (Backend → API → Frontend → Deploy). Expect 12-14 hours total.

**Critical success factors:**
1. Test backend independently (Phase 1-2) before frontend
2. Use conda for RDKit (avoid pip binary issues)
3. Run SA in BackgroundTasks (don't block event loop)
4. Send incremental SSE progress (not full history)
5. Close SSE stream on completion/error

**This research provides complete architectural blueprint in ARCHITECTURE.md with:**
- System overview diagram
- Component responsibilities table
- Recommended monorepo structure
- 4 detailed architectural patterns with code examples
- Integration points (new vs modified components)
- Anti-patterns to avoid
- Sources (20+ verified URLs)

Roadmap creator can generate phase plans directly from ARCHITECTURE.md patterns and this SUMMARY.

---
*Backend Integration Research for: WebFaulon v2.0*
*Researched: 2026-02-15*
*Confidence: HIGH*
