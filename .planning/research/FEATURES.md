# Feature Landscape

**Domain:** FastAPI Backend Migration with RDKit for Molecular SA (Subsequent Milestone)
**Researched:** 2026-02-15
**Confidence:** MEDIUM-HIGH

## Context

This research covers ONLY the NEW features for v2.0 backend migration. The existing browser-based demo (v1.x) already provides:
- Formula input with validation and preset molecules
- SA parameter controls (kT, cooling schedule, steps, cycles)
- Start/pause/reset execution controls
- Live Wiener Index vs step chart
- 2D molecule rendering with canonical SMILES display
- Responsive classroom-ready UI
- GitHub Pages deployment

**Migration goal:** Same UX experience, but with Python RDKit backend enabling multi-component target functions beyond Wiener Index (future spectroscopic scoring).

## Table Stakes (Users Expect These)

Features users assume exist. Missing these = migration feels broken.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| REST API for SA configuration | Standard backend pattern; frontend needs to configure algorithm parameters | Low | POST /api/sa/configure with JSON body for formula, kT, cooling, steps, cycles. Return session ID. |
| SSE streaming of SA progress | Required to maintain "live chart updates" UX from v1.x; polling is 2010s pattern | Medium | GET /api/sa/stream/{session_id} returns text/event-stream with step updates. Frontend EventSource consumes this. |
| Start/pause/reset via API | Frontend controls depend on it; users expect same control as v1.x | Low | POST /api/sa/{session_id}/start, /pause, /reset endpoints. State machine validates transitions. |
| Current state query endpoint | Frontend needs to reconstruct state on page refresh or reconnection | Low | GET /api/sa/{session_id}/status returns current step, energy, temperature, best molecule. |
| Backend molecule validation | RDKit must validate formula and generated structures; v1.x does this client-side with RDKit.js | Medium | Use rdkit.Chem.MolFromSmiles() and rdkit.Chem.Descriptors for validation. Return 400 with clear error messages for invalid formulas. |
| 2D SVG molecule rendering | v1.x renders molecules; backend must provide SVG to replace RDKit.js WASM rendering | Medium | Use rdkit.Chem.Draw.MolDraw2DSVG to generate SVG. Return via API or embed in SSE events. |
| Session state persistence | Users expect to pause and resume without losing work | Medium | Store SAEngine state in-memory (dict) initially. Session expires after configurable timeout (30 min default). |
| CORS configuration | Frontend (GitHub Pages) and backend (separate origin) require CORS | Low | FastAPI CORSMiddleware with allowed_origins from env var. Critical for deployment. |
| Error handling with HTTP codes | REST API standard; frontend needs clear error semantics | Low | 400 (bad input), 404 (session not found), 409 (invalid state transition), 500 (internal error). JSON error responses. |
| Health check endpoint | Deployment standard; container orchestration and monitoring expect this | Low | GET /api/health returns 200 with {status: "ok", version: "2.0.0"}. |

## Differentiators (Competitive Advantage)

Features that set v2.0 apart. Not in v1.x, not required for migration, but enable future value.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Multi-component target function framework | MAIN VALUE of migration; enables spectroscopic scoring beyond Wiener Index | High | Pluggable component architecture. Each component (WienerIndexComponent, NMRPredictionComponent) implements score(mol) -> float. Weighted sum combines them. |
| Adjustable component weights via API | Users can tune "60% Wiener, 40% NMR prediction" to explore tradeoffs | Medium | POST /api/sa/configure accepts weights dict: {wiener: 0.6, nmr: 0.4}. Validates sum = 1.0. Enables multi-objective optimization experiments. |
| Component-wise score streaming | SSE events include breakdown: {total: 42.3, wiener: 28.1, nmr: 14.2} so users see contribution of each | Low-Medium | Extends SSE event format. Minimal complexity once framework exists. High educational value (shows which components dominate). |
| RDKit descriptor calculation API | Expose RDKit's 200+ descriptors (MolWt, LogP, TPSA) for custom scoring components | Low | GET /api/descriptors/{smiles} returns JSON with all descriptors. Enables users to explore molecular properties. |
| Backend displacement on native RDKit | Python RDKit Chem.RWMol allows native bond manipulation; v1.x uses custom graph | High | Replicate Faulon displacement (eqs 7-11) using RDKit's SetBondBondOrder(), AddBond(), RemoveBond(). Validate with Chem.SanitizeMol(). More chemically rigorous than adjacency matrix. |
| Canonical SMILES via RDKit | RDKit's canonicalization is authoritative; v1.x uses custom DFS algorithm | Low | Use Chem.MolToSmiles(mol, canonical=True). Replaces frontend SMILES generation. More reliable. |
| Async SA execution with task status | Non-blocking server; multiple users can run SA simultaneously without worker contention | Medium | Use FastAPI BackgroundTasks initially. Each session runs SA in background thread. Status updates via shared state. |
| Graceful cancellation | Users can stop long-running SA without orphaning resources | Medium | POST /api/sa/{session_id}/cancel sets stop flag. SA loop checks flag each iteration. Clean resource cleanup. |

## Anti-Features (Commonly Requested, Often Problematic)

Features that seem good but create problems for this migration.

| Feature | Why Requested | Why Problematic | Alternative |
|---------|---------------|-----------------|-------------|
| Celery/Redis for task queue | "Production scalability" mindset | Overengineering for educational demo. v1.x runs client-side; comparable scale is single-server. Adds operational complexity (Redis deployment, worker monitoring). | Use FastAPI BackgroundTasks. If scaling needed later, add Celery then. YAGNI principle applies. |
| WebSocket bidirectional communication | "More modern than SSE" | Adds complexity (connection management, reconnection logic, message framing) for zero benefit. SA is server→client streaming only. SSE has built-in reconnection. | SSE is perfect fit. WebSockets solve different problem (bidirectional real-time). |
| Database persistence for sessions | "Sessions should survive restarts" | Educational demo users don't expect multi-day sessions. Adds DB schema, migrations, queries for minimal value. | In-memory dict with TTL. Document that sessions expire after 30 min. Export results if important. |
| Authentication/authorization | "Secure the API" | Users are students in classroom, not customers with accounts. Auth adds login flow, token management, user DB. v1.x is public; v2.0 should be too. | Rate limiting if abuse occurs. Focus on education, not security theater. |
| GraphQL instead of REST | "More flexible queries" | SA workflow is simple CRUD + SSE. GraphQL adds complexity (schema, resolvers, client libraries) for zero query flexibility benefit. | REST with clear endpoints. Simple is better for educational codebase students might read. |
| Batch job API | "Run multiple SA configurations" | Changes UX from "watch algorithm work" to "submit job, check later". Defeats educational goal of observing SA decisions. | Single focused run per session. Manual rerun with different params teaches exploration. |
| Docker containerization | "Modern deployment" | GitHub Pages frontend calls backend → backend needs public URL. Docker adds layer but doesn't solve hosting. Can run FastAPI on Heroku/Railway/Fly.io directly. | Deploy FastAPI directly. Add Dockerfile later if container orchestration needed. |
| Comprehensive logging framework | "Production observability" | Over-instrumentation for educational demo. Logs should debug dev issues, not monitor production SLA. | Python logging with INFO level. Structlog if complexity grows. Don't prematurely optimize for observability. |

## Feature Dependencies

```
REST API for SA configuration
    └──requires──> Backend molecule validation (validate formula before accepting)
    └──enables──> Session state persistence (config stored per session)
    └──enables──> Start/pause/reset via API (session must exist)

Start/pause/reset via API
    └──requires──> Async SA execution (non-blocking execution)
    └──requires──> Session state persistence (state between control commands)
    └──enables──> SSE streaming (must have running SA to stream)

SSE streaming of SA progress
    └──requires──> Async SA execution (background execution while streaming)
    └──enhances──> Component-wise score streaming (extends event format)
    └──maintains──> Live chart updates UX from v1.x (critical migration goal)

Multi-component target function framework
    └──requires──> Backend displacement on native RDKit (need RDKit mol objects, not adjacency matrix)
    └──enables──> Adjustable component weights via API (framework must exist to weight)
    └──enables──> Component-wise score streaming (breakdown requires components)
    └──future──> NMR prediction scoring (not in scope, but architecture enables)

Backend displacement on native RDKit
    └──requires──> RDKit Python installation
    └──conflicts──> Custom graph adjacency matrix approach (migration replaces frontend graph with RDKit mol)
    └──enables──> Canonical SMILES via RDKit (RDKit mol is source of truth)
    └──enables──> 2D SVG molecule rendering (RDKit Draw APIs)

2D SVG molecule rendering
    └──requires──> Backend displacement on native RDKit (RDKit mol object to render)
    └──replaces──> RDKit.js WASM rendering in v1.x (frontend consumes SVG from API)

CORS configuration
    └──enables──> Frontend-backend communication (GitHub Pages → FastAPI)
    └──required for──> ALL API endpoints (no CORS = API unusable from frontend)
```

### Dependency Notes

- **REST API requires validation:** Cannot accept invalid formulas; backend must validate before storing session
- **SSE requires async execution:** Streaming events while SA runs requires non-blocking execution model
- **Multi-component framework requires RDKit mols:** Adjacency matrix cannot support RDKit descriptor APIs or future spectroscopic scoring
- **Backend displacement conflicts with frontend graph:** Migration means RDKit mol is source of truth; cannot mix with adjacency matrix approach
- **CORS is foundational:** Without CORS, frontend cannot call backend; must configure before any testing
- **Component framework enables future:** NMR prediction not in v2.0, but architecture must support pluggable components
- **Session persistence scope:** In-memory dict sufficient for educational use; DB persistence is anti-feature unless multi-day sessions required

## MVP Definition

### Launch With (v2.0 Backend Migration)

Minimum viable migration — feature parity with v1.x plus multi-component architecture foundation.

**Feature Parity with v1.x (Table Stakes):**
- [ ] **REST API for SA configuration** — POST /api/sa/configure with formula, SA params, component weights
- [ ] **Backend molecule validation** — RDKit validates formula and structures
- [ ] **SSE streaming of SA progress** — Maintains live chart updates UX
- [ ] **Start/pause/reset via API** — POST /api/sa/{session_id}/{start|pause|reset}
- [ ] **Current state query endpoint** — GET /api/sa/{session_id}/status
- [ ] **2D SVG molecule rendering** — Backend replaces RDKit.js rendering
- [ ] **Session state persistence** — In-memory dict with 30-min TTL
- [ ] **CORS configuration** — Allow GitHub Pages origin
- [ ] **Error handling with HTTP codes** — 400/404/409/500 with JSON errors
- [ ] **Health check endpoint** — GET /api/health

**New Capabilities (Migration Value):**
- [ ] **Multi-component target function framework** — Pluggable components with weighted sum
- [ ] **Backend displacement on native RDKit** — Faulon displacement on RDKit.Chem.RWMol
- [ ] **Canonical SMILES via RDKit** — Authoritative canonicalization

### Add After Migration Validation (v2.1)

Features to add once backend migration is stable and frontend integrated.

- [ ] **Adjustable component weights via API** — Trigger: Users want to explore component tradeoffs
- [ ] **Component-wise score streaming** — Trigger: Educational value of seeing score breakdown
- [ ] **RDKit descriptor calculation API** — Trigger: Users want to explore molecular properties
- [ ] **Graceful cancellation** — Trigger: Long-running SA (>10 min) needs early termination
- [ ] **Async SA execution with task status** — Trigger: Multiple users or researchers want concurrent sessions

### Future Consideration (v2.x+)

Features to defer until multi-component framework is validated with real usage.

- [ ] **NMR prediction scoring component** — Why defer: Requires NMRShiftDB or ML model integration; high complexity; validates framework value first
- [ ] **Celery/Redis task queue** — Why defer: Overengineering unless concurrent usage proves BackgroundTasks insufficient
- [ ] **Database persistence for sessions** — Why defer: In-memory sufficient for educational use; add if multi-day research sessions emerge
- [ ] **Custom scoring component upload API** — Why defer: Research tool territory; framework must prove valuable first
- [ ] **Batch configuration comparison** — Why defer: Single-run focus for education; batch is research workflow
- [ ] **Docker containerization** — Why defer: Deploy directly to Heroku/Railway first; containerize if orchestration needed
- [ ] **Advanced rate limiting** — Why defer: Add if abuse occurs; don't prematurely optimize

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| REST API for SA configuration | HIGH | LOW | P1 |
| SSE streaming of SA progress | HIGH | MEDIUM | P1 |
| Start/pause/reset via API | HIGH | LOW | P1 |
| Backend molecule validation | HIGH | MEDIUM | P1 |
| 2D SVG molecule rendering | HIGH | MEDIUM | P1 |
| Multi-component framework | HIGH | HIGH | P1 |
| Backend RDKit displacement | HIGH | HIGH | P1 |
| Canonical SMILES via RDKit | MEDIUM | LOW | P1 |
| Session state persistence | MEDIUM | MEDIUM | P1 |
| CORS configuration | HIGH | LOW | P1 |
| Error handling | MEDIUM | LOW | P1 |
| Health check endpoint | LOW | LOW | P1 |
| Adjustable component weights | HIGH | MEDIUM | P2 |
| Component-wise score streaming | MEDIUM | LOW | P2 |
| RDKit descriptor API | MEDIUM | LOW | P2 |
| Graceful cancellation | MEDIUM | MEDIUM | P2 |
| Async execution with status | MEDIUM | MEDIUM | P2 |
| NMR prediction component | HIGH | VERY HIGH | P3 |
| Celery/Redis | LOW | HIGH | P3 |
| Database persistence | LOW | HIGH | P3 |
| Custom component upload | MEDIUM | VERY HIGH | P3 |
| Batch comparison | MEDIUM | HIGH | P3 |

**Priority key:**
- P1: Must have for v2.0 migration launch (feature parity + multi-component foundation)
- P2: Should have, add after migration validated (enhanced multi-component exploration)
- P3: Nice to have, future consideration (research-grade features)

**Prioritization rationale:**
- **P1 features** ensure migration doesn't regress UX (SSE streaming maintains live updates) while establishing multi-component architecture
- **P2 features** enhance multi-component exploration (component weights, score breakdown) after foundation validated
- **P3 features** expand into research territory (NMR prediction, custom components) after educational value proven

## Comparison: v1.x Browser-Based vs v2.0 Backend

| Feature | v1.x (RDKit.js WASM) | v2.0 (FastAPI + RDKit) | Migration Approach |
|---------|----------------------|------------------------|-------------------|
| Molecule validation | RDKit.js client-side | RDKit Python backend | Same library, different platform — high confidence |
| 2D rendering | RDKit.js canvas | Backend SVG via API | Replace WASM rendering with SVG endpoint |
| SMILES canonicalization | Custom DFS algorithm | RDKit Chem.MolToSmiles | Use authoritative RDKit implementation |
| SA progress updates | Web Worker postMessage | SSE streaming | SSE maintains real-time UX without client threading |
| Displacement operation | Custom adjacency matrix | RDKit RWMol bond manipulation | Higher fidelity; use RDKit's bond APIs |
| Target function | Wiener Index only | Multi-component weighted sum | Core migration value; pluggable architecture |
| State management | Web Worker in-memory | Backend session dict | Server-side state enables multi-component persistence |
| Deployment | GitHub Pages static | FastAPI server + GitHub Pages frontend | Frontend stays on Pages; backend needs hosting |
| Cost | Free (static hosting) | $5-10/mo (Railway/Heroku) | Consider serverless (Vercel/Cloudflare) if cost-sensitive |

**Key migration risks:**
- **SSE streaming latency:** Web Worker has zero network latency; SSE adds RTT. Mitigation: Event batching (10 updates/sec not 1000/sec).
- **CORS complexity:** Local dev needs CORS proxy or backend CORS config. Mitigation: vite proxy config for dev, env-based CORS for prod.
- **RDKit displacement fidelity:** Adjacency matrix is simple; RDKit RWMol has sanitization failures. Mitigation: Extensive test suite replicating v1.x displacement edge cases.
- **Session state explosion:** In-memory dict grows unbounded if no cleanup. Mitigation: TTL-based expiration (30 min) and max session limit.

## "Same UX" Definition for Migration

Migration succeeds if users CANNOT TELL the difference except for URL change. Specifically:

**Must maintain:**
1. Live chart updates (SSE must feel as responsive as Web Worker postMessage)
2. Instant start/pause/reset (API roundtrip <100ms on reasonable connection)
3. Molecule rendering quality (SVG must match RDKit.js canvas quality)
4. Progress update frequency (10-60 updates/sec depending on SA step rate)
5. Error messages (formula validation errors must be equally clear)

**Can change (acceptable differences):**
1. Loading time (backend SSE connection setup may add 200-500ms initial delay)
2. SMILES format (RDKit Python may canonicalize differently than custom DFS — this is improvement)
3. Molecule rendering style (SVG vs canvas different but equivalent)

**Must NOT regress:**
1. Mobile responsiveness (SVG must render well on mobile)
2. Offline usage (v1.x works offline; v2.0 requires backend connection — document this)
3. Control responsiveness (pause must stop within 1 SA iteration, not delayed by network)

## Domain-Specific Patterns for Backend Migration

### FastAPI + SSE for Algorithm Streaming (2026 Pattern)

**Standard approach:**
- SSE endpoint returns `text/event-stream` with keep-alive
- FastAPI `StreamingResponse` with async generator yields events
- Frontend `EventSource` consumes stream with auto-reconnect
- Event format: `data: {json}\n\n` per SSE spec

**Performance considerations:**
- Batch events (accumulate 10-100ms worth, then yield)
- Use `asyncio.sleep(0)` to yield control in tight loop
- Close SSE stream on session end (prevents resource leak)

**Error handling:**
- SSE errors sent as `event: error` messages
- Frontend EventSource.onerror triggers reconnection
- 409 Conflict if session not in RUNNING state

### Multi-Component Scoring Architecture (Plugin Pattern)

**Standard approach (scikit-learn inspired):**
- Base class `ScoringComponent` with `score(mol: Chem.Mol) -> float`
- Concrete components: `WienerIndexComponent`, `NMRPredictionComponent`
- Registry pattern: `ComponentRegistry.register("wiener", WienerIndexComponent)`
- Weighted sum: `total = sum(w * component.score(mol) for component, w in zip(components, weights))`

**Configuration pattern:**
```python
{
  "components": [
    {"name": "wiener", "weight": 0.6},
    {"name": "nmr", "weight": 0.4}
  ]
}
```

**Validation:**
- Weights sum to 1.0 (normalize or reject)
- Component names exist in registry
- Each component validates mol compatibility

### RDKit Molecule State Management

**Pattern:**
- Store `Chem.Mol` as SMILES string in session (serializable)
- Reconstruct `Chem.Mol` from SMILES when needed (stateless)
- Never pickle `Chem.Mol` (fragile, version-dependent)
- Use `Chem.MolToSmiles()` and `Chem.MolFromSmiles()` for serialization

**Displacement pattern:**
```python
mol = Chem.RWMol(current_mol)  # Make editable copy
mol.GetBondBetweenAtoms(i, j).SetBondType(Chem.BondType.DOUBLE)  # Modify
Chem.SanitizeMol(mol)  # Validate chemistry
return mol.GetMol()  # Convert back to ROMol
```

**Sanitization:**
- ALWAYS call `Chem.SanitizeMol()` after bond changes
- Catch `Chem.KekulizeException` for aromatic failures
- Reject displacement if sanitization fails (invalid chemistry)

### FastAPI State Machine for SA Sessions

**State transitions:**
```
IDLE → CONFIGURED (POST /configure)
CONFIGURED → RUNNING (POST /start)
RUNNING → PAUSED (POST /pause)
PAUSED → RUNNING (POST /start)
RUNNING → STOPPED (POST /reset or completion)
STOPPED → CONFIGURED (POST /configure)
any → ERROR (on exception)
```

**Validation:**
- 409 Conflict if invalid transition (e.g., pause when IDLE)
- Include allowed_actions in status response for client UX
- State stored in session dict: `sessions[id]['state']`

## Educational Context Insights for Backend Migration

**From research findings:**

1. **SSE maintains real-time learning** (FastAPI SSE research): Server-sent events preserve "watch algorithm work" pedagogy without WebSocket complexity
   - **Implication:** SSE is architecturally correct for one-way streaming; maintain 10-60 updates/sec from v1.x

2. **Multi-component framework teaches tradeoffs** (Multi-objective optimization research): Students learn by adjusting weights between competing objectives
   - **Implication:** Component weight API is educationally valuable; priority P2 after base framework validated

3. **Backend doesn't diminish transparency** (Computational chemistry education): Students benefit from reading Python RDKit code as much as JavaScript
   - **Implication:** Backend migration doesn't harm "students read source" value; Python may be more readable

4. **Deployment complexity is real** (Backend migration considerations): v1.x GitHub Pages = $0/mo; v2.0 needs FastAPI hosting
   - **Implication:** Document hosting options (Railway, Heroku, Fly.io); consider serverless (Vercel, Cloudflare Workers)

5. **Network latency changes UX** (SSE streaming latency): Web Worker postMessage ~0ms; SSE ~10-100ms RTT
   - **Implication:** Event batching (10 updates/sec) mitigates latency while maintaining "real-time" feel

## Sources

**FastAPI SSE Streaming:**
- [Implementing Server-Sent Events (SSE) with FastAPI: Real-Time Updates Made Simple | Medium](https://mahdijafaridev.medium.com/implementing-server-sent-events-sse-with-fastapi-real-time-updates-made-simple-6492f8bfc154)
- [FastAPI and SSE: How to Build Streamable MCP Servers | Aubergine](https://www.aubergine.co/insights/a-guide-to-building-streamable-mcp-servers-with-fastapi-and-sse)
- [How to use Server-Sent Events with FastAPI and React | Softgrade](https://www.softgrade.org/sse-with-fastapi-react-langgraph/)
- [Real-Time Notifications in Python: Using SSE with FastAPI | Medium](https://medium.com/@inandelibas/real-time-notifications-in-python-using-sse-with-fastapi-1c8c54746eb7)
- [FastAPI + SSE for LLM Tokens: Smooth Streaming without WebSockets | Medium](https://medium.com/@hadiyolworld007/fastapi-sse-for-llm-tokens-smooth-streaming-without-websockets-001ead4b5e53)

**RDKit Python API:**
- [Getting Started with the RDKit in Python — RDKit 2025.09.5 documentation](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [The RDKit Documentation — RDKit 2025.09.5](https://rdkit.org/docs/)
- [Manipulation of molecules with RDKit in Python for AI models | Medium](https://zoehlerbz.medium.com/manipulation-of-molecules-with-rdkit-in-python-for-ai-models-8023f1e677c7)
- [RDKit Python Tutorial for Chemists (With Examples) – Runcell Blog](https://www.runcell.dev/blog/rdkit-chemists-tutorial)

**FastAPI Architecture:**
- [Modern FastAPI Architecture Patterns for Scalable Production Systems | Medium](https://medium.com/algomart/modern-fastapi-architecture-patterns-for-scalable-production-systems-41a87b165a8b)
- [Building a Real-time Dashboard with FastAPI and Svelte | TestDriven.io](https://testdriven.io/blog/fastapi-svelte/)
- [Data Visualization using FastAPI and EasyCharts | Medium](https://medium.com/analytics-vidhya/data-visualization-using-fastapi-and-easycharts-493eda3b1f3d)

**Multi-Component Scoring:**
- [GitHub - laurentg/pyscoring: A generic scoring algorithm](https://github.com/laurentg/pyscoring)
- [How to Develop a Weighted Average Ensemble With Python | MachineLearningMastery](https://machinelearningmastery.com/weighted-average-ensemble-with-python/)
- [Weighted Scoring Model: Step-by-Step Implementation Guide | Product School](https://productschool.com/blog/product-fundamentals/weighted-scoring-model)

**FastAPI Background Tasks vs Celery:**
- [Celery and Background Tasks. Using FastAPI with long running tasks | Medium](https://medium.com/@hitorunajp/celery-and-background-tasks-aebb234cae5d)
- [Fastapi background tasks vs Celery | Medium](https://haseeb987.medium.com/fastapi-background-tasks-vs-celery-56d2394bf30d)
- [How to Implement Background Tasks in FastAPI | OneUptime](https://oneuptime.com/blog/post/2026-02-02-fastapi-background-tasks/view)
- [Background Tasks in FastAPI | FastAPI Tutorial](https://www.fastapitutorial.com/blog/fastapi-background-tasks/)
- [The Complete Guide to Background Processing with FastAPI × Celery/Redis | Greeden Blog](https://blog.greeden.me/en/2026/01/27/the-complete-guide-to-background-processing-with-fastapi-x-celery-redishow-to-separate-heavy-work-from-your-api-to-keep-services-stable/)

**REST API State Management:**
- [Designing a True REST State Machine | Nordic APIs](https://nordicapis.com/designing-a-true-rest-state-machine/)
- [How to design state machines for microservices | Red Hat Developer](https://developers.redhat.com/articles/2021/11/23/how-design-state-machines-microservices)
- [Designing Complex REST API Structure With State Machine | Medium](https://medium.com/geekculture/designing-complex-rest-api-structure-with-state-machine-cc2b8804d467)

**RDKit SVG Rendering:**
- [rdkit.Chem.Draw package — RDKit 2025.09.5 documentation](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html)
- [Visualization and Rendering | rdkit/rdkit | DeepWiki](https://deepwiki.com/rdkit/rdkit/4-2d3d-functionality)
- [Transforming a RDKit depiction function into a Web API with Flask @ hectormartinez.dev](https://hectormartinez.dev/posts/rdkit-depiction-web-api-flask/)

**FastAPI CORS:**
- [CORS (Cross-Origin Resource Sharing) - FastAPI](https://fastapi.tiangolo.com/tutorial/cors/)
- [Blocked by CORS in FastAPI? Here's How to Fix It | David Muraya](https://davidmuraya.com/blog/fastapi-cors-configuration/)
- [Configuring CORS in FastAPI | GeeksforGeeks](https://www.geeksforgeeks.org/python/configuring-cors-in-fastapi/)
- [Demystifying CORS in FastAPI & React: A Practical Guide | Medium](https://vinaysit.wordpress.com/2024/11/07/demystifying-cors-in-fastapi-react-a-practical-guide-%F0%9F%8C%90%F0%9F%9A%80/)

**NMR Prediction (Future Spectroscopic Scoring):**
- [NMR-Solver: Automated Structure Elucidation via Large-Scale Spectral Matching | arXiv](https://arxiv.org/html/2509.00640v1)
- [Molecular search by NMR spectrum | Scientific Reports](https://www.nature.com/articles/s41598-021-00488-z)
- [Rapid prediction of NMR spectral properties with quantified uncertainty | Journal of Cheminformatics](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0374-3)
- [Accurate Prediction of 1H NMR Chemical Shifts Using Machine Learning | PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11123270/)

**FastAPI State Management:**
- [Database Session Management | fastapi-practices | DeepWiki](https://deepwiki.com/fastapi-practices/fastapi_best_architecture/7.6-database-session-management)
- [Getting Started - FastAPI Sessions](https://jordanisaacs.github.io/fastapi-sessions/guide/getting_started/)
- [Mastering Session Management in FastAPI | Medium](https://medium.com/@rameshkannanyt0078/mastering-session-management-in-fastapi-a-human-friendly-guide-2b59701afa2e)

**RDKit Molecule Validation:**
- [rdkit.Chem.MolStandardize.rdMolStandardize module — RDKit 2025.09.4](https://www.rdkit.org/docs/source/rdkit.Chem.MolStandardize.rdMolStandardize.html)
- [MolVS: Molecule Validation and Standardization — MolVS 0.1.1](https://molvs.readthedocs.io/en/latest/)

---

**Research confidence: MEDIUM-HIGH**

**High confidence areas:**
- FastAPI SSE patterns (current 2026 sources, official docs)
- RDKit Python API (official documentation, version 2025.09.5)
- CORS configuration (official FastAPI docs)
- Background tasks vs Celery decision framework (multiple consistent sources)

**Medium confidence areas:**
- Multi-component scoring architecture patterns (inferred from scikit-learn patterns, needs domain validation)
- Migration UX equivalence (SSE latency vs Web Worker postMessage not directly compared in sources)
- Session state persistence approach (in-memory dict vs database tradeoff context-dependent)

**Low confidence areas:**
- NMR prediction integration complexity (future feature, limited integration examples found)
- Optimal SSE event batching rate (needs performance testing with actual SA step rates)
- Faulon displacement on RDKit RWMol fidelity (no existing implementation found, needs validation)

**Validation recommendations:**
- Prototype SSE streaming with simulated SA to measure latency vs Web Worker baseline
- Implement Faulon displacement on RDKit RWMol and validate against v1.x adjacency matrix test suite
- Survey chemistry instructors about multi-component scoring educational value
- Benchmark session state TTL and max session limits with realistic classroom usage (30 students)

---
*Feature research for: WebFaulon v2.0 Backend Migration*
*Researched: 2026-02-15*
*Researcher: Claude Opus 4.6*
