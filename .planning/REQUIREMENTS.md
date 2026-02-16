# Requirements: WebFaulon

**Defined:** 2026-02-14
**Core Value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.

## v1.0 Requirements (Complete)

All v1.0 requirements delivered. See MILESTONES.md for archive.

### Input

- [x] **INP-01**: User can enter any molecular formula (e.g., C6H14, C8H10) to define the constitutional isomer space
- [x] **INP-02**: Formula is validated before SA starts (invalid formulas rejected with clear error message)
- [x] **INP-03**: Preset example molecules (3-5 curated) available for quick start

### Algorithm

- [x] **ALG-01**: SA implements Faulon's displacement: pick 4 atoms, redistribute bond orders per equations 7-11
- [x] **ALG-02**: Initial structure generated deterministically from the molecular formula
- [x] **ALG-03**: Wiener Index computed as cost function
- [x] **ALG-04**: User can toggle between maximizing and minimizing the Wiener Index
- [x] **ALG-05**: SA accepts/rejects moves using Metropolis criterion
- [x] **ALG-06**: Graph connectivity validated after every SA move
- [x] **ALG-07**: Algorithm runs entirely in-browser via Web Worker

### Controls

- [x] **CTRL-01**: User can set initial temperature (kT)
- [x] **CTRL-02**: User can select cooling schedule from f0 through f32 family
- [x] **CTRL-03**: User can set number of steps per cycle
- [x] **CTRL-04**: User can set number of cycles
- [x] **CTRL-05**: User can play, pause, and reset the SA execution

### Visualization

- [x] **VIZ-01**: Live-updating chart shows Wiener Index vs step number
- [x] **VIZ-02**: Current best isomer rendered as 2D chemical structure
- [x] **VIZ-03**: Best Wiener Index score displayed alongside the structure
- [x] **VIZ-04**: Real-time algorithm state display (step, temperature, current/best score)

### User Experience

- [x] **UX-01**: Clean design suitable for classroom projection
- [x] **UX-02**: Mobile responsive layout

## v2.0 Requirements

Requirements for Python backend migration. Each maps to roadmap phases.

### Backend Foundation (BACK)

- [ ] **BACK-01**: FastAPI application serves REST endpoints with automatic OpenAPI docs
- [ ] **BACK-02**: CORS middleware allows frontend origin (GitHub Pages + localhost dev)
- [ ] **BACK-03**: Health check endpoint returns server status
- [ ] **BACK-04**: API returns structured JSON error responses with HTTP codes (400/404/409/500)

### RDKit Molecular Operations (MOL)

- [ ] **MOL-01**: Backend validates molecular formulas using RDKit before accepting SA configuration
- [ ] **MOL-02**: Faulon displacement (eqs 7-11) operates on RDKit RWMol with SanitizeMol validation
- [ ] **MOL-03**: Backend generates 2D SVG molecule depictions via RDKit MolDraw2DSVG
- [ ] **MOL-04**: Backend produces canonical SMILES via RDKit MolToSmiles

### SA Engine (SA)

- [ ] **SA-01**: Python SAEngine implements Metropolis criterion with configurable cooling schedules
- [ ] **SA-02**: Initial molecular structure generated deterministically from formula using RDKit
- [ ] **SA-03**: SA execution runs in background (non-blocking) — server remains responsive during optimization
- [ ] **SA-04**: Session state persisted in-memory with configurable TTL

### API & Streaming (API)

- [ ] **API-01**: POST endpoint accepts SA configuration (formula, parameters, component weights) and returns session ID
- [ ] **API-02**: SSE endpoint streams live SA progress events (step, temperature, energy, best molecule SVG)
- [ ] **API-03**: Control endpoints for start, pause, and reset SA execution
- [ ] **API-04**: Status endpoint returns current SA state for client reconnection

### Target Function (TGT)

- [ ] **TGT-01**: Multi-component target function framework with pluggable scoring components
- [ ] **TGT-02**: Wiener Index implemented as first scoring component
- [ ] **TGT-03**: Each component contributes weighted score to overall objective

### Frontend Migration (FE)

- [ ] **FE-01**: Frontend communicates with backend via REST API + EventSource (no Web Worker)
- [ ] **FE-02**: Live chart updates from SSE events maintain v1 UX responsiveness
- [ ] **FE-03**: Molecule rendering displays backend-generated SVG (no RDKit.js WASM)
- [ ] **FE-04**: Same UX preserved: formula input, presets, parameter controls, start/pause/reset

### Deployment (DEP)

- [ ] **DEP-01**: Backend deployable to cloud platform (Railway/Render/Fly.io)
- [ ] **DEP-02**: Frontend configurable for production backend URL
- [ ] **DEP-03**: Production CORS locked to deployment domain

## v2.1+ Requirements

Deferred to future release. Tracked but not in current roadmap.

### Enhanced Multi-Component

- **TGT-04**: Adjustable component weights configurable via API
- **TGT-05**: Component-wise score streaming in SSE events (breakdown per component)
- **TGT-06**: RDKit descriptor calculation API exposing 200+ molecular descriptors

### Operations

- **API-05**: Graceful cancellation endpoint for long-running SA
- **SA-05**: Concurrent session support via async task status

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| Spectroscopic target functions (NMR prediction) | Future milestone — requires NMR prediction engine |
| Celery/Redis task queue | YAGNI — BackgroundTasks sufficient for educational demo |
| Database persistence | In-memory dict sufficient — sessions are ephemeral classroom use |
| Authentication/authorization | Students in classroom, not customers with accounts |
| WebSocket communication | SSE is correct for one-way streaming |
| GraphQL | REST with clear endpoints is simpler for educational codebase |
| Docker containerization at launch | Deploy directly to Railway/Render first |
| 3D molecular visualization | Wiener Index is about graph topology, not 3D geometry |
| Batch processing | Students learn by watching individual runs |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

### v1.0 (Complete)

| Requirement | Phase | Status |
|-------------|-------|--------|
| INP-01 | Phase 2 | Done |
| INP-02 | Phase 2 | Done |
| INP-03 | Phase 2 | Done |
| ALG-01 | Phase 1 | Done |
| ALG-02 | Phase 1 | Done |
| ALG-03 | Phase 1 | Done |
| ALG-04 | Phase 1 | Done |
| ALG-05 | Phase 1 | Done |
| ALG-06 | Phase 1 | Done |
| ALG-07 | Phase 2 | Done |
| CTRL-01 | Phase 2 | Done |
| CTRL-02 | Phase 2 | Done |
| CTRL-03 | Phase 2 | Done |
| CTRL-04 | Phase 2 | Done |
| CTRL-05 | Phase 2 | Done |
| VIZ-01 | Phase 3 | Done |
| VIZ-02 | Phase 3 | Done |
| VIZ-03 | Phase 3 | Done |
| VIZ-04 | Phase 3 | Done |
| UX-01 | Phase 3 | Done |
| UX-02 | Phase 3 | Done |

### v2.0 (In Progress)

| Requirement | Phase | Status |
|-------------|-------|--------|
| BACK-01 | Phase 4 | Pending |
| BACK-02 | Phase 4 | Pending |
| BACK-03 | Phase 4 | Pending |
| BACK-04 | Phase 4 | Pending |
| MOL-01 | Phase 4 | Pending |
| MOL-02 | Phase 4 | Pending |
| MOL-03 | Phase 5 | Pending |
| MOL-04 | Phase 4 | Pending |
| SA-01 | Phase 4 | Pending |
| SA-02 | Phase 4 | Pending |
| SA-03 | Phase 5 | Pending |
| SA-04 | Phase 5 | Pending |
| API-01 | Phase 5 | Pending |
| API-02 | Phase 5 | Pending |
| API-03 | Phase 5 | Pending |
| API-04 | Phase 5 | Pending |
| TGT-01 | Phase 7 | Pending |
| TGT-02 | Phase 4 | Pending |
| TGT-03 | Phase 7 | Pending |
| FE-01 | Phase 6 | Pending |
| FE-02 | Phase 6 | Pending |
| FE-03 | Phase 6 | Pending |
| FE-04 | Phase 6 | Pending |
| DEP-01 | Phase 8 | Pending |
| DEP-02 | Phase 8 | Pending |
| DEP-03 | Phase 8 | Pending |

**Coverage:**
- v1.0 requirements: 21 total -- all Done
- v2.0 requirements: 26 total
- Mapped to phases: 26
- Unmapped: 0

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-15 after v2.0 roadmap created*
