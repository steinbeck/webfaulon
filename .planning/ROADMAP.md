# Roadmap: WebFaulon

## Milestones

- v1.0 Interactive SA Demo (Browser-Only) -- Phases 1-3.1 (shipped 2026-02-15)
- v2.0 Python Backend Architecture -- Phases 4-8 (in progress)

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

<details>
<summary>v1.0 Interactive SA Demo (Phases 1-3.1) -- SHIPPED 2026-02-15</summary>

- [x] **Phase 1: Molecular Graph & SA Core** - Chemical correctness foundation
- [x] **Phase 2: Browser Integration & Controls** - Web Workers, user input, execution controls
- [x] **Phase 3: Visualization & UX** - Charts, structure rendering, classroom-ready design
- [x] **Phase 3.1: README & GitHub Pages** - Documentation and public deployment (INSERTED)

### Phase 1: Molecular Graph & SA Core
**Goal**: Core algorithm produces chemically valid molecular structures through simulated annealing
**Depends on**: Nothing (first phase)
**Requirements**: ALG-01, ALG-02, ALG-03, ALG-04, ALG-05, ALG-06
**Plans:** 4 plans

Plans:
- [x] 01-01-PLAN.md -- Project scaffold, MolGraph data structure, and Wiener Index
- [x] 01-02-PLAN.md -- Initial structure generation from molecular formula (TDD)
- [x] 01-03-PLAN.md -- Faulon displacement equations 7-11 and seeded PRNG (TDD)
- [x] 01-04-PLAN.md -- SA engine with Metropolis criterion and cooling schedules (TDD)

### Phase 2: Browser Integration & Controls
**Goal**: Users can configure and execute SA algorithm in browser without UI freezing
**Depends on**: Phase 1
**Requirements**: ALG-07, CTRL-01, CTRL-02, CTRL-03, CTRL-04, CTRL-05, INP-01, INP-02, INP-03
**Plans:** 3 plans

Plans:
- [x] 02-01-PLAN.md -- SAEngine step-by-step execution API refactor (TDD)
- [x] 02-02-PLAN.md -- Web Worker + Comlink integration and Vite project scaffold
- [x] 02-03-PLAN.md -- Alpine.js UI controls, formula validation, presets, and RDKit.js WASM

### Phase 3: Visualization & UX
**Goal**: Students see optimization happening in real time with clear, classroom-ready visuals
**Depends on**: Phase 2
**Requirements**: VIZ-01, VIZ-02, VIZ-03, VIZ-04, UX-01, UX-02
**Plans:** 3 plans

Plans:
- [x] 03-01-PLAN.md -- SMILES generation for MolGraph + worker bestSMILES propagation (TDD)
- [x] 03-02-PLAN.md -- Chart.js live chart module + RDKit.js molecule renderer module
- [x] 03-03-PLAN.md -- UI integration, responsive redesign, and visual verification checkpoint

### Phase 3.1: README & GitHub Pages (INSERTED)
**Goal**: App is publicly accessible at https://steinbeck.github.io/webfaulon/ with comprehensive documentation
**Depends on**: Phase 3
**Plans:** 2 plans

Plans:
- [x] 03.1-01-PLAN.md -- Vite base path config, RDKit CDN fix, and GitHub Actions deployment workflow
- [x] 03.1-02-PLAN.md -- Comprehensive README with project docs, demo link, and academic citation

</details>

### v2.0 Python Backend Architecture (Phases 4-8)

**Milestone Goal:** Re-architect from browser-only to FastAPI backend with Python RDKit, enabling full cheminformatics capabilities and a multi-component target function framework.

- [x] **Phase 4: Backend Core & RDKit Foundation** - FastAPI app with Python SA engine operating on native RDKit molecules
- [x] **Phase 5: API Layer & SSE Streaming** - REST endpoints, session management, and live progress streaming
- [x] **Phase 6: Frontend Integration** - Rewire frontend from Web Worker to backend API + SSE
- [ ] **Phase 7: Multi-Component Target Function** - Pluggable scoring framework with weighted components
- [ ] **Phase 8: Deployment & Production** - Cloud deployment with production configuration

## Phase Details

### Phase 4: Backend Core & RDKit Foundation
**Goal**: Python backend performs all molecular operations (validation, displacement, scoring, SMILES) on native RDKit molecules with a working SA engine
**Depends on**: Phase 3.1 (v1.0 complete)
**Requirements**: BACK-01, BACK-02, BACK-03, BACK-04, MOL-01, MOL-02, MOL-04, SA-01, SA-02, TGT-02
**Success Criteria** (what must be TRUE):
  1. FastAPI server starts, `/docs` shows auto-generated OpenAPI documentation, and `/health` returns server status
  2. Submitting an invalid molecular formula to the API returns a structured JSON error with HTTP 400
  3. Python SA engine runs a complete optimization (500 steps x 4 cycles) on C6H14, producing valid SMILES at every step
  4. Faulon displacement on RDKit RWMol preserves valence rules and molecular connectivity (validated by SanitizeMol after every move)
  5. Wiener Index computed via RDKit distance matrix matches v1.0 values for the same molecules
**Plans**: 4 plans

Plans:
- [x] 04-01-PLAN.md -- Poetry project scaffold, FastAPI app with CORS, health endpoint, and error handling
- [x] 04-02-PLAN.md -- Port SeededRandom, formula parser, and cooling schedule with TDD
- [x] 04-03-PLAN.md -- RDKit MoleculeGraph wrapper, Wiener Index, displacement, and initial structure with TDD
- [x] 04-04-PLAN.md -- Python SAEngine with Metropolis criterion and full integration tests (TDD)

### Phase 5: API Layer & SSE Streaming
**Goal**: Backend exposes a complete REST + SSE API for configuring, controlling, and streaming SA optimization in real time
**Depends on**: Phase 4
**Requirements**: API-01, API-02, API-03, API-04, SA-03, SA-04, MOL-03
**Success Criteria** (what must be TRUE):
  1. POST to configure endpoint returns a session ID; GET to stream endpoint delivers live SSE events with step number, temperature, energy, and best molecule SVG
  2. Start, pause, and reset commands via API control a running SA session (pause actually halts iteration, reset clears state)
  3. Server remains responsive to `/health` requests while an SA optimization is running in the background
  4. Reconnecting to the status endpoint after a disconnect returns the current SA state (step, best score, best molecule)
  5. Idle sessions are automatically cleaned up after TTL expiration
**Plans**: 3 plans

Plans:
- [x] 05-01-PLAN.md -- SessionManager + SVG renderer services with TDD
- [x] 05-02-PLAN.md -- REST endpoints (configure, control, status) + main.py wiring + API tests
- [x] 05-03-PLAN.md -- SSE stream endpoint + end-to-end integration tests

### Phase 6: Frontend Integration
**Goal**: Frontend communicates entirely through backend API and SSE, delivering the same UX as v1.0 without Web Workers or RDKit.js WASM
**Depends on**: Phase 5
**Requirements**: FE-01, FE-02, FE-03, FE-04
**Success Criteria** (what must be TRUE):
  1. User can enter a formula, select a preset, adjust SA parameters, and click Start -- the full v1.0 workflow works end-to-end via backend
  2. Live chart updates from SSE events with the same responsiveness as v1.0 Web Worker updates
  3. Best molecule renders as SVG received from the backend (no RDKit.js WASM loaded in the browser)
  4. Start, pause, and reset buttons work correctly during an active SA run
**Plans**: 2 plans

Plans:
- [x] 06-01-PLAN.md -- API client infrastructure (types, REST client, SSE wrapper) + SVG molecule renderer + Vite proxy
- [x] 06-02-PLAN.md -- Rewire app.ts from Web Worker to backend API + SSE, update index.html, human verification

### Phase 06.1: Fix presets and displacement efficiency for unsaturated molecule demo (INSERTED)

**Goal:** [Urgent work - to be planned]
**Depends on:** Phase 6
**Plans:** 0 plans

Plans:
- [ ] TBD (run /gsd:plan-phase 06.1 to break down)

### Phase 7: Multi-Component Target Function
**Goal**: SA optimization supports multiple pluggable scoring components with user-configurable weights
**Depends on**: Phase 5 (API contract for weights), Phase 6 (end-to-end verification)
**Requirements**: TGT-01, TGT-03
**Success Criteria** (what must be TRUE):
  1. SA engine evaluates a weighted sum of multiple scoring components (Wiener Index as default, framework accepts additional components)
  2. Each component contributes an independently computed score that is combined via configurable weights into the overall objective
**Plans**: 2 plans

Plans:
- [ ] 07-01-PLAN.md -- Scoring component framework: Protocol, WienerIndexComponent, MolecularWeightComponent, ComponentRegistry (TDD)
- [ ] 07-02-PLAN.md -- SA engine integration: extend SAParams with component_weights, refactor SAEngine to use weighted sum (TDD)

### Phase 8: Deployment & Production
**Goal**: Backend is deployed to a cloud platform and frontend is configured to use the production API
**Depends on**: Phase 6 (end-to-end functional)
**Requirements**: DEP-01, DEP-02, DEP-03
**Success Criteria** (what must be TRUE):
  1. Backend is running on a cloud platform (Railway/Render/Fly.io) and responds to health check requests from the internet
  2. Frontend at the production URL (GitHub Pages) successfully connects to the deployed backend and runs an SA optimization
  3. CORS is locked to the production domain -- requests from unauthorized origins are rejected
**Plans**: TBD

Plans:
- [ ] 08-01: TBD

## Progress

**Execution Order:**
Phases execute in numeric order: 4 -> 5 -> 6 -> 7 -> 8

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 1. Molecular Graph & SA Core | v1.0 | 4/4 | Complete | 2026-02-14 |
| 2. Browser Integration & Controls | v1.0 | 3/3 | Complete | 2026-02-15 |
| 3. Visualization & UX | v1.0 | 3/3 | Complete | 2026-02-15 |
| 3.1. README & GitHub Pages | v1.0 | 2/2 | Complete | 2026-02-15 |
| 4. Backend Core & RDKit Foundation | v2.0 | 4/4 | Complete | 2026-02-16 |
| 5. API Layer & SSE Streaming | v2.0 | 3/3 | Complete | 2026-02-16 |
| 6. Frontend Integration | v2.0 | 2/2 | Complete | 2026-02-16 |
| 7. Multi-Component Target Function | v2.0 | 0/? | Not started | - |
| 8. Deployment & Production | v2.0 | 0/? | Not started | - |

---

*Roadmap created: 2026-02-14*
*v2.0 phases added: 2026-02-15*
*Phase 4 planned: 2026-02-16*
*Phase 4 executed: 2026-02-16*
*Phase 5 planned: 2026-02-16*
*Phase 5 executed: 2026-02-16*
*Phase 6 planned: 2026-02-16*
*Phase 6 executed: 2026-02-16*
