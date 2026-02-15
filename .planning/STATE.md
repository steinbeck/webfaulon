# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.
**Current focus:** Phase 3 - Visualization & UX Polish (Phases 1-2 complete)

## Current Position

Phase: 3 of 3 (Visualization & UX)
Plan: 2 of 3 in current phase
Status: In Progress
Last activity: 2026-02-15 — Completed plan 03-02 (Chart & Molecule Renderer Modules)

Progress: [███████░░░] 67% (Phase 3: 2/3 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 8
- Average duration: 16.3 min
- Total execution time: 2.17 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-molecular-graph-sa-core | 4 | 28 min | 7.0 min |
| 02-browser-integration-controls | 3 | 96 min | 32.0 min |
| 03-visualization-ux | 1 | 2 min | 2.0 min |

| Phase | Plan | Duration | Tasks | Files |
|-------|------|----------|-------|-------|
| 01-molecular-graph-sa-core | 01 | 4 min | 3 | 8 |
| 01-molecular-graph-sa-core | 02 | 4 min | 3 | 4 |
| 01-molecular-graph-sa-core | 03 | 8 min | 3 | 4 |
| 01-molecular-graph-sa-core | 04 | 10 min | 3 | 4 |
| 02-browser-integration-controls | 01 | 3 min | 2 | 3 |
| 02-browser-integration-controls | 02 | 3 min | 3 | 10 |
| 02-browser-integration-controls | 03 | 90 min | 4 | 10 |
| 03-visualization-ux | 02 | 2 min | 2 | 4 |

## Accumulated Context

### Decisions

Decisions are logged in PROJECT.md Key Decisions table.
Recent decisions affecting current work:

- Initialization: RDKit.js chosen for 2D rendering (full cheminformatics in browser)
- Initialization: In-browser only architecture (no backend, simplest classroom deployment)
- Initialization: Wiener Index as sole cost function for v1 (matches paper's primary test case)
- Initialization: Modern SPA with Vite (good DX, tree-shaking, easy dependencies)
- [Phase 01-01]: Mutable MolGraph with explicit validation (SA displacement pattern)
- [Phase 01-01]: Array-based BFS queue for O(1) dequeue in Wiener computation
- [Phase 01-02]: Linear chain with iterative bond upgrades for unsaturation (deterministic, guarantees valid structures)
- [Phase 01-02]: Regex-based formula parser with known element validation (simple, robust, extensible)
- [Phase 01-03]: Mulberry32 PRNG for seeded randomness (fast, deterministic, no dependencies)
- [Phase 01-03]: Return null for invalid displacements rather than retry internally (SA engine controls retry logic)
- [Phase 01-04]: Class-based SAEngine with private state management (clean API, encapsulates SA state)
- [Phase 01-04]: Extract isBetter() helper method in refactor phase (eliminates duplication, single source of truth)
- [Phase 02-01]: Refactored run() to delegate to init()+step()+getResult() for single code path
- [Phase 02-01]: getState() returns snapshot (not live reference) for safe cross-thread sharing
- [Phase 02-01]: Error throwing on invalid method call sequences (step before init, getResult before complete)
- [Phase 02-02]: Singleton SAWorker instance exposed via Comlink (simpler lifecycle than class exposure)
- [Phase 02-02]: Moved SAParams and SAResult to shared types.ts for cross-module type safety
- [Phase 02-02]: Periodic yield every 100 steps ensures pause/resume messages can be processed
- [Phase 02-03]: Two-stage formula validation (regex format + chemical validity via parseFormula + HDI)
- [Phase 02-03]: 5 preset molecules selected to demonstrate different SA behaviors (size and unsaturation variety)
- [Phase 02-03]: Alpine.js for reactive UI (lightweight, suitable for classroom projection)
- [Phase 02-03]: RDKit.js loaded via CDN script tag (library exposes global, not ES module)
- [Phase 02-03]: Module-level worker references outside Alpine reactive scope (prevents Proxy-of-Proxy conflicts with Comlink)
- [Phase 03-02]: Chart.js tree-shaken imports (not auto bundle) for optimal bundle size
- [Phase 03-02]: Module-level chart instance storage outside Alpine reactive scope prevents Proxy conflicts
- [Phase 03-02]: Decimation plugin with min-max algorithm preserves SA peaks/valleys
- [Phase 03-02]: Animation disabled and parsing:false for real-time performance
- [Phase 03-02]: aspectRatio 2.5 for wide step timeline visualization
- [Phase 03-02]: Redundant SMILES re-render prevention via _lastRenderedSMILES tracking
- [Phase 03-02]: WASM memory cleanup guaranteed via try/finally with mol.delete()

### Pending Todos

None yet.

### Blockers/Concerns

**Phase 1 planning:**
- ~~Research gap: Faulon displacement equations (eqs 7-11) details needed from original 1996 paper~~ RESOLVED: Equations verified from Faulon 1996 paper page 733, implemented in 01-03
- ~~Decision needed: Initial structure generation algorithm (deterministic vs stochastic approach)~~ RESOLVED: Linear chain approach implemented in 01-02

**Phase 3 planning:**
- Chart decimation threshold needs performance testing (500 vs 1000 vs 2000 sample points)
- Mobile breakpoints need validation on actual classroom tablets/Chromebooks

## Session Continuity

Last session: 2026-02-15 (Phase 3 in progress)
Stopped at: Completed 03-02-PLAN.md — Chart & Molecule Renderer Modules. Phase 3: 2/3 plans complete. Note: Plan 03-01 may be incomplete (uncommitted toSMILES implementation, failing test).
Resume file: None
