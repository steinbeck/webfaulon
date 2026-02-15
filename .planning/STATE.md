# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.
**Current focus:** Phase 2 - Browser Integration & Controls (Phase 1 complete)

## Current Position

Phase: 2 of 3 (Browser Integration & Controls)
Plan: 1 of 3 in current phase
Status: In progress
Last activity: 2026-02-15 — Completed plan 02-01 (Step-by-step SA execution API)

Progress: [███░░░░░░░] 33% (Phase 2: 1/3 plans complete)

## Performance Metrics

**Velocity:**
- Total plans completed: 5
- Average duration: 6.2 min
- Total execution time: 0.52 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-molecular-graph-sa-core | 4 | 28 min | 7.0 min |
| 02-browser-integration-controls | 1 | 3 min | 3.0 min |

| Phase | Plan | Duration | Tasks | Files |
|-------|------|----------|-------|-------|
| 01-molecular-graph-sa-core | 01 | 4 min | 3 | 8 |
| 01-molecular-graph-sa-core | 02 | 4 min | 3 | 4 |
| 01-molecular-graph-sa-core | 03 | 8 min | 3 | 4 |
| 01-molecular-graph-sa-core | 04 | 10 min | 3 | 4 |
| 02-browser-integration-controls | 01 | 3 min | 2 | 3 |

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

Last session: 2026-02-15 (phase 2 in progress)
Stopped at: Completed 02-01-PLAN.md — Step-by-step SA execution API. Ready for plan 02-02 (Web Worker integration).
Resume file: None
