# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.
**Current focus:** Phase 1 - Molecular Graph & SA Core

## Current Position

Phase: 1 of 3 (Molecular Graph & SA Core)
Plan: 1 of 4 in current phase
Status: Executing
Last activity: 2026-02-14 — Completed plan 01-01 with MolGraph and Wiener Index implementation

Progress: [██░░░░░░░░] 25% (1/4 plans complete in Phase 1)

## Performance Metrics

**Velocity:**
- Total plans completed: 1
- Average duration: 4 min
- Total execution time: 0.07 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-molecular-graph-sa-core | 1 | 4 min | 4 min |

| Phase | Plan | Duration | Tasks | Files |
|-------|------|----------|-------|-------|
| 01-molecular-graph-sa-core | 01 | 4 min | 3 | 8 |

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

### Pending Todos

None yet.

### Blockers/Concerns

**Phase 1 planning:**
- Research gap: Faulon displacement equations (eqs 7-11) details needed from original 1996 paper or reference implementation
- Decision needed: Initial structure generation algorithm (deterministic vs stochastic approach)

**Phase 3 planning:**
- Chart decimation threshold needs performance testing (500 vs 1000 vs 2000 sample points)
- Mobile breakpoints need validation on actual classroom tablets/Chromebooks

## Session Continuity

Last session: 2026-02-14 (plan execution)
Stopped at: Completed 01-01-PLAN.md with 47 tests passing
Resume file: None
