# Project State

## Project Reference

See: .planning/PROJECT.md (updated 2026-02-14)

**Core value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.
**Current focus:** Phase 1 - Molecular Graph & SA Core

## Current Position

Phase: 1 of 3 (Molecular Graph & SA Core)
Plan: 3 of 4 in current phase
Status: Executing
Last activity: 2026-02-14 — Completed plan 01-03 with SeededRandom PRNG and Faulon displacement

Progress: [███████░░░] 75% (3/4 plans complete in Phase 1)

## Performance Metrics

**Velocity:**
- Total plans completed: 3
- Average duration: 6.7 min
- Total execution time: 0.33 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 01-molecular-graph-sa-core | 3 | 20 min | 6.7 min |

| Phase | Plan | Duration | Tasks | Files |
|-------|------|----------|-------|-------|
| 01-molecular-graph-sa-core | 01 | 4 min | 3 | 8 |
| 01-molecular-graph-sa-core | 02 | 4 min | 3 | 4 |
| 01-molecular-graph-sa-core | 03 | 8 min | 3 | 4 |

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

Last session: 2026-02-14 (plan execution)
Stopped at: Completed 01-03-PLAN.md with 104 tests passing (29 new: 14 PRNG + 15 displacement)
Resume file: None
