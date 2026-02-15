# Roadmap: WebFaulon

## Overview

WebFaulon delivers a browser-based educational tool that brings Faulon's 1996 simulated annealing algorithm to life. Students enter molecular formulas and watch the algorithm explore constitutional isomer space in real time, with live optimization charts and 2D chemical structures. The three-phase journey progresses from core algorithm correctness through browser integration to polished visualization, ensuring chemical validity before adding visual feedback.

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

- [x] **Phase 1: Molecular Graph & SA Core** - Chemical correctness foundation
- [x] **Phase 2: Browser Integration & Controls** - Web Workers, user input, execution controls
- [ ] **Phase 3: Visualization & UX** - Charts, structure rendering, classroom-ready design

## Phase Details

### Phase 1: Molecular Graph & SA Core
**Goal**: Core algorithm produces chemically valid molecular structures through simulated annealing
**Depends on**: Nothing (first phase)
**Requirements**: ALG-01, ALG-02, ALG-03, ALG-04, ALG-05, ALG-06
**Success Criteria** (what must be TRUE):
  1. Algorithm can generate a valid initial molecular structure from any reasonable formula (C6H14, C8H10)
  2. SA displacement operations preserve chemical valence rules (no carbon with 5 bonds)
  3. Molecular graph remains connected after every SA move (no disconnected fragments)
  4. Wiener Index computes correctly in under 5ms for 50-atom molecules
  5. SA accepts/rejects moves according to Metropolis criterion with configurable max/min optimization
**Plans:** 4 plans

Plans:
- [x] 01-01-PLAN.md — Project scaffold, MolGraph data structure, and Wiener Index
- [x] 01-02-PLAN.md — Initial structure generation from molecular formula (TDD)
- [x] 01-03-PLAN.md — Faulon displacement equations 7-11 and seeded PRNG (TDD)
- [x] 01-04-PLAN.md — SA engine with Metropolis criterion and cooling schedules (TDD)

### Phase 2: Browser Integration & Controls
**Goal**: Users can configure and execute SA algorithm in browser without UI freezing
**Depends on**: Phase 1
**Requirements**: ALG-07, CTRL-01, CTRL-02, CTRL-03, CTRL-04, CTRL-05, INP-01, INP-02, INP-03
**Success Criteria** (what must be TRUE):
  1. User can enter molecular formula and see validation feedback (accept C6H14, reject C6H99)
  2. User can select from 3-5 preset example molecules without knowing chemistry
  3. User can adjust SA parameters (initial kT, cooling schedule, steps, cycles) via UI controls
  4. User can start, pause, and reset SA execution at any time
  5. Browser remains responsive during 10,000-step SA execution (Web Worker isolation working)
  6. RDKit.js WASM loads successfully in worker and processes molecular graphs
**Plans:** 3 plans

Plans:
- [x] 02-01-PLAN.md — SAEngine step-by-step execution API refactor (TDD)
- [x] 02-02-PLAN.md — Web Worker + Comlink integration and Vite project scaffold
- [x] 02-03-PLAN.md — Alpine.js UI controls, formula validation, presets, and RDKit.js WASM

### Phase 3: Visualization & UX
**Goal**: Students see optimization happening in real time with clear, classroom-ready visuals
**Depends on**: Phase 2
**Requirements**: VIZ-01, VIZ-02, VIZ-03, VIZ-04, UX-01, UX-02
**Success Criteria** (what must be TRUE):
  1. Live chart displays Wiener Index vs step number, updating smoothly during SA execution (like Figure 4 in paper)
  2. Current best molecular structure renders as 2D chemical drawing using RDKit.js
  3. Algorithm state (current step, temperature, current/best Wiener Index) visible throughout execution
  4. UI is readable when projected (clean design, good contrast, large fonts)
  5. Application works on tablets and phones (mobile responsive layout for BYOD classrooms)
**Plans:** 3 plans

Plans:
- [ ] 03-01-PLAN.md — SMILES generation for MolGraph + worker bestSMILES propagation (TDD)
- [ ] 03-02-PLAN.md — Chart.js live chart module + RDKit.js molecule renderer module
- [ ] 03-03-PLAN.md — UI integration, responsive redesign, and visual verification checkpoint

## Progress

**Execution Order:**
Phases execute in numeric order: 1 → 2 → 3

| Phase | Plans Complete | Status | Completed |
|-------|----------------|--------|-----------|
| 1. Molecular Graph & SA Core | 4/4 | Complete | 2026-02-14 |
| 2. Browser Integration & Controls | 3/3 | Complete | 2026-02-15 |
| 3. Visualization & UX | 0/3 | Not started | - |

---

*Roadmap created: 2026-02-14*
*Last updated: 2026-02-15 after phase 2 execution complete*
