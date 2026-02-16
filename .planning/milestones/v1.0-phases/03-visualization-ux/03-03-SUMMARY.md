---
phase: 03-visualization-ux
plan: 03
subsystem: ui
tags: [alpine, chart.js, rdkit, responsive, molblock, visualization]

# Dependency graph
requires:
  - phase: 03-visualization-ux/03-01
    provides: toMolBlock() on MolGraph, bestMolBlock propagation through worker pipeline
  - phase: 03-visualization-ux/03-02
    provides: Chart.js live chart module, RDKit.js molecule renderer module
provides:
  - Fully integrated visualization UI with live chart, 2D molecule rendering, and responsive layout
  - Canonical SMILES derived from RDKit.js (not custom generation)
  - Classroom-ready design with responsive breakpoints
affects: [project-complete, milestone-1]

# Tech tracking
tech-stack:
  patterns:
    - MOL block pipeline: MolGraph → V2000 MOL block → RDKit.js parse → set_new_coords() → draw_to_canvas()
    - RDKit-derived canonical SMILES displayed in UI (not custom SMILES generation)
    - Alpine bestSMILES state separate from worker bestMolBlock data

key-files:
  modified:
    - src/core/MolGraph.ts (toSMILES → toMolBlock)
    - src/core/types.ts (bestSMILES → bestMolBlock)
    - src/worker/types.ts (bestSMILES → bestMolBlock)
    - src/core/SAEngine.ts (calls toMolBlock)
    - src/worker/sa-worker.ts (sends bestMolBlock)
    - src/ui/molecule-renderer.ts (accepts MOL block, set_new_coords, derives SMILES)
    - src/ui/app.ts (bestMolBlock pipeline, Alpine bestSMILES state)
    - index.html (SMILES display binds to Alpine bestSMILES)
    - src/ui/presets.ts (added C10H16)
    - src/core/__tests__/MolGraph.test.ts (toSMILES → toMolBlock tests)

key-decisions:
  - "Replaced toSMILES() with toMolBlock() — MOL blocks are trivially correct from adjacency matrix, no ring-closure logic needed"
  - "RDKit.js set_new_coords() required before draw_to_canvas() — all-zero MOL block coordinates cause WASM crash"
  - "Pipeline sends MOL block from worker, main thread uses RDKit to parse, render, and derive canonical SMILES"
  - "Added C10H16 (monoterpene isomers) to preset dropdown per user request"

patterns-established:
  - "MOL block rendering: Always call mol.set_new_coords() before mol.draw_to_canvas() when MOL block has zero coordinates"

# Metrics
duration: ~90min (across 2 sessions, including MOL block debugging)
completed: 2026-02-15
---

# Phase 03 Plan 03: UI Integration & Responsive Redesign Summary

**Wired Chart.js and molecule renderer into Alpine.js UI, redesigned layout with responsive CSS, fixed MOL block rendering pipeline**

## Performance

- **Duration:** ~90 min (across 2 sessions)
- **Tasks:** 3 (2 auto + 1 human verification checkpoint)
- **Files modified:** 10
- **Commits:** 5 (25d9e32, 4f5dba3, 78f893e, 9f43627, ac3f6ab)

## Accomplishments

- Chart.js live chart wired into Alpine progress callback, updates every 10 steps
- RDKit.js molecule renderer integrated, re-renders only when best energy improves
- Responsive layout: 2-column chart+molecule on desktop, single column on mobile
- Algorithm state panel shows step, temperature, current/best Wiener Index
- Replaced broken toSMILES() with toMolBlock() + RDKit.js set_new_coords() for reliable rendering
- Canonical SMILES derived by RDKit displayed below structure
- Added C10H16 (monoterpene isomers) preset
- All 159 tests pass, TypeScript clean

## Key Bug Fix

**Root cause of WASM crash:** `mol.draw_to_canvas()` throws a WASM exception when atom coordinates are all zeros (degenerate geometry). The MOL block from MolGraph stores topology only — no 2D coordinates. The fix: call `mol.set_new_coords()` before `mol.draw_to_canvas()` to let RDKit generate a 2D layout.

**Why toSMILES was replaced:** The DFS-based SMILES generation was fundamentally broken for complex ring systems (misattributed parent/child relationships, produced invalid structures). MOL block generation is trivially correct from the adjacency matrix — just list atoms and bonds from the upper triangle, no ring-closure logic needed.

## Human Verification Results

All 6 tests passed:
1. Chart visualization — smooth live updates for multiple molecules
2. Molecule rendering — linear chains, polycyclic rings, fused systems all correct
3. Algorithm state — all 4 stats visible and updating during execution
4. Responsive layout — single column at 375px, side-by-side at 1400px+
5. Classroom projection — large fonts, good contrast, no horizontal scroll
6. Full cycle — Start/Reset cycle works correctly

## Deviations from Plan

### Major: Replaced SMILES pipeline with MOL block pipeline

- **Found during:** Task 3 human verification (molecule rendering not working)
- **Issue:** toSMILES() generated invalid SMILES for complex ring systems; RDKit.js draw_to_canvas() crashed on zero-coordinate MOL blocks
- **Fix:** Replaced toSMILES() with toMolBlock() across entire pipeline (MolGraph → types → worker → renderer). Added set_new_coords() call in renderer.
- **Impact:** 10 files changed, but architecture is cleaner — MOL blocks are trivially correct, SMILES comes from RDKit (the authority)
- **Commits:** 78f893e (attempted SMILES fix), 9f43627 (paused for debugging), ac3f6ab (final MOL block fix)

---
*Phase: 03-visualization-ux*
*Completed: 2026-02-15*
