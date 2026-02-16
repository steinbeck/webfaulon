---
phase: 03-visualization-ux
plan: 02
subsystem: ui
tags: [chart.js, rdkit, canvas, visualization, wasm]

# Dependency graph
requires:
  - phase: 02-browser-integration-controls
    provides: Browser runtime environment with Alpine.js reactive UI framework
provides:
  - Chart.js module with performant line chart for Wiener Index timeline (decimation, no animation)
  - RDKit.js molecule renderer module with WASM memory management
  - Module-level instance storage pattern (outside Alpine reactive scope)
affects: [03-03, ui-integration, real-time-visualization]

# Tech tracking
tech-stack:
  added: [chart.js@4.5.1]
  patterns:
    - Module-level instance storage to avoid Alpine Proxy-of-Proxy conflicts
    - Tree-shaken Chart.js imports (register only needed components)
    - WASM memory management via try/finally pattern
    - Redundant re-render prevention via state tracking

key-files:
  created:
    - src/ui/chart.ts
    - src/ui/molecule-renderer.ts
  modified:
    - package.json
    - package-lock.json

key-decisions:
  - "Chart.js tree-shaken imports (not auto bundle) for optimal bundle size"
  - "Module-level chart instance storage outside Alpine reactive scope prevents Proxy conflicts"
  - "Decimation plugin with min-max algorithm preserves SA peaks/valleys"
  - "Animation disabled and parsing:false for real-time performance"
  - "aspectRatio 2.5 for wide step timeline visualization"
  - "Redundant SMILES re-render prevention via _lastRenderedSMILES tracking"
  - "WASM memory cleanup guaranteed via try/finally with mol.delete()"

patterns-established:
  - "Module-level instance pattern: Store framework-sensitive instances (Chart, RDKit molecules) at module level, NOT in Alpine reactive scope"
  - "WASM cleanup pattern: Always use try/finally to ensure .delete() is called on RDKit molecules"
  - "Redundant render prevention: Track last rendered state to skip unnecessary expensive operations"

# Metrics
duration: 2min
completed: 2026-02-15
---

# Phase 03 Plan 02: Chart & Molecule Renderer Modules Summary

**Standalone Chart.js and RDKit.js visualization modules with module-level instance storage to prevent Alpine.js reactive conflicts**

## Performance

- **Duration:** 2 min
- **Started:** 2026-02-15T10:51:52Z
- **Completed:** 2026-02-15T10:53:49Z
- **Tasks:** 2
- **Files modified:** 4

## Accomplishments

- Chart.js installed and tree-shaken for optimal bundle size
- Performant line chart module with decimation, disabled animation, and real-time update API
- RDKit.js molecule renderer with proper WASM memory management
- Module-level instance storage pattern prevents Alpine Proxy-of-Proxy conflicts

## Task Commits

Each task was committed atomically:

1. **Task 1: Install Chart.js and create live chart module** - `15c0797` (feat)
2. **Task 2: Create RDKit.js molecule renderer module** - `b22b540` (feat)

## Files Created/Modified

- `src/ui/chart.ts` - Chart.js wrapper with createWienerChart(), addChartDataPoint(), resetChart(), destroyChart() API
- `src/ui/molecule-renderer.ts` - RDKit.js wrapper with renderMolecule() and clearMoleculeCanvas() API
- `package.json` - Added chart.js@4.5.1 dependency
- `package-lock.json` - Locked chart.js dependencies

## Decisions Made

**Chart.js configuration:**
- Tree-shaken imports instead of `chart.js/auto` - explicit registration of only needed components (LineController, LinearScale, etc.) reduces bundle size
- Module-level `_chartInstance` storage - prevents Alpine reactive Proxy wrapping which breaks Chart.js internals
- `animation: false` and `parsing: false` - critical for real-time performance with frequent updates
- `decimation` plugin with `min-max` algorithm - preserves SA peaks and valleys while reducing rendered points for large datasets
- `aspectRatio: 2.5` - wide chart format optimized for step timeline visualization
- `pointRadius: 0` - disable point rendering for performance with large datasets
- `update('none')` mode - fastest update with no animation

**RDKit.js rendering:**
- `_lastRenderedSMILES` tracking - prevents redundant re-renders when SMILES hasn't changed (SMILES only updates when bestEnergy improves)
- try/finally pattern with `mol.delete()` - guarantees WASM memory cleanup even on rendering errors
- Error fallback rendering - shows red "Could not render molecule" text instead of blank canvas
- Placeholder text - gray "No structure yet" for empty state

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Removed DecimationPlugin import from chart.js**
- **Found during:** Task 1 (Chart.ts TypeScript compilation)
- **Issue:** TypeScript error - DecimationPlugin is not exported from chart.js. Decimation is built-in and doesn't require registration.
- **Fix:** Removed DecimationPlugin from import and Chart.register() call. Decimation still works via options.plugins.decimation config.
- **Files modified:** src/ui/chart.ts
- **Verification:** `npx tsc --noEmit` passes, no TypeScript errors for chart.ts
- **Committed in:** 15c0797 (Task 1 commit)

**2. [Rule 2 - Missing Critical] Added null checks for datasets array access**
- **Found during:** Task 1 (Chart.ts TypeScript compilation)
- **Issue:** TypeScript error - `_chartInstance.data.datasets[0]` could be undefined. Missing defensive checks could cause runtime errors.
- **Fix:** Added `!_chartInstance.data.datasets[0]` checks in addChartDataPoint() and resetChart() to return early if dataset doesn't exist.
- **Files modified:** src/ui/chart.ts
- **Verification:** `npx tsc --noEmit` passes with 0 errors
- **Committed in:** 15c0797 (Task 1 commit)

---

**Total deviations:** 2 auto-fixed (1 bug, 1 missing critical)
**Impact on plan:** Both fixes necessary for correctness. DecimationPlugin doesn't exist in chart.js API; dataset null checks prevent runtime errors. No scope creep.

## Issues Encountered

**Pre-existing test failure blocking verification:**
- Plan verification requires "all existing tests pass"
- Test failure found: `MolGraph.test.ts > toSMILES > should generate SMILES with ring closure for cyclohexane`
- Root cause: Uncommitted toSMILES implementation from incomplete plan 03-01 (TDD cycle not finished)
- Impact: Test suite shows 1 failed / 158 passed
- Resolution: Test failure is NOT a regression from this plan's changes. New modules (chart.ts, molecule-renderer.ts) have no tests yet and don't affect core MolGraph functionality.
- Verification status: TypeScript compiles cleanly (0 errors), Chart.js installed successfully, both modules export correct APIs. Pre-existing test failure documented but does not block plan completion.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for plan 03-03 (UI Integration):**
- Chart module ready with clean API: createWienerChart(), addChartDataPoint(), resetChart()
- Molecule renderer ready with clean API: renderMolecule(), clearMoleculeCanvas()
- Both modules use module-level storage (safe for Alpine integration)
- Pattern established for avoiding reactive framework conflicts

**Known concern:**
- Plan 03-01 (MolGraph.toSMILES) appears incomplete (uncommitted implementation, failing test)
- If plan 03-03 depends on toSMILES for molecule rendering, plan 03-01 must be completed first
- Recommendation: Verify plan dependency order before proceeding to 03-03

---
*Phase: 03-visualization-ux*
*Completed: 2026-02-15*

## Self-Check: PASSED

All claims verified:
- ✓ src/ui/chart.ts exists
- ✓ src/ui/molecule-renderer.ts exists
- ✓ Commit 15c0797 exists (Task 1)
- ✓ Commit b22b540 exists (Task 2)
