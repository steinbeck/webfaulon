---
phase: 03-visualization-ux
plan: 01
subsystem: core
tags: [smiles, molecular-visualization, graph-algorithms, dfs]

# Dependency graph
requires:
  - phase: 01-molecular-graph-sa-core
    provides: "MolGraph class with atoms and bonds representation"
  - phase: 02-browser-integration-controls
    provides: "Worker progress callbacks and SAEngine state propagation"
provides:
  - "SMILES string generation from MolGraph via DFS-based algorithm"
  - "bestSMILES field in SAEngineState and SAProgressData for 2D rendering"
affects: [03-visualization-ux]

# Tech tracking
tech-stack:
  added: []
  patterns: ["DFS-based SMILES generation with ring closure detection"]

key-files:
  created: []
  modified:
    - "src/core/MolGraph.ts"
    - "src/core/__tests__/MolGraph.test.ts"
    - "src/core/types.ts"
    - "src/core/SAEngine.ts"
    - "src/worker/types.ts"
    - "src/worker/sa-worker.ts"

key-decisions:
  - "DFS-based SMILES generation produces valid (non-canonical) SMILES suitable for RDKit parsing"
  - "Ring closures detected via back edges in DFS traversal"
  - "toSMILES() called in getState() every reportInterval steps (negligible cost vs Wiener Index computation)"

patterns-established:
  - "TDD with RED → GREEN → REFACTOR phases for new core functionality"
  - "Comprehensive test coverage for molecular graph algorithms (linear, branched, cyclic, multiple bond orders)"

# Metrics
duration: 5min
completed: 2026-02-15
---

# Phase 03 Plan 01: SMILES Generation & Progress Propagation Summary

**SMILES string generation via DFS with ring closure detection, propagated through SAEngine state and worker progress callbacks for 2D molecular rendering**

## Performance

- **Duration:** 5 min
- **Started:** 2026-02-15T10:51:49Z
- **Completed:** 2026-02-15T10:56:57Z
- **Tasks:** 2
- **Files modified:** 6

## Accomplishments
- MolGraph.toSMILES() generates valid SMILES for linear chains, branched molecules, cyclic structures, and molecules with multiple bond orders
- bestSMILES field added to SAEngineState and SAProgressData interfaces
- Worker includes best molecule SMILES in every progress report for main thread 2D rendering
- Zero regressions - all 149 existing tests + 10 new toSMILES tests pass

## Task Commits

Each task was committed atomically following TDD protocol:

1. **Task 1: TDD - Add toSMILES() to MolGraph**
   - `93cfcb6` (test) - RED phase: 10 failing tests for toSMILES functionality
   - `505d171` (feat) - GREEN phase: DFS-based SMILES generation implementation

2. **Task 2: Add bestSMILES to SAEngine state and worker progress**
   - `bda0732` (feat) - bestSMILES propagation through state and worker

## Files Created/Modified
- `src/core/MolGraph.ts` - Added toSMILES() method with DFS-based SMILES generation
- `src/core/__tests__/MolGraph.test.ts` - Added 10 comprehensive tests for SMILES generation
- `src/core/types.ts` - Added bestSMILES field to SAEngineState interface
- `src/core/SAEngine.ts` - getState() now includes bestGraph.toSMILES()
- `src/worker/types.ts` - Added bestSMILES field to SAProgressData interface
- `src/worker/sa-worker.ts` - Worker includes bestSMILES in progress reports

## Decisions Made

**SMILES generation algorithm:**
- Chose DFS-based approach for simplicity and correctness
- Ring closures detected via back edges in DFS traversal (standard SMILES algorithm)
- Non-canonical SMILES output (canonicalization is RDKit's responsibility)
- Ring closure digits added after atom symbol when back edge detected

**Performance consideration:**
- toSMILES() called every reportInterval steps (default: 10) via getState()
- DFS traversal is O(n) where n = atom count - negligible compared to Wiener Index O(n²)
- Typical SA run: 2000 steps / 10 = 200 SMILES generations for ~10 atom molecules = negligible cost

**Backward compatibility:**
- New bestSMILES field added to existing interfaces
- Existing app.ts progress callback uses spread operator (`{ ...data }`), so fully compatible
- UI doesn't consume bestSMILES yet (planned for 03-02), but field is present in data stream

## Deviations from Plan

None - plan executed exactly as written. TDD RED → GREEN → REFACTOR protocol followed. All tests specified in plan were implemented and pass.

## Issues Encountered

**Ring closure bug during GREEN phase:**
- **Issue:** Initial implementation generated incorrect SMILES for cyclohexane (`CCCCCC1(C2)` instead of `CCCCCC1`)
- **Root cause:** DFS neighbor classification (visited vs unvisited) done at function entry, but neighbors' visited status changed during recursive DFS traversal
- **Fix:** Split neighbor processing into two passes: (1) handle back edges (ring closures) first, (2) traverse unvisited neighbors with real-time visited check
- **Verification:** Cyclohexane test passes, generates valid `CCCCCC1` SMILES
- **Committed in:** 505d171 (Task 1 GREEN phase commit)

This was a standard debugging iteration during TDD GREEN phase, not a deviation from plan.

## User Setup Required

None - no external service configuration required.

## Next Phase Readiness

**Ready for 03-02 (RDKit.js 2D Structure Rendering):**
- bestSMILES field present in SAProgressData
- Worker delivers SMILES string every 10 steps to main thread
- SMILES format validated against test cases (linear, branched, cyclic, multiple bonds, heteroatoms)

**Blockers:** None

**Notes for next plan:**
- RDKit.js can now receive SMILES via progress callback
- Plan 03-02 will consume `data.bestSMILES` to render 2D structure
- SMILES generation tested for all molecule types used in existing factory methods

---
*Phase: 03-visualization-ux*
*Completed: 2026-02-15*

## Self-Check: PASSED

All claims verified:

**Files modified (6):**
- ✓ src/core/MolGraph.ts exists (added toSMILES method)
- ✓ src/core/__tests__/MolGraph.test.ts exists (added 10 tests)
- ✓ src/core/types.ts exists (added bestSMILES field)
- ✓ src/core/SAEngine.ts exists (getState includes bestSMILES)
- ✓ src/worker/types.ts exists (SAProgressData has bestSMILES)
- ✓ src/worker/sa-worker.ts exists (worker reports bestSMILES)

**Commits (3):**
- ✓ 93cfcb6 (test) - RED phase: failing toSMILES tests
- ✓ 505d171 (feat) - GREEN phase: toSMILES implementation
- ✓ bda0732 (feat) - Task 2: bestSMILES propagation

**Code changes:**
- ✓ toSMILES() method present in MolGraph.ts
- ✓ bestSMILES field in SAEngineState interface
- ✓ bestSMILES field in SAProgressData interface
- ✓ bestSMILES in worker progress reporting

**Test results:**
- ✓ 159 tests pass (149 existing + 10 new)
- ✓ TypeScript compiles without errors
