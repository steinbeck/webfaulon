---
phase: 01-molecular-graph-sa-core
plan: 01
subsystem: core
tags: [foundation, data-structures, algorithms, testing]
dependency_graph:
  requires: []
  provides: [MolGraph, Wiener-Index, project-scaffold]
  affects: [all-future-plans]
tech_stack:
  added: [Vite, TypeScript, Vitest]
  patterns: [BFS, adjacency-matrix, factory-methods, TDD]
key_files:
  created:
    - package.json
    - tsconfig.json
    - vitest.config.ts
    - src/core/types.ts
    - src/core/MolGraph.ts
    - src/core/wiener.ts
    - src/core/__tests__/MolGraph.test.ts
    - src/core/__tests__/wiener.test.ts
  modified: []
decisions:
  - choice: "Mutable bond matrix with explicit validation methods"
    rationale: "SA displacement needs to mutate graphs, but validation should be explicit to allow testing invalid intermediate states"
    alternatives: "Immutable with copy-on-write (rejected: performance overhead for SA)"
  - choice: "Constructor does not auto-validate connectivity/valences"
    rationale: "Allows SA engine to create intermediate states without throwing exceptions"
    impact: "Users must call validate() explicitly"
  - choice: "Array-based BFS queue with head index"
    rationale: "O(1) dequeue operations vs shift() being O(n)"
    impact: "Performance critical for n^2 Wiener computation"
metrics:
  duration_minutes: 4
  tasks_completed: 3
  tests_written: 47
  files_created: 8
  commits: 3
  completed_date: "2026-02-14"
---

# Phase 1 Plan 01: Project Scaffold and Core Graph Summary

TypeScript project with Vite+Vitest, MolGraph adjacency matrix data structure, and BFS-based Wiener Index calculation verified against Faulon paper values.

## Objective Achievement

Set up the foundational TypeScript project infrastructure and implemented the two core components that all future plans depend on: the MolGraph molecular representation and the Wiener Index cost function for simulated annealing optimization.

**Why this matters:** Every component in Phase 1 (displacement, initial structure generation, SA engine) operates on MolGraph. The Wiener Index is the objective function being optimized. Getting these right with comprehensive tests prevents cascading errors in later development.

## Tasks Completed

### Task 1: Initialize TypeScript Project with Vite and Vitest
**Commit:** `72129e9`

Created modern TypeScript development environment with:
- Vite for fast builds and development
- Vitest for unit testing with node environment (algorithmic code, no DOM)
- Strict TypeScript configuration with path aliases (`@/*` -> `src/*`)
- Core type definitions: `Atom`, `STANDARD_VALENCE`, `SAStepResult`

**Key files:**
- `package.json` - Test scripts and dependencies
- `tsconfig.json` - Strict mode, ES2020 target, path aliases
- `vitest.config.ts` - Node environment, test patterns
- `src/core/types.ts` - Foundation types for entire codebase

### Task 2: Implement MolGraph Class
**Commit:** `7f78f99`

Adjacency matrix-based molecular graph with:
- **Constructor:** Validates matrix dimensions, computes implicitH
- **Getters:** Readonly access to atoms, bonds, bond order sums
- **Validation:** `isConnected()` via BFS, `hasValidValences()` checks implicitH >= 0
- **Mutation:** `setBond(i, j, order)` with automatic implicitH recomputation
- **Clone:** Deep copy for SA displacement (mutate clone, validate, accept/reject)
- **Factory methods:** `createLinearAlkane(n)`, `createCyclohexane()`, `createBranched(pattern)`

**Design decisions:**
- Mutable bonds (SA needs it), but validation is explicit (not in constructor)
- implicitH eagerly computed on setBond (catches valence errors immediately)
- Readonly getters return copies (prevent accidental mutation)

**Tests:** 29 passing tests covering construction, factories, connectivity, valences, setBond, clone independence

### Task 3: Implement BFS-based Wiener Index
**Commit:** `31c2ae5`

O(n^2) Wiener Index calculation via breadth-first search:
- For each atom i, BFS to compute distances to all atoms j > i
- Sum all pairwise distances (topological descriptor)
- Throw error on disconnected graphs
- Array-based queue with head pointer for O(1) dequeue

**Verified against known values:**
- n-pentane (C5): 20 (matches formula n(n^2-1)/6)
- n-hexane (C6): 35
- Cyclohexane: 27
- Isobutane: 9
- Neopentane: 16

**Performance:** <1ms for 50-atom molecules (requirement was <5ms)

**Tests:** 18 passing tests including formula verification, manual distance counting, performance validation

## Verification Results

All success criteria met:
- ✅ 47 tests passing (29 MolGraph + 18 Wiener)
- ✅ Zero TypeScript errors
- ✅ MolGraph validates connectivity and valences correctly
- ✅ Wiener Index matches all known values from literature
- ✅ Performance <5ms for 50-atom molecules (actual: <1ms)

## Deviations from Plan

None - plan executed exactly as written. All three tasks completed with no blocking issues, no architectural changes needed, no unexpected bugs requiring fixes.

## Technical Notes

**MolGraph Design Pattern:**
The separation of construction from validation is critical for SA. The displacement operation needs to:
1. Clone the graph
2. Mutate bonds (setBond)
3. Validate the result
4. Accept or reject based on Metropolis criterion

If the constructor threw on invalid graphs, we couldn't create intermediate states for testing. The explicit `validate()` method gives the SA engine fine-grained control.

**Wiener Index Implementation:**
The BFS implementation uses an array-based queue with a head index instead of `Array.shift()` because shift is O(n). For the n^2 algorithm (running BFS n times), this optimization matters. Measured performance is well under the 5ms requirement.

**Test Coverage:**
Both components have comprehensive test suites. The Wiener tests include manual distance verification to ensure the BFS algorithm is correct, not just matching numbers. Factory methods make test setup trivial for future plans.

## Impact on Future Work

**Immediate dependencies (Plan 02):**
- Displacement operations will use `MolGraph.clone()` and `setBond()`
- Equations 7-11 from Faulon paper will operate on the bond matrix
- Validation methods will check resulting structures

**Phase 1 continuation:**
- Plan 03 (initial structure) will use factory methods as starting points
- Plan 04 (SA engine) will use Wiener Index as cost function
- All plans depend on this foundation

**Technical debt:** None. Code is production-ready with full test coverage.

## Next Steps

Proceed to **Plan 02: Displacement Operations (Equations 7-11)**. The MolGraph and Wiener Index foundation is solid and ready for the next layer.

---

## Self-Check: PASSED

**Created files verified:**
```
FOUND: package.json
FOUND: tsconfig.json
FOUND: vitest.config.ts
FOUND: src/core/types.ts
FOUND: src/core/MolGraph.ts
FOUND: src/core/wiener.ts
FOUND: src/core/__tests__/MolGraph.test.ts
FOUND: src/core/__tests__/wiener.test.ts
```

**Commits verified:**
```
FOUND: 72129e9 (Task 1: Project initialization)
FOUND: 7f78f99 (Task 2: MolGraph implementation)
FOUND: 31c2ae5 (Task 3: Wiener Index implementation)
```

**Test execution:**
```
47 tests passing
0 tests failing
TypeScript: 0 errors
```

All claims in this summary are verified against actual project state.
