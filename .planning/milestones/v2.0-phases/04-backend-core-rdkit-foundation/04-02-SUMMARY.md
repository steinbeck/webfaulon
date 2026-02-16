---
phase: 04-backend-core-rdkit-foundation
plan: 02
subsystem: backend-core
tags: [tdd, port, utilities, prng, formula-parser, cooling-schedule]
dependency-graph:
  requires: []
  provides:
    - backend/app/core/random.py (SeededRandom PRNG)
    - backend/app/core/formula_parser.py (Formula parsing + HDI)
    - backend/app/core/cooling.py (Faulon cooling schedule)
  affects:
    - Plan 03 (Displacement will use SeededRandom)
    - Plan 04 (SAEngine will use SeededRandom + cooling)
tech-stack:
  added:
    - pytest (testing framework)
    - Python regex (formula parsing)
  patterns:
    - TDD with RED-GREEN-REFACTOR cycles
    - Cross-language determinism (Python/TypeScript)
    - 32-bit arithmetic emulation in Python
key-files:
  created:
    - backend/pyproject.toml
    - backend/app/__init__.py
    - backend/app/core/__init__.py
    - backend/app/core/random.py
    - backend/app/core/formula_parser.py
    - backend/app/core/cooling.py
    - backend/tests/__init__.py
    - backend/tests/test_random.py
    - backend/tests/test_formula_parser.py
    - backend/tests/test_cooling.py
  modified: []
decisions:
  - decision: Use pip for pytest instead of Poetry
    rationale: Poetry not available in environment, pip works for simple dependencies
    impact: Can switch to Poetry when Plan 01 creates full setup
  - decision: Emulate JavaScript Math.imul with signed 32-bit conversion
    rationale: Python has arbitrary precision, need exact match with TypeScript Mulberry32
    impact: Achieved perfect cross-language determinism
  - decision: Use regex for formula parsing instead of custom parser
    rationale: Simple, readable, matches TypeScript implementation pattern
    impact: Clean implementation, all edge cases handled
metrics:
  duration: 280 seconds (4 minutes)
  completed: 2026-02-16T08:59:26Z
  tasks: 2
  tests: 41 (all passing)
  commits: 5
---

# Phase 04 Plan 02: Pure Python Utilities (TDD Port) Summary

**One-liner:** Port SeededRandom (Mulberry32), formula parser with HDI, and Faulon cooling schedule from TypeScript to Python with TDD, achieving cross-language determinism.

## Objective Achieved

Successfully ported three pure-Python utility modules from TypeScript to Python using strict TDD methodology. All modules produce identical output to their TypeScript counterparts, verified by comprehensive test suites ported from the original vitest tests.

**Output:** Three tested Python modules with 41 passing tests, ready for use by displacement and SAEngine implementations.

## Tasks Completed

### Task 1: Port SeededRandom with TDD
**Duration:** ~150 seconds | **Commits:** 9c1b91d (RED), 3db90df (GREEN)

**RED Phase:**
- Created backend directory structure (pyproject.toml, __init__.py files)
- Installed pytest as test runner
- Ported all 14 test cases from random.test.ts to test_random.py
- Tests covered: determinism, next() range, nextInt() inclusive range, selectNDistinct()
- Verified tests failed (no implementation)

**GREEN Phase:**
- Implemented SeededRandom class with Mulberry32 algorithm
- Critical challenge: Python arbitrary-precision integers vs JavaScript 32-bit
- Solution: Emulate Math.imul with signed 32-bit conversion and masking
- Implemented `_imul32()` helper to match JavaScript's Math.imul exactly
- All 14 tests passed
- **Cross-language verification:** SeededRandom(42) produces IDENTICAL 10-value sequences in Python and TypeScript

**REFACTOR Phase:**
- No refactoring needed - code clean with type hints and documentation

**Verification:**
```
cd backend && python3 -m pytest tests/test_random.py -v
14/14 tests passed
```

### Task 2: Port formula parser and cooling schedule with TDD
**Duration:** ~130 seconds | **Commits:** 3d86dce (RED), 3f9c5bb (GREEN)

**Formula Parser (RED-GREEN):**

**RED Phase:**
- Ported 15 test cases from formulaParser.test.ts to test_formula_parser.py
- Tests covered: simple alkanes, heteroatoms, multi-digit counts, HDI computation
- Verified tests failed (no implementation)

**GREEN Phase:**
- Implemented parse_formula() with regex pattern `([A-Z][a-z]?)(\d*)`
- Handles implicit count of 1, double-digit counts, order independence
- Supports all known elements: H, C, N, O, S, P, F, Cl, Br, I
- Implemented compute_hdi() with formula: (2C + 2 + N - H - Halogens) / 2
- Properly handles halogens as hydrogen equivalents
- All 15 tests passed

**Cooling Schedule (RED-GREEN):**

**RED Phase:**
- Ported 12 test cases from cooling.test.ts to test_cooling.py
- Tests covered: k=0 (constant), k=1 (linear), k=8 (fast decay), k=32 (very fast)
- Verified tests failed (no implementation)

**GREEN Phase:**
- Implemented compute_temperature() with Faulon formula: T = T0 - k * T0 * t / delta_t
- Proper clamping to MIN_TEMPERATURE = 0.01 (avoids Metropolis division by zero)
- All 12 tests passed

**REFACTOR Phase:**
- No refactoring needed - implementations clean and match TypeScript

**Verification:**
```
cd backend && python3 -m pytest tests/test_formula_parser.py tests/test_cooling.py -v
27/27 tests passed
```

## Verification Results

All verification criteria from plan met:

1. `cd backend && pytest tests/test_random.py -v` -- 14/14 passed
2. `cd backend && pytest tests/test_formula_parser.py -v` -- 15/15 passed
3. `cd backend && pytest tests/test_cooling.py -v` -- 12/12 passed
4. SeededRandom(42) produces same sequence as TypeScript new SeededRandom(42) -- VERIFIED
5. parse_formula("C6H14") returns {"C": 6, "H": 14} -- VERIFIED
6. compute_hdi({"C": 6, "H": 6}) returns 4 -- VERIFIED
7. compute_temperature(125, 1000, 100, 8) returns 0.01 -- VERIFIED

**Total:** 41/41 tests passing

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 3 - Blocking] Poetry not available in environment**
- **Found during:** Task 1 initialization
- **Issue:** Plan assumed Poetry would be available, but `poetry` command not found
- **Fix:** Used pip to install pytest directly, created minimal pyproject.toml manually
- **Files modified:** backend/pyproject.toml (created manually)
- **Commit:** 9c1b91d
- **Rationale:** Unblocks execution. Can switch to full Poetry setup when Plan 01 completes.

No other deviations. Plan executed exactly as written after initial Poetry workaround.

## Key Decisions Made

1. **32-bit arithmetic emulation strategy**
   - Chose signed 32-bit conversion approach for Math.imul
   - Alternative considered: numpy uint32 arrays (rejected as overkill)
   - Result: Perfect cross-language determinism achieved

2. **Test framework setup**
   - Used pytest directly via pip instead of Poetry
   - Alternative: Wait for Plan 01 to complete (rejected - wave 1 parallel)
   - Result: Can proceed independently, easy migration to Poetry later

## Technical Insights

1. **Mulberry32 cross-language porting**
   - JavaScript's `>>> 0` converts to unsigned 32-bit by masking high bits
   - JavaScript's `Math.imul` treats operands as signed 32-bit, returns signed result
   - Python solution: Check for sign bit (0x80000000), subtract 0x100000000 if set

2. **Formula parsing edge cases**
   - Empty element matches at end of regex - filter with `if not element: continue`
   - Order independence automatic with dict accumulation
   - HDI formula: Oxygen/Sulfur don't affect (divalent), halogens count as H

3. **Cooling schedule behavior**
   - k=8 reaches min at 1/8 of total steps (paper's recommended schedule)
   - MIN_TEMPERATURE = 0.01 prevents Metropolis division by zero
   - Formula is linear reduction, clamped to minimum

## Success Criteria Met

- [x] All three Python modules are exact behavioral ports of TypeScript originals
- [x] Test suites ported from TypeScript all pass (41/41)
- [x] SeededRandom determinism verified across languages (identical 10-value sequences)
- [x] Formula parser handles all element types including halogens and multi-digit counts
- [x] Cooling schedule matches Faulon paper formulas for all k values (0, 1, 8, 32)

## Artifacts Produced

**Core Modules (3 files, 313 LOC):**
- `backend/app/core/random.py` - SeededRandom class with Mulberry32 PRNG
- `backend/app/core/formula_parser.py` - parse_formula() and compute_hdi()
- `backend/app/core/cooling.py` - compute_temperature() with Faulon formula

**Test Suite (3 files, 379 LOC):**
- `backend/tests/test_random.py` - 14 tests for PRNG determinism and behavior
- `backend/tests/test_formula_parser.py` - 15 tests for parsing and HDI
- `backend/tests/test_cooling.py` - 12 tests for cooling schedules

**Infrastructure:**
- `backend/pyproject.toml` - Minimal project configuration
- `backend/app/__init__.py`, `backend/app/core/__init__.py`, `backend/tests/__init__.py` - Package structure

## Dependencies Created

**Requires:** None (wave 1 parallel with Plan 01)

**Provides:**
- SeededRandom for displacement (Plan 03)
- SeededRandom + cooling for SAEngine (Plan 04)
- Formula parser for initial structure generation (Plan 04)

**Affects:**
- Plan 03 can now use SeededRandom for reproducible atom selection
- Plan 04 can now use all three utilities for SA implementation

## Next Steps

1. Plan 01 will create full backend FastAPI structure - merge this minimal setup
2. Plan 03 can use SeededRandom for displacement implementation
3. Plan 04 can use all three modules for SAEngine

## Self-Check: PASSED

**Created files verified:**
- FOUND: backend/pyproject.toml
- FOUND: backend/app/__init__.py
- FOUND: backend/app/core/__init__.py
- FOUND: backend/app/core/random.py
- FOUND: backend/app/core/formula_parser.py
- FOUND: backend/app/core/cooling.py
- FOUND: backend/tests/__init__.py
- FOUND: backend/tests/test_random.py
- FOUND: backend/tests/test_formula_parser.py
- FOUND: backend/tests/test_cooling.py

**Commits verified:**
- FOUND: 9c1b91d (test: RED phase SeededRandom)
- FOUND: 3db90df (feat: GREEN phase SeededRandom)
- FOUND: 3d86dce (test: RED phase formula parser + cooling)
- FOUND: 3f9c5bb (feat: GREEN phase formula parser + cooling)

**Test results verified:**
- All 41 tests passing
- Cross-language determinism verified for SeededRandom(42)
- Formula parser correctly handles all edge cases
- Cooling schedule matches TypeScript values for all k parameters
