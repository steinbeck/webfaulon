---
phase: 04-backend-core-rdkit-foundation
plan: 03
subsystem: backend-core
tags: [tdd, port, rdkit, molecules, displacement, critical]
dependency-graph:
  requires:
    - Plan 02 (SeededRandom, formula_parser, cooling)
  provides:
    - backend/app/core/molecule.py (MoleculeGraph RDKit wrapper)
    - backend/app/core/displacement.py (Faulon displacement)
    - backend/app/core/wiener.py (Wiener Index)
    - backend/app/core/initial_structure.py (Initial structure generation)
    - backend/app/utils/rdkit_helpers.py (Bond type conversion)
  affects:
    - Plan 04 (SAEngine will use all molecular operations)
tech-stack:
  added:
    - RDKit 2025.09.5 (Python chemistry library)
    - RDKit GetDistanceMatrix (Wiener Index computation)
    - RDKit SanitizeMol (valence validation)
  patterns:
    - RDKit Mol wrapper with adjacency-matrix-like API
    - TDD with RED-GREEN-REFACTOR cycles
    - Automatic sanitization after bond modifications
    - Factory methods for test molecule construction
key-files:
  created:
    - backend/app/utils/rdkit_helpers.py
    - backend/app/core/molecule.py
    - backend/app/core/wiener.py
    - backend/app/core/displacement.py
    - backend/app/core/initial_structure.py
    - backend/tests/test_molecule.py
    - backend/tests/test_wiener.py
    - backend/tests/test_displacement.py
    - backend/tests/test_initial_structure.py
  modified: []
decisions:
  - decision: Wrap RDKit RWMol instead of pure adjacency matrix
    rationale: RDKit provides built-in valence validation, sanitization, and distance matrix computation - don't hand-roll chemistry
    impact: Automatic valence checking, simpler implementation, leverages battle-tested chemistry library
  - decision: Call SanitizeMol after every bond modification
    rationale: RDKit requires sanitization to recompute implicit hydrogens and validate chemistry
    impact: Bond modifications may fail if intermediate states invalid, but guarantees chemical correctness
  - decision: Use RDKit GetDistanceMatrix for Wiener Index
    rationale: Built-in BFS implementation, already optimized, handles disconnected graphs
    impact: Single-line Wiener computation instead of 80 lines of custom BFS
metrics:
  duration: 347 seconds (5 minutes 47 seconds)
  completed: 2026-02-16T09:09:18Z
  tasks: 2
  tests: 59 (all passing)
  commits: 2
---

# Phase 04 Plan 03: RDKit Molecular Operations (TDD Port) Summary

**One-liner:** Port RDKit-dependent molecular modules from TypeScript to Python with TDD - MoleculeGraph wrapper, Faulon displacement, Wiener Index, and initial structure generation - validating TypeScript-to-RDKit translation against v1.0 test values.

## Objective Achieved

Successfully ported the four critical RDKit-dependent molecular modules from TypeScript to Python using strict TDD methodology. All Wiener Index values match v1.0 exactly. The 500-displacement stress test produces zero invalid molecules, confirming that Faulon displacement on RDKit RWMol maintains chemical correctness.

**Output:** Four tested Python modules with 59 passing tests, providing the molecular foundation for SAEngine implementation.

## Tasks Completed

### Task 1: RDKit helpers, MoleculeGraph wrapper, and Wiener Index (TDD)
**Duration:** ~180 seconds | **Commits:** 931c82a

**Prerequisite:**
- Installed RDKit 2025.09.5 via `poetry run pip install rdkit`
- Verified installation: `from rdkit import Chem; print(Chem.MolToSmiles(Chem.MolFromSmiles('CCO')))` → "CCO"

**RDKit Helpers (backend/app/utils/rdkit_helpers.py):**
- `order_to_bond_type(order: int) -> Chem.BondType`: Map Faulon integer (0,1,2,3) to RDKit BondType enum
- `bond_type_to_order(bond_type: Chem.BondType) -> int`: Reverse mapping
- Explicit dict mapping (NOT direct cast) per research guidance

**MoleculeGraph Wrapper (RED-GREEN-REFACTOR):**

**RED Phase:**
- Created `backend/tests/test_molecule.py` with 18 tests ported from `MolGraph.test.ts`
- Tests covered: factory methods (linear alkane, cyclohexane, branched), connectivity, valences, bond manipulation, cloning, SMILES generation, implicit H counts
- Verified tests failed (no implementation)

**GREEN Phase:**
- Created `backend/app/core/molecule.py`: RDKit RWMol wrapper with adjacency-matrix-like API
- Key design: Unlike TypeScript's pure adjacency matrix, this wraps RDKit's native Mol object
- Properties: `get_atom_count()`, `get_bond_order()`, `get_bond_order_sum()`, `get_implicit_h()`, `get_atom_element()`
- Mutation: `set_bond()` - handles bond addition/removal/modification with automatic SanitizeMol
- Validation: `is_connected()` via RDKit GetMolFrags, `has_valid_valences()` via SanitizeMol try/catch
- Utilities: `clone()`, `to_smiles()`, `get_mol()`
- Factory methods: `create_linear_alkane(n)`, `create_cyclohexane()`, `create_branched(pattern)`
- CRITICAL: RDKit computes implicit H automatically after sanitization via `GetNumImplicitHs()`
- All 18 tests passed

**Wiener Index (RED-GREEN-REFACTOR):**

**RED Phase:**
- Created `backend/tests/test_wiener.py` with 13 tests ported from `wiener.test.ts`
- Tests covered: known values for linear alkanes (methane=0, ethane=1, propane=4, butane=10, pentane=20, hexane=35), cyclic (cyclohexane=27), branched (isobutane=9, neopentane=16), formula verification
- Verified tests failed (no implementation)

**GREEN Phase:**
- Created `backend/app/core/wiener.py`: Single-function implementation using RDKit GetDistanceMatrix
- Algorithm: `Chem.GetDistanceMatrix(mol)` → sum upper triangle → Wiener Index
- Disconnected graph detection: distance > 1000 (effectively infinite)
- All 13 tests passed

**Verification:**
```bash
cd backend && poetry run pytest tests/test_molecule.py tests/test_wiener.py -v
31/31 tests passed
```

**Wiener Index Validation (matches v1.0 exactly):**
- hexane=35 ✓
- cyclohexane=27 ✓
- isobutane=9 ✓
- neopentane=16 ✓

### Task 2: Displacement and initial structure generation (TDD)
**Duration:** ~167 seconds | **Commits:** 9d420bd

**Faulon Displacement (RED-GREEN-REFACTOR):**

**RED Phase:**
- Created `backend/tests/test_displacement.py` with 15 test classes ported from `displacement.test.ts`
- Tests covered: basic validation (< 4 atoms, atom count/type preservation, no mutation), bond order constraints (range [0,3], existing double bonds), connectivity/valence validation, reproducibility (same seed), stress test (500 displacements), cyclohexane, bond conservation equations
- Verified tests failed (no implementation)

**GREEN Phase:**
- Created `backend/app/core/displacement.py`: Port of `attemptDisplacement` and `computeDisplacementBonds`
- Constants: `MAX_BOND_ORDER = 3`, `MIN_ATOMS_FOR_DISPLACEMENT = 4`
- `compute_displacement_bonds()`: Faulon equations 7-11 implementation
  - Equation 10: `b11_min = max(0, a11-a22, a11-a12, a11+a12-3, a11+a21-3)`
  - Equation 11: `b11_max = min(3, a11+a12, a11+a21, a11-a22+3)`
  - Returns None if no valid range
  - Computes b12, b21, b22 via equations 7-9
  - Safety check: all bond orders in [0, 3]
- `attempt_displacement()`:
  - Select 4 distinct atoms via `rng.select_n_distinct(4, atom_count)`
  - Read current bond orders (a11, a12, a21, a22)
  - Compute new bonds via equations
  - Clone graph and apply 4 bond changes
  - Each `set_bond()` call may raise ValueError if SanitizeMol fails
  - Validate: connected AND valid valences
  - Return new graph or None
- All 15 test classes passed

**Initial Structure Generation (RED-GREEN-REFACTOR):**

**RED Phase:**
- Created `backend/tests/test_initial_structure.py` with 13 tests ported from `initialStructure.test.ts`
- Tests covered: saturated alkanes (hexane, butane, methane), unsaturation (HDI=1, HDI=4), heteroatoms (oxygen, nitrogen), simple unsaturated (ethene, ethyne), error handling (empty, impossible), determinism, edge cases (H2)
- Verified tests failed (no implementation)

**GREEN Phase:**
- Created `backend/app/core/initial_structure.py`: Port of `generateInitialStructure`
- Algorithm:
  1. Parse formula via `parse_formula()` and `compute_hdi()`
  2. Create RWMol with heavy atoms in deterministic order: C, N, O, S, P, F, Cl, Br, I
  3. Connect atoms in linear chain with single bonds
  4. Add unsaturation: upgrade bonds iteratively until HDI satisfied
  5. Sanitize and wrap in MoleculeGraph
  6. Verify implicit H count matches formula
- Edge cases: no heavy atoms (H2), single atom (CH4)
- Unsaturation strategy: Try upgrading each bond (1→2 or 2→3), check valences, upgrade if valid
- All 13 tests passed

**Verification:**
```bash
cd backend && poetry run pytest tests/test_displacement.py tests/test_initial_structure.py -v
28/28 tests passed
```

## Verification Results

All verification criteria from plan met:

1. `cd backend && poetry run pytest tests/test_molecule.py -v` -- 18/18 passed
2. `cd backend && poetry run pytest tests/test_wiener.py -v` -- 13/13 passed, Wiener values match: hexane=35, cyclohexane=27, isobutane=9, neopentane=16
3. `cd backend && poetry run pytest tests/test_displacement.py -v` -- 15/15 passed, 500-displacement stress test passes
4. `cd backend && poetry run pytest tests/test_initial_structure.py -v` -- 13/13 passed, all formulas produce valid structures
5. `cd backend && poetry run python -c "from app.core.molecule import MoleculeGraph; m = MoleculeGraph.create_linear_alkane(6); print(m.to_smiles())"` -- prints "CCCCCC" (valid canonical SMILES)

**Total:** 59/59 tests passing

## Deviations from Plan

No deviations. Plan executed exactly as written. RDKit installation completed successfully. All test values match TypeScript reference values.

## Key Decisions Made

1. **RDKit Mol wrapper instead of pure adjacency matrix**
   - TypeScript used pure adjacency matrix: `bonds: number[][]`
   - Python wraps RDKit RWMol: `self._mol = RWMol(Chem.Mol(mol))`
   - Rationale: RDKit provides built-in valence validation, sanitization, distance matrix computation
   - Result: Simpler implementation, automatic implicit H computation, leverages battle-tested chemistry library

2. **SanitizeMol after every bond modification**
   - TypeScript setBond: no validation, just matrix update
   - Python set_bond: calls `Chem.SanitizeMol(self._mol)` after every change
   - Rationale: RDKit requires sanitization to recompute implicit H and validate chemistry
   - Impact: Bond modifications may fail if intermediate states invalid, but guarantees chemical correctness

3. **RDKit GetDistanceMatrix for Wiener Index**
   - TypeScript: 80 lines of custom BFS implementation
   - Python: Single line `dist_matrix = Chem.GetDistanceMatrix(mol)`
   - Rationale: Don't hand-roll graph algorithms when library provides them
   - Result: Cleaner code, faster execution, identical results

## Technical Insights

1. **RDKit implicit hydrogen handling**
   - TypeScript: Manual calculation `implicitH = standardValence - bondOrderSum`
   - Python: RDKit computes automatically via `GetNumImplicitHs()` after sanitization
   - Key: Must call SanitizeMol first, then implicit H is correct

2. **RDKit bond manipulation**
   - `GetBondBetweenAtoms(i, j)` returns None if no bond
   - To remove bond: `RemoveBond(i, j)` (NOT set to UNSPECIFIED)
   - To add bond: `AddBond(i, j, bond_type)`
   - To update bond: `bond.SetBondType(bond_type)`
   - After ANY modification: MUST call SanitizeMol

3. **Displacement intermediate state handling**
   - Challenge: RDKit validates after EACH bond modification
   - Solution: Apply all 4 bond changes, let SanitizeMol fail gracefully
   - If any set_bond raises ValueError, catch and return None
   - Final validation: connected AND valid valences

4. **Initial structure unsaturation strategy**
   - Linear chain construction: Always valid starting point
   - Upgrade bonds one at a time: Check valences before upgrading
   - Each upgrade adds 1 HDI: double bond = +1, triple bond = +2
   - Stop when HDI satisfied or no more upgrades possible

## Success Criteria Met

- [x] MoleculeGraph wraps RDKit Mol with clean API matching TypeScript MolGraph interface
- [x] Faulon displacement on RDKit RWMol preserves valence rules (SanitizeMol validates)
- [x] Wiener Index via GetDistanceMatrix matches all v1.0 reference values
- [x] Initial structure from formula is connected, valid, and deterministic
- [x] Canonical SMILES produced via RDKit MolToSmiles
- [x] All ported test suites pass (59/59)

## Artifacts Produced

**Core Modules (5 files, 687 LOC):**
- `backend/app/utils/rdkit_helpers.py` - Bond type/order conversion (59 LOC)
- `backend/app/core/molecule.py` - MoleculeGraph RDKit wrapper (307 LOC)
- `backend/app/core/wiener.py` - Wiener Index via GetDistanceMatrix (43 LOC)
- `backend/app/core/displacement.py` - Faulon displacement equations (157 LOC)
- `backend/app/core/initial_structure.py` - Initial structure generation (221 LOC)

**Test Suite (4 files, 810 LOC):**
- `backend/tests/test_molecule.py` - 18 tests for MoleculeGraph (200 LOC)
- `backend/tests/test_wiener.py` - 13 tests for Wiener Index (141 LOC)
- `backend/tests/test_displacement.py` - 15 test classes for displacement (260 LOC)
- `backend/tests/test_initial_structure.py` - 13 tests for initial structure (209 LOC)

## Dependencies Created

**Requires:**
- Plan 02: SeededRandom (for displacement atom selection)
- Plan 02: formula_parser (for initial structure generation)
- RDKit 2025.09.5 (Python chemistry library)

**Provides:**
- MoleculeGraph: RDKit wrapper with adjacency-matrix-like API
- Faulon displacement: Core SA mutation operation
- Wiener Index: Topological descriptor for SA cost function
- Initial structure: Starting point for SA optimization

**Affects:**
- Plan 04: SAEngine can now use all molecular operations for SA implementation

## Critical Validation

**Wiener Index Validation (v1.0 reference values):**
- Linear hexane: 35 ✓ (matches v1.0)
- Cyclohexane: 27 ✓ (matches v1.0)
- Isobutane: 9 ✓ (matches v1.0)
- Neopentane: 16 ✓ (matches v1.0)

**Displacement Stress Test:**
- 500 consecutive displacements on hexane
- Result: Zero invalid molecules ✓
- All results connected AND valid valences ✓
- Confirms: TypeScript-to-RDKit translation is chemically correct

## Next Steps

1. Plan 04 will implement SAEngine using these molecular operations
2. SAEngine will use: MoleculeGraph (structure), attempt_displacement (mutation), compute_wiener_index (cost), generate_initial_structure (initialization), SeededRandom (reproducibility), compute_temperature (cooling)

## Self-Check: PASSED

**Created files verified:**
- FOUND: backend/app/utils/rdkit_helpers.py
- FOUND: backend/app/core/molecule.py
- FOUND: backend/app/core/wiener.py
- FOUND: backend/app/core/displacement.py
- FOUND: backend/app/core/initial_structure.py
- FOUND: backend/tests/test_molecule.py
- FOUND: backend/tests/test_wiener.py
- FOUND: backend/tests/test_displacement.py
- FOUND: backend/tests/test_initial_structure.py

**Commits verified:**
- FOUND: 931c82a (feat: implement RDKit molecule wrapper and Wiener Index)
- FOUND: 9d420bd (feat: implement Faulon displacement and initial structure generation)

**Test results verified:**
- All 59 tests passing
- Wiener Index values match v1.0: hexane=35, cyclohexane=27, isobutane=9, neopentane=16
- 500-displacement stress test produces zero invalid molecules
- Initial structure generation produces valid connected molecules with correct hydrogen counts
- SMILES generation works: hexane → "CCCCCC"
