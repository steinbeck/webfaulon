---
phase: 07-multi-component-target-function
plan: 02
subsystem: sa-engine
tags: [tdd, multi-component, weighted-scoring, pydantic]

dependency_graph:
  requires:
    - phase: 07-01
      provides: [ScoringComponent protocol, WienerIndexComponent, LogPComponent, ComponentRegistry]
  provides:
    - SAEngine._compute_energy() weighted sum implementation
    - SAParams.component_weights field with validation
    - Multi-component SA optimization capability
    - Backward compatible default weights
  affects: [sa-api, frontend-ui]

tech_stack:
  added: []
  patterns: [Weighted sum scoring, fail-fast validation, zero-weight optimization skip]

key_files:
  created: []
  modified:
    - backend/app/models/sa_params.py
    - backend/app/core/sa_engine.py
    - backend/app/core/scoring/logp.py
    - backend/tests/test_scoring_integration.py

key_decisions:
  - "Default component_weights={'wiener_index': 1.0} ensures perfect backward compatibility"
  - "Zero-weight components are skipped in _compute_energy() for efficiency"
  - "Component validation happens at SAEngine.__init__() for fail-fast errors"
  - "disconnected_moves counted separately in accounting (existing behavior preserved)"
  - "LogP requires UpdatePropertyCache() + SanitizeMol() before MolLogP() call"

patterns_established:
  - "Field validators on Pydantic models for complex validation logic"
  - "Registry validation at constructor time with helpful error messages"
  - "Weighted sum pattern extensible to future scoring components"

metrics:
  duration: 186s
  completed: 2026-02-16T16:26:26Z
---

# Phase 07 Plan 02: Multi-Component Target Function Summary

**SAEngine refactored for pluggable multi-component scoring with weighted sums, replacing hardcoded Wiener Index while maintaining full backward compatibility via default weights**

## Performance

- **Duration:** 3 min 6 sec (186s)
- **Started:** 2026-02-16T16:23:20Z
- **Completed:** 2026-02-16T16:26:26Z
- **Tasks:** 2 (TDD: RED → GREEN)
- **Files modified:** 4
- **Tests added:** 14 integration tests

## Accomplishments

- SAEngine now uses weighted sum of pluggable scoring components instead of hardcoded `compute_wiener_index()`
- SAParams extended with `component_weights: Dict[str, float]` field with comprehensive validation
- Default `component_weights={'wiener_index': 1.0}` produces identical results to original engine (backward compatible)
- Multi-component optimization (e.g., `{'wiener_index': 0.5, 'logp': 0.5}`) works correctly
- All 238 tests pass including 59 existing SA engine tests (zero regressions)

## Task Commits

Each task was committed atomically following TDD workflow:

1. **Task 1: RED - Write failing tests** - `23c694b` (test)
   - 14 integration tests across 4 test classes
   - SAParams validation tests
   - SAEngine weighted energy tests
   - Component registry validation tests
   - Backward compatibility tests

2. **Task 2: GREEN - Implement multi-component scoring** - `af0a286` (feat)
   - Extended SAParams with component_weights field and validator
   - Refactored SAEngine to use _compute_energy() with registry
   - Fixed LogP component to sanitize molecules before MolLogP()
   - Fixed test accounting assertions to include disconnected_moves
   - All 238 tests pass

## Files Created/Modified

- `backend/app/models/sa_params.py` - Added `component_weights: Dict[str, float]` field with validator (no negatives, no all-zero, non-empty)
- `backend/app/core/sa_engine.py` - Added `_compute_energy()` method, replaced hardcoded `compute_wiener_index()` calls, added registry validation in `__init__()`
- `backend/app/core/scoring/logp.py` - Added `UpdatePropertyCache()` and `SanitizeMol()` before `MolLogP()` computation (fixes RDKit implicit valence requirement)
- `backend/tests/test_scoring_integration.py` - 14 integration tests for multi-component SA engine

## Decisions Made

1. **Default weights for backward compatibility**: `component_weights={'wiener_index': 1.0}` ensures existing code without component_weights produces identical results to the original hardcoded engine.

2. **Zero-weight optimization**: Components with `weight == 0.0` are skipped in `_compute_energy()` to avoid unnecessary computation.

3. **Fail-fast validation**: Component name validation happens at `SAEngine.__init__()` rather than during execution, providing immediate feedback with helpful error messages listing available components.

4. **disconnected_moves accounting**: Preserved existing behavior where disconnected_moves are counted separately from invalid_moves (not summed into rejected_moves). Tests updated to reflect this.

5. **LogP sanitization requirement**: RDKit's `MolLogP()` requires implicit valences to be calculated, so LogPComponent now calls `UpdatePropertyCache()` and `SanitizeMol()` before computing LogP.

## Deviations from Plan

### Auto-fixed Issues

**1. [Rule 1 - Bug] Fixed LogP component to sanitize molecules before MolLogP()**
- **Found during:** Task 2 (GREEN phase test run)
- **Issue:** RDKit's `MolLogP()` raised "Pre-condition Violation: getNumImplicitHs() called without preceding call to calcImplicitValence()" when computing LogP during SA run
- **Root cause:** `MoleculeGraph.get_mol()` returns internal RWMol which may not have implicit valences calculated
- **Fix:** Added `mol.UpdatePropertyCache(strict=False)` and `Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)` before `Descriptors.MolLogP(mol)` call
- **Files modified:** `backend/app/core/scoring/logp.py`
- **Verification:** All tests pass including multi-component SA runs with LogP
- **Committed in:** af0a286 (Task 2 commit)

**2. [Rule 1 - Bug] Fixed test accounting assertions to include disconnected_moves**
- **Found during:** Task 2 (GREEN phase test run)
- **Issue:** Three tests failed with assertion error: `accepted_moves + rejected_moves + invalid_moves != total_steps`
- **Root cause:** Tests assumed `accepted + rejected + invalid == total`, but existing SAEngine counts `disconnected_moves` separately
- **Fix:** Updated three test assertions to `accepted + rejected + invalid + disconnected == total`
- **Files modified:** `backend/tests/test_scoring_integration.py`
- **Verification:** All 238 tests pass
- **Committed in:** af0a286 (Task 2 commit)

---

**Total deviations:** 2 auto-fixed (2 bugs)
**Impact on plan:** Both fixes necessary for correctness. LogP sanitization is essential for RDKit descriptor computation. Accounting fix preserves existing engine behavior. No scope creep.

## Implementation Details

### SAParams Extension

```python
component_weights: Dict[str, float] = Field(
    default_factory=lambda: {"wiener_index": 1.0},
    description="Weights for scoring components"
)

@field_validator("component_weights")
@classmethod
def validate_weights(cls, v: Dict[str, float]) -> Dict[str, float]:
    if not v:
        raise ValueError("At least one component weight must be specified")
    for name, weight in v.items():
        if weight < 0:
            raise ValueError(f"Component '{name}' has negative weight: {weight}")
    if all(w == 0.0 for w in v.values()):
        raise ValueError("At least one component weight must be non-zero")
    return v
```

### SAEngine Changes

**Registry validation in `__init__()`:**
```python
registry = get_registry()
available = set(registry.list_components())
requested = set(params.component_weights.keys())
unknown = requested - available
if unknown:
    raise ValueError(
        f"Unknown scoring component(s): {unknown}. "
        f"Available: {sorted(available)}"
    )
```

**Weighted energy computation:**
```python
def _compute_energy(self, mol_graph: MoleculeGraph) -> float:
    total_energy = 0.0
    for component_name, weight in self._component_weights.items():
        if weight == 0.0:
            continue
        component = self._registry.get(component_name)
        score = component.compute(mol_graph)
        total_energy += weight * score
    return total_energy
```

**Replaced two calls:**
- `init()`: Changed `self._current_energy = compute_wiener_index(self._current_graph)` to `self._current_energy = self._compute_energy(self._current_graph)`
- `_iterate()`: Changed `proposed_energy = compute_wiener_index(proposed_graph)` to `proposed_energy = self._compute_energy(proposed_graph)`

### Test Coverage

14 integration tests verify:
- **SAParams validation**: default weights, custom weights, negative rejection, all-zero rejection, empty rejection
- **Weighted energy**: weight=1.0 produces 35, weight=2.0 produces 70, two components sum correctly, zero-weight components excluded
- **Registry validation**: unknown components rejected with helpful message listing available components
- **Backward compatibility**: existing code patterns work unchanged, full SA runs complete successfully

All 238 tests pass including:
- 59 existing SA engine tests (unchanged)
- 17 scoring component tests
- 13 Wiener index tests
- 14 new integration tests

## Verification Results

All success criteria met:

1. ✅ SAParams has `component_weights: Dict[str, float]` field with default `{"wiener_index": 1.0}`
2. ✅ SAParams validates: no negative weights, no all-zero weights, no empty weights dict
3. ✅ SAEngine.__init__() validates component_weights against registry (fail-fast)
4. ✅ SAEngine._compute_energy() returns weighted sum of component scores
5. ✅ `compute_wiener_index` import removed from sa_engine.py
6. ✅ All 59 existing SA engine tests pass without modification (backward compatibility)
7. ✅ All 14 new integration tests pass
8. ✅ Full test suite passes (238 tests, zero failures)

## Technical Notes

- **Backward compatibility mechanism**: `default_factory=lambda: {"wiener_index": 1.0}` ensures existing code that doesn't pass `component_weights` gets the original Wiener-only behavior
- **Zero-weight optimization**: Skipping components with weight 0.0 avoids unnecessary computation, especially useful when experimenting with component selection
- **Fail-fast philosophy**: Component name validation at initialization (not during iteration) provides immediate feedback with clear error messages
- **RDKit requirements**: LogP computation requires sanitized molecules with calculated implicit valences - this is now handled automatically in LogPComponent
- **Extensibility**: Adding new scoring components requires zero changes to SAEngine - just register in registry and reference in component_weights

## Next Phase Readiness

Phase 7 (Multi-Component Target Function) is now complete:
- ✅ Plan 01: Scoring component framework (protocol, components, registry)
- ✅ Plan 02: Multi-component SA engine integration

**Ready for Phase 8** (if exists) or any future work requiring:
- Multi-objective optimization with weighted components
- Custom scoring functions (e.g., TPSA, rotatable bonds, ring count)
- User-configurable target functions via API

**No blockers or concerns.** System is production-ready for multi-component optimization.

## Self-Check: PASSED

All modified files exist:
- ✓ backend/app/models/sa_params.py
- ✓ backend/app/core/sa_engine.py
- ✓ backend/app/core/scoring/logp.py
- ✓ backend/tests/test_scoring_integration.py

All commits exist:
- ✓ 23c694b (Task 1: RED)
- ✓ af0a286 (Task 2: GREEN)

All tests pass:
- ✓ 238/238 tests passed
- ✓ Zero regressions in existing test suite
- ✓ Full backward compatibility verified

---
*Phase: 07-multi-component-target-function*
*Plan: 02*
*Completed: 2026-02-16*
