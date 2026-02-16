---
phase: 07-multi-component-target-function
plan: 01
subsystem: scoring-framework
tags: [tdd, protocol, components, registry]
dependency_graph:
  requires: [backend/app/core/wiener.py, backend/app/core/molecule.py]
  provides: [ScoringComponent protocol, WienerIndexComponent, LogPComponent, ComponentRegistry]
  affects: []
tech_stack:
  added: [typing.Protocol with @runtime_checkable]
  patterns: [Protocol-based polymorphism, singleton registry pattern]
key_files:
  created:
    - backend/app/core/scoring/protocol.py
    - backend/app/core/scoring/wiener.py
    - backend/app/core/scoring/logp.py
    - backend/app/core/scoring/registry.py
    - backend/tests/test_scoring_components.py
  modified:
    - backend/app/core/scoring/__init__.py
decisions:
  - decision: "Replace MolecularWeightComponent with LogPComponent"
    rationale: "Molecular weight is constant across all constitutional isomers of a given formula, making it useless for SA optimization. LogP (Wildman-Crippen partition coefficient) varies with atom connectivity and is practically relevant for drug-likeness."
    impact: "Foundation for multi-component optimization now includes a property that actually varies during SA"
metrics:
  duration: 147s
  completed: 2026-02-16T16:20:49Z
---

# Phase 07 Plan 01: Scoring Component Framework Summary

**One-liner:** Protocol-based scoring framework with Wiener Index and LogP components, validated registry with duplicate detection and helpful error messages.

## Tasks Completed

| Task | Type | Commit | Description |
|------|------|--------|-------------|
| 1 | RED | 1969cc0 | Failing tests for components and registry (17 tests) |
| 2 | GREEN | 2a77f8c | Implement protocol, components, registry - all tests pass |

## What Was Built

### ScoringComponent Protocol
- Structural interface using `typing.Protocol` with `@runtime_checkable`
- Requires: `name: str` attribute and `compute(mol_graph) -> float` method
- Enables duck-typed polymorphism without inheritance

### WienerIndexComponent
- Wraps existing `compute_wiener_index()` algorithm into component class
- Identical implementation: RDKit GetDistanceMatrix + upper triangle sum
- Backward compatibility verified: component returns same values as standalone function
- Handles edge cases: single atom (0), disconnected graphs (ValueError)

### LogPComponent
- Wildman-Crippen LogP via `RDKit Descriptors.MolLogP()`
- Varies across constitutional isomers (unlike molecular weight)
- Positive for hydrocarbons, suitable for drug-likeness optimization
- **Design decision:** Replaced MolecularWeight because MW is constant for all isomers

### ComponentRegistry
- Dictionary-based storage with validation
- `register(component)`: Rejects duplicates with ValueError
- `get(name)`: Raises KeyError with helpful message listing available components
- `list_components()`: Returns all registered names
- Singleton pattern: `get_registry()` returns pre-populated global instance

### Test Coverage
- 17 tests across 4 test classes
- Protocol compliance verification (hasattr, callable checks)
- Known values for multiple test molecules
- Backward compatibility test (component vs standalone function)
- Registry validation (duplicates, unknown components, error messages)
- Edge case: LogP differs between linear and branched isomers

## Deviations from Plan

### Auto-handled Issues (Rule 3)

**1. Pre-existing implementation files with wrong component**
- **Found during:** Task 1 verification
- **Issue:** Planner had created protocol.py, wiener.py, registry.py with MolecularWeightComponent instead of LogPComponent
- **Fix:** Created logp.py, updated registry.py imports and registration to use LogPComponent
- **Files modified:** backend/app/core/scoring/logp.py (created), backend/app/core/scoring/registry.py (imports updated)
- **Commit:** 2a77f8c (incorporated into Task 2)
- **Rationale:** Plan explicitly stated LogP replaced MolecularWeight; pre-existing files blocked test execution

## Verification Results

All success criteria met:

1. ✅ `python -m pytest tests/test_scoring_components.py -v` -- 17/17 tests pass
2. ✅ `python -m pytest tests/test_wiener.py -v` -- 13/13 tests pass (backward compat)
3. ✅ `python -m pytest tests/test_sa_engine.py -v` -- 59/59 tests pass (no regressions)
4. ✅ `WienerIndexComponent().compute(graph)` returns identical values to `compute_wiener_index(graph)` for all test molecules
5. ✅ `ComponentRegistry` validates names and raises clear errors with available component list

## Technical Notes

- **Protocol pattern:** `@runtime_checkable` enables `isinstance(obj, ScoringComponent)` checks at runtime, useful for debugging
- **LogP choice:** Wildman-Crippen method (RDKit default) balances accuracy and speed; varies with connectivity unlike topological descriptors
- **No REFACTOR phase:** Code was clean on first implementation, no cleanup needed
- **Singleton registry:** Module-level `_registry` instance pre-populated with built-in components; `get_registry()` returns reference

## Next Steps

This framework enables Plan 02: Multi-component target function refactoring of SA engine. The registry pattern allows dynamic component selection and future extensibility (e.g., add TPSA, rotatable bonds, ring count without modifying core engine).

## Self-Check: PASSED

All created files exist:
- ✓ backend/app/core/scoring/protocol.py
- ✓ backend/app/core/scoring/wiener.py
- ✓ backend/app/core/scoring/logp.py
- ✓ backend/app/core/scoring/registry.py
- ✓ backend/tests/test_scoring_components.py

All commits exist:
- ✓ 1969cc0 (Task 1: RED)
- ✓ 2a77f8c (Task 2: GREEN)
