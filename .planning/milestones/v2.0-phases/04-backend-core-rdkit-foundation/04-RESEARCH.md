# Phase 4: Backend Core & RDKit Foundation - Research

**Researched:** 2026-02-16
**Domain:** FastAPI + Python RDKit backend core for molecular SA optimization
**Confidence:** HIGH

## Summary

Phase 4 establishes the Python backend foundation by porting the TypeScript SA engine to Python, integrating native RDKit for molecular operations, and creating the FastAPI application skeleton. This phase is backend-only—no frontend changes, no API endpoints yet, no SSE streaming. The goal is to validate that Python RDKit can faithfully execute Faulon displacement on native molecule objects while maintaining chemical correctness.

The critical technical challenge is translating TypeScript's adjacency matrix `MolGraph` approach to Python RDKit's `RWMol` (editable molecule) API. RDKit requires explicit sanitization after bond modifications, uses a non-obvious BondType enum, and has C++ memory management implications. The existing v1.0 test suite (8 test files, 100+ unit tests) provides ground truth for validation: Python implementation must match TypeScript results for Wiener Index, displacement validity, and SMILES generation.

**Primary recommendation:** Port SAEngine, MolGraph (as RDKit wrapper), displacement, and Wiener Index first. Write pytest unit tests that mirror the TypeScript test suite. Validate Python results match TypeScript results before proceeding to Phase 5 (API layer).

## Standard Stack

### Core

| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| FastAPI | 0.129.0+ | ASGI web framework | Industry standard for async Python APIs in 2026. Native async/await, automatic OpenAPI docs, 5-50x faster than Flask. Requires Python 3.10+. |
| Python RDKit | 2025.09.5+ | Molecular informatics toolkit | Authoritative open-source cheminformatics library. Provides `Chem.Mol`, `Chem.RWMol`, `GetDistanceMatrix()` for Wiener Index, `MolDraw2DSVG` for rendering. |
| Pydantic | v2.x | Data validation & serialization | Required by FastAPI 0.128.0+ (v1 removed). 5-50x faster than v1 due to Rust core. Runtime type checking for request/response models. |
| Uvicorn | 0.34.0+ | ASGI server | Lightning-fast ASGI server for running FastAPI. Single-process dev, multi-worker production. Use `uvicorn[standard]` for HTTP/2. |

### Supporting

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| pytest | 8.x | Testing framework | Python standard. Use with `pytest-asyncio` for async tests. Required for validating Python SA matches TypeScript results. |
| pytest-asyncio | 0.24.x | Async test support | Enables `async def` test functions. Set `asyncio_default_fixture_loop_scope = "function"` in `pyproject.toml`. |
| pytest-cov | 6.x | Code coverage | Built on coverage.py 7.13.4. Run `pytest --cov=app --cov-report=html` for reports. |
| Ruff | 0.8.x | Linting & formatting | Fast Rust-based linter/formatter. Replaces Black + isort + Flake8. Configure in `pyproject.toml`. |
| mypy | 1.14.x | Static type checking | Verify type hints match Pydantic models. Configure strict mode. Catches type errors before runtime. |
| Poetry | 1.8.x | Dependency management | 2026 Python standard. Uses `pyproject.toml`, generates `poetry.lock` for reproducible builds. Run `poetry config virtualenvs.in-project true`. |

### Alternatives Considered

| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| FastAPI | Flask + Flask-RESTX | Flask has larger ecosystem but no native async. FastAPI 5-50x faster for async workloads. |
| Poetry | pip + requirements.txt | Poetry has learning curve but superior dependency management (lockfiles, prod/dev separation). |
| RDKit conda | RDKit pip wheels | Conda better for binary dependencies (recommended). Pip faster for CI/CD if available. |
| pytest | unittest | pytest is 2026 standard. unittest is stdlib but more verbose. Use pytest unless stdlib-only constraint. |

**Installation:**

```bash
# Backend setup with Poetry
cd backend/
poetry init --name webfaulon-backend --python "^3.10"

# Core dependencies
poetry add fastapi uvicorn[standard] pydantic pydantic-settings

# RDKit (RECOMMENDED: conda for development, pip for CI/CD)
# Option 1: Conda
conda install -c conda-forge rdkit python=3.12

# Option 2: Pip (if conda unavailable)
poetry add rdkit

# Dev dependencies
poetry add --group dev pytest pytest-asyncio pytest-cov ruff mypy

# Install all
poetry install
```

## Architecture Patterns

### Recommended Project Structure

```
backend/
├── app/
│   ├── __init__.py
│   ├── main.py              # FastAPI app, CORS, startup/shutdown
│   ├── config.py            # Pydantic settings from .env
│   ├── models/              # Pydantic request/response models
│   │   ├── __init__.py
│   │   ├── sa_params.py     # SAParams, SAResult
│   │   └── progress.py      # ProgressEvent (for Phase 5)
│   ├── core/                # Core SA implementation
│   │   ├── __init__.py
│   │   ├── sa_engine.py     # Python port of SAEngine.ts
│   │   ├── molecule.py      # RDKit Mol wrapper (port of MolGraph.ts)
│   │   ├── displacement.py  # Faulon displacement (port of displacement.ts)
│   │   ├── cooling.py       # Temperature schedules
│   │   ├── wiener.py        # Wiener Index via RDKit
│   │   ├── formula_parser.py # Parse "C6H14" to atom counts
│   │   └── random.py        # Seeded RNG (port of random.ts)
│   └── utils/
│       ├── __init__.py
│       └── rdkit_helpers.py # RDKit BondType mapping, sanitization utilities
├── tests/
│   ├── __init__.py
│   ├── conftest.py          # Pytest fixtures
│   ├── test_molecule.py     # RDKit Mol wrapper tests (port of MolGraph.test.ts)
│   ├── test_displacement.py # Faulon displacement tests
│   ├── test_wiener.py       # Wiener Index tests
│   ├── test_sa_engine.py    # SAEngine tests
│   ├── test_formula_parser.py
│   ├── test_cooling.py
│   └── test_random.py
├── pyproject.toml           # Poetry config, dependencies, tool config
├── poetry.lock              # Locked dependencies
├── .env                     # Local environment variables (gitignored)
└── README.md
```

### Pattern 1: RDKit Mol Wrapper (Python equivalent of TypeScript MolGraph)

**What:** Encapsulate RDKit `Chem.RWMol` with adjacency-matrix-like API to minimize porting friction.

**When to use:** When porting TypeScript code that expects `MolGraph` interface.

**Example:**

```python
# app/core/molecule.py
from rdkit import Chem
from rdkit.Chem import rdMolOps
from typing import List, Tuple

class MoleculeGraph:
    """RDKit wrapper providing adjacency-matrix-like API."""

    def __init__(self, mol: Chem.Mol):
        """Initialize from RDKit Mol (immutable)."""
        self.mol = Chem.Mol(mol)  # Deep copy

    def get_atom_count(self) -> int:
        return self.mol.GetNumAtoms()

    def get_bond_order(self, i: int, j: int) -> int:
        """Get bond order between atoms i and j (0 if no bond)."""
        bond = self.mol.GetBondBetweenAtoms(i, j)
        if bond is None:
            return 0
        return int(bond.GetBondTypeAsDouble())

    def set_bond(self, i: int, j: int, order: int) -> None:
        """Set bond order between atoms i and j (modifies molecule)."""
        # Convert to RWMol for editing
        rw_mol = Chem.RWMol(self.mol)

        bond = rw_mol.GetBondBetweenAtoms(i, j)
        bond_type = order_to_bond_type(order)

        if order == 0:
            # Remove bond
            if bond is not None:
                rw_mol.RemoveBond(i, j)
        elif bond is None:
            # Add new bond
            rw_mol.AddBond(i, j, bond_type)
        else:
            # Modify existing bond
            bond.SetBondType(bond_type)

        # CRITICAL: Sanitize after modification
        try:
            Chem.SanitizeMol(rw_mol)
            self.mol = rw_mol.GetMol()
        except Exception as e:
            raise ValueError(f"Invalid molecule after bond modification: {e}")

    def is_connected(self) -> bool:
        """Check if molecule is a single connected component."""
        frags = Chem.GetMolFrags(self.mol)
        return len(frags) == 1

    def has_valid_valences(self) -> bool:
        """Check if all atoms have valid valences."""
        try:
            Chem.SanitizeMol(self.mol)
            return True
        except Exception:
            return False

    def clone(self) -> 'MoleculeGraph':
        """Return deep copy."""
        return MoleculeGraph(self.mol)

    def to_smiles(self) -> str:
        """Generate canonical SMILES."""
        return Chem.MolToSmiles(self.mol)
```

### Pattern 2: Explicit BondType Mapping

**What:** Map Faulon integer bond orders (0=none, 1=single, 2=double, 3=triple) to RDKit BondType enum.

**When to use:** Every time setting bond orders from Faulon displacement.

**Example:**

```python
# app/utils/rdkit_helpers.py
from rdkit.Chem import BondType

def order_to_bond_type(order: int) -> BondType:
    """
    Convert Faulon integer bond order to RDKit BondType enum.

    Faulon: 0=none, 1=single, 2=double, 3=triple
    RDKit: UNSPECIFIED=0, SINGLE=1, DOUBLE=2, TRIPLE=3, AROMATIC=12
    """
    mapping = {
        0: BondType.UNSPECIFIED,
        1: BondType.SINGLE,
        2: BondType.DOUBLE,
        3: BondType.TRIPLE,
    }
    if order not in mapping:
        raise ValueError(f"Invalid bond order: {order}")
    return mapping[order]

def bond_type_to_order(bond_type: BondType) -> int:
    """Convert RDKit BondType to Faulon integer bond order."""
    mapping = {
        BondType.UNSPECIFIED: 0,
        BondType.SINGLE: 1,
        BondType.DOUBLE: 2,
        BondType.TRIPLE: 3,
    }
    if bond_type not in mapping:
        raise ValueError(f"Unsupported bond type: {bond_type}")
    return mapping[bond_type]
```

### Pattern 3: Mandatory Sanitization After Bond Modification

**What:** Always call `Chem.SanitizeMol()` after modifying bond orders. Wrap in try/except to reject invalid moves.

**When to use:** After every `SetBondType()`, `AddBond()`, or `RemoveBond()` operation.

**Example:**

```python
# app/core/displacement.py
def apply_displacement(
    mol: Chem.Mol,
    x1: int, y1: int, x2: int, y2: int,
    b11: int, b12: int, b21: int, b22: int
) -> Chem.Mol | None:
    """
    Apply Faulon displacement bond orders to RDKit molecule.
    Returns new molecule if valid, None if invalid.
    """
    rw_mol = Chem.RWMol(mol)

    # Apply new bond orders
    set_bond_order(rw_mol, x1, y1, b11)
    set_bond_order(rw_mol, y1, y2, b12)
    set_bond_order(rw_mol, x1, x2, b21)
    set_bond_order(rw_mol, x2, y2, b22)

    # CRITICAL: Validate with sanitization
    try:
        new_mol = rw_mol.GetMol()
        Chem.SanitizeMol(new_mol)

        # Additional check: must be connected
        frags = Chem.GetMolFrags(new_mol)
        if len(frags) != 1:
            return None

        return new_mol
    except Exception:
        # Invalid chemistry (wrong valence, etc.)
        return None
```

### Anti-Patterns to Avoid

- **Skipping sanitization after displacement:** RDKit does NOT auto-validate. Produces silently invalid molecules.
- **Storing Mol objects in history:** Memory leaks due to C++ object retention. Store SMILES strings instead.
- **Direct integer to BondType:** `BondType` enum is not zero-indexed. Use explicit mapping function.
- **Running SA in async event loop:** SA is CPU-bound. Will block FastAPI. Use BackgroundTasks or ProcessPoolExecutor (Phase 5).

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Wiener Index computation | Nested loops over adjacency matrix | `Chem.GetDistanceMatrix(mol)` + sum upper triangle | RDKit uses optimized graph algorithms, handles edge cases (disconnected fragments, aromatic rings). |
| SMILES generation | DFS with ring closure tracking | `Chem.MolToSmiles(mol)` | RDKit produces canonical SMILES with correct aromaticity perception, stereochemistry. Custom DFS is unreliable (v1.0 lesson). |
| Molecular formula parsing | Regex parsing | Use existing formula parser logic or `Chem.rdMolDescriptors.CalcMolFormula()` | Edge cases: capitalization, implicit hydrogen, isotopes. |
| Bond validation | Manual valence checking | `Chem.SanitizeMol()` | RDKit checks valence, aromaticity, kekulization, hybridization. Hand-rolled validation misses edge cases. |
| Initial structure generation | Random bond assignment | Port TypeScript linear chain logic or use RDKit's `AllChem.EmbedMolecule()` | Linear chain is deterministic and matches v1.0 behavior. Random graphs may fail sanitization. |

**Key insight:** RDKit is the authoritative cheminformatics library. When in doubt, delegate to RDKit rather than implementing custom graph algorithms.

## Common Pitfalls

### Pitfall 1: Missing Sanitization After RDKit Bond Order Changes

**What goes wrong:**

After performing Faulon displacement (modifying bond orders on 4 atoms using `SetBondType`), failing to call `Chem.SanitizeMol()` results in molecules with incorrect valence states, broken aromaticity perception, and invalid hydrogen counts. RDKit will generate incorrect SMILES strings and may crash on subsequent operations.

**Why it happens:**

Developers assume RDKit auto-validates like TypeScript `MolGraph` (which recomputes `implicitH` in `setBond()`). In RDKit, bond modifications create an intermediate "dirty" state—valence, aromaticity, and hybridization are NOT automatically updated.

**How to avoid:**

```python
from rdkit import Chem

# After any bond order modification
rw_mol = Chem.RWMol(mol)
rw_mol.GetBondBetweenAtoms(x1, y1).SetBondType(Chem.BondType.DOUBLE)

# CRITICAL: Sanitize to recompute valence, aromaticity, etc.
try:
    Chem.SanitizeMol(rw_mol)
    new_mol = rw_mol.GetMol()
except Exception as e:
    # Displacement created invalid molecule (wrong valence)
    return None  # Reject this move
```

**Warning signs:**

- `Explicit valence for atom #X Y, Z, is greater than permitted` errors
- SMILES strings differ from expected (e.g., `C1CC1` instead of `C=CC`)
- Test failures with message "Sanitization error"
- Kekulization errors on aromatic systems

**Validation:**

Unit test asserts SMILES correctness after displacement matches TypeScript reference values.

### Pitfall 2: Incorrect RDKit BondType Enum Mapping

**What goes wrong:**

After Faulon displacement calculates new bond orders as integers (0, 1, 2, 3), naively passing them to RDKit's `BondType` results in wrong bond types. Single bonds become aromatic, double bonds become single, etc.

**Why it happens:**

RDKit's `BondType` enum is NOT zero-indexed integers:
- `BondType.UNSPECIFIED = 0`
- `BondType.SINGLE = 1`
- `BondType.DOUBLE = 2`
- `BondType.TRIPLE = 3`
- `BondType.AROMATIC = 12`

But Faulon displacement uses:
- `0` = no bond
- `1` = single bond
- `2` = double bond
- `3` = triple bond

**How to avoid:**

```python
from rdkit.Chem import BondType

def faulon_order_to_rdkit_bondtype(order: int) -> BondType:
    """Convert Faulon integer bond order to RDKit BondType enum."""
    mapping = {
        0: BondType.UNSPECIFIED,  # Or remove bond entirely
        1: BondType.SINGLE,
        2: BondType.DOUBLE,
        3: BondType.TRIPLE,
    }
    if order not in mapping:
        raise ValueError(f"Invalid bond order: {order}")
    return mapping[order]
```

**Warning signs:**

- SMILES strings don't match expected structures (e.g., `C1=CC1` instead of `C1CC1`)
- Aromatic bonds appear in saturated alkanes
- Bond order validation failures

**Validation:**

Test all bond types (no bond, single, double, triple) produce expected SMILES.

### Pitfall 3: RDKit Molecule Memory Leaks in Long-Running SA

**What goes wrong:**

During thousands of SA iterations, memory usage grows unbounded. Python's GC doesn't immediately release C++ RDKit objects.

**Why it happens:**

RDKit's Python bindings use C++ objects with complex lifetimes. Storing `Mol` objects in lists (e.g., `history: List[Mol]`) creates circular references or delayed GC.

**How to avoid:**

```python
# BAD: Stores thousands of Mol objects
class SAEngine:
    def __init__(self):
        self.history: List[Mol] = []

    def step(self):
        proposed_mol = self.current_mol.clone()
        self.history.append(proposed_mol)  # Memory leak

# GOOD: Store SMILES strings (lightweight, serializable)
class SAEngine:
    def __init__(self):
        self.history: List[dict] = []

    def step(self):
        self.history.append({
            "step": self.step_num,
            "energy": self.current_energy,
            "smiles": Chem.MolToSmiles(self.current_mol)  # String, not Mol
        })
```

**Warning signs:**

- Server memory grows with each completed job (doesn't return to baseline)
- Memory profiling shows `rdkit.Chem` objects dominating

**Validation:**

Assert memory returns to baseline after `run()` completes.

### Pitfall 4: CPU-Bound SA Loop Blocking Async FastAPI (DEFERRED TO PHASE 5)

**What goes wrong:**

Running the SA loop (thousands of `engine.step()` calls) inside an `async def` endpoint blocks the entire event loop. FastAPI becomes unresponsive.

**Why it happens:**

FastAPI's async is I/O concurrency, not parallelism. SA loop is CPU-bound (RDKit operations).

**How to avoid (Phase 5):**

Run SA in BackgroundTasks or ProcessPoolExecutor. Not applicable to Phase 4 (no API endpoints yet).

## Code Examples

Verified patterns from official sources:

### Wiener Index via RDKit

```python
# Source: https://www.rdkit.org/docs/Cookbook.html
from rdkit import Chem

def compute_wiener_index(mol: Chem.Mol) -> float:
    """
    Compute Wiener index (sum of pairwise distances).

    RDKit GetDistanceMatrix returns NxN matrix of shortest path lengths.
    """
    dist_matrix = Chem.GetDistanceMatrix(mol)
    n = mol.GetNumAtoms()

    # Sum upper triangle (avoid double-counting)
    wiener = sum(
        dist_matrix[i, j]
        for i in range(n)
        for j in range(i + 1, n)
    )

    return wiener
```

### Initial Structure from Formula

```python
# Port of TypeScript initialStructure.ts
from rdkit import Chem
from rdkit.Chem import AllChem

def generate_initial_structure(formula: str) -> Chem.Mol:
    """
    Create initial linear chain structure from molecular formula.

    Matches v1.0 deterministic behavior (not random).
    """
    # Parse formula to atom counts
    atom_counts = parse_formula(formula)

    # Create RWMol
    mol = Chem.RWMol()

    # Add heavy atoms in linear chain
    for element, count in atom_counts.items():
        if element == 'H':
            continue  # RDKit handles implicit hydrogens
        for _ in range(count):
            mol.AddAtom(Chem.Atom(element))

    # Add single bonds in linear chain
    for i in range(mol.GetNumAtoms() - 1):
        mol.AddBond(i, i + 1, Chem.BondType.SINGLE)

    # Sanitize and return
    result = mol.GetMol()
    Chem.SanitizeMol(result)

    return result
```

### Faulon Displacement on RDKit Mol

```python
# Port of TypeScript displacement.ts
from rdkit import Chem
from typing import Optional

def attempt_displacement(
    mol: Chem.Mol,
    rng: SeededRandom
) -> Optional[Chem.Mol]:
    """
    Attempt Faulon displacement (equations 7-11).

    Returns new molecule if valid, None if invalid.
    """
    atom_count = mol.GetNumAtoms()

    if atom_count < 4:
        return None

    # Select 4 distinct atoms
    x1, y1, x2, y2 = rng.select_n_distinct(4, atom_count)

    # Read current bond orders
    a11 = get_bond_order(mol, x1, y1)
    a12 = get_bond_order(mol, y1, y2)
    a21 = get_bond_order(mol, x1, x2)
    a22 = get_bond_order(mol, x2, y2)

    # Compute new bond orders (equations 7-11)
    new_bonds = compute_displacement_bonds(a11, a12, a21, a22, rng)

    if new_bonds is None:
        return None

    # Apply displacement
    return apply_displacement(mol, x1, y1, x2, y2, **new_bonds)

def get_bond_order(mol: Chem.Mol, i: int, j: int) -> int:
    """Get bond order between atoms i and j (0 if no bond)."""
    bond = mol.GetBondBetweenAtoms(i, j)
    if bond is None:
        return 0
    return int(bond.GetBondTypeAsDouble())
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| TypeScript adjacency matrix | Python RDKit native Mol objects | Phase 4 (2026) | Native RDKit enables full descriptor library (200+ properties), correct SMILES generation, production-ready chemistry. |
| Custom SMILES generation (DFS) | `Chem.MolToSmiles()` | Phase 4 (2026) | Eliminates custom graph traversal code (100+ LOC), produces canonical SMILES with correct aromaticity. |
| Web Worker isolation | Backend process isolation | Phase 4-5 (2026) | Enables multi-component scoring, future NMR prediction, deployment to cloud. |
| RDKit.js WASM (2025.3.4) | Python RDKit 2025.09.5 | Phase 4 (2026) | Latest RDKit features, better performance, no WASM memory limits. |
| Pydantic v1 | Pydantic v2 | FastAPI 0.128.0 (2024) | 5-50x faster validation due to Rust core. v1 removed from FastAPI. |

**Deprecated/outdated:**

- TypeScript `MolGraph` class: Replaced by Python `MoleculeGraph` wrapper around `Chem.RWMol`
- Custom SMILES generation: Replaced by `Chem.MolToSmiles()`
- RDKit.js WASM: Removed from frontend in Phase 6

## Open Questions

1. **Faulon displacement on RDKit RWMol validation**
   - What we know: TypeScript implementation uses adjacency matrix. Python needs to use RDKit RWMol API.
   - What's unclear: Does RDKit bond manipulation preserve Faulon's 4-atom displacement semantics? Will sanitization reject valid moves?
   - Recommendation: Port displacement.ts to Python with RDKit API, run against full TypeScript test suite (displacement.test.ts), compare results. If mismatch, investigate sanitization flags or fall back to adjacency matrix with RDKit only for rendering.

2. **Memory management for SA history**
   - What we know: Storing Mol objects causes memory leaks. Storing SMILES strings is safe.
   - What's unclear: Does v2.0 need full history for replay? Or only best molecule?
   - Recommendation: Phase 4 stores only best molecule (SMILES + Mol object). History storage deferred to Phase 5 (API layer decides what to persist).

3. **Formula parser edge cases**
   - What we know: TypeScript formula parser handles simple cases (C6H14, C10H22).
   - What's unclear: Does it handle isotopes, charged species, implicit multipliers?
   - Recommendation: Port TypeScript formula parser as-is. Document limitations. RDKit's `CalcMolFormula()` available for validation but not primary input method.

## Sources

### Primary (HIGH confidence)

**Python RDKit:**
- [RDKit Getting Started](https://www.rdkit.org/docs/GettingStartedInPython.html) — Python API
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html) — Wiener Index examples
- [RDKit Installation](https://www.rdkit.org/docs/Install.html) — Conda vs pip
- [RDKit rdchem module](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html) — BondType enum
- [RDKit Sanitization](https://greglandrum.github.io/rdkit-blog/posts/2025-06-27-sanitization-and-file-parsing.html) — Sanitization options and edge cases

**FastAPI:**
- [FastAPI Release Notes](https://fastapi.tiangolo.com/release-notes/) — Version 0.129.0, Python 3.10+ requirement
- [FastAPI Testing](https://fastapi.tiangolo.com/tutorial/testing/) — pytest patterns

**Pydantic:**
- [Pydantic v2 Docs](https://docs.pydantic.dev/latest/) — Rust-based validation, migration from v1

### Secondary (MEDIUM confidence)

**RDKit Pitfalls:**
- [Memory Leakage with Molecule Objects — GitHub #3239](https://github.com/rdkit/rdkit/issues/3239) — Mol object retention
- [Explicit Valence Error — GitHub Discussion #8181](https://github.com/rdkit/rdkit/discussions/8181) — Sanitization failures
- [Sanitization issues on identical molecules — GitHub #8156](https://github.com/rdkit/rdkit/discussions/8156) — Edge cases

**Python Testing:**
- [pytest-asyncio documentation](https://pytest-asyncio.readthedocs.io/en/latest/) — Async test patterns

**Poetry:**
- [Poetry Documentation](https://python-poetry.org/docs/) — Dependency management

### Tertiary (LOW confidence, validate during implementation)

- Faulon displacement on RDKit RWMol (no existing implementation found — needs validation)
- Memory profiling for C++ object retention (general guidance, not RDKit-specific)

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH — FastAPI 0.129.0, RDKit 2025.09.5, Pydantic v2, pytest all verified via official docs
- Architecture patterns: HIGH — RDKit wrapper, BondType mapping, sanitization documented in official RDKit guides
- Pitfalls: HIGH — All pitfalls sourced from RDKit official docs, GitHub issues, and v2.0 project research

**Research date:** 2026-02-16
**Valid until:** 60 days (stable domain, mature libraries)
