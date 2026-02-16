# Phase 7: Multi-Component Target Function - Research

**Researched:** 2026-02-16
**Domain:** Plugin architectures, weighted scoring systems, simulated annealing objective functions
**Confidence:** MEDIUM-HIGH

## Summary

Phase 7 refactors the SA engine to support pluggable scoring components with configurable weights. Currently, the engine hardcodes `compute_wiener_index()` as the sole scoring function (line 226 in sa_engine.py). The goal is to create a framework where multiple scoring components (Wiener Index, molecular weight, LogP, etc.) can contribute weighted scores to the overall objective function.

This is a well-understood software pattern problem. The Strategy pattern provides the architectural foundation for pluggable components, while weighted sum methods are standard in multi-criteria optimization. RDKit's QED implementation demonstrates production-grade weighted scoring with tuples of weights corresponding to individual components. The key challenge is balancing extensibility with YAGNI principles: the framework must support multiple components without over-engineering for the v2.1+ features (adjustable weights via API, component-wise streaming) that are explicitly deferred.

For an educational codebase targeting classroom demonstrations, the design should be simple and transparent. Python's Protocol typing offers duck-typed interfaces without requiring formal inheritance, making component addition straightforward. The weighted sum calculation is trivial (sum of weight * score), but careful consideration is needed for optimization mode handling (MINIMIZE vs MAXIMIZE) and component normalization.

**Primary recommendation:** Use Protocol-based scoring components with a simple registry pattern. Each component implements `compute(mol_graph: MoleculeGraph) -> float` and has an associated weight. The SA engine replaces the hardcoded `compute_wiener_index()` call with a loop over registered components, computes the weighted sum, and uses that as the energy value. Start with two components (Wiener Index with weight 1.0, molecular weight with weight 0.0) to validate the framework without changing behavior, then allow weight configuration in future phases.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| Python typing.Protocol | 3.11+ | Structural subtyping for scoring components | Allows pluggable components without inheritance coupling; type-safe duck typing |
| Pydantic v2 | ^2.x | Configuration models with validation | Already used for SAParams; extend to include component weights |
| RDKit | ^2023.x | Descriptor calculation for components | Already in stack; provides 200+ molecular descriptors for future components |

### Supporting
| Tool | Purpose | When to Use |
|------|---------|-------------|
| dataclasses | Lightweight component metadata | If Pydantic is overkill for component definitions |
| dict[str, ScoringComponent] | Component registry | Simple lookup by name; no need for importlib plugin discovery in v2.0 |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Protocol | ABC (Abstract Base Class) | ABC requires formal inheritance; Protocol allows independent components. Use Protocol unless runtime validation is critical. |
| dict registry | Entry points (importlib.metadata) | Entry points enable third-party plugins but add complexity. YAGNI for classroom tool. |
| Weighted sum | Pareto front (MOSA) | Multi-objective SA approximates Pareto optimal solutions. Deferred to v3.0+; weighted sum sufficient for v2.0 requirements. |

**Installation:**
No new dependencies required. Uses existing Python 3.11+ typing module and RDKit.

## Architecture Patterns

### Recommended Project Structure
```
backend/app/core/
├── sa_engine.py              # SAEngine refactored to use component registry
├── scoring/
│   ├── __init__.py           # Exports ScoringComponent Protocol
│   ├── protocol.py           # Protocol definition + base utilities
│   ├── wiener.py             # WienerIndexComponent (refactored from wiener.py)
│   ├── molecular_weight.py   # MolecularWeightComponent (example second component)
│   └── registry.py           # ComponentRegistry (dict wrapper with validation)
backend/app/models/
├── sa_params.py              # Extended with component_weights: dict[str, float]
```

### Pattern 1: Protocol-Based Scoring Components

**What:** Define a structural interface for scoring components using typing.Protocol. Components implement `compute(mol_graph: MoleculeGraph) -> float` without formal inheritance.

**When to use:** When you need pluggable behavior without coupling to a base class. Ideal for educational codebases where transparency matters.

**Example:**
```python
# backend/app/core/scoring/protocol.py
from typing import Protocol
from app.core.molecule import MoleculeGraph

class ScoringComponent(Protocol):
    """Structural interface for scoring components.

    Any class implementing compute(mol_graph) -> float
    satisfies this protocol without explicit inheritance.
    """

    name: str  # Human-readable component name

    def compute(self, mol_graph: MoleculeGraph) -> float:
        """Compute score for the given molecular graph.

        Args:
            mol_graph: Molecule to score

        Returns:
            Component score (higher or lower is better depends on optimization mode)
        """
        ...

# backend/app/core/scoring/wiener.py
from app.core.molecule import MoleculeGraph
from rdkit import Chem

class WienerIndexComponent:
    """Wiener Index scoring component (sum of pairwise distances)."""

    name = "wiener_index"

    def compute(self, mol_graph: MoleculeGraph) -> float:
        mol = mol_graph.get_mol()
        n = mol.GetNumAtoms()

        if n <= 1:
            return 0

        dist_matrix = Chem.GetDistanceMatrix(mol)

        # Check for disconnected graph
        for i in range(n):
            for j in range(n):
                if dist_matrix[i][j] > 1000:
                    raise ValueError('Graph is disconnected')

        # Sum upper triangle
        wiener = sum(
            dist_matrix[i][j]
            for i in range(n)
            for j in range(i + 1, n)
        )

        return wiener
```

**Benefits:**
- No inheritance coupling (components are independent)
- Type checker validates interface compliance
- Easy to add new components (just implement the protocol)
- Transparent for educational use (no metaclass magic)

**Source:** [Real Python - Python Protocols](https://realpython.com/python-protocol/)

### Pattern 2: Weighted Sum Objective Function

**What:** Combine multiple scoring components using a weighted sum: `energy = Σ(weight_i × score_i)`. Each component contributes independently, weights control relative importance.

**When to use:** When you have multiple criteria and a single scalar objective (as opposed to Pareto optimization with multiple objectives).

**Example:**
```python
# backend/app/core/sa_engine.py (refactored _iterate method)

def _compute_energy(self, mol_graph: MoleculeGraph) -> float:
    """Compute weighted sum of scoring component contributions.

    Args:
        mol_graph: Molecule to score

    Returns:
        Total energy (weighted sum of all components)
    """
    total_energy = 0.0

    for component_name, weight in self._component_weights.items():
        if weight == 0.0:
            continue  # Skip components with zero weight (optimization)

        component = self._component_registry[component_name]
        score = component.compute(mol_graph)
        total_energy += weight * score

    return total_energy
```

**Optimization mode handling:**
The weighted sum is agnostic to optimization mode. The Metropolis criterion in `_iterate()` already handles MINIMIZE vs MAXIMIZE by computing delta_e appropriately (lines 232-234). No changes needed.

**Component normalization:**
For v2.0, components are NOT normalized. Users choose weights that make sense for their problem. For example:
- Wiener Index: weight 1.0 (default, values ~20-100 for small molecules)
- Molecular Weight: weight 0.0 (disabled by default, values ~100-500 Da)

Future phases may add automatic normalization (z-scores, min-max scaling) if users request it.

**Source:** [GeeksforGeeks - Weighted Sum Method](https://www.geeksforgeeks.org/dsa/weighted-sum-method-multi-criteria-decision-making/)

### Pattern 3: Component Registry with Validation

**What:** Centralized registry mapping component names to instances. Validates that requested components exist and weights are sensible.

**When to use:** When you have a small, fixed set of components (3-10) and don't need dynamic plugin loading.

**Example:**
```python
# backend/app/core/scoring/registry.py
from typing import Dict
from app.core.scoring.protocol import ScoringComponent
from app.core.scoring.wiener import WienerIndexComponent
from app.core.scoring.molecular_weight import MolecularWeightComponent

class ComponentRegistry:
    """Registry of available scoring components."""

    def __init__(self):
        self._components: Dict[str, ScoringComponent] = {}

        # Register built-in components
        self.register(WienerIndexComponent())
        self.register(MolecularWeightComponent())

    def register(self, component: ScoringComponent) -> None:
        """Register a scoring component.

        Args:
            component: Component to register

        Raises:
            ValueError: If component name is already registered
        """
        if component.name in self._components:
            raise ValueError(f"Component '{component.name}' already registered")

        self._components[component.name] = component

    def get(self, name: str) -> ScoringComponent:
        """Get component by name.

        Args:
            name: Component name

        Returns:
            Component instance

        Raises:
            KeyError: If component not found
        """
        if name not in self._components:
            available = ', '.join(self._components.keys())
            raise KeyError(
                f"Unknown component '{name}'. Available: {available}"
            )

        return self._components[name]

    def list_components(self) -> list[str]:
        """List all registered component names."""
        return list(self._components.keys())

# Global registry instance
_registry = ComponentRegistry()

def get_registry() -> ComponentRegistry:
    """Get the global component registry."""
    return _registry
```

**Why dict over importlib entry points:**
- YAGNI: No third-party plugins needed for classroom tool
- Simplicity: Explicit registration is easier to debug
- Control: Registry lives in application code, not package metadata

**Source:** [Building a Plugin Architecture with Python](https://mwax911.medium.com/building-a-plugin-architecture-with-python-7b4ab39ad4fc)

### Pattern 4: Configuration with Pydantic

**What:** Extend SAParams to include `component_weights: dict[str, float]` with validation and defaults.

**When to use:** When you need user-configurable weights with type safety and validation.

**Example:**
```python
# backend/app/models/sa_params.py (extended)
from pydantic import BaseModel, Field, field_validator
from typing import Dict

class SAParams(BaseModel):
    """Simulated annealing configuration parameters."""

    formula: str = Field(
        ..., pattern=r"^([A-Z][a-z]?\d*)+$", examples=["C6H14"]
    )
    initial_temp: float = Field(
        100.0, gt=0, description="Initial temperature kT_0"
    )
    cooling_schedule_k: float = Field(
        8.0, ge=0, description="Cooling rate parameter k"
    )
    steps_per_cycle: int = Field(
        500, gt=0, description="Steps per temperature cycle"
    )
    num_cycles: int = Field(4, gt=0, description="Number of cooling cycles")
    optimization_mode: Literal["MINIMIZE", "MAXIMIZE"] = "MINIMIZE"
    seed: int = Field(42, description="Random seed for reproducibility")

    # NEW: Component weights configuration
    component_weights: Dict[str, float] = Field(
        default_factory=lambda: {"wiener_index": 1.0},
        description="Weights for scoring components (sum need not equal 1.0)"
    )

    @field_validator("component_weights")
    @classmethod
    def validate_weights(cls, v: Dict[str, float]) -> Dict[str, float]:
        """Validate component weights are non-negative."""
        if not v:
            raise ValueError("At least one component weight must be specified")

        for name, weight in v.items():
            if weight < 0:
                raise ValueError(f"Component '{name}' has negative weight: {weight}")

        # Check at least one non-zero weight
        if all(w == 0.0 for w in v.values()):
            raise ValueError("At least one component weight must be non-zero")

        return v
```

**Default behavior (backward compatible):**
If `component_weights` is not specified, defaults to `{"wiener_index": 1.0}`, which replicates current behavior exactly.

**Future extensibility:**
When TGT-04 (adjustable weights via API) is implemented, the frontend can send:
```json
{
  "formula": "C6H14",
  "component_weights": {
    "wiener_index": 0.7,
    "molecular_weight": 0.3
  }
}
```

No backend changes needed beyond what's in Phase 7.

**Source:** Existing pattern in `backend/app/models/sa_params.py` (lines 11-28)

### Anti-Patterns to Avoid

**1. Premature Normalization**
- **Don't:** Normalize all component scores to [0,1] or z-scores in v2.0
- **Why:** Adds complexity without clear user benefit. Users can achieve equivalent results by adjusting weights.
- **When to add:** If users complain that components have incompatible scales (e.g., Wiener Index ~50, LogP ~2) and can't find sensible weights.

**2. Over-Abstraction with Factory Classes**
- **Don't:** Create ComponentFactory, ComponentBuilder, AbstractComponentProvider when you have 2-3 components
- **Why:** "Developers build plugin systems for a single plugin and write factory classes when there's only one implementation" (YAGNI violation)
- **Keep it simple:** Direct instantiation in registry is fine for <10 components

**Source:** [CodeOpinion - Do You Really Need That Abstraction?](https://codeopinion.com/do-you-really-need-that-abstraction-or-generic-code-yagni/)

**3. Runtime Component Loading from Config Files**
- **Don't:** Load component names from JSON/YAML and use importlib to dynamically import modules
- **Why:** Security risk (arbitrary code execution), debugging nightmare, YAGNI for classroom tool
- **Keep it simple:** Hardcode component registration in `registry.py`

**4. Component State Mutation**
- **Don't:** Let components maintain state between `compute()` calls
- **Why:** SA engine may call `compute()` thousands of times; stateful components leak memory and break reproducibility
- **Enforce:** `compute()` should be a pure function (same input → same output, no side effects)

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Molecular descriptors | Custom descriptor calculators | RDKit Chem.Descriptors module | RDKit has 200+ descriptors battle-tested in production. Example: `Descriptors.MolWt(mol)`, `Descriptors.MolLogP(mol)` |
| Distance matrix | Custom BFS/Floyd-Warshall | RDKit `Chem.GetDistanceMatrix()` | Already used in wiener.py (line 30); handles edge cases correctly |
| Weighted average | Custom sum loop | NumPy `np.average(scores, weights=weights)` | Handles edge cases (zero weights, NaN) and is 10-100x faster for large arrays |
| Configuration validation | Manual if/else checks | Pydantic field validators | Already in stack; automatic type coercion, clear error messages |

**Key insight:** RDKit is comprehensive. Before writing any molecular property calculation, check `rdkit.Chem.Descriptors` and `rdkit.Chem.rdMolDescriptors` modules. The only descriptor we hand-roll is Wiener Index, and only because it's the pedagogical focus.

**Source:** [RDKit Descriptors Documentation](https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html)

## Common Pitfalls

### Pitfall 1: Inconsistent Optimization Direction

**What goes wrong:** Component A (Wiener Index) should be minimized, but component B (LogP drug-likeness) should be maximized. Naive weighted sum breaks.

**Why it happens:** Different descriptors have different "good" directions. Wiener Index low = compact molecule, but LogP moderate (~2-3) = drug-like.

**How to avoid:**
- **Option 1 (v2.0):** All components use same optimization mode. User picks MINIMIZE or MAXIMIZE and components are chosen accordingly. Wiener Index + molecular weight both minimize.
- **Option 2 (future):** Each component declares `optimization_direction: Literal["MINIMIZE", "MAXIMIZE"]`. Weighted sum inverts MAXIMIZE components: `if component.direction == "MAXIMIZE": score = -score`.
- **Recommendation for Phase 7:** Use Option 1. Keep it simple. Defer per-component directions to v2.1+.

**Warning signs:**
- Best energy oscillates wildly instead of converging
- SA acceptance ratio is ~100% (all moves accepted) or ~0% (all rejected)
- User says "I wanted to maximize LogP but it's getting worse"

**Source:** Multi-objective optimization literature (weighted sum assumes aligned objectives)

### Pitfall 2: Component Score Explosion

**What goes wrong:** Wiener Index ~50, but molecular weight ~200. Weight vector `{wiener: 0.5, mw: 0.5}` means energy dominated by MW (50% of 50 = 25 vs 50% of 200 = 100).

**Why it happens:** Components have different scales. Weighted sum is scale-sensitive.

**How to avoid:**
- **Document clearly:** Explain in API docs that weights are not percentages; they're multiplicative factors.
- **Provide guidance:** Include example weight combinations in documentation. "For Wiener + MW, try `{wiener: 1.0, mw: 0.01}` to balance scales."
- **Future fix:** Add optional normalization in v2.1+ if users struggle with scaling.

**Warning signs:**
- User complains "changing weight from 0.5 to 0.6 had no effect"
- One component always dominates (check via component-wise score logging)

**Source:** Standard multi-criteria decision making (MCDM) literature on weighted sum limitations

### Pitfall 3: Zero-Weight Component Overhead

**What goes wrong:** User sets `{"wiener_index": 1.0, "molecular_weight": 0.0}` but SA engine still computes molecular weight 2000 times per run, wasting CPU.

**Why it happens:** Naive loop over all components without checking weights.

**How to avoid:**
```python
def _compute_energy(self, mol_graph: MoleculeGraph) -> float:
    total_energy = 0.0

    for component_name, weight in self._component_weights.items():
        if weight == 0.0:
            continue  # Skip zero-weight components

        component = self._component_registry[component_name]
        score = component.compute(mol_graph)
        total_energy += weight * score

    return total_energy
```

**Warning signs:**
- SA runs slower when users add components with zero weight
- Profiling shows time spent in unused components

**Source:** Standard optimization practice (don't compute what you won't use)

### Pitfall 4: Component Registry Desynchronization

**What goes wrong:** User requests `{"wiener_index": 1.0, "logp": 0.5}` but LogP component isn't registered. Error occurs deep in SA engine, not at request validation.

**Why it happens:** No validation that requested component names exist in registry.

**How to avoid:**
```python
# In SAEngine.__init__()
def __init__(self, params: SAParams):
    self._params = params
    self._component_registry = get_registry()

    # Validate component weights reference existing components
    available = set(self._component_registry.list_components())
    requested = set(params.component_weights.keys())
    unknown = requested - available

    if unknown:
        raise ValueError(
            f"Unknown scoring components: {unknown}. "
            f"Available: {available}"
        )

    self._component_weights = params.component_weights
```

**Warning signs:**
- Confusing errors like "KeyError: 'logp'" during SA execution
- User has to read source code to discover available components

**Source:** Fail-fast validation principle

## Code Examples

Verified patterns from research:

### Example 1: Complete Component Implementation

```python
# backend/app/core/scoring/molecular_weight.py
from app.core.molecule import MoleculeGraph
from rdkit.Chem import Descriptors

class MolecularWeightComponent:
    """Molecular weight scoring component (lower MW = lower energy)."""

    name = "molecular_weight"

    def compute(self, mol_graph: MoleculeGraph) -> float:
        """Compute molecular weight in Daltons.

        Args:
            mol_graph: Molecule to score

        Returns:
            Molecular weight (Da)
        """
        mol = mol_graph.get_mol()
        return Descriptors.MolWt(mol)
```

**This is a complete, working component.** No inheritance, no boilerplate, just implement the protocol.

### Example 2: SA Engine Refactored

```python
# backend/app/core/sa_engine.py (lines 200-226 refactored)

def _iterate(self, temperature: float, step_number: int) -> None:
    """Execute a single SA iteration.

    Disconnected displacements are discarded (counted as invalid).
    Only connected structures proceed to Metropolis scoring.

    Args:
        temperature: Current temperature (kT)
        step_number: Current step number (1-indexed for display)
    """
    # Attempt displacement
    proposed_graph = attempt_displacement(self._current_graph, self._rng)

    if proposed_graph is None:
        self._invalid_moves += 1
        self._record_step(step_number, temperature, False)
        return

    if not proposed_graph.is_connected():
        self._disconnected_moves += 1
        self._record_step(step_number, temperature, False)
        return

    # CHANGED: Compute energy using multi-component framework
    proposed_energy = self._compute_energy(proposed_graph)

    # Compute energy delta (rest unchanged)
    if self._params.optimization_mode == "MINIMIZE":
        delta_e = proposed_energy - self._current_energy
    else:
        delta_e = self._current_energy - proposed_energy

    # ... rest of method unchanged
```

**Only one line changes:** `compute_wiener_index(proposed_graph)` becomes `self._compute_energy(proposed_graph)`.

### Example 3: Testing Components

```python
# backend/tests/test_scoring_components.py
import pytest
from app.core.molecule import MoleculeGraph
from app.core.scoring.wiener import WienerIndexComponent
from app.core.scoring.molecular_weight import MolecularWeightComponent

def test_wiener_component_linear_chain():
    """Wiener index of linear chain: n(n²-1)/6."""
    # Create linear C4 chain
    mol_graph = MoleculeGraph.from_smiles("CCCC")

    component = WienerIndexComponent()
    wiener = component.compute(mol_graph)

    # Linear chain of n=4: 4*(16-1)/6 = 10
    assert wiener == 10.0

def test_molecular_weight_component():
    """Molecular weight of C6H14 = 86.18 Da."""
    mol_graph = MoleculeGraph.from_smiles("CCCCCC")

    component = MolecularWeightComponent()
    mw = component.compute(mol_graph)

    # C6H14: 6*12.01 + 14*1.008 = 86.178
    assert 86.0 < mw < 87.0

def test_components_satisfy_protocol():
    """Components satisfy ScoringComponent protocol."""
    from app.core.scoring.protocol import ScoringComponent

    # Type checker will validate at static analysis time
    # Runtime validation requires @runtime_checkable decorator
    wiener = WienerIndexComponent()
    mw = MolecularWeightComponent()

    assert hasattr(wiener, 'name')
    assert hasattr(wiener, 'compute')
    assert callable(wiener.compute)
```

**Test coverage should include:**
- Individual component correctness (known molecules, known scores)
- Protocol compliance (has name, has compute)
- Registry operations (register, get, list, duplicate names)
- Weighted sum calculation (multiple components, zero weights)
- Validation (unknown components, negative weights, all-zero weights)

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Single hardcoded objective | Pluggable multi-component framework | Phase 7 (TGT-01) | Enables future descriptor exploration without SA engine changes |
| Weights hardcoded to 1.0 | User-configurable via SAParams | Phase 7 (TGT-03) | Users can balance multiple criteria |
| ABC inheritance for plugins | Protocol structural typing | Python 3.8+ (PEP 544) | Simpler component implementation, no inheritance required |
| Pareto frontiers (MOSA) | Weighted sum scalarization | v2.0 choice | Weighted sum sufficient for educational use; MOSA deferred to v3.0+ |

**Deprecated/outdated:**
- **ABC-only plugin systems:** Protocols are now standard for duck-typed interfaces (Python 3.8+, mainstream since 3.10+)
- **Entry points for internal plugins:** Overkill for applications (vs libraries). Use entry points when distributing third-party plugin APIs, not for internal component discovery.
- **scipy.optimize.anneal:** Deprecated since SciPy 0.14 (2014). Use custom SA implementations or third-party libraries like simanneal.

**Source:** [Python Typing Spec - Protocols](https://typing.python.org/en/latest/spec/protocol.html)

## Open Questions

### Question 1: Component Normalization Strategy

**What we know:**
- Weighted sum is scale-sensitive (Wiener ~50, MW ~200)
- RDKit QED uses ADS (Applicability Domain Score) functions to map raw values to [0,1]
- Users may struggle finding weights that balance components

**What's unclear:**
- Will classroom users understand raw scores + weights, or do they need normalized scores?
- If normalization is needed, should it be z-scores (mean=0, std=1) or min-max [0,1]?

**Recommendation:**
- **v2.0:** No normalization. Document scale differences clearly. Include example weight combinations.
- **Future:** Add optional normalization if users request it. Easiest is min-max scaling based on typical molecular ranges (e.g., C4-C20 alkanes have Wiener 10-500).

### Question 2: Component Metadata (Units, Descriptions, Ranges)

**What we know:**
- Each component computes a score (float)
- Components have a name (string)

**What's unclear:**
- Should components declare units? (e.g., "Daltons", "dimensionless")
- Should components declare typical ranges? (e.g., "Wiener: 10-500 for C4-C20")
- Should components have descriptions for UI display?

**Recommendation:**
- **v2.0:** Keep protocol minimal (`name: str`, `compute() -> float`). Add metadata later if frontend needs it.
- **Future (TGT-05):** When component-wise score streaming is implemented, add metadata for UI display:
  ```python
  class ScoringComponent(Protocol):
      name: str
      display_name: str  # "Wiener Index"
      description: str   # "Sum of all pairwise distances"
      units: str         # "dimensionless"
  ```

### Question 3: Component Caching

**What we know:**
- SA engine evaluates 2000 molecules per run (4 cycles × 500 steps)
- Each evaluation computes all components
- Some components may be expensive (RDKit 3D descriptors)

**What's unclear:**
- Should we cache component scores for previously-seen SMILES strings?
- Is cache overhead (hashing, lookup) worth it for cheap components like MW?

**Recommendation:**
- **v2.0:** No caching. All current components are O(n²) or faster, negligible overhead.
- **Future:** If 3D descriptors are added (require conformer generation, ~10-100ms), implement LRU cache keyed by SMILES.

## Sources

### Primary (HIGH confidence)
- [RDKit QED Module](https://www.rdkit.org/docs/source/rdkit.Chem.QED.html) - Weighted scoring implementation patterns
- [RDKit Descriptors](https://www.rdkit.org/docs/source/rdkit.Chem.Descriptors.html) - Available molecular descriptors
- [Real Python - Python Protocols](https://realpython.com/python-protocol/) - Protocol vs ABC tradeoffs
- [Python Typing Spec - Protocols](https://typing.python.org/en/latest/spec/protocol.html) - Official Protocol specification
- Existing codebase: `backend/app/core/sa_engine.py`, `backend/app/models/sa_params.py`

### Secondary (MEDIUM confidence)
- [Strategy Pattern in Python](https://refactoring.guru/design-patterns/strategy/python/example) - Pluggable algorithm pattern
- [GeeksforGeeks - Weighted Sum Method](https://www.geeksforgeeks.org/dsa/weighted-sum-method-multi-criteria-decision-making/) - Multi-criteria decision making
- [Building a Plugin Architecture with Python](https://mwax911.medium.com/building-a-plugin-architecture-with-python-7b4ab39ad4fc) - Registry pattern examples
- [MOSA GitHub](https://github.com/rgaveiga/mosa) - Multi-objective SA reference (Pareto approach)

### Tertiary (LOW confidence)
- [CodeOpinion - YAGNI](https://codeopinion.com/do-you-really-need-that-abstraction-or-generic-code-yagni/) - Avoiding premature abstraction (general advice, not Python-specific)
- WebSearch results on weighted average implementations (multiple sources, general patterns)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - Protocol typing is well-established (Python 3.8+), RDKit descriptors are documented
- Architecture: MEDIUM-HIGH - Patterns are standard (Strategy, weighted sum) but application to SA engine is novel
- Pitfalls: MEDIUM - Identified from general MCDM literature + engineering experience, not domain-specific SA pitfalls
- Don't hand-roll: HIGH - RDKit descriptor library is comprehensive and well-documented

**Research date:** 2026-02-16
**Valid until:** ~60 days (stable domain; Python typing patterns unlikely to change; RDKit releases quarterly but descriptor APIs are stable)

**Primary gaps:**
- No direct examples of Protocol-based SA component systems in literature (novel application)
- Uncertain whether classroom users will need normalization (requires user testing)
- Component metadata requirements unclear until frontend integration (Phase 6+ feedback needed)

**Overall confidence: MEDIUM-HIGH**
- Architecture patterns are solid and well-researched
- Implementation is straightforward (Protocol + weighted sum)
- Uncertainty around UX decisions (normalization, metadata) but those are deferred to future phases
