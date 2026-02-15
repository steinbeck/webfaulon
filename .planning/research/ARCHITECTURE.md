# Architecture Research: FastAPI + Python RDKit Backend Integration

**Domain:** FastAPI Backend with Python RDKit for Molecular Optimization
**Researched:** 2026-02-15
**Confidence:** HIGH

## Executive Summary

This research covers the integration of a FastAPI Python backend with RDKit computational chemistry library into the existing browser-based WebFaulon SA demo. The architecture replaces the Web Worker pattern with a client-server model using Server-Sent Events (SSE) for real-time progress streaming.

**Key architectural shift:** Web Worker + RDKit.js WASM → FastAPI REST API + SSE + Python RDKit

## System Overview

### v2.0 Architecture with Backend

```
┌──────────────────────────────────────────────────────────────────────────┐
│                     FRONTEND (Vite + TypeScript)                          │
├──────────────────────────────────────────────────────────────────────────┤
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌─────────────┐  │
│  │   Alpine.js  │  │   Chart.js   │  │   Molecule   │  │  Parameter  │  │
│  │   UI State   │  │  Wiener Plot │  │   Display    │  │   Controls  │  │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘  └──────┬──────┘  │
│         │                 │                 │                 │          │
│         └─────────────────┴─────────────────┴─────────────────┘          │
│                                    │                                      │
│                     ┌──────────────▼───────────────┐                      │
│                     │   API Client (fetch/axios)   │                      │
│                     │   - POST /api/sa/optimize    │                      │
│                     │   - EventSource /api/sa/...  │                      │
│                     └──────────────┬───────────────┘                      │
│                                    │                                      │
└────────────────────────────────────┼──────────────────────────────────────┘
                                     │
                              HTTP + SSE
                                     │
┌────────────────────────────────────▼──────────────────────────────────────┐
│                     BACKEND (FastAPI + Python)                            │
├───────────────────────────────────────────────────────────────────────────┤
│                        API LAYER (FastAPI Router)                         │
│  ┌────────────────────────────────────────────────────────────────────┐   │
│  │  POST /api/sa/optimize (params) → 202 + job_id                     │   │
│  │  GET  /api/sa/stream/{job_id}  → SSE progress stream               │   │
│  │  GET  /api/sa/status/{job_id}  → Job status/result                 │   │
│  └────────────────────────────────┬───────────────────────────────────┘   │
│                                   │                                       │
├───────────────────────────────────┼───────────────────────────────────────┤
│                       SERVICE LAYER (Business Logic)                      │
│  ┌────────────────────────────────▼───────────────────────────────────┐   │
│  │  SAService                                                          │   │
│  │  - run_optimization(params, callback)                               │   │
│  │  - get_job_status(job_id)                                           │   │
│  │  - stream_progress(job_id)                                          │   │
│  └────────────────────────────────┬───────────────────────────────────┘   │
│                                   │                                       │
│  ┌────────────────────────────────▼───────────────────────────────────┐   │
│  │  TargetFunctionService (Multi-Component Framework)                  │   │
│  │  - evaluate(mol, components, weights)                               │   │
│  │  - register_component(name, fn)                                     │   │
│  └────────────────────────────────┬───────────────────────────────────┘   │
│                                   │                                       │
├───────────────────────────────────┼───────────────────────────────────────┤
│                       CORE LAYER (SA + RDKit)                             │
│  ┌────────────────────────────────▼───────────────────────────────────┐   │
│  │  SAEngine (Python port of TypeScript SAEngine)                      │   │
│  │  - init(formula, params)                                            │   │
│  │  - step() → SAStepResult                                            │   │
│  │  - get_state() → SAEngineState                                      │   │
│  └────────────────────────────────┬───────────────────────────────────┘   │
│                                   │                                       │
│  ┌────────────────────────────────▼───────────────────────────────────┐   │
│  │  MoleculeService (RDKit wrapper)                                    │   │
│  │  - create_from_formula(formula) → rdkit.Mol                         │   │
│  │  - apply_displacement(mol, atom_i, atom_j) → rdkit.Mol             │   │
│  │  - compute_wiener_index(mol) → float                                │   │
│  │  - generate_svg(mol, width, height) → str                           │   │
│  │  - mol_to_smiles(mol) → str                                         │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                            │
├───────────────────────────────────────────────────────────────────────────┤
│                       STATE LAYER (Job Management)                        │
│  ┌─────────────────────────────────────────────────────────────────────┐  │
│  │  In-Memory Job Store (asyncio.Queue per job_id)                     │  │
│  │  - job_state: Dict[str, JobState]                                   │  │
│  │  - progress_queues: Dict[str, asyncio.Queue[ProgressEvent]]         │  │
│  └─────────────────────────────────────────────────────────────────────┘  │
│                                                                            │
└───────────────────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

### Frontend Components (Modified from v1.0)

| Component | v1.0 Responsibility | v2.0 Changes |
|-----------|---------------------|--------------|
| **Alpine.js UI State** | Worker communication via Comlink | HTTP client + EventSource SSE listener |
| **API Client** | N/A (was Comlink) | REST API calls + SSE stream management |
| **Chart.js** | Worker progress callbacks | SSE event handlers for progress data |
| **Molecule Display** | RDKit.js WASM canvas rendering | Backend SVG rendering (img src="data:image/svg+xml;base64,...") |

### Backend Components (NEW)

| Component | Responsibility | Implementation |
|-----------|----------------|----------------|
| **FastAPI Router** | HTTP endpoints for SA optimization | APIRouter with async handlers |
| **SAService** | Orchestrate SA execution, manage job lifecycle | Async service with background tasks |
| **TargetFunctionService** | Multi-component scoring framework | Plugin registry + weighted composition |
| **SAEngine** | Python port of TypeScript SAEngine | Stateful class with step() method |
| **MoleculeService** | RDKit operations wrapper | Stateless utility functions |
| **Job Store** | Track active/completed jobs, route progress events | In-memory Dict + asyncio.Queue |
| **SSE Streamer** | Stream progress events to client | async generator yielding `data: {...}\n\n` |

## Recommended Project Structure (Monorepo)

```
webfaulon/                           # Monorepo root
├── frontend/                        # Vite frontend (moved from src/)
│   ├── src/
│   │   ├── core/                    # Types only (shared with backend via OpenAPI)
│   │   │   └── types.ts             # Frontend-specific types
│   │   ├── ui/                      # Alpine.js components
│   │   │   ├── app.ts               # Main Alpine component
│   │   │   ├── chart.ts             # Chart.js integration
│   │   │   └── presets.ts           # Molecule presets
│   │   ├── services/                # NEW: API client layer
│   │   │   ├── api-client.ts        # HTTP client (fetch wrapper)
│   │   │   └── sse-client.ts        # EventSource wrapper
│   │   └── main.ts                  # Vite entry point
│   ├── index.html
│   ├── package.json
│   ├── tsconfig.json
│   └── vite.config.ts
│
├── backend/                         # NEW: FastAPI backend
│   ├── app/
│   │   ├── __init__.py
│   │   ├── main.py                  # FastAPI app + CORS + startup
│   │   ├── api/                     # API routes
│   │   │   ├── __init__.py
│   │   │   └── sa.py                # /api/sa/* endpoints
│   │   ├── services/                # Business logic
│   │   │   ├── __init__.py
│   │   │   ├── sa_service.py        # SA orchestration
│   │   │   ├── molecule_service.py  # RDKit wrapper
│   │   │   └── target_fn_service.py # Multi-component scoring
│   │   ├── core/                    # Core SA implementation
│   │   │   ├── __init__.py
│   │   │   ├── sa_engine.py         # Python port of SAEngine.ts
│   │   │   ├── displacement.py      # Faulon displacement (eqs 7-11)
│   │   │   ├── cooling.py           # Temperature schedules
│   │   │   └── wiener.py            # Wiener index via RDKit
│   │   ├── models/                  # Pydantic models
│   │   │   ├── __init__.py
│   │   │   ├── sa_params.py         # SAParams, SAResult
│   │   │   └── progress.py          # ProgressEvent, JobStatus
│   │   ├── state/                   # Job state management
│   │   │   ├── __init__.py
│   │   │   └── job_store.py         # In-memory job tracking
│   │   └── utils/                   # Utilities
│   │       ├── __init__.py
│   │       └── rdkit_helpers.py     # RDKit initialization/setup
│   ├── tests/                       # pytest tests
│   │   ├── test_sa_engine.py
│   │   ├── test_molecule_service.py
│   │   └── test_api.py
│   ├── requirements.txt             # fastapi, uvicorn, rdkit, pydantic
│   ├── pyproject.toml               # Python project metadata (optional)
│   └── README.md
│
├── shared/                          # NEW: Shared types (optional)
│   └── schema.yaml                  # OpenAPI schema for type generation
│
├── .planning/                       # GSD planning (unchanged)
├── package.json                     # Root workspace config (pnpm/npm workspaces)
├── README.md                        # Project-level README
└── docker-compose.yml               # Dev environment (optional)
```

### Structure Rationale

- **Monorepo layout:** Keeps frontend/backend in sync, simplifies development workflow, enables shared types via OpenAPI codegen
- **Separate `frontend/` and `backend/`:** Clear boundaries, independent build systems (Vite vs Python), independent deployments
- **Backend layered architecture:** API → Service → Core → Models (standard FastAPI pattern)
- **Job store in-memory:** Sufficient for demo/classroom use; can swap for Redis/DB later
- **No shared code between TS/Python:** Use OpenAPI schema as source of truth for types (generate TS types from Pydantic models)

## Data Flow Patterns

### Pattern 1: SA Optimization Request Flow

**What:** Client initiates SA optimization, receives progress updates via SSE, displays final result

**Request Flow:**
```
[Alpine.js UI]
    │ POST /api/sa/optimize {formula, params}
    ▼
[FastAPI Router]
    │ Create job_id, spawn background task
    ▼
[SAService.run_optimization()]
    │ Initialize SAEngine, run step() loop
    │ Emit progress to asyncio.Queue
    ▼
[Job Store]
    │ Store progress events in queue
    │ Store final result in job_state

PARALLEL:

[Alpine.js UI]
    │ GET /api/sa/stream/{job_id} (EventSource)
    ▼
[FastAPI SSE Endpoint]
    │ async for event in job_store.progress_queues[job_id]
    │ yield f"data: {json.dumps(event)}\n\n"
    ▼
[EventSource onmessage]
    │ Update chart, molecule display
```

**Code Example:**
```python
# backend/app/api/sa.py
from fastapi import APIRouter, BackgroundTasks
from fastapi.responses import StreamingResponse
from sse_starlette.sse import EventSourceResponse

router = APIRouter(prefix="/api/sa")

@router.post("/optimize")
async def optimize(params: SAParams, background_tasks: BackgroundTasks):
    job_id = generate_job_id()
    job_store.create_job(job_id, params)
    background_tasks.add_task(sa_service.run_optimization, job_id, params)
    return {"job_id": job_id, "status": "started"}

@router.get("/stream/{job_id}")
async def stream_progress(job_id: str):
    async def event_generator():
        queue = job_store.get_progress_queue(job_id)
        while True:
            event = await queue.get()
            if event["type"] == "complete":
                yield {"data": json.dumps(event)}
                break
            yield {"data": json.dumps(event)}

    return EventSourceResponse(event_generator())
```

```typescript
// frontend/src/services/sse-client.ts
export async function startOptimization(params: SAParams, onProgress: (data: ProgressEvent) => void) {
  const { job_id } = await fetch('/api/sa/optimize', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(params)
  }).then(r => r.json());

  const eventSource = new EventSource(`/api/sa/stream/${job_id}`);
  eventSource.onmessage = (e) => {
    const event = JSON.parse(e.data);
    onProgress(event);
    if (event.type === 'complete') {
      eventSource.close();
    }
  };

  return { job_id, eventSource };
}
```

**Trade-offs:**
- ✅ Simpler than WebSockets (HTTP-based, auto-reconnect)
- ✅ Scales well with async Python (thousands of concurrent SSE streams)
- ✅ Built-in browser EventSource API
- ⚠️ Unidirectional (server → client only; need separate POST for pause/cancel)
- ⚠️ HTTP/1.1 connection limit (6 per domain; use HTTP/2 in production)

### Pattern 2: Multi-Component Target Function

**What:** Pluggable scoring framework supporting weighted combination of multiple molecular properties

**Architecture:**
```python
# backend/app/services/target_fn_service.py
from typing import Callable, Dict
from rdkit.Chem import Mol

ComponentFn = Callable[[Mol], float]

class TargetFunctionService:
    def __init__(self):
        self._components: Dict[str, ComponentFn] = {}

    def register(self, name: str, fn: ComponentFn):
        """Register a scoring component"""
        self._components[name] = fn

    def evaluate(self, mol: Mol, weights: Dict[str, float]) -> float:
        """Compute weighted sum of components"""
        score = 0.0
        for name, weight in weights.items():
            if name not in self._components:
                raise ValueError(f"Unknown component: {name}")
            score += weight * self._components[name](mol)
        return score

# Register components at startup
target_fn_service = TargetFunctionService()
target_fn_service.register("wiener_index", molecule_service.compute_wiener_index)
target_fn_service.register("logp", lambda mol: Descriptors.MolLogP(mol))
target_fn_service.register("molecular_weight", lambda mol: Descriptors.MolWt(mol))
```

**API Usage:**
```json
POST /api/sa/optimize
{
  "formula": "C6H14",
  "target_components": {
    "wiener_index": 1.0,
    "logp": 0.0
  },
  "optimization_mode": "MINIMIZE",
  ...
}
```

**Trade-offs:**
- ✅ Extensible: add new components without changing core SA logic
- ✅ Composable: combine multiple properties with weights
- ✅ Testable: each component is a pure function
- ⚠️ Performance: weighted sum computed on every SA step (minimize component count)
- ⚠️ Weight normalization: client must ensure weights sum to desired total

### Pattern 3: RDKit Molecule Lifecycle

**What:** Manage RDKit Mol objects through SA iterations (creation, mutation, validation, rendering)

**Lifecycle:**
```python
# backend/app/services/molecule_service.py
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Draw, AllChem

class MoleculeService:
    @staticmethod
    def create_from_formula(formula: str) -> Chem.Mol:
        """Create initial linear structure from molecular formula"""
        # Parse formula (e.g., "C6H14" → 6 carbons, 14 hydrogens)
        atom_counts = parse_formula(formula)

        # Create linear chain of heavy atoms
        mol = Chem.RWMol()
        for elem, count in atom_counts.items():
            if elem == 'H':
                continue  # Implicit hydrogens
            for _ in range(count):
                mol.AddAtom(Chem.Atom(elem))

        # Add single bonds in linear chain
        for i in range(mol.GetNumAtoms() - 1):
            mol.AddBond(i, i + 1, Chem.BondType.SINGLE)

        # Sanitize and return
        mol = mol.GetMol()
        Chem.SanitizeMol(mol)
        return mol

    @staticmethod
    def apply_displacement(mol: Chem.Mol, atom_i: int, atom_j: int) -> Chem.Mol:
        """Apply Faulon displacement (eqs 7-11): modify bond between atoms i, j"""
        mol_rw = Chem.RWMol(mol)  # Make editable copy

        bond = mol_rw.GetBondBetweenAtoms(atom_i, atom_j)
        if bond is None:
            # No bond → try to add one (if valid)
            if can_add_bond(mol_rw, atom_i, atom_j):
                mol_rw.AddBond(atom_i, atom_j, Chem.BondType.SINGLE)
        else:
            # Bond exists → cycle order (1 → 2 → 3 → 0)
            current_order = bond.GetBondTypeAsDouble()
            new_order = cycle_bond_order(current_order)
            if new_order == 0:
                mol_rw.RemoveBond(atom_i, atom_j)
            else:
                bond.SetBondType(order_to_bond_type(new_order))

        # Validate and return
        new_mol = mol_rw.GetMol()
        try:
            Chem.SanitizeMol(new_mol)
            if is_connected(new_mol):
                return new_mol
        except:
            pass

        return None  # Invalid displacement

    @staticmethod
    def compute_wiener_index(mol: Chem.Mol) -> float:
        """Compute Wiener index (sum of pairwise distances)"""
        dist_matrix = Chem.GetDistanceMatrix(mol)
        # Sum upper triangle (avoid double-counting)
        n = mol.GetNumAtoms()
        return sum(dist_matrix[i, j] for i in range(n) for j in range(i+1, n))

    @staticmethod
    def generate_svg(mol: Chem.Mol, width: int = 300, height: int = 300) -> str:
        """Generate SVG depiction of molecule"""
        AllChem.Compute2DCoords(mol)  # Generate 2D coordinates
        drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    @staticmethod
    def mol_to_smiles(mol: Chem.Mol) -> str:
        """Convert RDKit Mol to canonical SMILES"""
        return Chem.MolToSmiles(mol)
```

**Integration with SAEngine:**
```python
# backend/app/core/sa_engine.py
class SAEngine:
    def __init__(self, params: SAParams, molecule_service: MoleculeService, target_fn_service: TargetFunctionService):
        self.params = params
        self.mol_service = molecule_service
        self.target_fn = target_fn_service
        self.current_mol = molecule_service.create_from_formula(params.formula)
        self.current_energy = target_fn_service.evaluate(self.current_mol, params.target_components)
        ...

    def step(self):
        # Select random atom pair
        i, j = random_atom_pair(self.current_mol)

        # Attempt displacement
        proposed_mol = self.mol_service.apply_displacement(self.current_mol, i, j)
        if proposed_mol is None:
            self.invalid_moves += 1
            return

        # Evaluate new energy
        proposed_energy = self.target_fn.evaluate(proposed_mol, self.params.target_components)

        # Metropolis acceptance
        delta_e = proposed_energy - self.current_energy  # Assuming MINIMIZE
        if self.metropolis_accept(delta_e, self.temperature):
            self.current_mol = proposed_mol
            self.current_energy = proposed_energy
            ...
```

**Trade-offs:**
- ✅ RDKit handles valence, aromaticity, sanitization automatically
- ✅ Rich descriptor library (logP, MW, TPSA, etc.)
- ✅ High-quality 2D/3D coordinate generation
- ⚠️ RDKit is synchronous (blocks event loop; use `asyncio.to_thread()` if needed)
- ⚠️ Memory: RDKit Mol objects are not trivial; avoid keeping full history in memory

## Integration Points

### New Components vs Modified Components

| Component | Status | Changes |
|-----------|--------|---------|
| **Alpine.js app.ts** | MODIFIED | Replace Comlink worker calls with `fetch()` + `EventSource` |
| **chart.ts** | MODIFIED | Accept progress data from SSE instead of worker callbacks |
| **molecule-renderer.ts** | MODIFIED | Render SVG from backend instead of canvas rendering with RDKit.js |
| **presets.ts** | UNCHANGED | Same preset molecules |
| **MolGraph.ts** | REMOVED | Replaced by Python `MoleculeService` (RDKit `Chem.Mol`) |
| **SAEngine.ts** | REMOVED | Replaced by Python `SAEngine` |
| **sa-worker.ts** | REMOVED | Replaced by FastAPI backend |
| **API Client** | NEW | HTTP client wrapper for `/api/sa/*` endpoints |
| **SSE Client** | NEW | EventSource wrapper with reconnection logic |
| **FastAPI app** | NEW | Backend application with routing, CORS, startup |
| **SAService** | NEW | SA orchestration with background tasks |
| **MoleculeService** | NEW | RDKit wrapper (Python equivalent of MolGraph) |
| **TargetFunctionService** | NEW | Multi-component scoring framework |
| **Job Store** | NEW | In-memory job state and progress queue management |

### External Service Integration

| Service | Integration Pattern | Notes |
|---------|---------------------|-------|
| **RDKit (Python)** | Direct import (`from rdkit import Chem`) | Pre-installed in backend environment (conda or pip) |
| **FastAPI** | ASGI server (uvicorn) | `uvicorn app.main:app --reload` for dev |
| **Vite Dev Server** | Proxy to backend in dev mode | `vite.config.ts`: `proxy: { '/api': 'http://localhost:8000' }` |
| **GitHub Pages** | Static frontend only | Backend deployed separately (e.g., Render, Railway, Fly.io) |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| **Frontend ↔ Backend** | HTTP (REST) + SSE | CORS enabled for dev (localhost:5173) and production (GitHub Pages) |
| **API Router ↔ Service** | Direct function calls (dependency injection) | FastAPI `Depends()` for service injection |
| **Service ↔ Core** | Direct function calls | Services are thin orchestration layer, Core contains SA logic |
| **SAEngine ↔ MoleculeService** | Composition (service passed to constructor) | Engine delegates all RDKit operations to MoleculeService |

## Architectural Patterns

### Pattern 1: Server-Sent Events (SSE) for Progress Streaming

**What:** Unidirectional HTTP streaming from server to client for real-time progress updates

**When to use:**
- Server needs to push frequent updates to client (every N steps)
- Client does not need to send data back during stream (pause/cancel via separate POST)
- Simpler than WebSockets (no bidirectional protocol handshake)

**Trade-offs:**
- ✅ Built-in browser API (`EventSource`)
- ✅ Automatic reconnection on disconnect
- ✅ Works over HTTP (no WebSocket firewall issues)
- ✅ Scales well with async Python (`asyncio.Queue` + generators)
- ⚠️ Unidirectional (need separate endpoint for client → server commands)
- ⚠️ HTTP/1.1 connection limit (6 per domain; use HTTP/2)
- ⚠️ No binary data (text-only; base64 encode if needed)

**Example:**
```python
# Backend: SSE endpoint
from sse_starlette.sse import EventSourceResponse

@router.get("/stream/{job_id}")
async def stream_progress(job_id: str):
    async def event_generator():
        queue = job_store.get_progress_queue(job_id)
        try:
            while True:
                event = await queue.get()
                yield {
                    "event": event["type"],  # "progress" | "complete" | "error"
                    "data": json.dumps(event["data"])
                }
                if event["type"] in ("complete", "error"):
                    break
        except asyncio.CancelledError:
            # Client disconnected
            pass

    return EventSourceResponse(event_generator())

# Frontend: EventSource client
const eventSource = new EventSource(`/api/sa/stream/${jobId}`);
eventSource.addEventListener('progress', (e) => {
  const data = JSON.parse(e.data);
  updateChart(data.step, data.bestEnergy);
  updateMolecule(data.bestSvg);
});
eventSource.addEventListener('complete', (e) => {
  const result = JSON.parse(e.data);
  displayFinalResult(result);
  eventSource.close();
});
```

### Pattern 2: Dependency Injection with FastAPI `Depends()`

**What:** Service instances injected into route handlers via FastAPI's DI system

**When to use:**
- Services need shared state (e.g., job store, RDKit session)
- Want to mock services in tests
- Enforce layered architecture (API → Service → Core)

**Trade-offs:**
- ✅ Testable (easy to mock dependencies)
- ✅ Explicit dependencies (clear from function signature)
- ✅ Lazy initialization (services created on-demand)
- ✅ Scoped lifetimes (request, application, singleton)
- ⚠️ Learning curve (FastAPI DI syntax)
- ⚠️ Runtime overhead (minimal, but not zero)

**Example:**
```python
# backend/app/api/sa.py
from fastapi import Depends
from app.services.sa_service import SAService, get_sa_service
from app.state.job_store import JobStore, get_job_store

@router.post("/optimize")
async def optimize(
    params: SAParams,
    background_tasks: BackgroundTasks,
    sa_service: SAService = Depends(get_sa_service),
    job_store: JobStore = Depends(get_job_store)
):
    job_id = generate_job_id()
    job_store.create_job(job_id, params)
    background_tasks.add_task(sa_service.run_optimization, job_id, params)
    return {"job_id": job_id}

# backend/app/services/sa_service.py
_sa_service_instance = None

def get_sa_service() -> SAService:
    global _sa_service_instance
    if _sa_service_instance is None:
        _sa_service_instance = SAService(
            molecule_service=get_molecule_service(),
            target_fn_service=get_target_fn_service()
        )
    return _sa_service_instance
```

### Pattern 3: Background Tasks with `BackgroundTasks`

**What:** Run long-running SA optimization in background without blocking HTTP response

**When to use:**
- Operation takes longer than reasonable HTTP timeout (~30s)
- Client needs immediate response (job_id) to start polling/streaming
- Want to return HTTP 202 Accepted (operation started, not complete)

**Trade-offs:**
- ✅ Non-blocking (HTTP response returns immediately)
- ✅ Built into FastAPI (no external task queue needed for simple cases)
- ✅ Runs in same process (no serialization overhead)
- ⚠️ Not persistent (lost on server restart; use Celery for production)
- ⚠️ No retry/failure handling (add manually or use Celery)
- ⚠️ Single server only (does not scale horizontally; use Celery for multi-worker)

**Example:**
```python
@router.post("/optimize")
async def optimize(params: SAParams, background_tasks: BackgroundTasks):
    job_id = generate_job_id()
    background_tasks.add_task(sa_service.run_optimization, job_id, params)
    return {"job_id": job_id, "status": "started"}

# backend/app/services/sa_service.py
async def run_optimization(self, job_id: str, params: SAParams):
    try:
        engine = SAEngine(params, self.molecule_service, self.target_fn_service)
        engine.init()

        queue = job_store.get_progress_queue(job_id)

        while not engine.is_complete():
            engine.step()

            # Emit progress every 10 steps
            if engine.step_number % 10 == 0:
                state = engine.get_state()
                await queue.put({
                    "type": "progress",
                    "data": {
                        "step": state.step,
                        "bestEnergy": state.best_energy,
                        "bestSvg": self.molecule_service.generate_svg(state.best_mol)
                    }
                })

        # Emit completion
        result = engine.get_result()
        await queue.put({
            "type": "complete",
            "data": {
                "bestEnergy": result.best_energy,
                "bestSmiles": self.molecule_service.mol_to_smiles(result.best_mol)
            }
        })
    except Exception as e:
        await queue.put({"type": "error", "data": {"message": str(e)}})
```

### Pattern 4: Pydantic Models for Type Safety

**What:** Use Pydantic models for request/response validation and OpenAPI schema generation

**When to use:**
- Always (FastAPI best practice)
- Shared types between frontend/backend (generate TS types from Pydantic)
- Input validation with helpful error messages

**Trade-offs:**
- ✅ Automatic validation (type checking, range checks, regex patterns)
- ✅ OpenAPI schema generation (auto-docs at `/docs`)
- ✅ Type hints for IDE autocomplete
- ✅ JSON serialization/deserialization
- ⚠️ Runtime overhead (validation on every request; minimal but not zero)
- ⚠️ Learning curve (Pydantic v2 syntax)

**Example:**
```python
# backend/app/models/sa_params.py
from pydantic import BaseModel, Field
from typing import Dict

class SAParams(BaseModel):
    formula: str = Field(..., pattern=r"^([A-Z][a-z]?\d*)+$", example="C6H14")
    initial_temp: float = Field(100.0, gt=0, description="Initial temperature (kT)")
    cooling_schedule_k: int = Field(8, ge=0, le=32, description="Cooling schedule index")
    steps_per_cycle: int = Field(500, gt=0)
    num_cycles: int = Field(4, gt=0)
    optimization_mode: str = Field("MINIMIZE", pattern="^(MINIMIZE|MAXIMIZE)$")
    target_components: Dict[str, float] = Field(
        {"wiener_index": 1.0},
        description="Weighted scoring components"
    )

class SAResult(BaseModel):
    best_energy: float
    best_smiles: str
    best_svg: str
    final_energy: float
    initial_energy: float
    total_steps: int
    accepted_moves: int
    rejected_moves: int
    invalid_moves: int
    acceptance_ratio: float

# Generate TypeScript types:
# npm install -g openapi-typescript
# openapi-typescript http://localhost:8000/openapi.json -o frontend/src/types/api.ts
```

## Development Workflow

### Local Development Setup

1. **Backend (Terminal 1):**
   ```bash
   cd backend
   python -m venv venv
   source venv/bin/activate  # or `venv\Scripts\activate` on Windows
   pip install -r requirements.txt
   uvicorn app.main:app --reload --port 8000
   ```

2. **Frontend (Terminal 2):**
   ```bash
   cd frontend
   npm install
   npm run dev  # Vite dev server on http://localhost:5173
   ```

3. **Vite Proxy Configuration:**
   ```typescript
   // frontend/vite.config.ts
   export default defineConfig({
     server: {
       proxy: {
         '/api': {
           target: 'http://localhost:8000',
           changeOrigin: true
         }
       }
     }
   })
   ```

4. **CORS Configuration:**
   ```python
   # backend/app/main.py
   from fastapi.middleware.cors import CORSMiddleware

   app.add_middleware(
       CORSMiddleware,
       allow_origins=[
           "http://localhost:5173",  # Vite dev server
           "https://steinbeck.github.io"  # GitHub Pages
       ],
       allow_credentials=True,
       allow_methods=["*"],
       allow_headers=["*"]
   )
   ```

### Build Order (Dependency-Aware)

**Phase 1: Core Backend (No Frontend Changes)**
1. Project structure setup (create `backend/` directory)
2. Pydantic models (`SAParams`, `SAResult`, `ProgressEvent`)
3. `MoleculeService` (RDKit wrapper)
4. `SAEngine` (Python port of TypeScript version)
5. Unit tests for core components

**Phase 2: API Layer**
6. FastAPI app setup (CORS, startup/shutdown)
7. Job store (in-memory state management)
8. `SAService` (orchestration + background tasks)
9. API routes (`POST /optimize`, `GET /stream/{job_id}`)
10. Integration tests for API

**Phase 3: Frontend Integration**
11. API client service (`services/api-client.ts`)
12. SSE client service (`services/sse-client.ts`)
13. Modify `app.ts` to use API client instead of Web Worker
14. Modify `chart.ts` to consume SSE events
15. Modify `molecule-renderer.ts` to display SVG from backend

**Phase 4: Multi-Component Target Function**
16. `TargetFunctionService` (plugin registry)
17. Register Wiener Index component
18. Add additional components (logP, MW, etc.)
19. Update frontend to configure component weights

**Phase 5: Polish & Deploy**
20. Error handling (invalid formulas, job not found, etc.)
21. Loading states and progress indicators
22. Backend deployment (Render, Railway, Fly.io)
23. Update frontend to point to production backend URL

## Scaling Considerations

| Scale | Architecture | Notes |
|-------|--------------|-------|
| **0-10 concurrent users** | Single FastAPI server + in-memory job store | Current design sufficient |
| **10-100 concurrent users** | Add Redis for job store + multiple FastAPI workers | Replace `Dict[str, JobState]` with Redis |
| **100+ concurrent users** | Celery for background tasks + Redis + load balancer | Replace `BackgroundTasks` with Celery |

### First Bottleneck: In-Memory Job Store

**What breaks:** Multiple FastAPI workers don't share memory; job started on worker 1 can't be streamed from worker 2

**Fix:** Replace in-memory `Dict` with Redis
```python
# backend/app/state/job_store.py (Redis version)
import redis.asyncio as redis
import json

class JobStore:
    def __init__(self):
        self.redis = redis.from_url("redis://localhost")

    async def create_job(self, job_id: str, params: SAParams):
        await self.redis.set(f"job:{job_id}", json.dumps(params.dict()))

    async def push_progress(self, job_id: str, event: dict):
        await self.redis.lpush(f"progress:{job_id}", json.dumps(event))

    async def stream_progress(self, job_id: str):
        while True:
            event_json = await self.redis.brpop(f"progress:{job_id}", timeout=1)
            if event_json:
                yield json.loads(event_json[1])
```

### Second Bottleneck: Background Tasks in HTTP Process

**What breaks:** Long-running SA jobs block FastAPI workers (limited concurrency)

**Fix:** Offload to Celery workers
```python
# backend/app/tasks.py (Celery)
from celery import Celery

celery_app = Celery('webfaulon', broker='redis://localhost')

@celery_app.task
def run_optimization(job_id: str, params: dict):
    sa_service.run_optimization(job_id, SAParams(**params))

# backend/app/api/sa.py
@router.post("/optimize")
async def optimize(params: SAParams):
    job_id = generate_job_id()
    run_optimization.delay(job_id, params.dict())
    return {"job_id": job_id}
```

## Anti-Patterns

### Anti-Pattern 1: Blocking RDKit Operations in Async Handlers

**What people do:** Call synchronous RDKit functions directly in `async def` route handlers

**Why it's wrong:** Blocks the event loop, prevents other requests from being processed, reduces concurrency

**Do this instead:**
```python
# ❌ BAD: Blocks event loop
@router.get("/molecule/{smiles}")
async def get_molecule(smiles: str):
    mol = Chem.MolFromSmiles(smiles)  # Synchronous RDKit call
    svg = molecule_service.generate_svg(mol)  # Synchronous
    return {"svg": svg}

# ✅ GOOD: Run in thread pool
from asyncio import to_thread

@router.get("/molecule/{smiles}")
async def get_molecule(smiles: str):
    mol = await to_thread(Chem.MolFromSmiles, smiles)
    svg = await to_thread(molecule_service.generate_svg, mol)
    return {"svg": svg}

# ✅ BETTER: For SA engine, run entire step() loop in background task
background_tasks.add_task(sa_service.run_optimization, job_id, params)
```

### Anti-Pattern 2: Returning Full Molecule History in API Response

**What people do:** Include entire SA history (`history: SAStepResult[]`) in `/stream/{job_id}` or `/status/{job_id}` response

**Why it's wrong:**
- JSON payload grows linearly with steps (500 steps/cycle × 4 cycles = 2000 objects)
- Each object includes MOL block or SVG (kilobytes per step)
- SSE streams send duplicate data (step 100 includes steps 1-100)

**Do this instead:**
```python
# ❌ BAD: Send full history on every progress event
{
  "step": 100,
  "history": [
    {"step": 1, "energy": 50, "molBlock": "..."},
    {"step": 2, "energy": 49, "molBlock": "..."},
    ...
    {"step": 100, "energy": 30, "molBlock": "..."}
  ]
}

# ✅ GOOD: Send only current step data
{
  "step": 100,
  "currentEnergy": 30,
  "bestEnergy": 25,
  "bestSvg": "...",  # Only best molecule, not all
  "temperature": 8.5
}

# ✅ Client aggregates history locally
const history = [];
eventSource.onmessage = (e) => {
  const data = JSON.parse(e.data);
  history.push({step: data.step, energy: data.bestEnergy});
  updateChart(history);
};
```

### Anti-Pattern 3: Not Closing SSE Streams

**What people do:** Forget to send completion event or close EventSource on client

**Why it's wrong:**
- Server-side generator runs forever, holding resources
- Client keeps connection open, counts against HTTP/1.1 limit (6 per domain)
- Memory leak (progress queue grows indefinitely)

**Do this instead:**
```python
# ✅ Server: Always send completion event
async def event_generator():
    try:
        while True:
            event = await queue.get()
            yield {"data": json.dumps(event)}
            if event["type"] in ("complete", "error"):
                break  # Exit generator
    except asyncio.CancelledError:
        pass  # Client disconnected

# ✅ Client: Close on completion
eventSource.addEventListener('complete', () => {
  eventSource.close();
});
eventSource.addEventListener('error', () => {
  eventSource.close();
});
```

### Anti-Pattern 4: Using `*` for CORS `allow_origins` with `allow_credentials=True`

**What people do:**
```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True  # ⚠️ Browser rejects this!
)
```

**Why it's wrong:** Browsers reject wildcard origins when credentials are enabled (security violation)

**Do this instead:**
```python
# ✅ Explicit origins list
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:5173",
        "https://steinbeck.github.io"
    ],
    allow_credentials=True
)

# ✅ Or disable credentials if not needed
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=False
)
```

## Sources

### FastAPI + SSE Streaming
- [Implementing Server-Sent Events (SSE) with FastAPI](https://mahdijafaridev.medium.com/implementing-server-sent-events-sse-with-fastapi-real-time-updates-made-simple-6492f8bfc154)
- [How to use Server-Sent Events with FastAPI and React](https://www.softgrade.org/sse-with-fastapi-react-langgraph/)
- [FastAPI + SSE for LLM Tokens: Smooth Streaming without WebSockets](https://medium.com/@hadiyolworld007/fastapi-sse-for-llm-tokens-smooth-streaming-without-websockets-001ead4b5e53)
- [Streaming Response In FastAPI](https://medium.com/@ab.hassanein/streaming-responses-in-fastapi-d6a3397a4b7b)
- [Custom Response - HTML, Stream, File, others - FastAPI](https://fastapi.tiangolo.com/advanced/custom-response/)

### Python RDKit
- [Getting Started with the RDKit in Python — The RDKit 2025.09.5 documentation](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [RDKit Cookbook — The RDKit 2025.09.5 documentation](https://www.rdkit.org/docs/Cookbook.html)
- [rdkit.Chem.Draw package — The RDKit 2025.09.5 documentation](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.html)
- [rdkit.Chem.GraphDescriptors module — The RDKit 2025.09.4 documentation](https://www.rdkit.org/docs/source/rdkit.Chem.GraphDescriptors.html)
- [Revisiting a Classic Cheminformatics Paper: The Wiener Index](https://bertiewooster.github.io/2023/03/10/Revisiting-a-Classic-Cheminformatics-Paper-The-Wiener-Index.html)

### FastAPI Architecture & Patterns
- [Mastering Dependency Injection in FastAPI](https://medium.com/@azizmarzouki/mastering-dependency-injection-in-fastapi-clean-scalable-and-testable-apis-5f78099c3362)
- [How to Use Dependency Injection in FastAPI](https://oneuptime.com/blog/post/2026-02-02-fastapi-dependency-injection/view)
- [FastAPI Complete Guide: Build High-Performance Python APIs in 2026](https://devtoolbox.dedyn.io/blog/fastapi-complete-guide)
- [Combining FastAPI Dependency Injection with Service and Repository Layers](https://blog.dotcs.me/posts/fastapi-dependency-injection-x-layers)
- [The Service Layer Pattern - Marc Puig - Notes](https://mpuig.github.io/Notes/fastapi_basics/04.service_layer_pattern/)

### Monorepo Structure
- [How do you typically structure your project if it includes both frontend and fastapi?](https://github.com/fastapi/fastapi/discussions/4344)
- [Embedding a React Frontend Inside a FastAPI Python Package (in a Monorepo)](https://medium.com/@asafshakarzy/embedding-a-react-frontend-inside-a-fastapi-python-package-in-a-monorepo-c00f99e90471)
- [Generating API clients in monorepos with FastAPI & Next.js](https://www.vintasoftware.com/blog/nextjs-fastapi-monorepo)
- [A vertical monorepo architecture for FastAPI client-server codebases](https://sqr-075.lsst.io/)

### CORS & Development Setup
- [CORS (Cross-Origin Resource Sharing) - FastAPI](https://fastapi.tiangolo.com/tutorial/cors/)
- [FastAPI: Configuring CORS for Python's ASGI Framework](https://www.stackhawk.com/blog/configuring-cors-in-fastapi/)
- [Setup Proxy in React Vite for CORs issues](https://medium.com/@aparna1002v/setup-proxy-in-react-vite-for-cors-issues-167c6a1eb569)
- [Blocked by CORS in FastAPI? Here's How to Fix It](https://davidmuraya.com/blog/fastapi-cors-configuration/)

---
*Architecture research for: FastAPI + Python RDKit Backend Integration*
*Researched: 2026-02-15*
*Confidence: HIGH (verified with official docs + recent 2026 sources)*
