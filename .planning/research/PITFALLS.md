# Pitfalls Research

**Domain:** Browser-to-Python Backend Migration (RDKit + FastAPI + SSE)
**Researched:** 2026-02-15
**Confidence:** HIGH

---

## Critical Pitfalls

### Pitfall 1: Missing Sanitization After RDKit Bond Order Changes

**What goes wrong:**

After performing Faulon displacement (modifying bond orders on 4 atoms using `SetBondOrder`), failing to call `Chem.SanitizeMol()` results in molecules with incorrect valence states, broken aromaticity perception, and invalid hydrogen counts. RDKit will generate incorrect SMILES strings and may crash on subsequent operations.

**Why it happens:**

Developers assume RDKit auto-validates like the TypeScript `MolGraph` class (which recomputes `implicitH` in `setBond()`). In RDKit, bond modifications create an intermediate "dirty" state—valence, aromaticity, and hybridization are NOT automatically updated. The molecule appears valid but is chemically inconsistent.

**How to avoid:**

```python
from rdkit import Chem
from rdkit.Chem import rdMolOps

# After any bond order modification
mol = Chem.RWMol(...)
mol.GetBondBetweenAtoms(x1, y1).SetBondType(Chem.BondType.DOUBLE)
mol.GetBondBetweenAtoms(y1, y2).SetBondType(Chem.BondType.SINGLE)
mol.GetBondBetweenAtoms(x1, x2).SetBondType(Chem.BondType.SINGLE)
mol.GetBondBetweenAtoms(x2, y2).SetBondType(Chem.BondType.DOUBLE)

# CRITICAL: Sanitize to recompute valence, aromaticity, etc.
try:
    Chem.SanitizeMol(mol)
except Exception as e:
    # Displacement created invalid molecule (wrong valence)
    return None  # Reject this move

# Now safe to use mol.GetMol() or generate SMILES
```

**Alternative: Partial Sanitization** (if full sanitization fails on valid structures):

```python
# More granular control
from rdkit.Chem import SanitizeFlags

flags = SanitizeFlags.SANITIZE_ALL ^ SanitizeFlags.SANITIZE_KEKULIZE
Chem.SanitizeMol(mol, sanitizeOps=flags)
```

**Warning signs:**

- `Explicit valence for atom #X Y, Z, is greater than permitted` errors
- SMILES strings differ from expected (e.g., `C1CC1` instead of `C=CC`)
- Atom valence queries return stale values
- Kekulization errors on aromatic systems
- Test failures with message "Sanitization error"

**Phase to address:**

**Phase 1: RDKit Backend Core** — Implement displacement with mandatory post-sanitization. Write unit tests that verify SMILES correctness after displacement (catches sanitization omissions).

---

### Pitfall 2: CPU-Bound SA Loop Blocking Async FastAPI

**What goes wrong:**

Running the SA loop (thousands of `engine.step()` calls) inside an `async def` endpoint with `await` on synchronous code blocks the entire event loop. FastAPI becomes unresponsive—health checks timeout, new requests queue, SSE connections drop. The server appears hung.

**Why it happens:**

FastAPI's async is I/O concurrency, not parallelism. The SA loop is CPU-bound (RDKit molecule operations, Wiener index computation). Running it in the main event loop is like running the TypeScript `SAEngine.run()` in the browser UI thread instead of a Web Worker.

Python's Global Interpreter Lock (GIL) prevents threads from providing true parallelism for CPU-bound tasks. Developers incorrectly assume `async def` makes code non-blocking, but synchronous operations still block.

**How to avoid:**

**Option A: ProcessPoolExecutor** (TRUE parallelism, recommended for CPU-bound SA):

```python
from fastapi import FastAPI
from concurrent.futures import ProcessPoolExecutor
import asyncio

app = FastAPI()
executor = ProcessPoolExecutor(max_workers=4)

@app.post("/optimize")
async def run_optimization(params: SAParams):
    # Offload to separate process (no GIL contention)
    loop = asyncio.get_event_loop()
    result = await loop.run_in_executor(
        executor,
        run_sa_in_process,  # Function must be picklable
        params
    )
    return result

def run_sa_in_process(params: SAParams):
    # This runs in a separate process
    engine = SAEngine(params)
    return engine.run()  # Blocks process, not event loop
```

**Option B: Background Task + Polling** (simpler, no SSE):

```python
from fastapi import BackgroundTasks

tasks = {}  # In production, use Redis/DB

@app.post("/optimize")
async def start_optimization(params: SAParams, background: BackgroundTasks):
    task_id = uuid4()
    background.add_task(run_sa_sync, task_id, params)
    return {"task_id": task_id}

def run_sa_sync(task_id, params):
    # Runs in thread pool (still has GIL, but doesn't block event loop)
    result = SAEngine(params).run()
    tasks[task_id] = result

@app.get("/status/{task_id}")
async def get_status(task_id: str):
    return tasks.get(task_id, {"status": "running"})
```

**DO NOT do this** (blocks event loop):

```python
@app.post("/optimize")
async def run_optimization(params: SAParams):
    engine = SAEngine(params)
    # BAD: Synchronous CPU work in async function
    for _ in range(1000):
        engine.step()  # Blocks entire server
    return engine.getResult()
```

**Warning signs:**

- FastAPI health check endpoint timeouts
- `curl localhost:8000/docs` hangs during optimization
- SSE connections disconnect after 30-60 seconds
- `uvicorn` logs show long request durations (>1s for simple endpoints)
- CPU pegged at 100% on single core while other cores idle

**Phase to address:**

**Phase 1: RDKit Backend Core** — Design SA engine to run in ProcessPoolExecutor from day 1. Phase 2 (SSE streaming) depends on this foundation.

---

### Pitfall 3: SSE Backpressure Causing Memory Exhaustion

**What goes wrong:**

SA generates progress updates every 10 steps (potentially 100+ events/second). If the client consumes SSE slowly (slow network, browser throttling, or client-side processing lag), the server's event queue grows unbounded. Memory usage spikes, the server OOMs, or SSE connections drop due to buffer overflow.

**Why it happens:**

FastAPI's `StreamingResponse` with async generators doesn't automatically throttle if the client can't keep up. The Web Worker version had implicit backpressure (main thread processes `onProgress` before worker continues), but SSE is fire-and-forget.

Developers assume SSE "just works" like WebSockets (which have flow control). SSE is one-way HTTP—no client ACKs. If you `yield` faster than the client reads, events queue in OS buffers, then in Python memory.

**How to avoid:**

**Adaptive Reporting Interval** (best for SA use case):

```python
async def sse_progress_stream(params: SAParams):
    engine = SAEngine(params)
    last_report_time = time.time()
    MIN_REPORT_INTERVAL = 0.1  # Max 10 updates/sec

    for step in range(total_steps):
        engine.step()

        now = time.time()
        if now - last_report_time >= MIN_REPORT_INTERVAL:
            yield f"data: {json.dumps(engine.getState())}\n\n"
            last_report_time = now
            await asyncio.sleep(0)  # Yield control to event loop

    # Always send final state
    yield f"data: {json.dumps(engine.getResult())}\n\n"
```

**Queue-Based Throttling** (if client processing time varies):

```python
from asyncio import Queue

async def sse_with_bounded_queue(params: SAParams):
    queue = Queue(maxsize=10)  # Backpressure if client falls behind

    async def producer():
        engine = SAEngine(params)
        for step in range(total_steps):
            engine.step()
            if step % 10 == 0:
                # Put blocks when queue full (backpressure)
                await queue.put(engine.getState())
        await queue.put(None)  # Sentinel

    async def consumer():
        while True:
            state = await queue.get()
            if state is None:
                break
            yield f"data: {json.dumps(state)}\n\n"

    asyncio.create_task(producer())
    async for event in consumer():
        yield event
```

**Monitor Event Queue Size**:

```python
import sys

async def sse_stream(params: SAParams):
    events_queued = 0
    MAX_QUEUE_SIZE = 100

    for state in generate_states():
        if events_queued > MAX_QUEUE_SIZE:
            # Client too slow, skip this update
            continue

        yield f"data: {json.dumps(state)}\n\n"
        events_queued = sys.getsizeof(...)  # Estimate pending data
```

**Warning signs:**

- Server memory usage grows linearly during optimization
- SSE connection drops mid-stream with no error
- Client receives bursts of events (20+ at once) instead of steady stream
- `uvicorn` logs show "send buffer full" warnings
- Slower clients (mobile, throttled networks) fail more often

**Phase to address:**

**Phase 2: SSE Streaming** — Implement adaptive reporting interval from start. Add queue monitoring in Phase 3 (Multi-component targets) when update frequency increases.

---

### Pitfall 4: SSE Client Disconnect Not Canceling SA Computation

**What goes wrong:**

User closes browser tab or navigates away. SSE connection drops, but the FastAPI endpoint continues running the SA loop for 10+ minutes, wasting CPU and blocking ProcessPoolExecutor workers. Server capacity degrades as zombie tasks accumulate.

**Why it happens:**

FastAPI's `StreamingResponse` doesn't auto-cancel the generator when the client disconnects. The async generator continues yielding events into the void. Unlike Web Workers (which have explicit termination), SSE disconnects are passive—the server doesn't know unless it actively checks.

**How to avoid:**

**Monitor `request.is_disconnected()`** (FastAPI 0.100+):

```python
from fastapi import Request
from fastapi.responses import StreamingResponse

@app.post("/optimize")
async def stream_optimization(params: SAParams, request: Request):
    async def event_generator():
        engine = SAEngine(params)

        for step in range(total_steps):
            # Check for disconnect every N steps
            if step % 100 == 0 and await request.is_disconnected():
                # Clean up and exit
                engine.cleanup()
                return

            engine.step()

            if step % 10 == 0:
                yield f"data: {json.dumps(engine.getState())}\n\n"

        yield f"data: {json.dumps(engine.getResult())}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")
```

**Using `anyio.TaskGroup` for Cancellation**:

```python
import anyio

@app.post("/optimize")
async def stream_optimization(params: SAParams, request: Request):
    async def event_generator():
        async def monitor_disconnect(cancel_scope):
            while True:
                if await request.is_disconnected():
                    cancel_scope.cancel()
                    return
                await asyncio.sleep(1)

        async def run_sa():
            engine = SAEngine(params)
            for step in range(total_steps):
                engine.step()
                if step % 10 == 0:
                    yield f"data: {json.dumps(engine.getState())}\n\n"
                await asyncio.sleep(0)  # Check for cancellation

        async with anyio.create_task_group() as tg:
            tg.start_soon(monitor_disconnect, tg.cancel_scope)
            async for event in run_sa():
                yield event

    return StreamingResponse(event_generator(), media_type="text/event-stream")
```

**For ProcessPoolExecutor** (harder—processes can't be cancelled):

```python
# Use timeout and shared state
from multiprocessing import Manager

manager = Manager()
cancellation_flags = manager.dict()

def run_sa_cancellable(task_id, params):
    engine = SAEngine(params)
    for step in range(total_steps):
        if cancellation_flags.get(task_id, False):
            return {"cancelled": True}
        engine.step()
    return engine.getResult()

@app.post("/optimize")
async def stream_optimization(params: SAParams, request: Request):
    task_id = str(uuid4())
    cancellation_flags[task_id] = False

    loop = asyncio.get_event_loop()
    future = loop.run_in_executor(executor, run_sa_cancellable, task_id, params)

    # Monitor disconnect
    while not future.done():
        if await request.is_disconnected():
            cancellation_flags[task_id] = True
            break
        await asyncio.sleep(0.5)

    result = await future
    return result
```

**Warning signs:**

- CPU usage stays high after closing browser tabs
- ProcessPoolExecutor workers show 100% usage but no active clients
- Server capacity degrades over time (workers stuck on abandoned tasks)
- Memory usage doesn't drop after clients disconnect
- `ps aux | grep python` shows more SA processes than active connections

**Phase to address:**

**Phase 2: SSE Streaming** — Implement disconnect detection. Test by killing browser mid-optimization and verifying server logs show cancellation.

---

### Pitfall 5: CORS Middleware Ordering Breaks Error Responses

**What goes wrong:**

During development, CORS is configured correctly, but when FastAPI raises exceptions (404, 422 validation errors, 500 internal errors), the browser shows "CORS policy: No 'Access-Control-Allow-Origin' header" instead of the actual error message. This hides bugs and makes debugging impossible.

**Why it happens:**

FastAPI middleware runs in order. If `CORSMiddleware` is added AFTER another middleware (e.g., authentication, error handling), and that middleware raises an exception, the response bypasses `CORSMiddleware`. The browser receives an error response WITHOUT CORS headers and blocks it.

This is especially problematic with the Vite dev server at `http://localhost:5173` calling FastAPI at `http://localhost:8000`.

**How to avoid:**

**ALWAYS add CORSMiddleware FIRST**:

```python
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

# CRITICAL: Add CORS middleware BEFORE any other middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173"],  # Vite dev server
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Other middleware AFTER CORS
app.add_middleware(SomeAuthMiddleware)
app.add_middleware(LoggingMiddleware)

@app.get("/optimize")
async def optimize():
    raise ValueError("Something broke")  # Will include CORS headers
```

**Development-Specific Origins**:

```python
import os

ALLOWED_ORIGINS = [
    "http://localhost:5173",  # Vite
    "http://127.0.0.1:5173",  # Alternative localhost
]

if os.getenv("ENV") == "production":
    ALLOWED_ORIGINS = ["https://webfaulon.example.com"]

app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    # ...
)
```

**Alternative: Vite Proxy** (avoids CORS entirely in dev):

```javascript
// vite.config.ts
export default {
  server: {
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
        rewrite: (path) => path.replace(/^\/api/, '')
      }
    }
  }
}

// Frontend calls /api/optimize instead of http://localhost:8000/optimize
// Browser sees same-origin request, no CORS needed
```

**Warning signs:**

- Browser console shows "CORS policy" errors instead of actual HTTP errors
- `curl` to FastAPI works, but browser requests fail
- Authentication errors return "No Access-Control-Allow-Origin"
- 422 validation errors are invisible in browser (show generic CORS error)
- Works in Postman/Insomnia but not in browser

**Phase to address:**

**Phase 1: RDKit Backend Core** — Configure CORS correctly from start. Add to README so it's not forgotten when adding auth middleware later.

---

### Pitfall 6: RDKit Molecule Memory Leaks in Long-Running SA

**What goes wrong:**

During thousands of SA iterations, memory usage grows unbounded. The server OOMs after processing 10-20 optimization jobs. Each job completes successfully, but memory isn't released.

**Why it happens:**

RDKit's Python bindings use C++ objects with complex lifetimes. When storing `Mol` objects in lists (e.g., `history: List[Mol]`), circular references or delayed garbage collection can prevent cleanup. Python's GC doesn't immediately release C++ resources.

Additionally, RDKit's `RWMol` and `EditableMol` may hold references to parent molecules, creating implicit retention.

**How to avoid:**

**Don't Store Full Molecules in History**:

```python
# BAD: Stores thousands of Mol objects
class SAEngine:
    def __init__(self):
        self.history: List[Mol] = []

    def step(self):
        proposed_mol = self.current_mol.clone()  # New C++ object
        self.history.append(proposed_mol)  # Leaks

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

**Explicit Cleanup for RWMol**:

```python
from rdkit import Chem

def attempt_displacement(mol: Mol) -> Mol:
    rw_mol = Chem.RWMol(mol)

    # Modify bonds
    rw_mol.GetBondBetweenAtoms(x1, y1).SetBondType(Chem.BondType.DOUBLE)

    try:
        Chem.SanitizeMol(rw_mol)
        result = rw_mol.GetMol()  # Convert to immutable Mol
        del rw_mol  # Explicit cleanup (helps CPython)
        return result
    except:
        del rw_mol
        return None
```

**Batch Processing with Cleanup**:

```python
def run_sa_batched(params: SAParams):
    engine = SAEngine(params)

    for batch_start in range(0, total_steps, 1000):
        for _ in range(1000):
            engine.step()

        # Periodic cleanup
        gc.collect()  # Force Python GC

        # Yield intermediate result
        yield engine.getState()

    # Final cleanup
    del engine
    gc.collect()
```

**Monitor Memory**:

```python
import tracemalloc
import gc

tracemalloc.start()

@app.post("/optimize")
async def optimize(params: SAParams):
    snapshot_before = tracemalloc.take_snapshot()

    result = run_sa(params)

    snapshot_after = tracemalloc.take_snapshot()
    top_stats = snapshot_after.compare_to(snapshot_before, 'lineno')

    # Log top memory increases
    for stat in top_stats[:10]:
        print(stat)

    return result
```

**Warning signs:**

- Server memory grows with each completed job (doesn't return to baseline)
- `docker stats` shows container memory at 80%+ after 10 jobs
- OOM kills in production logs
- Faster jobs (smaller molecules) trigger OOM sooner (more iterations/time)
- Memory profiling shows `rdkit.Chem` objects dominating

**Phase to address:**

**Phase 1: RDKit Backend Core** — Design history storage to use SMILES, not Mol objects. Add memory monitoring to tests (assert memory returns to baseline after `run()`).

---

### Pitfall 7: Incorrect RDKit BondType Enum Mapping

**What goes wrong:**

After Faulon displacement calculates new bond orders as integers (0, 1, 2, 3), naively mapping them to RDKit's `BondType` enum results in wrong bond types. Single bonds become aromatic, double bonds become single, etc. Molecules are chemically incorrect.

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

Directly passing integer to `SetBondType(bond_order)` causes misinterpretation.

**How to avoid:**

**Explicit Mapping Function**:

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

# Usage
bond = mol.GetBondBetweenAtoms(x1, y1)
bond.SetBondType(faulon_order_to_rdkit_bondtype(new_b11))
```

**Handle Bond Removal** (order 0):

```python
def apply_faulon_bonds(mol: RWMol, x1, y1, x2, y2, b11, b12, b21, b22):
    """Apply Faulon displacement bond orders to RDKit molecule."""

    # Helper to set or remove bond
    def set_bond_order(i, j, order):
        bond = mol.GetBondBetweenAtoms(i, j)

        if order == 0:
            if bond is not None:
                mol.RemoveBond(i, j)
        else:
            bond_type = faulon_order_to_rdkit_bondtype(order)
            if bond is None:
                mol.AddBond(i, j, bond_type)
            else:
                bond.SetBondType(bond_type)

    set_bond_order(x1, y1, b11)
    set_bond_order(y1, y2, b12)
    set_bond_order(x1, x2, b21)
    set_bond_order(x2, y2, b22)
```

**Validation Test**:

```python
def test_bond_type_mapping():
    """Ensure integer bond orders map correctly to RDKit BondType."""
    mol = Chem.MolFromSmiles("CC")
    bond = mol.GetBondBetweenAtoms(0, 1)

    # Original is single bond
    assert bond.GetBondType() == BondType.SINGLE

    # Change to double
    rw_mol = Chem.RWMol(mol)
    rw_mol.GetBondBetweenAtoms(0, 1).SetBondType(faulon_order_to_rdkit_bondtype(2))

    result_mol = rw_mol.GetMol()
    Chem.SanitizeMol(result_mol)

    # Verify SMILES reflects double bond
    assert Chem.MolToSmiles(result_mol) == "C=C"
```

**Warning signs:**

- SMILES strings don't match expected structures (e.g., `C1=CC1` instead of `C1CC1`)
- Aromatic bonds appear in saturated alkanes
- Bond order validation failures (`bond.GetBondTypeAsDouble()` returns unexpected values)
- Molecules with no bonds after displacement (all bonds accidentally removed)

**Phase to address:**

**Phase 1: RDKit Backend Core** — Implement explicit mapping function from day 1. Add SMILES validation tests for all bond types (no bond, single, double, triple).

---

### Pitfall 8: EventSource Auto-Reconnect Causing Duplicate Optimizations

**What goes wrong:**

User starts an optimization. Network hiccup causes SSE disconnect. Browser's `EventSource` auto-reconnects and sends another POST request, starting a SECOND optimization. Server now runs two identical jobs, doubling cost and confusing the user (two result sets).

**Why it happens:**

`EventSource` automatically reconnects when the connection drops (every 3 seconds by default). If the SSE endpoint is a POST request that STARTS the optimization (not just subscribes to it), each reconnect triggers a new job.

This differs from Web Workers, which have no implicit retry logic—a crashed worker stays dead until manually restarted.

**How to avoid:**

**Separate Job Creation from Streaming** (recommended):

```python
# 1. POST to create job (returns task_id)
@app.post("/optimize")
async def create_optimization(params: SAParams):
    task_id = str(uuid4())

    # Start background task
    background_tasks[task_id] = asyncio.create_task(run_sa(task_id, params))

    return {"task_id": task_id}

# 2. GET to stream progress via SSE
@app.get("/optimize/{task_id}/stream")
async def stream_progress(task_id: str, request: Request):
    async def event_generator():
        while True:
            state = await get_task_state(task_id)
            yield f"data: {json.dumps(state)}\n\n"

            if state["isComplete"]:
                break

            await asyncio.sleep(0.1)

    return StreamingResponse(event_generator(), media_type="text/event-stream")
```

**Client Side**:

```typescript
// 1. Start optimization
const response = await fetch('/api/optimize', {
  method: 'POST',
  body: JSON.stringify(params)
});
const { task_id } = await response.json();

// 2. Subscribe to progress (GET request, safe to reconnect)
const eventSource = new EventSource(`/api/optimize/${task_id}/stream`);
eventSource.onmessage = (event) => {
  const state = JSON.parse(event.data);
  updateUI(state);
};
```

**Idempotency with Client-Provided Task ID**:

```python
@app.post("/optimize/{client_task_id}")
async def create_optimization(client_task_id: str, params: SAParams):
    # Client generates UUID, ensures idempotency
    if client_task_id in background_tasks:
        return {"task_id": client_task_id, "status": "already_running"}

    background_tasks[client_task_id] = asyncio.create_task(run_sa(client_task_id, params))
    return {"task_id": client_task_id, "status": "created"}
```

**Last-Event-ID for Resume** (advanced):

```python
@app.get("/optimize/{task_id}/stream")
async def stream_progress(task_id: str, request: Request):
    last_event_id = request.headers.get("Last-Event-ID", "0")
    start_step = int(last_event_id)

    async def event_generator():
        async for state in get_states_from_step(task_id, start_step):
            yield f"id: {state['step']}\ndata: {json.dumps(state)}\n\n"

    return StreamingResponse(event_generator(), media_type="text/event-stream")
```

**Warning signs:**

- Server logs show multiple POST requests from same client for same params
- CPU usage doubles after brief network hiccup
- Client receives two complete result sets
- Database has duplicate job records with timestamps 3 seconds apart
- Users report "optimization finished twice"

**Phase to address:**

**Phase 2: SSE Streaming** — Separate job creation from streaming from start. Document EventSource behavior in API design.

---

### Pitfall 9: Frontend Assumes Synchronous Responses (Web Worker Mindset)

**What goes wrong:**

Frontend code ported from Web Worker uses `await worker.run()` pattern, expecting a single response. When migrated to REST API, developers use `fetch('/optimize')` and await the response, causing the browser to hang for 10+ minutes (timeout) while the server runs SA synchronously.

**Why it happens:**

Web Workers allow blocking calls in a separate thread—`await worker.run()` returns when complete. HTTP requests are fundamentally different: long-running requests tie up browser connections (6-per-domain limit) and timeout.

The correct pattern is fire-and-forget (POST to start) + polling/SSE (GET for updates), but this requires frontend architectural changes that developers underestimate.

**How to avoid:**

**Client-Side State Machine**:

```typescript
// Before (Web Worker pattern)
async function runOptimization(params: SAParams) {
  const worker = new SAWorker();
  worker.onProgress = (data) => updateUI(data);
  const result = await worker.run(params);  // Blocks in worker thread
  return result;
}

// After (REST API pattern)
async function runOptimization(params: SAParams) {
  // 1. Start job
  const { task_id } = await fetch('/api/optimize', {
    method: 'POST',
    body: JSON.stringify(params)
  }).then(r => r.json());

  // 2. Stream progress
  const eventSource = new EventSource(`/api/optimize/${task_id}/stream`);

  return new Promise((resolve) => {
    eventSource.onmessage = (event) => {
      const state = JSON.parse(event.data);
      updateUI(state);

      if (state.isComplete) {
        eventSource.close();
        resolve(state);
      }
    };

    eventSource.onerror = () => {
      eventSource.close();
      resolve(null);  // Handle error
    };
  });
}
```

**Shared TypeScript Types** (maintain type safety):

```typescript
// shared/types.ts (used by both frontend and backend)
export interface SAParams {
  formula: string;
  initialTemp: number;
  coolingScheduleK: number;
  stepsPerCycle: number;
  numCycles: number;
}

export interface SAProgressData {
  step: number;
  totalSteps: number;
  currentEnergy: number;
  bestEnergy: number;
  bestMolBlock: string;
  isComplete: boolean;
}

// Frontend uses directly
// Backend: FastAPI Pydantic models mirror these
```

**Adapter Pattern** (minimize frontend changes):

```typescript
// Create SAClient that mimics Web Worker API
class SAClient {
  onProgress: (data: SAProgressData) => void;

  async run(params: SAParams): Promise<SAResult> {
    const { task_id } = await this.startJob(params);
    return this.streamProgress(task_id);
  }

  private async startJob(params: SAParams) {
    return fetch('/api/optimize', {
      method: 'POST',
      body: JSON.stringify(params)
    }).then(r => r.json());
  }

  private async streamProgress(task_id: string): Promise<SAResult> {
    return new Promise((resolve, reject) => {
      const eventSource = new EventSource(`/api/optimize/${task_id}/stream`);

      eventSource.onmessage = (event) => {
        const state = JSON.parse(event.data);
        if (this.onProgress) this.onProgress(state);

        if (state.isComplete) {
          eventSource.close();
          resolve(state);
        }
      };

      eventSource.onerror = (err) => {
        eventSource.close();
        reject(err);
      };
    });
  }
}

// Frontend code barely changes
const client = new SAClient();
client.onProgress = (data) => updateChart(data);
const result = await client.run(params);
```

**Warning signs:**

- Browser "pending" requests that never complete
- Fetch requests timeout after 30-60 seconds
- Frontend uses single `await fetch()` call for entire optimization
- No SSE or polling code in frontend
- Browser dev tools show 1 active request for entire optimization duration

**Phase to address:**

**Phase 2: SSE Streaming** — Build `SAClient` adapter class to mimic Web Worker API. Update Alpine.js app to use new client. Test reconnection behavior.

---

### Pitfall 10: Plugin Architecture Prevents Pickling for ProcessPoolExecutor

**What goes wrong:**

Multi-component target function uses dependency injection and plugin pattern (classes with abstract base classes). When passing `TargetFunction` instance to `ProcessPoolExecutor`, Python fails with `PickleError: Can't pickle <class>: attribute lookup failed`.

**Why it happens:**

`ProcessPoolExecutor` serializes function arguments via `pickle` to send to worker processes. Lambda functions, nested classes, and instances with unpicklable dependencies (e.g., database connections, file handles) fail.

Plugin architectures often use dynamic class loading, decorators, or closures—all problematic for pickle.

**How to avoid:**

**Top-Level Functions + Config Objects**:

```python
# BAD: Class-based plugins can't pickle
class SpectralTargetFunction:
    def __init__(self, spectrum_db):
        self.db = spectrum_db  # Unpicklable

    def __call__(self, mol):
        return self.db.similarity(mol)

# ProcessPoolExecutor fails
executor.submit(run_sa, params, target_fn=SpectralTargetFunction(db))

# GOOD: Top-level function + serializable config
from dataclasses import dataclass

@dataclass
class TargetConfig:
    components: List[str]  # e.g., ["wiener", "spectral"]
    weights: List[float]
    spectral_data_path: str  # Not the DB connection itself

def evaluate_target(mol, config: TargetConfig):
    """Top-level function (picklable)."""
    total = 0.0

    if "wiener" in config.components:
        total += compute_wiener(mol) * config.weights[0]

    if "spectral" in config.components:
        # Load data inside worker process
        with open(config.spectral_data_path) as f:
            spectrum = json.load(f)
        total += compute_spectral_similarity(mol, spectrum) * config.weights[1]

    return total

# Usage
config = TargetConfig(
    components=["wiener", "spectral"],
    weights=[0.5, 0.5],
    spectral_data_path="/data/spectra.json"
)

executor.submit(run_sa, params, config)
```

**Plugin Registry Pattern**:

```python
# plugins.py (importable by worker processes)
TARGET_FUNCTIONS = {}

def register_target(name):
    def decorator(func):
        TARGET_FUNCTIONS[name] = func
        return func
    return decorator

@register_target("wiener")
def wiener_target(mol):
    return compute_wiener_index(mol)

@register_target("spectral")
def spectral_target(mol):
    # Load data inside function
    return compute_spectral_similarity(mol)

# Worker function
def run_sa_with_plugins(params, target_names, weights):
    functions = [TARGET_FUNCTIONS[name] for name in target_names]

    def combined_target(mol):
        return sum(f(mol) * w for f, w in zip(functions, weights))

    engine = SAEngine(params, target_fn=combined_target)
    return engine.run()

# Picklable arguments (strings, not functions)
executor.submit(run_sa_with_plugins, params, ["wiener", "spectral"], [0.5, 0.5])
```

**Alternative: ThreadPoolExecutor** (if CPU-bound isn't critical):

```python
# ThreadPoolExecutor doesn't pickle (shared memory)
from concurrent.futures import ThreadPoolExecutor

# Can use class instances directly
target_fn = SpectralTargetFunction(db)
executor = ThreadPoolExecutor(max_workers=4)
executor.submit(run_sa, params, target_fn)

# Downside: GIL limits parallelism for CPU-bound SA
```

**Warning signs:**

- `PickleError` or `AttributeError` when submitting to ProcessPoolExecutor
- Works with ThreadPoolExecutor but not ProcessPoolExecutor
- Plugin classes use `__init__` with database connections, file handles
- Errors mention `__getstate__`, `__setstate__`, or `Can't pickle`
- Works when running SA directly, fails in executor

**Phase to address:**

**Phase 3: Multi-component Target Function** — Design plugin system with picklability in mind from start. Use registry pattern or top-level functions + config objects.

---

## Technical Debt Patterns

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Skip sanitization after displacement | Faster iteration (no RDKit overhead) | Silent chemistry errors, invalid SMILES, production bugs | **Never** — sanitization is mandatory for correctness |
| Use `async def` without ProcessPoolExecutor | Simpler code (no executor boilerplate) | Server hangs under load, timeouts, unusable in production | **Only for prototyping** — replace before any demo |
| Store full `Mol` objects in history | Easier debugging (inspect molecules directly) | Memory leaks, OOM in production | **Only for short SA runs** (< 100 steps) |
| Single endpoint for POST + SSE | Less code (no job/stream separation) | Duplicate jobs on reconnect, no idempotency | **Only if network never fails** (i.e., never acceptable) |
| ThreadPoolExecutor instead of ProcessPool | Easier to debug, pickle not required | GIL limits CPU parallelism, slower SA | **Acceptable for demo** if SA runs < 10 seconds |
| CORS `allow_origins=["*"]` in dev | No hardcoded localhost ports | Security risk if deployed, accidental production exposure | **Acceptable in local dev** if env-gated |
| Synchronous blocking in FastAPI | Matches Web Worker mental model | Server becomes unresponsive | **Only for MVP with single user** — refactor before v1 |
| No disconnect detection in SSE | Simpler generator code | CPU waste on zombie tasks, capacity issues | **Only for short streams** (< 10 seconds) |

---

## Integration Gotchas

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| **RDKit ↔ Faulon** | Passing integer bond orders directly to `SetBondType()` | Use explicit mapping function (`faulon_order_to_rdkit_bondtype`) |
| **FastAPI ↔ ProcessPoolExecutor** | Passing lambda or class instances | Use top-level functions or registry pattern with string names |
| **SSE ↔ EventSource** | POST endpoint that starts job + streams | Separate POST (create job) and GET (stream progress) endpoints |
| **Vite ↔ FastAPI CORS** | Adding CORS middleware after other middleware | Add `CORSMiddleware` FIRST, before all other middleware |
| **Frontend ↔ SSE** | Assuming synchronous `fetch()` like Web Worker | Use `EventSource` + Promise wrapper to mimic Web Worker API |
| **RDKit sanitization ↔ bond changes** | Assuming valence auto-updates like TypeScript `MolGraph` | Call `Chem.SanitizeMol()` after EVERY bond modification |
| **Python GIL ↔ async** | Believing `async def` makes CPU-bound code non-blocking | Use `ProcessPoolExecutor` for CPU-bound, async for I/O-bound |

---

## Performance Traps

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| **SSE without backpressure throttling** | Server memory grows during streaming, disconnects | Adaptive reporting interval (max 10 events/sec) | Client on slow network (< 1 Mbps) or mobile |
| **SA in event loop** | All FastAPI endpoints timeout during optimization | Run SA in ProcessPoolExecutor, not `async def` | First optimization request |
| **RDKit Mol objects in history** | Memory doesn't return to baseline after SA | Store SMILES strings, not Mol objects | SA > 1000 steps or multiple jobs |
| **No SSE disconnect detection** | CPU usage stays high after client disconnect | Check `request.is_disconnected()` every 100 steps | User closes tab mid-optimization |
| **Synchronous RDKit in async endpoint** | Event loop blocked, 10+ second request times | Offload to thread pool via `run_in_executor` | Molecules > 20 atoms or complex operations |
| **EventSource auto-reconnect on POST** | Duplicate jobs after network hiccup | Separate POST (create) from GET (stream) | Network unstable (< 95% uptime) |

---

## Security Mistakes

| Mistake | Risk | Prevention |
|---------|------|------------|
| **CORS `allow_origins=["*"]` in production** | Any site can call API, data exfiltration | Use explicit origin list, env-based config |
| **No timeout on ProcessPoolExecutor tasks** | Resource exhaustion (workers stuck on infinite loops) | Set executor timeout (e.g., 10 minutes max) |
| **User-provided formula without validation** | Code injection via malicious SMILES parsing | Validate formula regex (`^[A-Z][a-z]?\d*$`) before RDKit |
| **SSE streams never close** | Memory leak, connection exhaustion | Always send terminating event, implement timeout |
| **No rate limiting on `/optimize` endpoint** | DoS via parallel job spam | Use FastAPI `SlowAPI` or Redis-based rate limiter |

---

## UX Pitfalls

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| **No feedback on SSE disconnect** | User doesn't know optimization stopped, waits indefinitely | Detect `eventSource.onerror`, show "Connection lost, retrying..." |
| **EventSource reconnects silently** | User confused by duplicate results, progress jumps | Use `Last-Event-ID` to resume, or disable auto-reconnect |
| **No progress during RDKit sanitization** | UI freezes for 2-3 seconds after each displacement | Send "sanitizing..." event before slow operations |
| **SSE errors show generic CORS message** | User can't debug actual error | Add CORS middleware first, log SSE errors server-side |
| **Long-running SA with no cancel button** | User stuck waiting, can't abort | Implement cancel endpoint, use `request.is_disconnected()` |
| **No indication if job already running** | User clicks "optimize" multiple times, spawns N jobs | Return `409 Conflict` if task_id already active |

---

## "Looks Done But Isn't" Checklist

- [ ] **RDKit displacement:** Verified `Chem.SanitizeMol()` called after bond modifications (check unit tests assert SMILES correctness)
- [ ] **FastAPI async:** Verified SA runs in ProcessPoolExecutor, not event loop (check `/health` endpoint responsive during optimization)
- [ ] **SSE backpressure:** Verified adaptive reporting interval implemented (check memory doesn't grow during streaming)
- [ ] **SSE disconnect:** Verified server cancels SA when client disconnects (check CPU drops after closing browser tab)
- [ ] **CORS middleware:** Verified `CORSMiddleware` added FIRST (check error responses include CORS headers in browser)
- [ ] **Memory cleanup:** Verified history stores SMILES not Mol objects (check memory returns to baseline after `run()`)
- [ ] **BondType mapping:** Verified explicit mapping function used (check SMILES match expected structures in tests)
- [ ] **EventSource reconnect:** Verified job creation separated from streaming (check only one job runs after reconnect)
- [ ] **Frontend adapter:** Verified `SAClient` mimics Web Worker API (check minimal frontend changes needed)
- [ ] **Plugin pickling:** Verified target functions use registry pattern (check ProcessPoolExecutor doesn't raise `PickleError`)

---

## Recovery Strategies

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| **Missing sanitization** | **MEDIUM** | Add `Chem.SanitizeMol()` after displacement, re-run tests, check SMILES correctness |
| **SA in event loop** | **MEDIUM** | Wrap SA in `run_in_executor(ProcessPoolExecutor, ...)`, update endpoint to async |
| **SSE memory leak** | **LOW** | Implement adaptive reporting interval, redeploy (no data loss) |
| **No disconnect detection** | **LOW** | Add `request.is_disconnected()` checks, redeploy |
| **CORS middleware order** | **LOW** | Move `CORSMiddleware` to first position, restart server |
| **Mol objects in history** | **MEDIUM** | Change storage to SMILES, may lose historical debugging data |
| **Wrong BondType enum** | **HIGH** | Add explicit mapping, re-validate all SA results (past results may be chemically wrong) |
| **Duplicate jobs on reconnect** | **MEDIUM** | Separate POST/GET endpoints, update frontend to use new API, requires client update |
| **Frontend synchronous fetch** | **HIGH** | Implement `SAClient` adapter, refactor all optimization UI code |
| **Unpicklable plugins** | **HIGH** | Redesign plugin system with registry pattern, may require architecture change |

---

## Pitfall-to-Phase Mapping

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| **Missing sanitization** | Phase 1: RDKit Backend Core | Unit test asserts SMILES correct after displacement |
| **SA in event loop** | Phase 1: RDKit Backend Core | `/health` endpoint responsive during optimization |
| **SSE backpressure** | Phase 2: SSE Streaming | Memory stable during 1000-step SA |
| **SSE disconnect** | Phase 2: SSE Streaming | CPU drops to baseline after killing browser tab |
| **CORS middleware order** | Phase 1: RDKit Backend Core | Browser shows 422 validation error, not CORS error |
| **RDKit memory leak** | Phase 1: RDKit Backend Core | Memory returns to baseline after SA completes |
| **BondType mapping** | Phase 1: RDKit Backend Core | SMILES for all bond types (0,1,2,3) match expected |
| **EventSource reconnect** | Phase 2: SSE Streaming | Network hiccup doesn't spawn duplicate jobs |
| **Frontend sync fetch** | Phase 2: SSE Streaming | `SAClient` adapter passes all Web Worker tests |
| **Plugin pickling** | Phase 3: Multi-component Targets | ProcessPoolExecutor doesn't raise `PickleError` |

---

## Sources

### RDKit Sanitization & Bond Manipulation
- [Molecular Operations and Sanitization | DeepWiki](https://deepwiki.com/rdkit/rdkit/2.2-atom-and-bond-representation)
- [Sanitization options and molecule parsing | RDKit blog](https://greglandrum.github.io/rdkit-blog/posts/2025-06-27-sanitization-and-file-parsing.html)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
- [Sanitization issues on identical molecules | GitHub Discussion #8156](https://github.com/rdkit/rdkit/discussions/8156)
- [Explicit Valence Error | GitHub Discussion #8181](https://github.com/rdkit/rdkit/discussions/8181)
- [FrequentlyAskedQuestions | RDKit Wiki](https://github.com/rdkit/rdkit/wiki/FrequentlyAskedQuestions)
- [rdkit.Chem.rdchem module](https://www.rdkit.org/docs/source/rdkit.Chem.rdchem.html)
- [Editing, merging, and replacing molecules in RDKit](https://asteeves.github.io/blog/2015/01/14/editing-in-rdkit/)
- [The RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html)

### FastAPI SSE & Async Patterns
- [Implementing Server-Sent Events (SSE) with FastAPI | Medium](https://mahdijafaridev.medium.com/implementing-server-sent-events-sse-with-fastapi-real-time-updates-made-simple-6492f8bfc154)
- [Using Asyncio queues in SSE implementation | Medium](https://medium.com/@Rachita_B/lookout-for-these-cryptids-while-working-with-server-sent-events-43afabb3a868)
- [FastAPI + SSE for LLM Tokens | Medium](https://medium.com/@hadiyolworld007/fastapi-sse-for-llm-tokens-smooth-streaming-without-websockets-001ead4b5e53)
- [Concurrency and async / await | FastAPI](https://fastapi.tiangolo.com/async/)
- [FastAPI endpoint with IO and CPU bound tasks | GitHub Issue #3725](https://github.com/fastapi/fastapi/issues/3725)
- [Managing Background Tasks in FastAPI | Leapcell](https://leapcell.io/blog/managing-background-tasks-and-long-running-operations-in-fastapi)

### SSE Client Disconnect & Reconnection
- [Stop streaming response when client disconnects | GitHub Discussion #7572](https://github.com/fastapi/fastapi/discussions/7572)
- [Stop Burning CPU on Dead FastAPI Streams](https://jasoncameron.dev/posts/fastapi-cancel-on-disconnect)
- [Using server-sent events | MDN](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events)
- [Server Sent Events | javascript.info](https://javascript.info/server-sent-events)
- [The Complete Guide to Server-Sent Events (SSE)](https://singhajit.com/server-sent-events-explained/)

### CORS & FastAPI
- [CORS (Cross-Origin Resource Sharing) | FastAPI](https://fastapi.tiangolo.com/tutorial/cors/)
- [How to Fix CORS in FastAPI | Medium](https://cleverzone.medium.com/how-to-fix-cors-in-fastapi-a9b1f597661b)
- [Blocked by CORS in FastAPI? Here's How to Fix It](https://davidmuraya.com/blog/fastapi-cors-configuration/)
- [FastAPI CORS Error | Sentry](https://sentry.io/answers/fastapi-error-no-access-control-allow-origin-header-is-present-on-the-requested-resource/)

### RDKit Memory Management
- [Memory leakage when storing molecule objects | GitHub Issue #3239](https://github.com/rdkit/rdkit/issues/3239)
- [Getting Started with the RDKit in Python](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [Process aborts during gc | GitHub Issue #579](https://github.com/rdkit/rdkit/issues/579)

### Python Plugin Architecture
- [How to Design and Implement a Plugin Architecture in Python](https://mathieularose.com/plugin-architecture-in-python)
- [Python Dependency Injection: Best Practices | ArjanCodes](https://arjancodes.com/blog/python-dependency-injection-best-practices/)
- [Developing Plugin Architecture with Pluggy | Medium](https://medium.com/@garzia.luke/developing-plugin-architecture-with-pluggy-8eb7bdba3303)
- [Dependency injection in Python | Snyk](https://snyk.io/blog/dependency-injection-python/)

### Monorepo Structure
- [Python Monorepo: Structure and Tooling | Tweag](https://www.tweag.io/blog/2023-04-04-python-monorepo-1/)
- [Best practices for managing frontend and backend in monorepo](https://graphite.com/guides/monorepo-frontend-backend-best-practices)
- [Building a Monorepo with Python | Earthly Blog](https://earthly.dev/blog/python-monorepo/)

### Frontend Migration Patterns
- [Modernizing a React Application | Medium](https://medium.com/@sriram_in/modernizing-a-react-application-a-phased-approach-to-backend-migration-and-frontend-refactoring-bf170caf79ef)
- [Navigating Frontend Migration | Medium](https://medium.com/syngenta-digitalblog/navigating-frontend-migration-strategies-for-refactoring-rewriting-and-embracing-microfrontends-331520cde2bb)

---

*Pitfalls research for: WebFaulon Browser-to-Python Backend Migration*
*Researched: 2026-02-15*
*Confidence: HIGH — findings based on official RDKit/FastAPI documentation, recent community discussions (2024-2026), and established patterns from production deployments*
