# Architecture Research

**Domain:** Browser-based Cheminformatics with Computationally Intensive Algorithms
**Researched:** 2026-02-14
**Confidence:** HIGH

## Standard Architecture

### System Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                          UI LAYER (Main Thread)                              │
├─────────────────────────────────────────────────────────────────────────────┤
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐    │
│  │   UI         │  │   Chart      │  │   Controls   │  │   Status     │    │
│  │  Container   │  │  Renderer    │  │              │  │   Display    │    │
│  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘  └──────┬───────┘    │
│         │                 │                 │                 │             │
│         └─────────────────┴─────────────────┴─────────────────┘             │
│                                    │                                         │
│                            ┌───────▼────────┐                                │
│                            │  State Manager │                                │
│                            │  (Application  │                                │
│                            │     State)     │                                │
│                            └───────┬────────┘                                │
│                                    │                                         │
├────────────────────────────────────┼─────────────────────────────────────────┤
│                       COMMUNICATION LAYER                                    │
│                                    │                                         │
│                      ┌─────────────▼──────────────┐                          │
│                      │  Message Queue/Dispatcher  │                          │
│                      │    (postMessage + rAF)     │                          │
│                      └─────────────┬──────────────┘                          │
│                                    │                                         │
├────────────────────────────────────┼─────────────────────────────────────────┤
│                      COMPUTATION LAYER (Web Worker)                          │
│                                    │                                         │
│  ┌─────────────────────────────────▼──────────────────────────────────┐     │
│  │                          Worker Thread                              │     │
│  │                                                                      │     │
│  │  ┌────────────────┐    ┌────────────────┐    ┌────────────────┐   │     │
│  │  │  SA Algorithm  │───▶│  RDKit.js      │───▶│  Graph Engine  │   │     │
│  │  │  Controller    │    │  (WASM)        │    │  (MolGraph)    │   │     │
│  │  └────────────────┘    └────────────────┘    └────────────────┘   │     │
│  │         │                      │                      │            │     │
│  │         │                      │                      │            │     │
│  │         ▼                      ▼                      ▼            │     │
│  │  ┌──────────────────────────────────────────────────────────┐     │     │
│  │  │          Worker State (Iteration Data)                   │     │     │
│  │  └──────────────────────────────────────────────────────────┘     │     │
│  │                                                                      │     │
│  └──────────────────────────────────────────────────────────────────────┘     │
│                                                                               │
└───────────────────────────────────────────────────────────────────────────────┘
```

### Component Responsibilities

| Component | Responsibility | Typical Implementation |
|-----------|----------------|------------------------|
| **UI Container** | DOM manipulation, user interactions, display updates | React/Vue component or vanilla JS with DOM refs |
| **Chart Renderer** | Real-time visualization using canvas/WebGL | Chart.js with decimation + rAF throttling |
| **Controls** | Start/stop/pause controls, parameter inputs | Event-driven UI components |
| **Status Display** | Current iteration, temperature, energy metrics | Real-time DOM updates from batched worker messages |
| **State Manager** | Application state (algorithm status, parameters, history) | Zustand/Jotai for client state or useState/useReducer |
| **Message Dispatcher** | Worker communication, message batching, rAF scheduling | Custom queue with requestAnimationFrame throttling |
| **SA Algorithm Controller** | Simulated annealing loop, temperature scheduling, acceptance criteria | Pure JS/TS implementation in worker |
| **RDKit.js (WASM)** | Molecular calculations (if needed for validation) | Async-initialized WASM module |
| **Graph Engine (MolGraph)** | Graph representation, bond matrix, valence checking, Wiener Index, BFS/DFS | TypeScript class with adjacency matrix operations |
| **Worker State** | Current molecule, energy history, iteration counters | In-memory data structures in worker scope |

## Recommended Project Structure

```
src/
├── components/           # UI components (if using React/Vue)
│   ├── ChartView/        # Real-time chart visualization
│   ├── Controls/         # Start/stop/parameter controls
│   └── MoleculeDisplay/  # Molecular structure display
├── core/                 # Shared logic (both main + worker)
│   ├── types.ts          # TypeScript interfaces/types
│   └── constants.ts      # Shared constants
├── worker/               # Web Worker code
│   ├── worker.ts         # Worker entry point
│   ├── SAEngine.ts       # Simulated annealing implementation
│   ├── MolGraph.ts       # Molecular graph data structure
│   └── rdkit-loader.ts   # RDKit.js initialization
├── services/             # Main thread services
│   ├── WorkerManager.ts  # Worker communication abstraction
│   └── MessageBatcher.ts # Message queue + batching
├── state/                # State management
│   └── store.ts          # Application state (Zustand/Jotai)
├── utils/                # Utilities
│   ├── rAFThrottle.ts    # requestAnimationFrame throttling
│   └── chartConfig.ts    # Chart.js configuration
└── main.ts               # Application entry point
```

### Structure Rationale

- **components/:** UI layer isolated from computation logic, easily testable and swappable
- **core/:** Shared types prevent duplication between main thread and worker
- **worker/:** All computation isolated in Web Worker to prevent UI blocking
- **services/:** Communication layer abstracts worker complexity from UI
- **state/:** Centralized state management for UI consistency
- **utils/:** Reusable helpers for performance optimization (throttling, batching)

## Architectural Patterns

### Pattern 1: Web Worker for Computation Isolation

**What:** Offload the SA algorithm to a dedicated Web Worker, keeping the main thread free for UI updates and user interactions.

**When to use:** Always for computationally intensive algorithms that would block the UI thread (>16ms per iteration at 60fps).

**Trade-offs:**
- **Pros:** Non-blocking UI, better responsiveness, prevents INP degradation
- **Cons:** Serialization overhead for data transfer, debugging complexity, no DOM access in worker

**Example:**
```typescript
// Main thread: WorkerManager.ts
export class WorkerManager {
  private worker: Worker;
  private messageQueue: Message[] = [];
  private batchSize = 10;

  constructor() {
    this.worker = new Worker(new URL('../worker/worker.ts', import.meta.url));
    this.worker.onmessage = this.handleMessage.bind(this);
    this.worker.onerror = this.handleError.bind(this);
  }

  start(params: SAParams) {
    this.worker.postMessage({ type: 'START', params });
  }

  stop() {
    this.worker.postMessage({ type: 'STOP' });
  }

  private handleMessage(e: MessageEvent) {
    const { type, data } = e.data;

    if (type === 'ITERATION_UPDATE') {
      this.messageQueue.push(data);
      if (this.messageQueue.length >= this.batchSize) {
        this.flushQueue();
      }
    }
  }

  private flushQueue() {
    // Batch processing for UI updates
    const batch = this.messageQueue.splice(0, this.batchSize);
    requestAnimationFrame(() => {
      // Update UI with batched data
      this.updateChart(batch);
    });
  }
}

// Worker thread: worker.ts
let saEngine: SAEngine | null = null;
let isRunning = false;

self.onmessage = (e: MessageEvent) => {
  const { type, params } = e.data;

  switch (type) {
    case 'START':
      saEngine = new SAEngine(params);
      isRunning = true;
      runLoop();
      break;
    case 'STOP':
      isRunning = false;
      break;
  }
};

function runLoop() {
  if (!isRunning || !saEngine) return;

  // Run one SA iteration
  const result = saEngine.iterate();

  // Send update to main thread
  self.postMessage({
    type: 'ITERATION_UPDATE',
    data: {
      iteration: result.iteration,
      energy: result.energy,
      temperature: result.temperature
    }
  });

  // Continue loop
  setTimeout(runLoop, 0); // Or use setInterval pattern
}
```

### Pattern 2: Message Batching + requestAnimationFrame Throttling

**What:** Batch multiple worker messages together and throttle UI updates to the browser's refresh rate (60fps) using requestAnimationFrame.

**When to use:** When worker sends many messages per second (hundreds to thousands) but UI only needs ~60 updates/sec.

**Trade-offs:**
- **Pros:** Prevents UI thrashing, reduces render overhead, smoother animations, better Core Web Vitals (INP)
- **Cons:** Slight latency between computation and display (negligible for visualization), increased complexity

**Example:**
```typescript
// services/MessageBatcher.ts
export class MessageBatcher {
  private queue: any[] = [];
  private rafId: number | null = null;
  private callback: (batch: any[]) => void;

  constructor(callback: (batch: any[]) => void) {
    this.callback = callback;
  }

  push(message: any) {
    this.queue.push(message);

    // Schedule flush if not already scheduled
    if (this.rafId === null) {
      this.rafId = requestAnimationFrame(() => this.flush());
    }
  }

  private flush() {
    if (this.queue.length > 0) {
      const batch = this.queue.splice(0);
      this.callback(batch);
    }
    this.rafId = null;
  }
}

// Usage in WorkerManager
this.batcher = new MessageBatcher((batch) => {
  // Update chart with last N points or decimated data
  const decimatedData = this.decimateIfNeeded(batch);
  this.updateChart(decimatedData);
});

worker.onmessage = (e) => {
  this.batcher.push(e.data);
};
```

### Pattern 3: Transferable Objects for Large Data

**What:** Use transferable objects (ArrayBuffer, TypedArray) for zero-copy data transfer between main thread and worker.

**When to use:** When transferring large datasets (adjacency matrices, molecular coordinates, large arrays) between threads.

**Trade-offs:**
- **Pros:** No serialization/deserialization overhead, instant transfer, lower memory usage
- **Cons:** Original buffer becomes neutered (unusable), requires careful ownership management

**Example:**
```typescript
// Main thread: Transferring adjacency matrix to worker
const adjacencyMatrix = new Uint8Array(1000 * 1000); // Large matrix
// ... populate matrix ...

// Transfer ownership to worker
worker.postMessage(
  { type: 'INIT_GRAPH', matrix: adjacencyMatrix.buffer },
  [adjacencyMatrix.buffer] // Transferable list
);
// adjacencyMatrix.buffer is now neutered/unusable on main thread

// Worker thread: Receiving transferred data
self.onmessage = (e) => {
  if (e.data.type === 'INIT_GRAPH') {
    const matrix = new Uint8Array(e.data.matrix);
    // Now worker owns this memory
  }
};
```

### Pattern 4: Async WASM Initialization Pattern

**What:** Initialize RDKit.js WASM module asynchronously in worker before starting computation.

**When to use:** Always when using WASM libraries that require initialization (RDKit.js, TensorFlow.js, etc.).

**Trade-offs:**
- **Pros:** Proper initialization, error handling, prevents race conditions
- **Cons:** Adds initialization latency, requires promise handling

**Example:**
```typescript
// worker/rdkit-loader.ts
import initRDKitModule from '@rdkit/rdkit';

let RDKit: any = null;

export async function initRDKit(): Promise<void> {
  if (RDKit) return; // Already initialized

  try {
    RDKit = await initRDKitModule({
      locateFile: (filename: string) => {
        return `/wasm/${filename}`;
      }
    });
    console.log(`RDKit version: ${RDKit.version()}`);
  } catch (error) {
    console.error('Failed to initialize RDKit:', error);
    throw error;
  }
}

export function getRDKit() {
  if (!RDKit) {
    throw new Error('RDKit not initialized. Call initRDKit() first.');
  }
  return RDKit;
}

// worker/worker.ts
import { initRDKit } from './rdkit-loader';

let initialized = false;

self.onmessage = async (e) => {
  if (!initialized) {
    self.postMessage({ type: 'STATUS', message: 'Initializing RDKit...' });
    await initRDKit();
    initialized = true;
    self.postMessage({ type: 'READY' });
  }

  // Handle messages after initialization
  // ...
};
```

### Pattern 5: Chart.js Decimation + Performance Optimization

**What:** Configure Chart.js for optimal performance with large, streaming datasets using decimation, disabled animations, and pre-prepared data.

**When to use:** When displaying thousands+ data points with real-time updates.

**Trade-offs:**
- **Pros:** Smooth 60fps updates even with large datasets, lower memory usage
- **Cons:** Some visual detail lost with decimation, animations disabled

**Example:**
```typescript
// utils/chartConfig.ts
export const chartConfig: ChartConfiguration = {
  type: 'line',
  data: {
    datasets: [{
      label: 'Energy',
      data: [],
      parsing: false, // Pre-formatted data
      normalized: true // Data already sorted
    }]
  },
  options: {
    animation: false, // Critical for performance
    responsive: true,
    maintainAspectRatio: false,
    elements: {
      line: {
        tension: 0 // Disable Bézier curves (faster)
      },
      point: {
        radius: 0 // Don't render points (faster)
      }
    },
    scales: {
      x: {
        type: 'linear',
        min: 0,
        max: 1000, // Fixed range avoids recalculation
        ticks: {
          sampleSize: 10 // Sample ticks for speed
        }
      },
      y: {
        min: 0,
        max: 100, // Fixed range avoids recalculation
        ticks: {
          sampleSize: 10
        }
      }
    },
    plugins: {
      decimation: {
        enabled: true,
        algorithm: 'lttb', // Largest Triangle Three Buckets
        samples: 500 // Max points to display
      }
    }
  }
};

// components/ChartView/index.ts
export class ChartView {
  private chart: Chart;
  private dataBuffer: Point[] = [];
  private maxPoints = 10000;

  updateData(newPoints: Point[]) {
    // Add to buffer
    this.dataBuffer.push(...newPoints);

    // Trim if exceeds max
    if (this.dataBuffer.length > this.maxPoints) {
      this.dataBuffer = this.dataBuffer.slice(-this.maxPoints);
    }

    // Update chart (Chart.js decimation handles rendering)
    this.chart.data.datasets[0].data = this.dataBuffer;
    this.chart.update('none'); // Skip animations
  }
}
```

### Pattern 6: State Machine for Algorithm Control

**What:** Use explicit state machine pattern for algorithm lifecycle (IDLE → INITIALIZING → RUNNING → PAUSED → STOPPED → ERROR).

**When to use:** Always for complex asynchronous operations with multiple states and transitions.

**Trade-offs:**
- **Pros:** Clear state transitions, easier debugging, prevents invalid operations
- **Cons:** More boilerplate code, state management overhead

**Example:**
```typescript
// worker/SAEngine.ts
enum SAState {
  IDLE = 'IDLE',
  INITIALIZING = 'INITIALIZING',
  RUNNING = 'RUNNING',
  PAUSED = 'PAUSED',
  STOPPED = 'STOPPED',
  ERROR = 'ERROR'
}

export class SAEngine {
  private state: SAState = SAState.IDLE;

  async initialize(params: SAParams) {
    if (this.state !== SAState.IDLE) {
      throw new Error(`Cannot initialize from state ${this.state}`);
    }

    this.state = SAState.INITIALIZING;
    // ... initialization logic ...
    this.state = SAState.RUNNING;
  }

  pause() {
    if (this.state !== SAState.RUNNING) {
      throw new Error(`Cannot pause from state ${this.state}`);
    }
    this.state = SAState.PAUSED;
  }

  resume() {
    if (this.state !== SAState.PAUSED) {
      throw new Error(`Cannot resume from state ${this.state}`);
    }
    this.state = SAState.RUNNING;
  }

  stop() {
    if (![SAState.RUNNING, SAState.PAUSED].includes(this.state)) {
      throw new Error(`Cannot stop from state ${this.state}`);
    }
    this.state = SAState.STOPPED;
  }

  iterate() {
    if (this.state !== SAState.RUNNING) {
      return null;
    }
    // ... SA iteration logic ...
  }
}
```

## Data Flow

### Request Flow

```
[User clicks START]
    │
    ▼
[Controls Component] → [State Manager] → [WorkerManager.start()]
    │                       │                   │
    │                       ▼                   ▼
    │               [Update UI state]    [postMessage to Worker]
    │                                            │
    │                                            ▼
    │                                  [Worker: Initialize SA]
    │                                            │
    │                                            ▼
    │                                  [Worker: Run SA loop]
    │                                            │
    │                          ┌─────────────────┴─────────────────┐
    │                          ▼                                   ▼
    │              [Worker: Calculate energy]        [Worker: Update graph]
    │                          │                                   │
    │                          └─────────────────┬─────────────────┘
    │                                            ▼
    │                              [Worker: postMessage with results]
    │                                            │
    │                                            ▼
    └──────────────────────┐         [MessageBatcher.push()]
                           │                     │
                           │                     ▼
                           │         [requestAnimationFrame]
                           │                     │
                           │                     ▼
                           │           [MessageBatcher.flush()]
                           │                     │
                           ▼                     ▼
                    [State Manager] ←────[Batched data]
                           │
                           ▼
            ┌──────────────┴──────────────┐
            ▼                              ▼
    [Chart Update]                [Status Display Update]
```

### State Management

```
[Application State Store]
    ├─ algorithmState: { status, iteration, temperature, energy }
    ├─ parameters: { initialTemp, coolingRate, maxIterations }
    ├─ history: { energyData[], acceptedMoves[], rejectedMoves[] }
    └─ ui: { isPaused, isRunning, error }
            │
            ▼
    [React Components / UI]
    (subscribe to state changes)
            │
            ▼
    [Render UI based on state]
```

### Key Data Flows

1. **Initialization Flow:**
   - User configures parameters → State Manager → Worker receives config → RDKit.js initializes (async) → Graph structure created → Ready signal to UI

2. **Iteration Flow (Hot Path):**
   - Worker runs SA iteration → Calculates energy → Posts message → MessageBatcher queues → rAF triggers → Batch processed → Chart updates

3. **Control Flow:**
   - User clicks pause/resume/stop → State Manager updates → Worker receives control message → Changes internal state → Acknowledges to UI

4. **Error Flow:**
   - Worker encounters error → Posts error message → State Manager sets error state → UI displays error → Algorithm stops

## Scaling Considerations

| Concern | Initial (Demo) | Medium (10K iterations) | Large (1M+ iterations) |
|---------|----------------|-------------------------|------------------------|
| **Data Storage** | In-memory arrays | In-memory with decimation | IndexedDB for history, keep window in memory |
| **Chart Rendering** | Chart.js with decimation | Same + aggressive decimation (500 points) | Consider WebGL (Plotly.js, ECharts GL) |
| **Message Frequency** | Every iteration | Batch every 10-100 iterations | Batch every 100-1000 iterations |
| **Worker Strategy** | Single dedicated worker | Same | Consider splitting (one for SA, one for validation) |
| **Memory Management** | Basic arrays | Circular buffers for history | Typed arrays + manual GC triggers |

### Scaling Priorities

1. **First bottleneck: Chart rendering with large datasets**
   - **What breaks:** Canvas redraws become expensive beyond ~10K points
   - **Fix:** Enable Chart.js decimation plugin with LTTB algorithm, limit visible points to 500-1000
   - **Alternative:** Switch to WebGL-based charting (Plotly.js, ECharts GL, LightningChart JS)

2. **Second bottleneck: Message passing overhead**
   - **What breaks:** Posting message every iteration (1000s/sec) causes serialization overhead
   - **Fix:** Batch messages (send every N iterations or every X ms), use transferable objects for large data
   - **Alternative:** Use SharedArrayBuffer for shared state (requires CORS headers)

3. **Third bottleneck: Memory from history accumulation**
   - **What breaks:** Storing millions of iteration records exhausts browser memory
   - **Fix:** Use circular buffer, store only last N iterations + sampled history
   - **Alternative:** Move history to IndexedDB, keep sliding window in memory

## Anti-Patterns

### Anti-Pattern 1: Running SA in Main Thread

**What people do:** Implement SA algorithm directly in React component or main thread code.

**Why it's wrong:**
- Blocks UI thread during each iteration (SA can take seconds/minutes)
- Causes janky scrolling, unresponsive controls
- Terrible Core Web Vitals (INP > 500ms)
- Browser may show "unresponsive script" warning

**Do this instead:**
Always use Web Worker for any computation taking >16ms. Even if initial performance seems OK, it will degrade with larger molecules or more iterations.

```typescript
// BAD: Main thread blocks
function runSA() {
  for (let i = 0; i < 10000; i++) {
    const newState = generateNeighbor(currentState);
    const energy = calculateEnergy(newState); // Expensive!
    // UI is frozen during this entire loop
  }
}

// GOOD: Worker handles computation
// worker.ts
self.onmessage = (e) => {
  if (e.data.type === 'START') {
    for (let i = 0; i < 10000; i++) {
      const newState = generateNeighbor(currentState);
      const energy = calculateEnergy(newState);

      if (i % 100 === 0) {
        self.postMessage({ iteration: i, energy }); // Periodic updates
      }
    }
  }
};
```

### Anti-Pattern 2: Sending Worker Messages Every Iteration Without Batching

**What people do:** Send postMessage for every SA iteration, causing thousands of messages per second.

**Why it's wrong:**
- Each postMessage has serialization overhead
- Main thread overwhelmed processing messages
- Chart updates 1000x/sec when display only needs 60fps
- Wastes CPU on redundant work

**Do this instead:**
Batch messages and throttle UI updates with requestAnimationFrame.

```typescript
// BAD: Flooding main thread
// worker.ts
for (let i = 0; i < 100000; i++) {
  self.postMessage({ iteration: i, energy }); // 100K messages!
}

// GOOD: Batched messages
// worker.ts
const batch = [];
for (let i = 0; i < 100000; i++) {
  batch.push({ iteration: i, energy });

  if (batch.length >= 100) {
    self.postMessage({ type: 'BATCH', data: batch });
    batch.length = 0;
  }
}
```

### Anti-Pattern 3: Initializing RDKit.js on Every Message

**What people do:** Call initRDKitModule() every time worker receives a message.

**Why it's wrong:**
- WASM initialization takes 100s of ms
- Creates multiple instances (memory leak)
- Async initialization races with computation

**Do this instead:**
Initialize once on worker startup, wait for ready signal before accepting work.

```typescript
// BAD: Reinitializing
self.onmessage = async (e) => {
  const RDKit = await initRDKitModule(); // Slow!
  // ... use RDKit ...
};

// GOOD: One-time initialization
let RDKit = null;
let isReady = false;

(async () => {
  RDKit = await initRDKitModule();
  isReady = true;
  self.postMessage({ type: 'READY' });
})();

self.onmessage = (e) => {
  if (!isReady) {
    self.postMessage({ type: 'ERROR', message: 'Not ready' });
    return;
  }
  // ... use RDKit ...
};
```

### Anti-Pattern 4: Not Using Chart.js Decimation for Large Datasets

**What people do:** Feed all 100K+ data points directly to Chart.js without decimation.

**Why it's wrong:**
- Canvas draws all points even if many pixels overlap
- Causes stuttering/lag on chart updates
- High memory usage from large data arrays
- Poor UX with slow panning/zooming

**Do this instead:**
Enable decimation plugin and limit visible points.

```typescript
// BAD: No decimation
const chartConfig = {
  type: 'line',
  data: {
    datasets: [{
      data: allHundredThousandPoints // Slow!
    }]
  }
};

// GOOD: With decimation
const chartConfig = {
  type: 'line',
  data: {
    datasets: [{
      data: dataPoints,
      parsing: false,
      normalized: true
    }]
  },
  options: {
    animation: false, // Critical!
    plugins: {
      decimation: {
        enabled: true,
        algorithm: 'lttb',
        samples: 500 // Only render 500 points
      }
    }
  }
};
```

### Anti-Pattern 5: Blocking Worker with Synchronous Operations

**What people do:** Use synchronous loops without yielding control in worker.

**Why it's wrong:**
- Worker becomes unresponsive to stop/pause messages
- Cannot interrupt long-running computations
- Poor UX (user clicks stop, nothing happens)

**Do this instead:**
Use async/await or setTimeout to yield control periodically.

```typescript
// BAD: Blocking loop
self.onmessage = (e) => {
  if (e.data.type === 'START') {
    for (let i = 0; i < 1000000; i++) {
      // ... expensive work ...
      // Cannot receive 'STOP' message during this loop!
    }
  }
};

// GOOD: Interruptible loop
let isRunning = false;

self.onmessage = (e) => {
  if (e.data.type === 'START') {
    isRunning = true;
    runLoop(0);
  } else if (e.data.type === 'STOP') {
    isRunning = false;
  }
};

async function runLoop(iteration) {
  if (!isRunning || iteration >= 1000000) return;

  // ... expensive work ...

  if (iteration % 1000 === 0) {
    // Yield control every 1000 iterations
    setTimeout(() => runLoop(iteration + 1), 0);
  } else {
    runLoop(iteration + 1);
  }
}
```

### Anti-Pattern 6: Copying Large Adjacency Matrices Between Threads

**What people do:** Send large typed arrays via postMessage without using transferables.

**Why it's wrong:**
- Data is copied (structured clone), wasting memory
- Slow serialization/deserialization for large arrays
- 2x memory usage (copy in each thread)

**Do this instead:**
Use transferable objects for zero-copy transfer.

```typescript
// BAD: Copying data
const matrix = new Float32Array(10000); // 40KB
worker.postMessage({ matrix }); // Copied, slow

// GOOD: Transferring ownership
const matrix = new Float32Array(10000);
worker.postMessage({ matrix }, [matrix.buffer]); // Zero-copy, fast
// matrix is now neutered on main thread
```

## Integration Points

### External Libraries

| Library | Integration Pattern | Notes |
|---------|---------------------|-------|
| **RDKit.js** | Async initialization in worker, singleton pattern | Initialize once, cache instance, use for validation/rendering |
| **Chart.js** | Canvas-based rendering on main thread | Configure with decimation, disable animations, use rAF throttling |
| **React/Vue** (optional) | Component-based UI, state management hooks | Use state manager (Zustand/Jotai) for worker communication |
| **Vite/Webpack** | Bundle with Web Worker support | Use `new Worker(new URL('./worker.ts', import.meta.url))` pattern |
| **TypeScript** | Shared types between main + worker | Define interfaces in `core/types.ts` |

### Internal Boundaries

| Boundary | Communication | Notes |
|----------|---------------|-------|
| **UI ↔ State Manager** | Direct function calls / hooks | Synchronous, React re-renders on state change |
| **State Manager ↔ WorkerManager** | Direct function calls | Synchronous method calls (start, stop, pause) |
| **WorkerManager ↔ Worker** | postMessage / onmessage | Asynchronous, structured cloning or transferables |
| **Worker ↔ RDKit.js** | Direct function calls (same thread) | Synchronous within worker |
| **State Manager ↔ Chart** | Direct function calls | Synchronous, update chart data array |
| **MessageBatcher ↔ Chart** | requestAnimationFrame callbacks | Throttled to 60fps |

## Build Order Implications

Based on dependencies and complexity, suggested implementation order:

### Phase 1: Core Infrastructure
1. **Molecular Graph Engine** - Build `MolGraph.ts` with adjacency matrix, valence checking, connectivity, Wiener Index
   - *Why first:* Foundation for all SA operations, can be tested independently
   - *Testing:* Unit tests for graph operations, no UI needed

2. **Web Worker Setup** - Create worker entry point, message handling skeleton
   - *Why second:* Test worker communication before complex logic
   - *Testing:* Echo messages back and forth, verify lifecycle

3. **Basic UI Shell** - Minimal HTML/CSS with start/stop buttons, status display
   - *Why third:* Visual feedback for testing worker integration
   - *Testing:* Manual testing of controls

### Phase 2: Algorithm Implementation
4. **SA Algorithm Core** - Implement SA loop, temperature schedule, acceptance criteria
   - *Why fourth:* Can test in worker with simple logging
   - *Testing:* Verify algorithm converges on known test cases

5. **Displacement Function** - Implement Faulon's 4-atom bond order modification
   - *Why fifth:* Core SA operation, depends on MolGraph
   - *Testing:* Unit tests for valence preservation, reversibility

6. **Energy Calculation** - Implement Wiener Index or other graph metric
   - *Why sixth:* Completes SA loop
   - *Testing:* Compare with reference values

### Phase 3: Communication & Optimization
7. **Message Batching** - Implement `MessageBatcher` with rAF throttling
   - *Why seventh:* Optimize worker→main communication
   - *Testing:* Verify batching reduces message frequency

8. **WorkerManager** - Abstract worker communication, lifecycle management
   - *Why eighth:* Clean API for UI layer
   - *Testing:* Test start/stop/pause/resume transitions

9. **State Management** - Set up Zustand/Jotai store for application state
   - *Why ninth:* Centralize state for UI consistency
   - *Testing:* Verify state updates propagate to UI

### Phase 4: Visualization
10. **Chart Integration** - Set up Chart.js with optimization (decimation, no animations)
    - *Why tenth:* Data is flowing, ready to visualize
    - *Testing:* Performance testing with 10K+ points

11. **Real-time Updates** - Connect batched worker messages to chart updates
    - *Why eleventh:* Complete the visualization loop
    - *Testing:* Verify smooth 60fps updates

### Phase 5: Polish & Advanced Features
12. **RDKit.js Integration** - Add async initialization, molecular rendering
    - *Why twelfth:* Optional enhancement, not critical path
    - *Testing:* Test WASM initialization, rendering

13. **UI Enhancement** - Add parameter controls, better styling, molecule display
    - *Why last:* Polish after core functionality works
    - *Testing:* User testing, visual QA

**Critical Dependencies:**
- MolGraph MUST be built before SA Algorithm
- Worker MUST be set up before SA Algorithm runs
- MessageBatcher SHOULD be implemented before chart integration
- State Management SHOULD be in place before complex UI

**Parallel Work Opportunities:**
- MolGraph can be developed in parallel with Worker setup
- UI shell can be developed in parallel with SA Algorithm
- Chart configuration can be prototyped in parallel with Message Batching

## Sources

**Web Workers Architecture:**
- [Using Web Workers - MDN](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Using_web_workers)
- [8 Web-Worker Patterns That Make Browser Apps Multicore](https://medium.com/@connect.hashblock/8-web-worker-patterns-that-make-browser-apps-multicore-e3f22e9f6f82)
- [Web Workers: Parallel Processing in the Browser](https://medium.com/@artemkhrenov/web-workers-parallel-processing-in-the-browser-e4c89e6cad77)

**WASM & Web Workers:**
- [Using WebAssembly with Web Workers - SitePen](https://www.sitepen.com/blog/using-webassembly-with-web-workers)
- [3W for In-Browser AI: WebLLM + WASM + WebWorkers](https://blog.mozilla.ai/3w-for-in-browser-ai-webllm-wasm-webworkers/)
- [Worker Communication – WebR](https://docs.r-wasm.org/webr/latest/communication.html)

**Real-time Visualization:**
- [Real-Time Dashboard Performance: WebGL vs Canvas Rendering](https://dev3lop.com/real-time-dashboard-performance-webgl-vs-canvas-rendering-benchmarks/)
- [SVG vs. Canvas vs. WebGL for Data Visualization](https://dev3lop.com/svg-vs-canvas-vs-webgl-rendering-choice-for-data-visualization/)
- [Chart.js Performance Documentation](https://www.chartjs.org/docs/latest/general/performance.html)

**Cheminformatics:**
- [RDKit.js GitHub Repository](https://github.com/rdkit/rdkit-js)
- [3Dmol.js: Molecular Visualization with WebGL](https://academic.oup.com/bioinformatics/article/31/8/1322/213186)

**Performance Patterns:**
- [When browsers throttle requestAnimationFrame](https://motion.dev/magazine/when-browsers-throttle-requestanimationframe)
- [Message Batching Patterns](https://medium.com/@pilovm/achieve-more-reactivity-with-web-workers-and-queues-dac461ec5f8e)
- [State Management in 2026](https://medium.com/@orami98/modern-state-management-in-vanilla-javascript-2026-patterns-and-beyond-ce00425f7ac5)

---
*Architecture research for: WebFaulon - Browser-based Cheminformatics SA Visualization*
*Researched: 2026-02-14*
