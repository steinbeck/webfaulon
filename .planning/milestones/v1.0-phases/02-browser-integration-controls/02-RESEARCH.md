# Phase 2: Browser Integration & Controls - Research

**Researched:** 2026-02-15
**Domain:** Web Workers, Browser UI Integration, WASM Cheminformatics
**Confidence:** HIGH

## Summary

Phase 2 integrates the Phase 1 SA engine into a browser environment with responsive UI controls. The core technical challenge is running computationally expensive SA iterations (10,000+ steps) without freezing the UI, achieved through Web Workers. The tech stack is already decided: Vite 6 + TypeScript 5 + Alpine.js 3 + Comlink 4 + RDKit.js WASM.

**Key findings:**
- Comlink abstracts Web Worker message passing into RPC-style async calls, making worker integration straightforward
- Vite 6 provides first-class Web Worker support with `?worker` import syntax and automatic TypeScript transpilation
- Alpine.js offers lightweight reactivity perfect for simple control panels (play/pause/reset, parameter inputs)
- RDKit.js WASM can load in workers but requires careful initialization pattern (locateFile option)
- Progress reporting from long-running SA requires callback pattern or periodic postMessage updates

**Primary recommendation:** Use Comlink to wrap SAEngine class in Web Worker, Alpine.js for reactive UI controls, and implement step-by-step progress callbacks for real-time visualization updates.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| [Comlink](https://github.com/GoogleChromeLabs/comlink) | 4.4.2 | Web Worker RPC abstraction | Official Google Chrome Labs library, ~1.1kB, industry standard for worker communication, removes postMessage complexity |
| [Alpine.js](https://alpinejs.dev/directives/model) | 3.x | Reactive UI framework | Lightweight (<15kB), no build step required, perfect for simple interactive UIs without React/Vue overhead |
| [Vite](https://vite.dev/guide/features) | 6.x | Build tool & dev server | First-class Web Worker support, native ESM in dev, automatic TypeScript transpilation for workers |
| [RDKit.js](https://github.com/rdkit/rdkit-js) | Latest | Cheminformatics WASM | Official JavaScript distribution of RDKit C++ library, full molecule rendering in browser |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| TypeScript | 5.6+ | Type safety | Worker message contracts, SAEngine API types, Alpine component interfaces |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| Comlink | Raw postMessage | More control but 10x boilerplate, manual serialization, no type safety |
| Alpine.js | React/Vue | Much heavier bundle, requires complex build setup, overkill for simple controls |
| Vite Workers | Webpack worker-loader | Older approach, less standard, more configuration needed |

**Installation:**
```bash
npm install comlink alpinejs @rdkit/rdkit
# Dev dependencies already present: vite, typescript
```

## Architecture Patterns

### Recommended Project Structure
```
src/
├── core/              # Phase 1 SA engine (already complete)
├── worker/            # Web Worker entry points
│   └── sa-worker.ts   # Comlink-exposed SAEngine wrapper
├── ui/                # Alpine.js components
│   ├── controls.ts    # Parameter inputs, play/pause/reset
│   └── validation.ts  # Formula validation logic
└── main.ts            # App initialization, worker setup
```

### Pattern 1: Comlink Worker Wrapper
**What:** Expose SAEngine class through Comlink for async RPC-style calls
**When to use:** Any class-based logic that needs to run off main thread
**Example:**
```typescript
// Source: https://github.com/GoogleChromeLabs/comlink (official README)
// worker/sa-worker.ts
import * as Comlink from 'comlink';
import { SAEngine } from '../core/SAEngine';

// Expose class for RPC calls
Comlink.expose(SAEngine);

// main.ts
import * as Comlink from 'comlink';
import type { SAEngine } from './core/SAEngine';

const worker = new Worker(new URL('./worker/sa-worker.ts', import.meta.url), {
  type: 'module'
});
const SAEngineWorker = Comlink.wrap<typeof SAEngine>(worker);

// Usage: instantiate and call methods
const engine = await new SAEngineWorker({
  formula: 'C6H14',
  initialTemp: 100,
  coolingScheduleK: 8,
  stepsPerCycle: 500,
  numCycles: 4,
  optimizationMode: 'MINIMIZE',
  seed: 42
});
const result = await engine.run();
```

### Pattern 2: Alpine.js Reactive Controls
**What:** Two-way bound form inputs with validation feedback
**When to use:** Parameter inputs, play/pause controls, formula entry
**Example:**
```typescript
// Source: https://alpinejs.dev/directives/model (official docs)
// ui/controls.ts
export function controlsComponent() {
  return {
    formula: 'C6H14',
    initialTemp: 100,
    coolingScheduleK: 8,
    stepsPerCycle: 500,
    numCycles: 4,
    isRunning: false,
    isPaused: false,

    // Validation
    get isValidFormula() {
      return /^([A-Z][a-z]?\d*)+$/.test(this.formula);
    },

    // Controls
    async start() {
      if (!this.isValidFormula) return;
      this.isRunning = true;
      // Call worker...
    },

    pause() {
      this.isPaused = true;
      // Signal worker...
    },

    reset() {
      this.isRunning = false;
      this.isPaused = false;
      // Reset worker...
    }
  };
}

// HTML
<div x-data="controlsComponent()">
  <input type="text" x-model="formula" :class="{ 'error': !isValidFormula }">
  <button @click="start()" :disabled="!isValidFormula || isRunning">Start</button>
  <button @click="pause()" :disabled="!isRunning || isPaused">Pause</button>
  <button @click="reset()">Reset</button>
</div>
```

### Pattern 3: Vite Web Worker Import
**What:** Use Vite's `?worker` suffix or `new Worker(new URL(...))` pattern
**When to use:** Importing worker files with TypeScript support
**Example:**
```typescript
// Source: https://vite.dev/guide/features (official Vite docs)
// Recommended approach (standards-compliant)
const worker = new Worker(new URL('./worker.ts', import.meta.url), {
  type: 'module'
});

// Alternative (Vite-specific)
import MyWorker from './worker.ts?worker';
const worker = new MyWorker();
```

### Pattern 4: Progressive SA Execution
**What:** Report progress during long-running SA optimization
**When to use:** Real-time chart updates, responsiveness feedback
**Example:**
```typescript
// Source: https://www.smashingmagazine.com/2020/10/tasks-react-app-web-workers/
// Callback-based progress (Comlink supports function proxying)
import * as Comlink from 'comlink';

// worker/sa-worker.ts
class SAEngineWithProgress {
  async runWithProgress(params, onProgress) {
    // onProgress is a Comlink.proxy callback
    for (let step = 0; step < totalSteps; step++) {
      this.iterate(temp, step);
      if (step % 10 === 0) {
        await onProgress({
          step,
          currentEnergy: this.currentEnergy,
          bestEnergy: this.bestEnergy
        });
      }
    }
    return this.getResult();
  }
}

// main.ts
const result = await engine.runWithProgress(
  params,
  Comlink.proxy((progress) => {
    console.log(`Step ${progress.step}: energy=${progress.currentEnergy}`);
    updateChart(progress);
  })
);
```

### Pattern 5: RDKit.js Initialization
**What:** Load WASM module asynchronously before first use
**When to use:** App startup or lazy-loaded when needed
**Example:**
```typescript
// Source: https://github.com/rdkit/rdkit-js (official README)
// Option 1: CDN (simpler for classroom deployment)
<script src="https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js"></script>
<script>
  window.initRDKitModule()
    .then((RDKit) => {
      console.log("RDKit version: " + RDKit.version());
      window.RDKit = RDKit;
    })
    .catch((error) => {
      console.error("Failed to load RDKit:", error);
    });
</script>

// Option 2: Local files (better control)
window.initRDKitModule({
  locateFile: () => '/assets/RDKit_minimal.wasm'
}).then((RDKit) => {
  // RDKit ready
});
```

### Anti-Patterns to Avoid
- **Blocking main thread:** Never run SAEngine.run() directly in main thread (UI will freeze)
- **Deep nesting in Alpine:** Keep x-data flat; split complex state into multiple components
- **Worker without cleanup:** Always call `worker.terminate()` when done to prevent memory leaks
- **Unvalidated formulas:** Never pass formula to SAEngine without validation (parseFormula checks required)
- **Synchronous worker calls:** All Comlink-wrapped methods are async; always await

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| Worker message passing | Custom postMessage handlers with manual serialization | Comlink | Type safety, automatic serialization, function proxying, error propagation |
| Reactive UI updates | Manual DOM manipulation or event listeners | Alpine.js x-model, x-bind | Declarative bindings, automatic updates, less code |
| Module bundling for workers | Custom Rollup config or bundling logic | Vite's built-in worker support | Automatic code splitting, HMR, TypeScript transpilation |
| Molecular formula parsing | Custom regex parser | Phase 1 parseFormula() | Handles HDI validation, edge cases (existing, tested code) |
| PRNG seeding | Math.random() or custom LCG | Phase 1 SeededRandom (Mulberry32) | Deterministic, reproducible SA runs (existing, tested code) |

**Key insight:** This phase is about integration, not invention. Leverage existing libraries (Comlink, Alpine) and Phase 1 core logic rather than rebuilding primitives. The value is in the glue code, not custom implementations.

## Common Pitfalls

### Pitfall 1: Worker Initialization Race Conditions
**What goes wrong:** Attempting to call worker methods before Comlink proxy is ready or before RDKit WASM loads
**Why it happens:** Async initialization (initRDKitModule, worker instantiation) not properly awaited
**How to avoid:**
- Wrap worker creation in async function, await Comlink.wrap()
- Disable UI controls until initialization complete
- Use Alpine.js reactive flag: `rdkitReady: false` → show loading state
**Warning signs:**
- "RDKit is not defined" errors
- Worker method calls throw "Cannot read property of undefined"
- Intermittent failures on page refresh

### Pitfall 2: Large Data Transfer Performance
**What goes wrong:** Sending entire MolGraph adjacency matrices (N×N arrays) via postMessage causes slowdowns
**Why it happens:** Structured cloning copies data instead of transferring ownership
**How to avoid:**
- Return only necessary data from worker (final graph, not intermediate states)
- For progress updates, send primitives (numbers, booleans) not complex objects
- Use Transferable objects (ArrayBuffer) if sending large typed arrays
- Batch progress updates (every 10 steps, not every step)
**Warning signs:**
- Progress updates cause UI jank
- Worker execution fast but result retrieval slow
- Memory usage grows during SA run

### Pitfall 3: Alpine.js Reactivity with Worker State
**What goes wrong:** UI doesn't update when worker sends progress updates
**Why it happens:** Worker messages arrive outside Alpine's reactivity system
**How to avoid:**
- Update Alpine state inside message handler explicitly
- Use `this.property = value` in Alpine component methods
- For frequent updates, use `Alpine.$nextTick()` to batch DOM updates
**Warning signs:**
- Console logs show progress updates but UI frozen
- Manual DOM inspection shows stale values
- UI updates only after worker completes

### Pitfall 4: Formula Validation Edge Cases
**What goes wrong:** Valid formulas rejected or invalid formulas accepted
**Why it happens:** Regex validation alone doesn't check chemical validity (valence, HDI)
**How to avoid:**
- Use two-stage validation:
  1. Regex for format (`/^([A-Z][a-z]?\d*)+$/`)
  2. Call `parseFormula()` to check parsability and HDI
- Show specific error messages (format error vs. chemical impossibility)
- Provide examples of valid formulas (C6H14, C8H10, etc.)
**Warning signs:**
- User enters "C6H99" and sees generic error
- Valid organometallic formulas rejected
- SA starts but crashes with "invalid graph"

### Pitfall 5: Memory Leaks from Unmanaged Workers
**What goes wrong:** Workers accumulate when user starts multiple SA runs without cleanup
**Why it happens:** No call to `worker.terminate()` when run completes or user navigates away
**How to avoid:**
- Store worker reference in Alpine component
- Call `worker.terminate()` in reset() method
- Use `window.addEventListener('beforeunload', cleanup)` for page exit
- Create new worker only when needed, reuse if possible
**Warning signs:**
- Browser memory usage grows with each SA run
- Multiple workers visible in DevTools
- Performance degrades over time

### Pitfall 6: TypeScript Types Lost Across Worker Boundary
**What goes wrong:** Type errors at runtime despite TypeScript compilation passing
**Why it happens:** Comlink.Remote<T> approximates types but doesn't preserve exact signatures
**How to avoid:**
- Define explicit interfaces for worker API (don't rely on inferred types)
- Use `Comlink.wrap<typeof Class>` for class constructors
- Validate worker responses at runtime (especially for user input)
- Test worker integration with real data, not just types
**Warning signs:**
- Runtime errors like "method is not a function"
- Unexpected undefined values from worker
- Type assertions needed everywhere (`as unknown as Type`)

## Code Examples

Verified patterns from official sources:

### Example 1: Comlink Worker Initialization
```typescript
// Source: https://github.com/GoogleChromeLabs/comlink
// main.ts - main thread
import * as Comlink from 'comlink';
import type { SAEngine, SAParams, SAResult } from './core/SAEngine';

// Create worker with Vite's URL import
const worker = new Worker(
  new URL('./worker/sa-worker.ts', import.meta.url),
  { type: 'module' }
);

// Wrap worker for typed RPC calls
const SAEngineWorker = Comlink.wrap<typeof SAEngine>(worker);

// Usage
const params: SAParams = {
  formula: 'C6H14',
  initialTemp: 100,
  coolingScheduleK: 8,
  stepsPerCycle: 500,
  numCycles: 4,
  optimizationMode: 'MINIMIZE',
  seed: 42
};

const engineInstance = await new SAEngineWorker(params);
const result: SAResult = await engineInstance.run();

// Cleanup
worker.terminate();
```

### Example 2: Worker-Side Comlink Exposure
```typescript
// Source: https://github.com/GoogleChromeLabs/comlink
// worker/sa-worker.ts - worker thread
import * as Comlink from 'comlink';
import { SAEngine } from '../core/SAEngine';

// Expose entire class for instantiation on main thread
Comlink.expose(SAEngine);

// Alternative: Expose singleton instance
// const engine = new SAEngine(params);
// Comlink.expose(engine);
```

### Example 3: Alpine.js Form Validation Component
```typescript
// Source: https://alpinejs.dev/directives/model
// ui/formula-input.ts
export function formulaInputComponent() {
  return {
    formula: 'C6H14',
    errorMessage: '',

    validateFormula() {
      // Stage 1: Format check
      const formatRegex = /^([A-Z][a-z]?\d*)+$/;
      if (!formatRegex.test(this.formula)) {
        this.errorMessage = 'Invalid format. Use format like C6H14.';
        return false;
      }

      // Stage 2: Chemical validity (use Phase 1 parser)
      try {
        parseFormula(this.formula);
        this.errorMessage = '';
        return true;
      } catch (e) {
        this.errorMessage = e.message;
        return false;
      }
    },

    get isValid() {
      return this.errorMessage === '';
    }
  };
}

// HTML usage
<div x-data="formulaInputComponent()">
  <input
    type="text"
    x-model="formula"
    @blur="validateFormula()"
    :class="{ 'border-red-500': !isValid }"
  >
  <p x-show="!isValid" x-text="errorMessage" class="text-red-500"></p>
</div>
```

### Example 4: Alpine.js Play/Pause/Reset Controls
```typescript
// Source: https://alpinejs.dev/essentials/state (official Alpine docs)
// ui/sa-controls.ts
export function saControlsComponent() {
  return {
    state: 'idle', // 'idle' | 'running' | 'paused' | 'complete'
    currentStep: 0,
    totalSteps: 2000,

    async start() {
      this.state = 'running';
      // Initialize and start worker...
    },

    pause() {
      this.state = 'paused';
      // Signal worker to pause...
    },

    resume() {
      this.state = 'running';
      // Signal worker to resume...
    },

    reset() {
      this.state = 'idle';
      this.currentStep = 0;
      // Terminate worker, clear results...
    },

    get progress() {
      return (this.currentStep / this.totalSteps) * 100;
    },

    get canStart() {
      return this.state === 'idle';
    },

    get canPause() {
      return this.state === 'running';
    },

    get canResume() {
      return this.state === 'paused';
    }
  };
}

// HTML usage
<div x-data="saControlsComponent()">
  <button @click="start()" :disabled="!canStart">Start</button>
  <button @click="pause()" :disabled="!canPause">Pause</button>
  <button @click="resume()" :disabled="!canResume">Resume</button>
  <button @click="reset()">Reset</button>

  <div class="progress-bar">
    <div :style="`width: ${progress}%`"></div>
  </div>

  <p>Step <span x-text="currentStep"></span> of <span x-text="totalSteps"></span></p>
</div>
```

### Example 5: Vite Worker Configuration
```typescript
// Source: https://vite.dev/config/worker-options
// vite.config.ts
import { defineConfig } from 'vite';

export default defineConfig({
  worker: {
    format: 'es', // Use ES modules in workers (modern browsers)
    plugins: () => [], // Worker-specific plugins if needed
  },
  build: {
    target: 'es2020', // Modern JS for workers
  }
});
```

### Example 6: Progress Reporting Pattern
```typescript
// Source: https://blog.logrocket.com/comlink-web-workers-match-made-in-heaven/
// worker/sa-worker.ts
import * as Comlink from 'comlink';
import { SAEngine, SAParams } from '../core/SAEngine';

class SAEngineWrapper {
  private engine: SAEngine | null = null;

  async runWithProgress(params: SAParams, onProgress: (data: any) => void) {
    this.engine = new SAEngine(params);

    // Override iterate to report progress
    const totalSteps = params.stepsPerCycle * params.numCycles;

    // Run SA with periodic progress callbacks
    for (let step = 0; step < totalSteps; step++) {
      // Single iteration logic here...

      // Report progress every 10 steps
      if (step % 10 === 0) {
        await onProgress({
          step,
          currentEnergy: this.engine['currentEnergy'],
          bestEnergy: this.engine['bestEnergy'],
          temperature: /* ... */
        });
      }
    }

    return this.engine['getResult']();
  }
}

Comlink.expose(new SAEngineWrapper());

// main.ts
const worker = new Worker(new URL('./worker/sa-worker.ts', import.meta.url), {
  type: 'module'
});
const wrapper = Comlink.wrap(worker);

await wrapper.runWithProgress(
  params,
  Comlink.proxy((progress) => {
    // Update Alpine.js state
    Alpine.store('sa').currentStep = progress.step;
    Alpine.store('sa').currentEnergy = progress.currentEnergy;
  })
);
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Web Workers with postMessage | Comlink RPC abstraction | ~2018 (Comlink v1) | 90% less boilerplate, type-safe worker calls |
| importScripts() in workers | ES modules (import/export) | Vite 3+ (2022) | Modern syntax, tree-shaking, better DX |
| jQuery DOM manipulation | Alpine.js reactive bindings | Alpine 3.x (2021) | Declarative UI, automatic updates, smaller bundle |
| worker-loader (Webpack) | Vite ?worker suffix or new Worker() | Vite 2+ (2021) | Standards-compliant, no loader config |
| Manual WASM loading | initRDKitModule() promise | RDKit.js current | Cleaner API, error handling, async initialization |

**Deprecated/outdated:**
- **Worker suffixes (?worker)**: Still works but `new Worker(new URL(...))` is now recommended as more standards-compliant
- **importScripts()**: Use ES module imports in workers (Vite transpiles for production)
- **Comlink.proxy() for all functions**: Only needed for callbacks; regular return values work fine
- **Alpine 2.x**: Use Alpine 3.x (better TypeScript support, smaller core)

## Open Questions

1. **RDKit.js in Web Worker**
   - What we know: RDKit.js uses WASM which can run in workers, community has done it (@iktos/rdkit-provider)
   - What's unclear: Best initialization pattern (load in main then transfer? load in worker directly?)
   - Recommendation: Load RDKit in main thread first (simpler), only move to worker if initialization blocks UI (unlikely for minimal build)

2. **Pause/Resume SA Execution**
   - What we know: Current SAEngine.run() is a blocking loop, no built-in pause mechanism
   - What's unclear: Best pattern for pause/resume (async generator? manual step-by-step control?)
   - Recommendation: Refactor run() to accept a step callback that can return 'pause' signal, or use async generator pattern for next phase

3. **Optimal Progress Reporting Frequency**
   - What we know: Too frequent → UI jank from postMessage overhead; too rare → poor UX
   - What's unclear: Ideal reporting interval (every 10 steps? every 100ms?)
   - Recommendation: Start with every 10 steps (20 updates per 200-step cycle), make configurable, test on target hardware

4. **Comlink Memory Management**
   - What we know: worker.terminate() needed to prevent leaks
   - What's unclear: Does Comlink.proxy() for callbacks need manual cleanup? Memory footprint of long-running workers?
   - Recommendation: Test memory usage over multiple SA runs, add worker pool pattern if needed (low priority for v1)

5. **TypeScript Worker Types**
   - What we know: Comlink.wrap<typeof Class> provides best-effort types
   - What's unclear: Will SAEngine's private methods break type inference? Need explicit interface?
   - Recommendation: Create ISAEngine interface with public API only, use for Comlink.wrap typing

## Sources

### Primary (HIGH confidence)
- [Comlink GitHub](https://github.com/GoogleChromeLabs/comlink) - Official repo, API documentation, TypeScript types
- [Alpine.js Official Docs - x-model](https://alpinejs.dev/directives/model) - Reactive binding documentation
- [Alpine.js Official Docs - State](https://alpinejs.dev/essentials/state) - State management patterns
- [Vite Official Docs - Features](https://vite.dev/guide/features) - Web Worker support documentation
- [Vite Official Docs - Worker Options](https://vite.dev/config/worker-options) - Worker configuration reference
- [RDKit.js GitHub](https://github.com/rdkit/rdkit-js) - Official JavaScript distribution, installation, initialization patterns

### Secondary (MEDIUM confidence)
- [MDN Web Worker postMessage](https://developer.mozilla.org/en-US/docs/Web/API/Worker/postMessage) - Official browser API reference
- [MDN Transferable Objects](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Transferable_objects) - Performance optimization guide
- [LogRocket: Comlink and Web Workers](https://blog.logrocket.com/comlink-web-workers-match-made-in-heaven/) - Tutorial verified against official docs
- [Smashing Magazine: Long-Running Tasks in React](https://www.smashingmagazine.com/2020/10/tasks-react-app-web-workers/) - Worker patterns (framework-agnostic principles)
- [David East: Comlink Article](https://davidea.st/articles/comlink-simple-web-worker/) - Google engineer, practical patterns

### Tertiary (LOW confidence - needs validation)
- Various Alpine.js community discussions on reactivity edge cases - general patterns match docs
- Molecular formula regex patterns - need to validate against chemistry rules (use Phase 1 parseFormula as source of truth)
- RDKit.js bundle size claims - need to verify with actual npm package inspection or build analysis

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All libraries verified from official sources, versions confirmed, ecosystem standard
- Architecture: HIGH - Patterns from official docs and verified tutorials, proven in production
- Pitfalls: MEDIUM-HIGH - Common issues documented in GitHub issues and community posts, some from experience reports
- Code examples: HIGH - All examples based on official documentation or verified tutorials
- Open questions: MEDIUM - Identified gaps are real but don't block implementation (can test/iterate)

**Research date:** 2026-02-15
**Valid until:** ~2026-04-15 (60 days - stack is mature/stable, slow-moving ecosystem)

**Notes:**
- No CONTEXT.md existed for this phase (user constraints section omitted)
- Tech stack was pre-decided in project initialization (Vite, TypeScript, Alpine.js, Comlink, RDKit.js)
- Phase 1 codebase is complete and tested (139 tests passing) - integration is the focus
- Classroom deployment requirement (no backend) confirmed feasible with this stack
