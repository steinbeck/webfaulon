# Phase 6: Frontend Integration - Research

**Researched:** 2026-02-16
**Domain:** Frontend-Backend Integration (Alpine.js + SSE + REST API)
**Confidence:** HIGH

## Summary

Phase 6 transforms the frontend from a self-contained Web Worker + WASM architecture to a thin client that communicates entirely through backend REST + SSE APIs. The core challenge is replacing Comlink-based worker communication with EventSource SSE streaming while maintaining the same responsive UX.

The research reveals a mature, well-understood technology stack. EventSource (SSE) provides automatic reconnection and is native to browsers. Alpine.js lacks explicit cleanup hooks but offers `destroy()` method patterns for third-party library teardown. Chart.js performance optimizations (animations disabled, `update('none')`) are well-documented and proven. SVG rendering via innerHTML is straightforward. The primary integration risks are EventSource memory leaks from missing cleanup and state synchronization between REST control calls and SSE stream events.

**Primary recommendation:** Use native EventSource with manual cleanup in Alpine's `destroy()` method. Proxy backend API through Vite dev server to avoid CORS. Keep Chart.js instance module-scoped (not Alpine-reactive). Display SVG via innerHTML. Create TypeScript types matching backend Pydantic models for type safety.

## Standard Stack

### Core
| Library | Version | Purpose | Why Standard |
|---------|---------|---------|--------------|
| EventSource (native) | Browser API | SSE streaming from backend | Native browser API, automatic reconnection, wide support |
| Alpine.js | 3.15.8 (existing) | UI reactivity and state management | Already in use, lightweight, reactive without virtual DOM |
| Chart.js | 4.5.1 (existing) | Live chart visualization | Already in use, proven real-time performance patterns |
| Fetch API (native) | Browser API | REST API calls (configure, control, status) | Native, TypeScript-friendly, standard HTTP client |

### Supporting
| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| Vite proxy | 6.0.11 (existing) | Dev server proxy to avoid CORS | Development only (proxy `/api` → `http://localhost:8000`) |
| TypeScript | 5.6.3 (existing) | Type-safe API contracts | Define types matching backend Pydantic models |

### Alternatives Considered
| Instead of | Could Use | Tradeoff |
|------------|-----------|----------|
| EventSource | fetch with ReadableStream | EventSource provides automatic reconnection, simpler API, SSE standard |
| Native fetch | Axios library | Fetch is native, sufficient for this use case, no dependencies needed |
| Alpine.js | React/Vue | Alpine already integrated, no need to rewrite entire frontend |

**Installation:**
```bash
# No new packages needed - using browser native APIs
# Backend already serves SSE via sse-starlette
```

## Architecture Patterns

### Recommended Project Structure
```
src/
├── ui/
│   ├── app.ts              # Alpine app component (replace worker with SSE)
│   ├── chart.ts            # Chart.js integration (keep module-scoped)
│   ├── molecule-renderer.ts # NEW: SVG renderer (replace RDKit canvas)
│   ├── validation.ts       # Formula validation (keep client-side)
│   └── presets.ts          # Molecule presets (keep)
├── api/
│   ├── client.ts           # NEW: Typed REST API client (fetch wrapper)
│   ├── types.ts            # NEW: TypeScript types (mirror backend Pydantic)
│   └── sse.ts              # NEW: SSE EventSource wrapper with cleanup
└── main.ts                 # Entry point (Alpine.start())
```

### Pattern 1: EventSource with Cleanup
**What:** Wrap EventSource in a class with explicit cleanup to prevent memory leaks
**When to use:** Always when using SSE in SPAs

**Example:**
```typescript
// src/api/sse.ts
class SSEConnection {
  private eventSource: EventSource | null = null;

  connect(sessionId: string, onProgress: (data: any) => void) {
    this.close(); // Cleanup any existing connection
    this.eventSource = new EventSource(`/api/sa/${sessionId}/stream`);

    this.eventSource.addEventListener('progress', (e) => {
      onProgress(JSON.parse(e.data));
    });

    this.eventSource.addEventListener('complete', (e) => {
      console.log('SA complete:', JSON.parse(e.data));
      this.close();
    });

    this.eventSource.addEventListener('error', (e) => {
      console.error('SSE error:', e);
      // EventSource auto-reconnects by default (3s delay)
    });
  }

  close() {
    if (this.eventSource) {
      this.eventSource.close();
      this.eventSource = null;
    }
  }
}
```
**Source:** [MDN EventSource documentation](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events)

### Pattern 2: Alpine Component with destroy() Method
**What:** Use Alpine's `destroy()` method for cleanup of third-party resources
**When to use:** When integrating libraries that need teardown (EventSource, Chart.js)

**Example:**
```typescript
// src/ui/app.ts
export function appComponent() {
  const sseConnection = new SSEConnection();

  return {
    // ... reactive data

    init() {
      // Initialize resources
    },

    destroy() {
      // Alpine calls this automatically on teardown
      sseConnection.close();
      // Chart.js cleanup handled in chart.ts module
    },

    async start() {
      sseConnection.connect(this.sessionId, (data) => {
        this.progress = data; // Alpine reactivity updates UI
        addChartDataPoint(data.step, data.best_energy);
      });
    }
  };
}
```
**Source:** [Alpine.js component cleanup discussion](https://github.com/alpinejs/alpine/discussions/2218)

### Pattern 3: Non-Reactive Chart.js Instance
**What:** Keep Chart.js instance module-scoped to avoid Alpine proxy wrapping
**When to use:** Always - Chart.js mutates DOM directly, conflicts with Alpine reactivity

**Example:**
```typescript
// src/ui/chart.ts (existing pattern - keep it)
let _chartInstance: Chart | null = null;

export function createWienerChart(canvas: HTMLCanvasElement): void {
  if (_chartInstance) _chartInstance.destroy();
  _chartInstance = new Chart(canvas, {
    type: 'line',
    data: { datasets: [{ data: [] }] },
    options: {
      animation: false, // CRITICAL for performance
      parsing: false,   // CRITICAL for performance
      // ... other options
    }
  });
}

export function addChartDataPoint(step: number, wienerIndex: number): void {
  if (!_chartInstance) return;
  _chartInstance.data.datasets[0].data.push({ x: step, y: wienerIndex });
  _chartInstance.update('none'); // 'none' = no animation, fastest
}
```
**Source:** [Chart.js Performance Docs](https://www.chartjs.org/docs/latest/general/performance.html), [Alpine + Chart.js integration](https://janostlund.com/2024-02-11/integrating-chartjs-with-alpine)

### Pattern 4: SVG Rendering via innerHTML
**What:** Display backend-generated SVG by setting innerHTML of a container
**When to use:** Replacing canvas-based molecule rendering with backend SVG

**Example:**
```typescript
// src/ui/molecule-renderer.ts
export function renderMoleculeSVG(svg: string, container: HTMLElement): void {
  if (!svg || !container) return;
  container.innerHTML = svg; // Browser parses SVG markup
}

export function clearMoleculeDisplay(container: HTMLElement): void {
  container.innerHTML = '<p class="hint">No structure yet</p>';
}
```
**Why innerHTML:** Browser interprets CSS rules properly, URL-encoded SVG is 20-30% smaller than Base64
**Source:** [SVG data URI optimization](https://www.svgbackgrounds.com/data-uris-are-a-wildly-underused-website-speed-optimization/)

### Pattern 5: Typed REST API Client
**What:** Create typed fetch wrapper matching backend Pydantic models
**When to use:** All REST API calls (configure, control, status)

**Example:**
```typescript
// src/api/types.ts (mirror backend/app/models/sa_params.py)
export interface SAParams {
  formula: string;
  initial_temp: number;
  cooling_schedule_k: number;
  steps_per_cycle: number;
  num_cycles: number;
  optimization_mode: 'MINIMIZE' | 'MAXIMIZE';
  seed: number;
}

export interface ConfigureResponse {
  session_id: string;
  status: string;
}

// src/api/client.ts
export class SAAPIClient {
  private baseURL = '/api/sa'; // Vite proxies to localhost:8000

  async configure(params: SAParams): Promise<ConfigureResponse> {
    const response = await fetch(`${this.baseURL}/configure`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(params),
    });

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || 'Configure failed');
    }

    return response.json();
  }

  async start(sessionId: string): Promise<void> {
    const response = await fetch(`${this.baseURL}/${sessionId}/start`, {
      method: 'POST',
    });
    if (!response.ok) throw new Error('Start failed');
  }

  // ... pause, reset methods
}
```
**Source:** [Fetch error handling with TypeScript](https://jessewarden.com/2025/02/error-handling-for-fetch-in-typescript.html)

### Pattern 6: Vite Dev Proxy Configuration
**What:** Proxy `/api` requests to backend during development to avoid CORS
**When to use:** Development only (production uses same-origin or CORS headers)

**Example:**
```typescript
// vite.config.ts
export default defineConfig({
  base: '/webfaulon/',
  server: {
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
        // rewrite not needed - backend already serves at /api
      }
    }
  },
  // ... rest of config
});
```
**Source:** [Vite proxy configuration guide](https://medium.com/@eric_abell/simplifying-api-proxies-in-vite-a-guide-to-vite-config-js-a5cc3a091a2f)

### Anti-Patterns to Avoid
- **EventSource without cleanup:** Causes memory leaks, zombie connections on page refresh. Always call `.close()` in Alpine `destroy()` method.
- **Chart.js in Alpine reactive data:** Causes proxy-of-proxy issues, infinite recursion. Keep Chart.js module-scoped.
- **Base64-encoded SVG:** 20-30% larger than URL-encoded SVG. Backend should send URL-encoded or raw SVG.
- **Animations enabled on Chart.js:** Redundant redraws kill performance. Always set `animation: false` for real-time data.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead | Why |
|---------|-------------|-------------|-----|
| SSE reconnection logic | Custom retry/backoff | Native EventSource | Automatic reconnection (3s default), retry header support, Last-Event-ID for deduplication |
| REST client with types | Manual fetch everywhere | Typed API client class | Centralized error handling, type safety, single source of truth |
| CORS during dev | Custom server setup | Vite proxy | Built-in, zero config, works with SSE and fetch |
| SVG sanitization | Manual XSS prevention | Backend RDKit output + innerHTML | RDKit generates safe SVG, browser sandboxes SVG by default |

**Key insight:** Browser APIs (EventSource, fetch) and Vite tooling already solve 90% of integration challenges. Custom solutions add complexity without benefit. The only custom code needed is the API client wrapper and SSE connection manager.

## Common Pitfalls

### Pitfall 1: EventSource Memory Leaks
**What goes wrong:** EventSource connections are not closed when user navigates away or component unmounts, causing memory leaks and zombie server connections.
**Why it happens:** Browser's `beforeunload` handling is inconsistent. Chrome delays closure up to 60s. Alpine has no automatic cleanup for EventSource.
**How to avoid:**
- Always store EventSource in a module-scoped variable or class instance
- Call `.close()` explicitly in Alpine's `destroy()` method
- Consider `window.addEventListener('beforeunload', () => eventSource.close())` as backup
**Warning signs:** Server logs show multiple active connections for same session. Browser DevTools shows multiple EventSource objects in memory.
**Source:** [EventSource cleanup best practices](https://tigerabrodi.blog/server-sent-events-a-practical-guide-for-the-real-world)

### Pitfall 2: Chart.js Alpine Proxy Conflict
**What goes wrong:** Storing Chart.js instance in Alpine reactive data (`this.chart = new Chart(...)`) causes "Maximum call stack size exceeded" or chart stops updating.
**Why it happens:** Alpine wraps data in a Proxy for reactivity. Chart.js is already a complex object with its own internal state. Proxy-wrapping it creates conflicts.
**How to avoid:**
- Keep Chart.js instance module-scoped (`let _chartInstance: Chart | null`)
- Export functions that operate on the module-scoped instance
- Never assign Chart.js instance to Alpine reactive properties
**Warning signs:** Console errors about stack overflow, chart stops rendering updates, Alpine devtools shows huge nested object.
**Source:** [Alpine + Chart.js integration patterns](https://janostlund.com/2024-02-11/integrating-chartjs-with-alpine)

### Pitfall 3: SSE Event Stream Buffering
**What goes wrong:** SSE events arrive in batches or with delays instead of real-time, making live visualization laggy.
**Why it happens:** Nginx/proxies buffer SSE responses by default. Backend needs to set `X-Accel-Buffering: no` header.
**How to avoid:**
- Backend already sets header in Phase 5 (`EventSourceResponse(..., headers={"X-Accel-Buffering": "no"})`)
- Verify in browser DevTools Network tab that header is present
- If using reverse proxy (nginx, Cloudflare), configure SSE pass-through
**Warning signs:** Events arrive in bursts instead of smoothly. Network tab shows delays between chunks.
**Source:** [SSE buffering issues](https://www.infoq.com/articles/reactive-notification-system-server-sent-events/)

### Pitfall 4: Fast Chart.js Updates Cause Crashes
**What goes wrong:** Rapidly updating Chart.js (e.g., every SSE event at 100 events/sec) causes browser freezes or crashes.
**Why it happens:** Even with `animation: false`, calling `chart.update()` triggers layout recalculation. Too many updates saturate the main thread.
**How to avoid:**
- Use `chart.update('none')` instead of `chart.update()` - bypasses some internal checks
- Batch updates: push data to array, update chart every N events (e.g., every 10 steps)
- Backend already throttles SSE to ~100 events/sec with `asyncio.sleep(0.01)`
- Consider `requestAnimationFrame()` batching if updates are still too fast
**Warning signs:** Browser becomes unresponsive during SA run. DevTools performance profiler shows chart update calls dominating CPU.
**Source:** [Chart.js performance optimization](https://www.chartjs.org/docs/latest/general/performance.html)

### Pitfall 5: Stale Session State After Network Interruption
**What goes wrong:** After network reconnection, frontend and backend session states diverge (e.g., frontend shows "running" but backend is paused).
**Why it happens:** EventSource auto-reconnects but doesn't re-sync session state. Frontend keeps old reactive data.
**How to avoid:**
- On EventSource `error` event, call `/api/sa/{id}/status` to re-sync state
- Use `Last-Event-ID` header for event deduplication (EventSource does this automatically)
- Reset frontend state on reconnection: `this.progress = null` then fetch status
**Warning signs:** Start/pause buttons become desynced. Chart stops updating but SSE connection is alive.
**Source:** [SSE reconnection patterns](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events)

### Pitfall 6: TypeScript Type Mismatch Between Frontend and Backend
**What goes wrong:** Frontend expects `initialTemp` but backend sends `initial_temp`. Runtime errors, silent data loss.
**Why it happens:** Python uses snake_case, TypeScript uses camelCase. Manual type definitions drift from backend Pydantic models.
**How to avoid:**
- Define TypeScript interfaces matching backend Pydantic models EXACTLY (use snake_case)
- Use a type generation tool OR manually copy from backend models and comment source file
- Add runtime validation in API client (throw on unexpected keys)
**Warning signs:** Data appears undefined in DevTools even though API returns 200. Chart shows NaN values.
**Source:** Backend models in `backend/app/models/sa_params.py`

## Code Examples

Verified patterns from official sources and existing codebase:

### EventSource Connection with Cleanup
```typescript
// src/api/sse.ts
export class SSEConnection {
  private eventSource: EventSource | null = null;
  private sessionId: string | null = null;

  connect(
    sessionId: string,
    handlers: {
      onProgress?: (data: any) => void;
      onComplete?: (data: any) => void;
      onError?: (event: Event) => void;
    }
  ): void {
    this.close(); // Cleanup previous connection
    this.sessionId = sessionId;
    this.eventSource = new EventSource(`/api/sa/${sessionId}/stream`);

    if (handlers.onProgress) {
      this.eventSource.addEventListener('progress', (e: MessageEvent) => {
        handlers.onProgress!(JSON.parse(e.data));
      });
    }

    if (handlers.onComplete) {
      this.eventSource.addEventListener('complete', (e: MessageEvent) => {
        handlers.onComplete!(JSON.parse(e.data));
        this.close(); // Auto-close on complete
      });
    }

    this.eventSource.onerror = (event) => {
      console.error('SSE connection error:', event);
      if (handlers.onError) handlers.onError(event);
      // EventSource auto-reconnects by default
    };
  }

  close(): void {
    if (this.eventSource) {
      this.eventSource.close();
      this.eventSource = null;
    }
  }

  isConnected(): boolean {
    return this.eventSource?.readyState === EventSource.OPEN;
  }
}
```

### Alpine App Component Integration
```typescript
// src/ui/app.ts
import { SAAPIClient } from '../api/client';
import { SSEConnection } from '../api/sse';
import { addChartDataPoint, resetChart } from './chart';
import { renderMoleculeSVG, clearMoleculeDisplay } from './molecule-renderer';

export function appComponent() {
  const apiClient = new SAAPIClient();
  const sseConnection = new SSEConnection();

  return {
    // State
    formula: '',
    sessionId: null as string | null,
    state: 'idle' as 'idle' | 'running' | 'paused' | 'complete',
    progress: null as any,

    // Alpine lifecycle
    init() {
      // Chart initialization already done in chart.ts via queueMicrotask
    },

    destroy() {
      // Cleanup on component unmount
      sseConnection.close();
      // Chart.js cleanup handled by chart.ts module
    },

    // Methods
    async start() {
      try {
        // 1. Configure session
        const response = await apiClient.configure({
          formula: this.formula,
          initial_temp: this.initialTemp,
          cooling_schedule_k: this.coolingScheduleK,
          steps_per_cycle: this.stepsPerCycle,
          num_cycles: this.numCycles,
          optimization_mode: this.optimizationMode,
          seed: Date.now(),
        });
        this.sessionId = response.session_id;

        // 2. Start session
        await apiClient.start(this.sessionId);

        // 3. Connect SSE for progress updates
        sseConnection.connect(this.sessionId, {
          onProgress: (data) => {
            this.progress = data;
            addChartDataPoint(data.step, data.best_energy);

            // Update molecule SVG
            const container = document.getElementById('molecule-display');
            if (container && data.best_svg) {
              renderMoleculeSVG(data.best_svg, container);
            }
          },
          onComplete: (data) => {
            this.state = 'complete';
            this.progress = data;
          },
          onError: (event) => {
            console.error('SSE error, re-syncing state');
            this.resyncState();
          },
        });

        this.state = 'running';
      } catch (error) {
        console.error('Start failed:', error);
        alert('Failed to start: ' + (error instanceof Error ? error.message : 'Unknown error'));
      }
    },

    async pause() {
      if (!this.sessionId) return;
      await apiClient.pause(this.sessionId);
      this.state = 'paused';
    },

    async resume() {
      if (!this.sessionId) return;
      await apiClient.resume(this.sessionId);
      this.state = 'running';
    },

    async reset() {
      if (!this.sessionId) return;
      sseConnection.close();
      await apiClient.reset(this.sessionId);
      resetChart();
      const container = document.getElementById('molecule-display');
      if (container) clearMoleculeDisplay(container);
      this.state = 'idle';
      this.progress = null;
    },

    async resyncState() {
      if (!this.sessionId) return;
      const status = await apiClient.getStatus(this.sessionId);
      this.state = status.session_state as any;
      this.progress = status;
    },
  };
}
```

### Molecule SVG Renderer (Replaces RDKit Canvas)
```typescript
// src/ui/molecule-renderer.ts
export function renderMoleculeSVG(svg: string, container: HTMLElement): void {
  if (!svg || !container) return;

  // Backend sends URL-encoded or raw SVG string
  // innerHTML parses and renders SVG correctly
  container.innerHTML = svg;
}

export function clearMoleculeDisplay(container: HTMLElement): void {
  if (!container) return;
  container.innerHTML = `
    <div style="display: flex; align-items: center; justify-content: center; height: 100%; color: #9ca3af;">
      <p>No structure yet</p>
    </div>
  `;
}
```

### TypeScript API Types (Mirror Backend)
```typescript
// src/api/types.ts
// Mirror backend/app/models/sa_params.py

export interface SAParams {
  formula: string;
  initial_temp: number;
  cooling_schedule_k: number;
  steps_per_cycle: number;
  num_cycles: number;
  optimization_mode: 'MINIMIZE' | 'MAXIMIZE';
  seed: number;
}

export interface ConfigureResponse {
  session_id: string;
  status: string;
}

export interface ControlResponse {
  session_id: string;
  status: string;
}

export interface StatusResponse {
  session_id: string;
  session_state: string;
  step: number;
  total_steps: number;
  cycle: number;
  current_energy: number;
  best_energy: number;
  best_smiles: string;
  best_svg: string;
  temperature: number;
  accepted_moves: number;
  rejected_moves: number;
  invalid_moves: number;
  is_complete: boolean;
}

// SSE event data types
export interface SSEProgressData {
  step: number;
  total_steps: number;
  cycle: number;
  temperature: number;
  current_energy: number;
  best_energy: number;
  best_smiles: string;
  best_svg: string;
  accepted_moves: number;
  rejected_moves: number;
  invalid_moves: number;
  is_complete: boolean;
}

export interface SSECompleteData {
  step: number;
  total_steps: number;
  best_energy: number;
  best_smiles: string;
  best_svg: string;
  accepted_moves: number;
  rejected_moves: number;
  invalid_moves: number;
}
```

### REST API Client
```typescript
// src/api/client.ts
import type { SAParams, ConfigureResponse, ControlResponse, StatusResponse } from './types';

export class SAAPIClient {
  private baseURL = '/api/sa';

  async configure(params: SAParams): Promise<ConfigureResponse> {
    const response = await fetch(`${this.baseURL}/configure`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(params),
    });

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || `HTTP ${response.status}: Configure failed`);
    }

    return response.json();
  }

  async start(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'start');
  }

  async pause(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'pause');
  }

  async resume(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'start'); // Resume = start (backend handles state)
  }

  async reset(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'reset');
  }

  async getStatus(sessionId: string): Promise<StatusResponse> {
    const response = await fetch(`${this.baseURL}/${sessionId}/status`);
    if (!response.ok) {
      throw new Error(`HTTP ${response.status}: Get status failed`);
    }
    return response.json();
  }

  private async controlRequest(sessionId: string, action: 'start' | 'pause' | 'reset'): Promise<ControlResponse> {
    const response = await fetch(`${this.baseURL}/${sessionId}/${action}`, {
      method: 'POST',
    });

    if (!response.ok) {
      const error = await response.json();
      throw new Error(error.error || `HTTP ${response.status}: ${action} failed`);
    }

    return response.json();
  }
}
```

## State of the Art

| Old Approach | Current Approach | When Changed | Impact |
|--------------|------------------|--------------|--------|
| Comlink + Web Worker | EventSource SSE | 2020s | SSE is simpler for unidirectional streaming, native browser API, no library needed |
| RDKit.js WASM in browser | Backend Python RDKit + SVG | 2025 (v2.0) | Eliminates 2.5MB WASM download, faster page load, simpler frontend |
| Callback-based fetch | Async/await fetch | 2017+ | Cleaner error handling, native TypeScript support |
| Chart.js 2.x animations | Chart.js 4.x with `animation: false` | 2022 | Performance improvement for real-time data, explicit opt-out |

**Deprecated/outdated:**
- **Comlink Web Workers for SA:** Replaced by backend execution. Workers still useful for CPU-heavy frontend tasks, but SA moved server-side.
- **RDKit.js WASM:** Replaced by backend RDKit. Still valid for offline apps or client-side chem tools, but not needed here.
- **XMLHttpRequest:** Replaced by fetch API. Still supported but deprecated pattern.

## Open Questions

1. **HTTP/2 for SSE connection pooling**
   - What we know: SSE limited to 6 concurrent connections on HTTP/1.1 per domain
   - What's unclear: Whether GitHub Pages supports HTTP/2 (needed for production SSE)
   - Recommendation: Test on GitHub Pages staging. If HTTP/1.1 only, document limitation (1 concurrent SA run per browser). Consider HTTP/2 CDN if needed.

2. **EventSource polyfill for legacy browsers**
   - What we know: EventSource supported in all modern browsers
   - What's unclear: Whether to support IE11 or older Safari versions
   - Recommendation: Skip polyfill unless user analytics show legacy traffic. Modern browsers only.

3. **SSE event ID usage for resume-from-step**
   - What we know: EventSource sends `Last-Event-ID` header on reconnect
   - What's unclear: Whether to implement server-side resume-from-step using event IDs
   - Recommendation: Not needed for v2.0 - SA runs are short (<30s). Full restart on reconnect is acceptable. Mark as future enhancement.

## Sources

### Primary (HIGH confidence)
- [MDN EventSource API](https://developer.mozilla.org/en-US/docs/Web/API/Server-sent_events/Using_server-sent_events) - Official browser API documentation
- [Alpine.js Reactivity](https://alpinejs.dev/advanced/reactivity) - Official framework docs
- [Chart.js Performance](https://www.chartjs.org/docs/latest/general/performance.html) - Official optimization guide
- Backend API: `backend/app/api/` (sa_stream.py, sa_configure.py, sa_control.py, sa_status.py)
- Backend models: `backend/app/models/sa_params.py`
- Existing frontend: `src/ui/app.ts`, `src/ui/chart.ts`

### Secondary (MEDIUM confidence)
- [Alpine.js + Chart.js Integration](https://janostlund.com/2024-02-11/integrating-chartjs-with-alpine) - Community pattern (2024)
- [SSE Practical Guide](https://tigerabrodi.blog/server-sent-events-a-practical-guide-for-the-real-world) - Community best practices
- [Vite Proxy Configuration](https://medium.com/@eric_abell/simpl simplifying-api-proxies-in-vite-a-guide-to-vite-config-js-a5cc3a091a2f) - Dev workflow pattern
- [Fetch TypeScript Patterns](https://jessewarden.com/2025/02/error-handling-for-fetch-in-typescript.html) - Type-safe error handling (2025)

### Tertiary (LOW confidence)
- [Alpine.js cleanup discussion](https://github.com/alpinejs/alpine/discussions/2218) - Community workarounds (no official API)
- [SVG data URI performance](https://www.svgbackgrounds.com/data-uris-are-a-wildly-underused-website-speed-optimization/) - General web perf (not Alpine-specific)

## Metadata

**Confidence breakdown:**
- Standard stack: HIGH - All technologies are mature, well-documented browser/framework APIs
- Architecture: HIGH - Patterns verified in existing codebase (chart.ts) and official docs (EventSource, Alpine)
- Pitfalls: HIGH - Documented in official sources (Chart.js perf, EventSource cleanup) and community discussions

**Research date:** 2026-02-16
**Valid until:** 2026-03-16 (30 days - stable tech stack, slow-moving browser APIs)
