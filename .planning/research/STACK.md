# Technology Stack

**Project:** WebFaulon
**Researched:** 2026-02-14
**Confidence:** MEDIUM

## Recommended Stack

### Core Framework

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **Vite** | ^6.x | Build tool and dev server | Industry standard for modern browser apps in 2026. Provides instant HMR, esbuild-powered TypeScript transpilation (20-30x faster than traditional tools), and optimal WASM integration via plugins. Zero-config dev experience with production-ready optimizations. |
| **TypeScript** | ^5.x | Type-safe development | Mandatory for maintainable scientific applications. Vite's esbuild integration provides instantaneous type stripping. Use strict mode for early error detection in graph algorithms. |
| **RDKit.js** | 2025.3.4-1.0.0 | Molecular graph manipulation and 2D rendering | Official JavaScript port of RDKit (WASM-based). Only production-ready cheminformatics library for browser. Handles SMILES parsing, graph operations, and native 2D structure rendering. Published 7 months ago (latest stable). |

### UI and Visualization

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **Apache ECharts** | ^5.x | Real-time SA progress charting | Canvas-based rendering handles live updates efficiently. Renders 10,000 points <100ms. Hybrid rendering (Canvas/SVG/WebGL) for performance at scale. Superior to Chart.js for IoT/scientific data streams. Well-documented streaming data patterns. |
| **Alpine.js** | ^3.x | Lightweight UI interactivity | 7.1KB framework perfect for parameter controls and UI state without SPA overhead. Declarative syntax in HTML reduces build complexity. Ideal for educational demos where students read source. Alternative: vanilla JS for minimal footprint. |

### Web Worker Infrastructure

| Technology | Version | Purpose | Why |
|------------|---------|---------|-----|
| **Comlink** | ^4.x | Web Worker RPC abstraction | Google Chrome Labs library (1.1KB) that eliminates postMessage boilerplate. Makes worker communication feel synchronous via ES6 Proxies. TypeScript-native with `Comlink.Remote<T>` types. Use with `vite-plugin-comlink` for seamless Vite integration. |

### Supporting Libraries

| Library | Version | Purpose | When to Use |
|---------|---------|---------|-------------|
| **vite-plugin-wasm** | ^3.x | WASM module loading | Required for RDKit.js WASM bundle. Handles `.wasm` file imports and ESM integration proposal limitations in Vite. |
| **vite-plugin-comlink** | ^5.x | Comlink integration | Simplifies worker setup with Vite. Auto-handles module resolution and HMR for workers. |

### Development Tools

| Tool | Purpose | Notes |
|------|---------|-------|
| **Vitest** | Unit and integration testing | Use `@vitest/web-worker` for worker tests. WASM testing has known issues; use globalSetup for WASM initialization. Browser mode recommended for RDKit.js integration tests. |
| **ESLint** | Code quality | TypeScript-aware rules. Catch algorithm bugs early. |
| **Prettier** | Code formatting | Consistent style for classroom code review. |

## Installation

```bash
# Core
npm install @rdkit/rdkit apache-echarts alpinejs comlink

# Build tooling
npm install -D vite typescript vite-plugin-wasm vite-plugin-comlink

# Dev dependencies
npm install -D vitest @vitest/web-worker eslint prettier @typescript-eslint/parser @typescript-eslint/eslint-plugin
```

## Alternatives Considered

| Category | Recommended | Alternative | Why Not |
|----------|-------------|-------------|---------|
| **Charting** | Apache ECharts | Chart.js | Chart.js degrades with large datasets. ECharts' canvas rendering is 10x faster for >1000 points. Chart.js has 2M+ weekly downloads (popular) but not optimized for scientific streaming data. |
| **Charting** | Apache ECharts | ApexCharts | ApexCharts performance degrades with multiple charts on same page. ECharts handles complex dashboards better. ApexCharts easier API but not needed for single SA chart. |
| **Charting** | Apache ECharts | LightningChart JS | LightningChart is GPU-accelerated (WebGL) overkill for SA visualization. Commercial license required. ECharts free and sufficient for classroom demo scale. |
| **UI Framework** | Alpine.js | React | React adds 40KB+ overhead unnecessary for parameter form + start button. Alpine's in-HTML syntax more educational (students see interactivity in source). No build complexity. |
| **UI Framework** | Alpine.js | Preact | Preact (3KB) lighter than React but still component-model overhead. Alpine's declarative attributes simpler for non-framework developers extending demo. |
| **UI Framework** | Alpine.js | Vanilla JS | Viable alternative. Choose Alpine if you want reactive data binding without manual DOM updates. Choose vanilla if targeting zero dependencies. |
| **Worker Communication** | Comlink | Manual postMessage | postMessage requires serialization boilerplate. Comlink's proxy pattern reduces worker code by ~60%. TypeScript types preserved across thread boundary. |

## What NOT to Use

| Avoid | Why | Use Instead |
|-------|-----|-------------|
| **CDN imports for RDKit.js** | WASM file requires proper CORS and path handling. npm package with Vite bundler ensures correct resolution. Unpkg CDN mentioned in docs but problematic for local dev. | `npm install @rdkit/rdkit` + `vite-plugin-wasm` |
| **Chart.js for SA visualization** | Canvas rendering struggles with live updates >1000 points. Documentation admits "performance issues with very large datasets." SA may generate 10K+ steps. | Apache ECharts with streaming data configuration |
| **Web Worker without transfer/Comlink** | Large molecular graphs (100+ atoms) will cause main thread jank if cloned via postMessage. Structured clone has overhead. | Comlink with transferable objects or shared array buffers for state |
| **navigator.hardwareConcurrency for worker pool** | Over-engineering for single SA run. Algorithm `(2 * cores) + 1` only relevant for parallelized SA (not in Faulon 1996). Single worker sufficient. | One dedicated worker for SA algorithm |
| **Service Workers** | Not needed for in-browser demo. Adds caching complexity without benefit. Confusion with Web Workers. | Web Workers for computation only |

## Stack Patterns by Variant

**If deploying to GitHub Pages (static hosting):**
- Use Vite's `base` option for correct asset paths
- Pre-load WASM in index.html to avoid CORS issues
- Configure `publicPath` for workers

**If adding backend persistence later:**
- Keep all SA computation client-side (educational transparency)
- Only use backend for saving/loading parameter presets
- WASM + Workers remain browser-only

**If targeting low-end devices (Chromebooks):**
- Reduce default SA steps to <5000
- Add progress throttling (update chart every 10th step, not every step)
- Consider switching to Chart.js if ECharts bundle size is prohibitive (test <100KB budget)

## Version Compatibility

| Package A | Compatible With | Notes |
|-----------|-----------------|-------|
| @rdkit/rdkit@2025.3.4-1.0.0 | Node.js >=14 | WASM requires browser with WebAssembly support (Chrome 57+, Firefox 52+, Safari 11+). No IE11. |
| vite-plugin-wasm@^3.x | Vite ^5.x or ^6.x | Check plugin compatibility with Vite major version. Breaking changes in Vite 6 around worker imports. |
| comlink@^4.x | vite-plugin-comlink@^5.x | Plugin abstracts Comlink setup. Version sync ensures proxy transfer works correctly. |
| apache-echarts@^5.x | TypeScript ^4.x+ | Types available via @types/echarts if using older TS. Modern TS has native support. |
| @vitest/web-worker | Vitest ^1.x or ^2.x | Requires Node 17+ for native structuredClone. Polyfill included for Node 16. WASM testing has known issues (see GitHub #4283, #6118). |

## Confidence Assessment

| Technology | Confidence | Source | Notes |
|------------|------------|--------|-------|
| RDKit.js | **HIGH** | npm package metadata, official GitHub releases | Version 2025.3.4-1.0.0 verified via npm-compare and libraries.io. Published 7 months ago (Jul 2025). Active maintenance. |
| Apache ECharts | **HIGH** | Multiple 2026 comparison articles, official docs | Consistently ranked top 3 for real-time data visualization. Canvas performance claims verified across sources. |
| Vite | **HIGH** | Official documentation, 2026 production setup guides | De facto standard for TypeScript browser apps in 2026. Multiple recent guides confirm continued dominance. |
| Comlink | **MEDIUM** | Google Chrome Labs GitHub, LogRocket articles | Actively maintained but niche library. Recent integrations (Next.js 15, 2025) confirm ongoing relevance. |
| Alpine.js | **MEDIUM** | Framework comparison articles 2026 | Lightweight category recommendation consistent across sources. Not as mainstream as React but established for progressive enhancement use cases. |
| Vitest Web Workers | **LOW** | GitHub issues, community discussions | `@vitest/web-worker` package exists but WASM compatibility has open issues. Recommend manual testing in browser for RDKit.js integration. |

## Sources

**RDKit.js:**
- [GitHub - rdkit/rdkit-js](https://github.com/rdkit/rdkit-js) — Official repository
- [@rdkit/rdkit on npm](https://libraries.io/npm/@rdkit%2Frdkit) — Version verification (2025.3.4-1.0.0)
- [RDKit.js Official Site](https://www.rdkitjs.com/) — Usage patterns (WebFetch blocked)

**Charting Libraries:**
- [JavaScript Chart Libraries in 2026 - Luzmo](https://www.luzmo.com/blog/javascript-chart-libraries) — ECharts performance analysis
- [6 Best JavaScript Charting Libraries for Dashboards in 2026](https://embeddable.com/blog/javascript-charting-libraries) — Real-time comparison
- [ApexCharts.js Comparison Table](https://apexcharts.com/javascript-charts-comparison/) — Feature matrix
- [npm-compare: chart.js vs echarts vs apexcharts](https://npm-compare.com/chart.js,echarts,apexcharts) — Package statistics

**Frontend Frameworks:**
- [Top 10 Best Front End Frameworks in 2026 - Imaginary Cloud](https://www.imaginarycloud.com/blog/best-frontend-frameworks) — Alpine.js positioning
- [Alpine.js: The Minimalist JavaScript Framework - Medium](https://medium.com/@zulfikarditya/alpine-js-the-minimalist-javascript-framework-for-modern-web-development-839382997988) — Size and approach
- [React vs Alpine.js Comparison](https://formations-developpeur-frontend.itgalaxy.io/react-vs-alpine-comparison/index.html) — Use case analysis

**Vite and TypeScript:**
- [Complete Guide to Setting Up React with TypeScript and Vite (2026) - Medium](https://medium.com/@robinviktorsson/complete-guide-to-setting-up-react-with-typescript-and-vite-2025-468f6556aaf2) — Best practices
- [Production-Ready React Project with TypeScript and Vite](https://oneuptime.com/blog/post/2026-01-08-react-typescript-vite-production-setup/view) — Production patterns
- [Vite Official Documentation](https://vite.dev/guide/) — Official source

**Web Workers:**
- [Using Web Workers - MDN](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Using_web_workers) — Official W3C API docs
- [GitHub - GoogleChromeLabs/comlink](https://github.com/GoogleChromeLabs/comlink) — Official Comlink repository
- [Comlink and Web Workers - LogRocket](https://blog.logrocket.com/comlink-web-workers-match-made-in-heaven/) — Integration patterns
- [Web Workers, Comlink, Vite - johnnyreilly](https://johnnyreilly.com/web-workers-comlink-vite-tanstack-query) — Vite setup

**Testing:**
- [@vitest/web-worker on npm](https://www.npmjs.com/package/@vitest/web-worker) — Official package
- [GitHub vitest/vitest #4283](https://github.com/vitest-dev/vitest/discussions/4283) — WASM testing issues
- [GitHub vitest/vitest #6118](https://github.com/vitest-dev/vitest/issues/6118) — Worker WASM path resolution

**WebAssembly:**
- [The State of WebAssembly – 2025 and 2026](https://platform.uno/blog/the-state-of-webassembly-2025-2026/) — Current state
- [Implementing WebAssembly for High-Performance Web Apps - Better Stack](https://betterstack.com/community/guides/scaling-nodejs/webassembly-web-apps/) — Best practices
- [vite-plugin-wasm on npm](https://www.npmjs.com/package/vite-plugin-wasm) — Vite WASM integration

---
*Stack research for: Browser-based cheminformatics SA visualization*
*Researched: 2026-02-14*
