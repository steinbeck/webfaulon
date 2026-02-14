# Project Research Summary

**Project:** WebFaulon
**Domain:** Browser-based Cheminformatics with Algorithm Visualization
**Researched:** 2026-02-14
**Confidence:** MEDIUM-HIGH

## Executive Summary

WebFaulon is an educational web application combining two specialized domains: chemistry education tools and algorithm visualization. The research reveals a clear consensus on architecture: Web Workers for computation isolation, WASM-based cheminformatics (RDKit.js), and real-time charting with performance optimization. This is a domain where both chemistry correctness and computational performance are equally critical.

The recommended approach follows established patterns from both domains: client-side computation for educational transparency, component isolation between UI and algorithm, and progressive enhancement starting with core correctness before visualization polish. The primary technical risk is balancing chemical validity constraints (valence rules, graph connectivity) with simulated annealing algorithm implementation complexity. The secondary risk is WASM integration in Web Workers, which has well-documented solutions but requires careful initialization handling.

Success requires implementing comprehensive validation at each SA step (connectivity + valence checking) BEFORE optimizing for performance or visualization. Research shows most failures stem from invalid molecular graph operations, not algorithm convergence issues. The architecture naturally mitigates UI performance risks through Web Worker isolation, but chemical correctness must be validated through extensive unit testing from day one.

## Key Findings

### Recommended Stack

Modern browser-based scientific computation in 2026 centers on Vite + TypeScript + WASM + Web Workers. The stack achieves instant dev experience (esbuild transpilation) while supporting computationally intensive cheminformatics through RDKit.js.

**Core technologies:**
- **Vite 6.x + TypeScript 5.x**: Industry standard for browser apps in 2026, zero-config dev with production-ready optimizations, native WASM plugin support
- **RDKit.js 2025.3.4-1.0.0**: Only production-ready cheminformatics library for browsers, handles SMILES parsing and molecular graph operations via WASM
- **Apache ECharts 5.x**: Canvas rendering handles 10K+ points <100ms, superior to Chart.js for scientific streaming data visualization
- **Comlink 4.x**: Eliminates Web Worker postMessage boilerplate, TypeScript-native with proxy-based RPC (1.1KB)
- **Alpine.js 3.x (or vanilla JS)**: Lightweight (7.1KB) for parameter controls without SPA overhead, educational transparency (students read source)

**Critical dependencies:**
- `vite-plugin-wasm` for WASM module loading
- `vite-plugin-comlink` for seamless worker integration
- Vitest for testing (note: WASM testing has known issues, use browser mode for RDKit.js tests)

**Alternatives rejected:** Chart.js (canvas rendering degrades >1000 points), React/Preact (unnecessary overhead for parameter forms), CDN imports for RDKit.js (CORS and path resolution issues).

### Expected Features

Research across algorithm visualization (VisuAlgo) and chemistry education (MolView, WebMO) reveals clear expectations. Users need both algorithm control (play/pause/step-by-step) and chemistry-specific visualization (2D structures, formula input).

**Must have (table stakes):**
- Molecular formula input with validation — standard chemistry tool starting point
- 2D structure visualization of current best — core educational requirement
- Play/Pause/Reset controls — algorithm visualization standard
- Real-time algorithm state (temperature, step, energy) — users need "what is algorithm doing now"
- Real-time optimization chart (Wiener Index vs step) — makes optimization visible not abstract
- SA parameter controls (kT, cooling, max steps) — students learn by adjusting parameters
- Preset example molecules — reduces barrier to entry (3-5 curated examples)
- Mobile responsive design — 2026 standard, students use personal devices

**Should have (competitive differentiators):**
- Step-by-step mode (manual advancement) — high educational value, aids comprehension
- Detailed acceptance log — shows probabilistic decisions, demystifies SA
- Configurable cooling schedules — lets students experiment with schedule impact
- Download results (structure image + optimization CSV) — enables use in assignments
- Visual molecule transitions — animate bond changes showing displacement operation

**Defer to v2+:**
- Temperature reheat (advanced SA technique, not in Faulon 1996 paper)
- Comparative runs (side-by-side parameter comparison, high development cost)
- Annotated molecular features (connects Wiener Index to topology, complex to implement well)

**Anti-features (explicitly avoid):**
- User accounts/login (adds complexity for one-off educational demo)
- 3D visualization (Wiener Index is about 2D graph topology, not 3D geometry)
- Database of molecules (scope creep, provide 5-10 curated presets instead)
- Backend computation (client-side for educational transparency)
- Multi-objective optimization (Faulon focuses on Wiener Index)

### Architecture Approach

The architecture follows a clear three-layer pattern: UI Layer (main thread), Communication Layer (message batching + rAF throttling), and Computation Layer (Web Worker). This separation prevents UI blocking during expensive SA iterations while maintaining responsive controls.

**Major components:**
1. **Web Worker (Computation Layer)** — Runs SA algorithm loop, RDKit.js WASM, molecular graph operations (MolGraph class with adjacency matrix, valence checking, connectivity via BFS, Wiener Index calculation)
2. **Main Thread (UI Layer)** — Chart rendering (ECharts), UI controls, status display, molecular structure display (from worker results)
3. **Communication Infrastructure** — WorkerManager abstracts worker lifecycle, MessageBatcher + requestAnimationFrame throttling (batch worker messages to 60fps updates), transferable objects for large adjacency matrices (zero-copy transfer)
4. **Molecular Graph Engine (MolGraph.ts)** — Core data structure with adjacency matrix operations, valence validation, connectivity checking (BFS), Wiener Index via BFS from each node (O(n²) for sparse graphs)
5. **SA Algorithm Controller** — Manages SA state machine (IDLE → INITIALIZING → RUNNING → PAUSED → STOPPED → ERROR), temperature scheduling, acceptance criteria

**Critical patterns identified:**
- **Async WASM initialization in worker** — RDKit.js must initialize INSIDE worker thread (cannot pass WASM across postMessage due to structured clone limitations)
- **Message batching + rAF throttling** — Worker may generate 1000s messages/sec, batch and throttle to 60fps for UI updates
- **BFS-based Wiener Index** — O(n × (n+m)) = O(n²) for sparse molecular graphs, NOT Floyd-Warshall O(n³)
- **Connectivity validation after every move** — BFS from arbitrary node to verify graph remains connected (reject moves creating disconnected fragments)
- **Chart.js decimation** — Enable LTTB algorithm, limit visible points to 500-1000, disable animations for performance

### Critical Pitfalls

Research reveals that chemical correctness failures dominate algorithm convergence issues. Most pitfalls stem from invalid graph operations or WASM integration mistakes, not SA parameter tuning.

1. **Incorrect bond order redistribution** — Faulon's displacement equations (eqs 7-11) must preserve valences; implementation errors produce chemically impossible molecules with invalid valences (carbon with 5 bonds). **Prevention:** Comprehensive unit tests for every move type, verify detailed balance (if A→B accepted, B→A must be possible), validate connectivity + valence + total bond order after each move.

2. **Graph connectivity loss after bond moves** — Removing single bond connecting fragments creates disconnected structures. **Prevention:** ALWAYS BFS connectivity check after every SA move, reject moves creating fragments, O(n+m) is fast enough.

3. **WASM initialization and memory in workers** — RDKit.js WASM cannot be imported ES6-style or passed via postMessage (DataCloneError). **Prevention:** Initialize RDKit.js INSIDE each worker via importScripts(), copy .wasm to public assets, pass serializable representations (adjacency lists, SMILES) between threads.

4. **Premature convergence from wrong cooling schedule** — Faulon describes 33 schedules (f0-f32); inappropriate choice causes local minimum trapping or excessive computation. **Prevention:** Implement multiple schedules, monitor acceptance ratio (should start ~80%, end ~5%), adaptive cooling if ratio drops too fast.

5. **Wiener Index performance bottleneck** — Naive Floyd-Warshall O(n³) freezes browser for >30 atoms. **Prevention:** Use BFS from each node O(n²) for sparse molecular graphs, cache when move rejected, profile separately (should be <10ms per call).

6. **Invalid initial structure generation** — Converting formula (C6H12) to valid graph non-trivial; random bond assignments violate valences. **Prevention:** Deterministic/stochastic generators, validate connectivity + valences before SA starts, consider SMILES input alternative.

7. **Implicit vs explicit hydrogen inconsistency** — Mixing representations causes valence errors and incorrect Wiener Index. **Prevention:** Choose implicit hydrogens (recommended), document clearly, exclude H from Wiener Index, validate bond_order_sum + implicit_H_count = standard_valence.

8. **Aromaticity and kekulization failures** — Aromatic systems require specific bond orders; mixing Kekule and aromatic forms causes RDKit.js errors. **Prevention:** Use Kekule structures (explicit double bonds) for constitutional isomer generation, let RDKit.js detect aromaticity AFTER generation.

## Implications for Roadmap

Based on research consensus, implementation must follow dependency order: molecular graph foundation → SA algorithm correctness → worker integration → visualization → polish. Attempting visualization before validating chemical correctness leads to hard-to-debug failures.

### Suggested Phase Structure

#### Phase 1: Molecular Graph Foundation & SA Core
**Rationale:** Chemistry correctness is the foundation. SA algorithm depends on valid graph operations. All downstream work (visualization, UI) assumes chemically valid structures. Research shows most failures stem from invalid graph operations, not algorithm issues.

**Delivers:**
- MolGraph.ts class (adjacency matrix, valence checking, connectivity via BFS)
- Wiener Index calculation (BFS-based, O(n²) for sparse graphs)
- Initial structure generation from molecular formula (with validation)
- SA algorithm core (temperature scheduling, acceptance criteria, state machine)
- Bond order displacement function (Faulon equations 7-11)
- Comprehensive unit test suite (100+ move validations, edge cases)

**Addresses from FEATURES.md:**
- Foundation for molecular formula input
- Foundation for SA parameter controls
- Core algorithm required for all visualization

**Avoids from PITFALLS.md:**
- Pitfall 1 (invalid bond redistribution) — extensive test suite
- Pitfall 2 (connectivity loss) — BFS check after every move
- Pitfall 5 (Wiener Index bottleneck) — BFS approach from start
- Pitfall 6 (invalid initial structure) — validated generation
- Pitfall 7 (hydrogen handling) — implicit H chosen at architecture level

**Research needed:** LOW — Graph algorithms and SA theory well-documented. Standard patterns.

**Verification:** Test suite passes, 1000-step SA runs produce no disconnected structures, Wiener calculation <5ms for 50-atom molecules.

---

#### Phase 2: Web Worker Integration & WASM
**Rationale:** Isolate computation from UI to prevent blocking. Must happen before UI work to establish communication patterns. WASM initialization in workers is well-documented but must be done carefully.

**Delivers:**
- Web Worker setup (worker entry point, message handling)
- RDKit.js WASM initialization in worker (async pattern, error handling)
- WorkerManager abstraction (lifecycle management, clean API for UI)
- MessageBatcher + rAF throttling (batch to 60fps)
- Worker state machine implementation
- Transferable object patterns for large adjacency matrices

**Addresses from FEATURES.md:**
- Non-blocking UI during SA execution
- Foundation for real-time state display
- Play/pause/reset controls infrastructure

**Avoids from PITFALLS.md:**
- Pitfall 3 (WASM initialization) — proper async init pattern in worker
- Anti-pattern 1 (running SA in main thread) — worker isolation
- Anti-pattern 2 (message flooding) — batching + throttling
- Anti-pattern 5 (blocking worker) — interruptible loop with control messages

**Research needed:** LOW — Web Worker + WASM integration has official MDN docs and established patterns. RDKit.js worker examples exist.

**Verification:** WASM loads in worker, error handling tested, 10K-step SA remains responsive to stop/pause commands.

---

#### Phase 3: Basic UI & Real-time Visualization
**Rationale:** With computation stable and isolated, add user-facing features. Chart performance critical (10K+ points), requires ECharts configuration and optimization from start.

**Delivers:**
- Minimal UI shell (HTML/CSS with controls)
- Parameter controls (kT initial, cooling rate, max steps with sliders)
- Play/pause/reset buttons
- Real-time algorithm state display (step, temperature, current/best Wiener Index)
- ECharts integration (with decimation, no animations, rAF throttling)
- Real-time optimization chart (Wiener Index vs step)
- Preset examples (3-5 molecules with pre-configured parameters)
- Mobile responsive layout

**Addresses from FEATURES.md:**
- Molecular formula text input
- Formula validation
- SA parameter controls (P1 table stakes)
- Play/pause/reset controls (P1 table stakes)
- Real-time optimization chart (P1 table stakes)
- Real-time state display (P1 table stakes)
- Preset examples (P1 table stakes)
- Mobile responsive (P1 table stakes)

**Avoids from PITFALLS.md:**
- Anti-pattern 4 (no decimation) — ECharts decimation configured from start
- UX pitfall (no progress indication) — real-time state updates
- UX pitfall (no cancel) — pause/stop controls

**Research needed:** LOW — ECharts configuration well-documented. Alpine.js patterns standard.

**Verification:** Smooth 60fps chart updates with 10K points, mobile responsive on tablets/phones, parameters affect SA behavior correctly.

---

#### Phase 4: 2D Structure Visualization
**Rationale:** Depends on RDKit.js integration (Phase 2) and core SA (Phase 1). Chemistry visualization completes the educational loop but isn't required for algorithm validation.

**Delivers:**
- 2D molecular structure rendering (RDKit.js native rendering)
- Display of current best structure (updates as SA improves)
- Conversion from internal graph to RDKit.js molecule
- Implicit hydrogen handling in rendering
- Kekulization error handling

**Addresses from FEATURES.md:**
- 2D structure visualization (P1 table stakes)
- Visual completion of optimization loop

**Avoids from PITFALLS.md:**
- Pitfall 7 (hydrogen handling) — implicit H set correctly for RDKit.js
- Pitfall 8 (kekulization failures) — error handling for aromatic structures

**Research needed:** MEDIUM — RDKit.js rendering API needs exploration. Canvas integration patterns.

**Verification:** Generated structures render correctly, aromatic test molecules (benzene, pyridine) convert successfully.

---

#### Phase 5: Enhanced Interaction & Export
**Rationale:** Polish features after core functionality validated. High educational value (step-by-step, acceptance log) but not critical path.

**Delivers:**
- Step-by-step mode (manual advancement one iteration at a time)
- Detailed acceptance log (proposed move, energy delta, acceptance decision)
- Download results (structure as PNG/SVG, optimization data as CSV)
- Configurable cooling schedules (exponential, linear, logarithmic)
- Optimization metrics dashboard (acceptance rate, diversity metrics)

**Addresses from FEATURES.md:**
- Step-by-step mode (P2 feature)
- Detailed acceptance log (P2 feature)
- Download results (P2 feature)
- Configurable cooling schedules (P2 feature)

**Research needed:** LOW — Standard patterns for file export and UI enhancements.

**Verification:** Step-by-step mode aids comprehension (user testing), export formats work in assignments.

---

### Phase Ordering Rationale

**Dependency-driven sequencing:**
- MolGraph must exist before SA algorithm (Phase 1 internal dependency)
- SA algorithm must work before worker integration (Phase 2 depends on Phase 1 correctness)
- Worker must work before UI (Phase 3 depends on Phase 2 communication layer)
- 2D rendering depends on RDKit.js in worker (Phase 4 depends on Phase 2 WASM setup)

**Risk mitigation through early validation:**
- Phase 1 addresses 5 of 8 critical pitfalls through comprehensive testing BEFORE integration
- Phase 2 addresses WASM integration complexity with clear error handling
- Phase 3 addresses performance (ECharts decimation) from initial configuration
- Phase 4 handles chemistry visualization edge cases (aromaticity, kekulization)

**Educational value delivery:**
- Phase 1: Core algorithm works (can log SA iterations to console for validation)
- Phase 2: Non-blocking execution (can see algorithm run without freezing browser)
- Phase 3: Visual feedback (students see optimization happening in real-time)
- Phase 4: Chemistry connection (students see molecular structures, not just numbers)
- Phase 5: Deep understanding (step-through and logs demystify SA decisions)

**Parallel work opportunities:**
- Phase 1 MolGraph can develop parallel to SA algorithm skeleton
- Phase 3 UI shell can prototype parallel to Phase 2 worker integration
- Phase 4 RDKit.js rendering research can happen during Phase 3 implementation

### Research Flags

**Phases needing deeper research during planning:**
- **Phase 4 (2D Structure Visualization):** RDKit.js rendering API has limited documentation, canvas integration patterns need exploration, kekulization edge cases require investigation
- **Phase 5 (Export functionality):** File format preferences (PNG vs SVG, CSV vs JSON) should survey chemistry instructors

**Phases with standard patterns (skip research-phase):**
- **Phase 1 (Molecular Graph & SA):** Graph algorithms textbook material, SA theory well-documented, numerous implementations exist
- **Phase 2 (Web Workers):** Official MDN documentation comprehensive, Comlink patterns established, WASM initialization examples plentiful
- **Phase 3 (UI & Charts):** ECharts documentation extensive, Alpine.js patterns standard, mobile responsive design established
- **Phase 5 (Enhanced Interaction):** Step-by-step UI patterns common in algorithm visualizations, file export browser APIs well-documented

## Confidence Assessment

| Area | Confidence | Notes |
|------|------------|-------|
| Stack | **HIGH** | Vite, RDKit.js, ECharts verified via npm metadata and official docs. Comlink from Google Chrome Labs. Active maintenance confirmed. |
| Features | **MEDIUM-HIGH** | Algorithm visualization patterns (VisuAlgo) and chemistry education tools (MolView, WebMO) show consistent expectations. SA-specific feature precedent from TSP demos. |
| Architecture | **HIGH** | Web Worker + WASM patterns extensively documented (MDN). Performance optimization patterns (rAF throttling, batching) established. Component boundaries clear. |
| Pitfalls | **MEDIUM** | WASM/Worker issues well-documented (MDN, GitHub issues). SA theory high confidence (academic consensus). Faulon-specific implementation LOW confidence (1996 paper not directly accessed, inferred from related work). |

**Overall confidence:** MEDIUM-HIGH

**High confidence areas:**
- Web Worker + WASM integration (official documentation, established patterns)
- Performance optimization (rAF throttling, decimation, BFS vs Floyd-Warshall)
- Technology choices (npm package verification, version compatibility)

**Medium confidence areas:**
- Feature prioritization (inferred from similar domains, needs user testing validation)
- Cooling schedule specifics (Faulon paper describes 33 schedules, transfer to JavaScript context may need tuning)
- Educational value of specific features (step-by-step, acceptance log assumed valuable from VisuAlgo research)

**Low confidence areas:**
- Faulon displacement equations (eqs 7-11) implementation details (1996 paper not directly accessed)
- Optimal cooling schedule parameters for browser context (may differ from 1996 Fortran implementation)

### Gaps to Address

**During Phase 1 planning:**
- **Faulon displacement equations:** Access original 1996 paper for equations 7-11 details, or reverse-engineer from reference implementations
- **Initial structure generation algorithm:** Choose between deterministic (MAYGEN-style) vs stochastic approach, benchmark validity rate

**During Phase 3 planning:**
- **Chart decimation threshold:** Determine optimal sample count (500 vs 1000 vs 2000) through performance testing on target devices
- **Mobile breakpoints:** Test preset parameters on actual classroom tablets/Chromebooks to validate performance

**During Phase 4 planning:**
- **RDKit.js rendering canvas size:** Determine optimal resolution for clarity vs performance
- **Kekulization failure UX:** Decide error handling strategy (silent fallback vs user notification)

**During Phase 5 planning:**
- **Export formats:** Survey chemistry instructors for PNG vs SVG preference, CSV vs JSON for data
- **Step-by-step performance:** Determine if manual stepping needs different optimization than continuous run

**Validation recommendations from research:**
- Conduct user testing with chemistry students to validate P1/P2 priority splits
- Benchmark cooling schedules against known test cases (all C6H12 isomers generation)
- Test visual molecule transitions with target users to assess if complexity justified by educational value

## Sources

### Primary (HIGH confidence)

**Official Documentation:**
- MDN Web Docs — Web Workers API, WASM integration, structured clone algorithm
- Vite Official Documentation — Build configuration, WASM plugins
- RDKit.js GitHub Repository — WASM usage patterns, API reference
- npm package metadata — Version verification (@rdkit/rdkit 2025.3.4-1.0.0, libraries.io)
- Apache ECharts Official Docs — Performance optimization, decimation configuration

**Academic/Authoritative:**
- Wikipedia: Simulated Annealing — Theory and cooling schedules
- Wikipedia: Wiener Index — Definition and graph theory foundation
- ACM/Springer — Algorithm visualization educational research
- PMC (PubMed Central) — Chemistry education pedagogy

### Secondary (MEDIUM confidence)

**Technical Articles (2026):**
- Luzmo, Embeddable — JavaScript chart library comparisons (2026)
- LogRocket, johnnyreilly — Comlink integration patterns
- Medium, Better Stack — WASM best practices and performance
- Dev3lop — WebGL vs Canvas rendering benchmarks

**Community Resources:**
- GitHub issues (Vitest #4283, #6118) — WASM testing known issues
- VisuAlgo.net — Algorithm visualization patterns (research-backed)
- MolView, WebMO — Chemistry education tool feature analysis
- TSP SA demos (vgarciasc, Inspiaaa) — SA visualization precedent

**Chemistry Domain:**
- Wiley Online Library — Molecular graphics roadmap, computational chemistry education
- ACS Journal of Chemical Education — Parameter tuning simulations (ICP-MS TuneSim)
- NFDI4Chem, ChemAxon — Cheminformatics access and tool patterns

### Tertiary (LOW confidence, needs validation)

**Faulon Algorithm Details:**
- Faulon Research Site (jfaulon.com) — Chemical structure generation background (1996 paper not directly accessed)
- Related isomer generators (MAYGEN, Surge) — Inferred patterns from similar approaches
- NextMove, Depth-First blog — Hydrogen handling and aromaticity (expert opinions, not original research)

**Inferred from Analogous Domains:**
- TSP SA cooling schedules — Transferred to molecular optimization context (may need tuning)
- Chart.js vs ECharts performance claims — Vendor comparisons, needs independent validation
- Alpine.js educational value — Assumed from progressive enhancement philosophy

---

*Research completed: 2026-02-14*
*Ready for roadmap: YES*

## Next Steps for Orchestrator

This synthesis is complete and ready for roadmap creation. The suggested five-phase structure provides clear dependency ordering and addresses all critical pitfalls from research.

**Recommended approach for roadmap:**
1. Adopt suggested phase structure as starting point
2. Skip `/gsd:research-phase` for Phases 1, 2, 3, 5 (standard patterns, well-documented)
3. Consider `/gsd:research-phase` for Phase 4 (RDKit.js rendering API exploration)
4. Flag Phase 1 for comprehensive test suite requirement (chemical correctness critical)
5. Flag Phase 2 for WASM initialization error handling focus

**Key research insights to carry forward:**
- Chemical correctness before optimization before polish
- Validate graph operations (connectivity + valence) after EVERY SA move
- Initialize RDKit.js INSIDE workers, never pass WASM via postMessage
- Use BFS for Wiener Index (O(n²)), not Floyd-Warshall (O(n³))
- ECharts with decimation configured from start, not Chart.js
