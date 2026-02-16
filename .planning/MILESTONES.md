# Milestones: WebFaulon

## v1.0: Interactive SA Demo (Browser-Only)

**Completed:** 2026-02-15
**Phases:** 4 (1, 2, 3, 3.1) | **Plans:** 12 | **Duration:** ~2.5 hours

**Delivered:**
- Core SA algorithm with Faulon displacement (eqs 7-11), Metropolis criterion, configurable cooling schedules (f0-f32)
- MolGraph data structure with Wiener Index computation, valence validation, connectivity checks
- Web Worker isolation for non-blocking SA execution
- Alpine.js UI with formula input/validation, 5 preset molecules, SA parameter controls, start/pause/reset
- Live Chart.js visualization (Wiener Index vs step number)
- RDKit.js WASM 2D molecule rendering from MOL blocks
- Responsive, classroom-ready layout (desktop + mobile)
- GitHub Pages deployment at https://steinbeck.github.io/webfaulon/
- Comprehensive README with academic citation

**Architecture:** Pure frontend (TypeScript/Vite, Web Worker, RDKit.js WASM, Alpine.js, Chart.js)

**Key technical decisions:**
- MOL block pipeline (not SMILES) for reliable molecule representation
- RDKit.js set_new_coords() required before draw_to_canvas() for zero-coord MOL blocks
- unpkg CDN for RDKit.js in production builds

**Last phase number:** 3.1

---
*Archive created: 2026-02-15*

## v2.0: Python Backend Architecture

**Completed:** 2026-02-16
**Phases:** 6 (4, 5, 6, 6.1, 7, 8) | **Plans:** 15 | **Duration:** ~9 hours (single day)
**Commits:** 52 | **Python:** 2,428 LOC app + 3,245 LOC tests | **238 tests passing**

**Delivered:**
- FastAPI backend with Python RDKit for all molecular operations (validation, displacement, scoring, SVG depiction)
- Faulon SA engine ported to Python using native RDKit RWMol molecules
- REST API (configure, control, status) + SSE streaming for live SA progress
- Session management with in-memory state and TTL expiration
- Multi-component target function framework with pluggable scoring (Wiener Index + LogP)
- Frontend rewired from Web Worker to backend API + SSE (no more RDKit.js WASM in browser)
- Unsaturated molecule support (benzene, naphthalene) via displacement fix
- Local deployment with CORS configuration for production

**Architecture:** FastAPI/Python RDKit backend + Vite/Alpine.js/Chart.js thin frontend

**Key technical decisions:**
- No SanitizeMol in set_bond() -- Faulon displacement preserves valence by construction
- AROMATIC bond type mapped to 1 as safety fallback for integer bond order compatibility
- No-op displacement retry (MAX_DISPLACEMENT_ATTEMPTS=10) for sparse graph efficiency
- Disconnected moves discarded and counted separately (not persisted)
- LogP replaces MW as second scoring component (MW constant across isomers)
- Protocol-based polymorphism for pluggable scoring components
- SSE generator runs SA steps inline (not BackgroundTasks) for pause/resume control
- Local deployment (not cloud PaaS) -- backend runs on user's machine

**Last phase number:** 8

---
*Archive created: 2026-02-16*

