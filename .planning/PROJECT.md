# WebFaulon

## What This Is

A web application implementing Faulon's 1996 simulated annealing algorithm for searching constitutional isomer space (Faulon, J. Chem. Inf. Comput. Sci. 1996, 36, 731-740). Users enter a molecular formula, configure SA parameters, and watch the algorithm explore isomer space in real time with live optimization charts and 2D structure rendering. Built on a FastAPI/Python RDKit backend with a lightweight web frontend.

## Core Value

Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.

## Current Milestone: v2.0 Python Backend Architecture

**Goal:** Re-architect from browser-only to FastAPI backend with Python RDKit, enabling full cheminformatics capabilities and a multi-component target function framework.

**Target features:**
- FastAPI backend with Python RDKit for all molecular operations
- Faulon SA algorithm ported to Python using native RDKit molecules
- Multi-component target function framework (pluggable, weighted) — Wiener Index as first component
- SSE streaming for live SA progress to frontend
- Backend SVG rendering (eliminates RDKit.js WASM dependency)
- Frontend becomes thin visualization layer (Chart.js + Alpine.js)

## Requirements

### Validated

- User can enter any molecular formula to define the constitutional space — v1.0
- SA algorithm implements Faulon's displacement (eqs 7-11) preserving valences — v1.0
- Initial structure generated deterministically from molecular formula — v1.0
- Wiener Index computed as cost function — v1.0
- User can toggle between maximizing and minimizing objective — v1.0
- SA accepts/rejects moves using Metropolis criterion — v1.0
- User can control SA parameters (kT, cooling schedule, steps, cycles) — v1.0
- Live-updating chart shows objective vs step number — v1.0
- Current best isomer rendered as 2D chemical structure — v1.0
- Best score displayed alongside the structure — v1.0
- UI is clean and suitable for classroom projection — v1.0

### Active

- [ ] Python RDKit backend handles all molecular operations natively
- [ ] Faulon displacement algorithm operates on RDKit molecules directly (no custom MolGraph)
- [ ] Multi-component target function framework with pluggable components and adjustable weights
- [ ] Wiener Index implemented as first target function component
- [ ] FastAPI serves REST endpoints for SA configuration and control
- [ ] SSE (Server-Sent Events) streams live SA progress to frontend
- [ ] Backend generates 2D SVG depictions via Python RDKit
- [ ] Frontend displays backend-rendered SVGs (no RDKit.js WASM)
- [ ] Same UX as v1 (formula input, presets, parameter controls, start/pause/reset, live chart, structure display)

### Out of Scope

- Spectroscopic target functions (1D NMR, HMBC) — future milestone
- NMR prediction engine integration — future milestone
- Native mobile app — web-only
- 3D conformational visualization — 2D structure only
- Database of known structures / isomorphism checking
- Multi-user / collaboration features

## Context

- Based on: Faulon, J.-L. "Stochastic Generator of Chemical Structure. 2. Using Simulated Annealing To Search the Space of Constitutional Isomers." J. Chem. Inf. Comput. Sci. 1996, 36, 731-740.
- The SA displacement (Table 2 in paper) works by selecting 4 atoms and redistributing bond orders between the 4 bonds connecting them, subject to valence constraints (eqs 1-11).
- Wiener Index = half the sum of all shortest-path distances between pairs of atoms in the molecular graph.
- The paper's best annealing schedule was f8 with initial kT=100, using 500 steps per cycle and 4-8 cycles.
- Cooling schedules: f_k: kT_t = kT_0 - k * kT_0 * t / delta_t (clamped to 0).
- v1 was browser-only (TypeScript, Web Worker, RDKit.js WASM). v2 moves computation to Python backend.
- Python RDKit provides full cheminformatics: molecular operations, descriptor computation, 2D depiction, and extensibility for future spectroscopic scoring.
- Target audience: chemistry students and researchers.
- The paper (PDF) is in `docs/` for reference.

## Constraints

- **Backend**: Python 3.11+ with FastAPI, RDKit installed via conda/pip
- **Frontend**: Vite + TypeScript + Alpine.js + Chart.js (keep existing frontend tooling)
- **Communication**: REST API for commands, SSE for streaming progress
- **Rendering**: Python RDKit generates SVG, frontend displays it
- **Deployment**: Docker or similar for backend; frontend still buildable as static assets
- **Browser**: Modern browsers (no WASM dependency for RDKit anymore)

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| RDKit.js for 2D rendering (v1) | Full cheminformatics in browser, proper depiction | v1 validated |
| In-browser only, no backend (v1) | Simplest deployment for classroom use | v1 validated |
| Wiener Index as sole cost function (v1) | Matches paper's primary test case | v1 validated |
| Modern SPA with Vite (v1) | Good DX, tree-shaking | v1 validated |
| Switch to Python RDKit backend (v2) | Full cheminformatics access, enables spectroscopic scoring in future | — Pending |
| FastAPI for backend (v2) | Async, fast, auto-generated OpenAPI docs, SSE support | — Pending |
| SSE for live SA updates (v2) | One-way stream, simpler than WebSocket, perfect for progress updates | — Pending |
| Backend SVG rendering (v2) | Eliminates RDKit.js WASM from frontend, single source of truth for molecule handling | — Pending |
| Multi-component target function (v2) | Extensible architecture for future spectroscopic components | — Pending |

---
*Last updated: 2026-02-15 after milestone v2.0 started*
