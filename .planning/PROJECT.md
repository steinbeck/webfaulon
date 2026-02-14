# WebFaulon

## What This Is

A web-based interactive demonstration of Faulon's 1996 simulated annealing algorithm for searching the space of constitutional isomers (Faulon, J. Chem. Inf. Comput. Sci. 1996, 36, 731-740). Students enter a molecular formula, configure SA parameters, and watch the algorithm explore constitutional isomer space in real time, optimizing the Wiener Index. The tool visualizes the SA progress as a live chart and renders the current best molecular structure as a 2D chemical drawing.

## Core Value

Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.

## Requirements

### Validated

(None yet — ship to validate)

### Active

- [ ] User can enter any molecular formula to define the constitutional space
- [ ] SA algorithm implements Faulon's displacement: pick 4 atoms (x1, y1, x2, y2), redistribute bond orders per eqs 7-11 while preserving valences
- [ ] Initial structure is generated deterministically from the molecular formula
- [ ] Wiener Index computed as cost function (half-sum of all shortest path distances)
- [ ] User can toggle between maximizing and minimizing the Wiener Index
- [ ] SA accepts/rejects moves using Metropolis criterion: accept if delta_e < 0 (improving), else accept with probability exp(-delta_e/kT)
- [ ] User can control SA parameters: initial kT, cooling schedule function (f0 through f32 family from paper), number of steps per cycle, number of cycles
- [ ] Live-updating chart shows Wiener Index vs step number as the algorithm runs (like Figure 4 in paper)
- [ ] Current best isomer rendered as 2D chemical structure using RDKit.js (WASM)
- [ ] Best score displayed alongside the structure
- [ ] Algorithm runs entirely in-browser (no backend)
- [ ] UI is clean and suitable for classroom projection/demonstration

### Out of Scope

- Native mobile app — web-only
- Other cost functions (log P, potential energy) — Wiener Index only for v1
- 3D conformational visualization — 2D structure only
- Database of known structures / isomorphism checking — not needed per paper's findings
- Multi-user / collaboration features — single-user demo
- Backend / server-side computation — everything in browser

## Context

- Based on: Faulon, J.-L. "Stochastic Generator of Chemical Structure. 2. Using Simulated Annealing To Search the Space of Constitutional Isomers." J. Chem. Inf. Comput. Sci. 1996, 36, 731-740.
- The SA displacement (Table 2 in paper) works by selecting 4 atoms and redistributing bond orders between the 4 bonds connecting them, subject to valence constraints (eqs 1-11).
- Wiener Index = half the sum of all shortest-path distances between pairs of atoms in the molecular graph. For n-paraffins, maximum Wiener = n(n^2-1)/6.
- The paper's best annealing schedule was f8 with initial kT=100, using 500 steps per cycle and 4-8 cycles.
- Cooling schedules: f_k: kT_t = kT_0 - k * kT_0 * t / delta_t (clamped to 0). f0=constant, f1=linear, f8=fast decay.
- RDKit.js provides WASM-based cheminformatics in the browser: structure generation from SMILES, 2D coordinate generation, SVG rendering.
- Target audience: chemistry students learning about computational chemistry and stochastic methods.
- The paper (PDF) is in `docs/` for reference.

## Constraints

- **Runtime**: All computation in browser (JavaScript/WASM) — no server required
- **Chemistry rendering**: RDKit.js for 2D structure depiction
- **Charting**: Lightweight chart library (e.g., Chart.js or similar) for live SA progress
- **Deployment**: Should work as a static site (can be served from any HTTP server or opened locally)
- **Browser**: Modern browsers with WASM support

## Key Decisions

| Decision | Rationale | Outcome |
|----------|-----------|---------|
| RDKit.js for 2D rendering | Full cheminformatics in browser, proper depiction, SMILES support | — Pending |
| In-browser only (no backend) | Simplest deployment for classroom use, no server setup | — Pending |
| Wiener Index as sole cost function for v1 | Matches paper's primary test case, easy to verify correctness | — Pending |
| Modern SPA with build tooling (Vite) | Good DX, tree-shaking, easy to add dependencies | — Pending |

---
*Last updated: 2026-02-14 after initialization*
