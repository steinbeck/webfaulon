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
