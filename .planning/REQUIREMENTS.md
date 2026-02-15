# Requirements: WebFaulon

**Defined:** 2026-02-14
**Core Value:** Students can see and interact with the SA algorithm exploring constitutional isomer space in real time — making the abstract algorithm from the paper tangible and intuitive.

## v1 Requirements

Requirements for initial release. Each maps to roadmap phases.

### Input

- [x] **INP-01**: User can enter any molecular formula (e.g., C6H14, C8H10) to define the constitutional isomer space
- [x] **INP-02**: Formula is validated before SA starts (invalid formulas rejected with clear error message)
- [x] **INP-03**: Preset example molecules (3-5 curated) available for quick start without needing to know valid formulas

### Algorithm

- [x] **ALG-01**: SA implements Faulon's displacement: pick 4 atoms (x1, y1, x2, y2), redistribute bond orders per equations 7-11 while preserving valences
- [x] **ALG-02**: Initial structure generated deterministically from the molecular formula (valid connectivity and valences)
- [x] **ALG-03**: Wiener Index computed as cost function (half-sum of all shortest path distances between atom pairs in heavy-atom graph)
- [x] **ALG-04**: User can toggle between maximizing and minimizing the Wiener Index
- [x] **ALG-05**: SA accepts/rejects moves using Metropolis criterion: accept if delta_e < 0 (improving), else accept with probability exp(-delta_e/kT)
- [x] **ALG-06**: Graph connectivity validated after every SA move (disconnected fragments rejected, molecule stays connected)
- [x] **ALG-07**: Algorithm runs entirely in-browser via Web Worker (no backend required, UI stays responsive during computation)

### Controls

- [x] **CTRL-01**: User can set initial temperature (kT) with sensible default (kT=100 per paper)
- [x] **CTRL-02**: User can select cooling schedule from f0 through f32 family as described in paper (kT_t = kT_0 - k * kT_0 * t / delta_t, clamped to 0)
- [x] **CTRL-03**: User can set number of steps per cycle (default 500 per paper)
- [x] **CTRL-04**: User can set number of cycles (default 4-8 per paper)
- [x] **CTRL-05**: User can play, pause, and reset the SA execution

### Visualization

- [ ] **VIZ-01**: Live-updating chart shows Wiener Index vs step number as the algorithm runs (like Figure 4 in paper)
- [ ] **VIZ-02**: Current best isomer rendered as 2D chemical structure using RDKit.js (WASM)
- [ ] **VIZ-03**: Best Wiener Index score displayed alongside the structure
- [ ] **VIZ-04**: Real-time algorithm state display (current step, temperature, current/best Wiener Index)

### User Experience

- [ ] **UX-01**: Clean design suitable for classroom projection/demonstration (readable fonts, good contrast)
- [ ] **UX-02**: Mobile responsive layout (works on tablets and phones for BYOD classroom use)

## v2 Requirements

Deferred to future release. Tracked but not in current roadmap.

### Enhanced Interaction

- **ENHI-01**: Step-by-step mode allowing manual advancement one SA iteration at a time
- **ENHI-02**: Detailed acceptance log showing each proposed move, energy delta, and acceptance decision with probability
- **ENHI-03**: Download results (structure as PNG/SVG image, optimization data as CSV)
- **ENHI-04**: Visual molecule transitions animating bond changes during displacement operation
- **ENHI-05**: Optimization metrics dashboard (acceptance rate over time, energy delta distribution)

## Out of Scope

Explicitly excluded. Documented to prevent scope creep.

| Feature | Reason |
|---------|--------|
| User accounts/login | Unnecessary complexity for educational demo; students won't create accounts for one-off use |
| 3D molecular visualization | Wiener Index is about graph topology (2D connectivity), not 3D geometry; adds no educational value |
| Database of molecules | Scope creep; 3-5 curated presets sufficient; link to PubChem if students want more |
| Backend/server computation | RDKit.js WASM proves complex cheminformatics runs in browser; no server simplifies deployment |
| Multi-objective optimization | Wiener Index is the sole metric in Faulon 1996 paper; adding objectives dilutes educational focus |
| Collaborative/sharing features | Complex infrastructure for minimal gain in solo learning tool |
| Temperature reheat | Advanced SA technique not in Faulon 1996 paper; adds conceptual complexity for students |
| Comparative runs | High development cost; students can run twice manually to compare |
| Annotated molecular features | High complexity to implement well; requires deep domain expertise |
| Batch processing | Students learn by watching individual runs; automation removes the learning experience |

## Traceability

Which phases cover which requirements. Updated during roadmap creation.

| Requirement | Phase | Status |
|-------------|-------|--------|
| INP-01 | Phase 2 | Done |
| INP-02 | Phase 2 | Done |
| INP-03 | Phase 2 | Done |
| ALG-01 | Phase 1 | Done |
| ALG-02 | Phase 1 | Done |
| ALG-03 | Phase 1 | Done |
| ALG-04 | Phase 1 | Done |
| ALG-05 | Phase 1 | Done |
| ALG-06 | Phase 1 | Done |
| ALG-07 | Phase 2 | Done |
| CTRL-01 | Phase 2 | Done |
| CTRL-02 | Phase 2 | Done |
| CTRL-03 | Phase 2 | Done |
| CTRL-04 | Phase 2 | Done |
| CTRL-05 | Phase 2 | Done |
| VIZ-01 | Phase 3 | Pending |
| VIZ-02 | Phase 3 | Pending |
| VIZ-03 | Phase 3 | Pending |
| VIZ-04 | Phase 3 | Pending |
| UX-01 | Phase 3 | Pending |
| UX-02 | Phase 3 | Pending |

**Coverage:**
- v1 requirements: 21 total
- Mapped to phases: 21
- Unmapped: 0

---
*Requirements defined: 2026-02-14*
*Last updated: 2026-02-15 after phase 2 execution complete*
