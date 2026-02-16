# Roadmap: WebFaulon

## Milestones

- v1.0 Interactive SA Demo (Browser-Only) -- Phases 1-3.1 (shipped 2026-02-15)
- v2.0 Python Backend Architecture -- Phases 4-8 (shipped 2026-02-16)

## Phases

**Phase Numbering:**
- Integer phases (1, 2, 3): Planned milestone work
- Decimal phases (2.1, 2.2): Urgent insertions (marked with INSERTED)

Decimal phases appear between their surrounding integers in numeric order.

<details>
<summary>v1.0 Interactive SA Demo (Phases 1-3.1) -- SHIPPED 2026-02-15</summary>

- [x] **Phase 1: Molecular Graph & SA Core** - Chemical correctness foundation
- [x] **Phase 2: Browser Integration & Controls** - Web Workers, user input, execution controls
- [x] **Phase 3: Visualization & UX** - Charts, structure rendering, classroom-ready design
- [x] **Phase 3.1: README & GitHub Pages** - Documentation and public deployment (INSERTED)

See `.planning/milestones/v1.0-ROADMAP.md` for full details.

</details>

<details>
<summary>v2.0 Python Backend Architecture (Phases 4-8) -- SHIPPED 2026-02-16</summary>

- [x] **Phase 4: Backend Core & RDKit Foundation** (4 plans) - FastAPI + Python SA engine on native RDKit molecules
- [x] **Phase 5: API Layer & SSE Streaming** (3 plans) - REST endpoints, session management, live streaming
- [x] **Phase 6: Frontend Integration** (2 plans) - Rewire frontend from Web Worker to backend API + SSE
- [x] **Phase 6.1: Fix unsaturated displacement** (2 plans) - AROMATIC fallback, no-op retry (INSERTED)
- [x] **Phase 7: Multi-Component Target Function** (2 plans) - Pluggable scoring with Wiener Index + LogP
- [x] **Phase 8: Deployment & Production** (2 plans) - Local deployment with CORS configuration

See `.planning/milestones/v2.0-ROADMAP.md` for full details.

</details>

## Progress

| Phase | Milestone | Plans | Status | Completed |
|-------|-----------|-------|--------|-----------|
| 1. Molecular Graph & SA Core | v1.0 | 4/4 | Complete | 2026-02-14 |
| 2. Browser Integration & Controls | v1.0 | 3/3 | Complete | 2026-02-15 |
| 3. Visualization & UX | v1.0 | 3/3 | Complete | 2026-02-15 |
| 3.1. README & GitHub Pages | v1.0 | 2/2 | Complete | 2026-02-15 |
| 4. Backend Core & RDKit Foundation | v2.0 | 4/4 | Complete | 2026-02-16 |
| 5. API Layer & SSE Streaming | v2.0 | 3/3 | Complete | 2026-02-16 |
| 6. Frontend Integration | v2.0 | 2/2 | Complete | 2026-02-16 |
| 6.1. Fix unsaturated displacement | v2.0 | 2/2 | Complete | 2026-02-16 |
| 7. Multi-Component Target Function | v2.0 | 2/2 | Complete | 2026-02-16 |
| 8. Deployment & Production | v2.0 | 2/2 | Complete | 2026-02-16 |

**Total: 10 phases, 27 plans, all complete across 2 milestones**

---

*Roadmap created: 2026-02-14*
*v1.0 shipped: 2026-02-15*
*v2.0 shipped: 2026-02-16*
