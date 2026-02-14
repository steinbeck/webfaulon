# Feature Landscape

**Domain:** Interactive Chemistry Education Web Apps (Simulated Annealing Molecular Optimization)
**Researched:** 2026-02-14
**Confidence:** MEDIUM

## Table Stakes

Features users expect. Missing = product feels incomplete.

| Feature | Why Expected | Complexity | Notes |
|---------|--------------|------------|-------|
| Molecular formula input | Standard starting point for chemistry tools; users expect simple text input (C6H12O6) | Low | MolView and similar tools use this as primary input method. Should validate formula before proceeding. |
| 2D structure visualization | Core output for chemistry education; students need to see molecular structure | Medium | RDKit.js handles this. ChemDoodle and MolView show this is expected in all chemistry web apps. |
| Play/Pause/Reset controls | Standard for all algorithm visualizations; VisuAlgo research shows users expect VCR-like controls | Low | Students need to control pacing for comprehension. |
| Real-time algorithm state display | Shows current temperature, step count, acceptance ratio; algorithm demos always show this | Medium | Users need to see "what is the algorithm doing now?" Similar to SA TSP demos. |
| Preset examples | Reduces barrier to entry; research shows presets help beginners understand parameter relationships | Low | 3-5 example molecules with pre-configured parameters. Students can click and run immediately. |
| Parameter controls for SA | kT (temperature), cooling schedule, max steps are the basic SA parameters users expect to tune | Medium | ICP-MS TuneSim research shows students learn by adjusting instrument parameters. Sliders with reasonable defaults. |
| Mobile responsiveness | 2026 standard; WebMO research shows BYOD paradigm is critical for chemistry education | Medium | Must work on tablets/phones. Students use personal devices in class. |
| Step-by-step mode | Allows manual advancement one step at a time; VisuAlgo research shows this aids comprehension | Medium | Critical for understanding stochastic acceptance decisions. |

## Differentiators

Features that set product apart. Not expected, but valued.

| Feature | Value Proposition | Complexity | Notes |
|---------|-------------------|------------|-------|
| Live optimization chart | Real-time graph showing Wiener Index vs step; makes optimization visible | Medium | Chart.js can handle streaming updates. VisuAlgo research: animation aids understanding. Differentiates from static demos. |
| Visual molecule transitions | Animate bond redistribution when SA accepts a move; shows the displacement operation | High | Visualizing the 4-atom displacement from Faulon paper. Highly educational but complex to implement smoothly. |
| Configurable cooling schedules | Exponential, linear, logarithmic options; lets students experiment with schedule impact | Medium | Most SA demos fix schedule. Allowing experimentation teaches schedule importance. Related to parameter tuning research. |
| Download results | Export final structure as image/SMILES/molfile; enables use in reports/assignments | Low-Medium | WebMO research shows export is valued. Students want to save work. PNG for structure, CSV for optimization data. |
| Comparative runs | Side-by-side comparison of multiple runs with different parameters | High | Shows stochastic nature and parameter sensitivity. Defer to v2 unless critical. |
| Detailed acceptance log | Shows each proposed move, energy delta, acceptance decision with probability | Medium | Educational value: demystifies probabilistic acceptance. Could be collapsible panel. |
| Temperature reheat | Increase temperature when stuck in local minimum; advanced SA feature | Medium | CADApps demo shows this. Teaches advanced optimization strategies. |
| Optimization metrics dashboard | Display current vs best Wiener Index, acceptance rate, diversity metrics | Medium | Helps students understand convergence and exploration/exploitation tradeoff. |
| Annotated molecular features | Highlight branching, rings, chains in structure; connects Wiener Index to topology | High | Deeply educational but complex. Shows WHY certain structures score better. |

## Anti-Features

Features to explicitly NOT build.

| Anti-Feature | Why Avoid | What to Do Instead |
|--------------|-----------|-------------------|
| User accounts/login | Adds complexity and privacy concerns for educational demo; students won't create accounts for one-off demos | Use browser localStorage to save recent runs. Focus on instant usability. |
| 3D molecular visualization | Wiener Index is about graph topology not 3D geometry; 3D doesn't add educational value here and adds significant complexity | Stick to clear 2D drawings. The algorithm operates on constitutional isomers (2D connectivity). |
| Database of molecules | Scope creep; students should focus on algorithm not browsing molecules | Provide 5-10 curated presets that demonstrate different phenomena. Link to PubChem if they want more. |
| Collaborative/sharing features | Complex infrastructure for minimal educational gain in solo learning tool | Export image/data for sharing via existing channels (email, LMS). |
| Backend computation | Overengineering; Faulon SA for small molecules runs fast in browser | Keep it client-side. RDKit.js WASM proves complex cheminformatics runs in browser. No server needed. |
| Social/gamification features | Wrong context; this is serious computational chemistry education not a game | Focus on clear explanations and parameter exploration. Educational value from understanding, not points. |
| Multi-objective optimization | Wiener Index is THE metric in Faulon 1996 paper; adding objectives dilutes focus | Stay true to the paper. If students want more, this teaches the foundation. |
| Batch processing | Students learn by watching individual runs; automation removes learning experience | Single focused run. The point is to WATCH the algorithm work. |

## Feature Dependencies

```
Molecular formula input
    └──requires──> Formula validation
    └──enables──> Initial structure generation
                      └──enables──> SA execution
                                       └──enables──> 2D visualization
                                       └──enables──> Optimization chart

SA execution
    └──requires──> Parameter controls (kT, cooling, steps)
    └──enables──> Real-time state display
    └──enables──> Step-by-step mode

Download results
    └──requires──> 2D visualization (to export image)
    └──requires──> Optimization chart (to export data)

Detailed acceptance log
    └──enhances──> Step-by-step mode
    └──enhances──> Real-time state display

Configurable cooling schedules
    └──enhances──> Parameter controls
    └──enables──> Temperature comparison experiments

Visual molecule transitions
    └──requires──> 2D visualization
    └──enhances──> Real-time algorithm understanding
    └──conflicts──> Fast animation speed (smooth animation needs slower pace)
```

### Dependency Notes

- **Molecular formula input requires validation:** Invalid formulas (C-1H4) must be caught before attempting structure generation
- **SA execution enables visualization:** All visual features depend on algorithm running and generating data
- **Download requires visualization components:** Can't export what isn't rendered
- **Visual transitions conflict with speed:** Smooth animation takes time; users watching transitions won't want 1000 steps/second
- **Preset examples enhance formula input:** Reduces friction for first-time users who don't know valid formulas
- **Step-by-step enhances detailed log:** Seeing acceptance decisions one at a time is more comprehensible than streaming log

## MVP Definition

### Launch With (v1)

Minimum viable product to validate that students can learn SA optimization from this tool.

- [x] **Molecular formula text input** — Core starting point; students enter what they want to explore
- [x] **Formula validation** — Prevents errors; ensures smooth experience
- [x] **Preset example molecules** — Reduces barrier to entry; 3 presets (small, medium, branched) let students click and learn
- [x] **Basic SA parameter controls** — kT initial, cooling rate, max steps with sliders and sensible defaults
- [x] **Play/Pause/Reset controls** — Students control pacing; critical for comprehension
- [x] **2D structure display of current best** — Shows what the algorithm found; core chemistry education requirement
- [x] **Real-time optimization chart** — Wiener Index vs step; makes "optimization" visible not abstract
- [x] **Real-time algorithm state** — Current step, temperature, current/best Wiener Index; shows algorithm internals
- [x] **Mobile responsive layout** — Students use phones/tablets; BYOD paradigm is standard in 2026

### Add After Validation (v1.x)

Features to add once students are using the tool and providing feedback.

- [ ] **Step-by-step mode** — Manual advancement; high educational value (trigger: students request finer control)
- [ ] **Detailed acceptance log** — Shows probabilistic decisions; demystifies SA (trigger: confusion about acceptance)
- [ ] **Configurable cooling schedules** — Exponential, linear, logarithmic; teaches schedule impact (trigger: advanced students want experimentation)
- [ ] **Download results** — Export structure image and optimization CSV; enables assignments (trigger: instructors request this)
- [ ] **Visual molecule transitions** — Animate bond changes; shows displacement operation (trigger: students don't understand how moves work)
- [ ] **Optimization metrics dashboard** — Acceptance rate, energy delta distribution; deeper insight (trigger: advanced use cases)

### Future Consideration (v2+)

Features to defer until product-market fit and user needs are clearer.

- [ ] **Temperature reheat** — Advanced SA technique; complex to explain (why defer: not in Faulon 1996 paper, adds conceptual complexity)
- [ ] **Comparative runs** — Side-by-side parameter comparison; high development cost (why defer: can run twice manually, unclear if valuable enough)
- [ ] **Annotated molecular features** — Highlight topology contributing to Wiener Index; deeply educational (why defer: high complexity, requires domain expertise to implement well)
- [ ] **Custom displacement rules** — Modify the 4-atom redistribution; research tool territory (why defer: changes from Faulon algorithm to generic framework)
- [ ] **Share/permalink URLs** — Encode state in URL for sharing; nice to have (why defer: export accomplishes similar goal, simpler)
- [ ] **Accessibility enhancements** — Screen reader support, keyboard navigation; important but specialized (why defer: validate core concept first, then invest in accessibility)

## Feature Prioritization Matrix

| Feature | User Value | Implementation Cost | Priority |
|---------|------------|---------------------|----------|
| Molecular formula input | HIGH | LOW | P1 |
| 2D structure visualization | HIGH | MEDIUM | P1 |
| Real-time optimization chart | HIGH | MEDIUM | P1 |
| Play/Pause/Reset controls | HIGH | LOW | P1 |
| Preset examples | HIGH | LOW | P1 |
| SA parameter controls | HIGH | MEDIUM | P1 |
| Real-time state display | HIGH | MEDIUM | P1 |
| Mobile responsive | HIGH | MEDIUM | P1 |
| Step-by-step mode | HIGH | MEDIUM | P2 |
| Detailed acceptance log | MEDIUM | MEDIUM | P2 |
| Download results | MEDIUM | LOW | P2 |
| Configurable cooling schedules | MEDIUM | MEDIUM | P2 |
| Visual molecule transitions | HIGH | HIGH | P2 |
| Optimization metrics dashboard | MEDIUM | MEDIUM | P2 |
| Temperature reheat | LOW | MEDIUM | P3 |
| Comparative runs | MEDIUM | HIGH | P3 |
| Annotated molecular features | MEDIUM | HIGH | P3 |
| Custom displacement rules | LOW | HIGH | P3 |
| Share/permalink URLs | LOW | MEDIUM | P3 |

**Priority key:**
- P1: Must have for launch (MVP)
- P2: Should have, add when MVP validated
- P3: Nice to have, future consideration

**Prioritization rationale:**
- **P1 features** form complete learning experience: enter molecule, configure SA, watch optimization, see results
- **P2 features** deepen understanding: step-through for comprehension, downloads for assignments, transitions for mechanism clarity
- **P3 features** expand scope: advanced techniques, comparative analysis, customization beyond Faulon paper

## Competitor/Comparable Analysis

| Feature | VisuAlgo (algorithm viz) | MolView (chemistry viz) | TSP SA Demos | Our Approach |
|---------|--------------------------|-------------------------|--------------|--------------|
| Play/Pause/Reset | Yes, standard controls | N/A (static) | Yes | Yes - table stakes |
| Step-by-step mode | Yes, critical feature | N/A | Yes in advanced demos | Yes - P2 after validation |
| Custom input | Yes, key differentiator | Yes, molecular formulas | Mixed - some demos only presets | Yes - P1 MVP |
| Preset examples | Yes, multiple algorithms | Yes, compound database | Yes, famous TSP instances | Yes - P1, 3-5 curated |
| Real-time metrics | Yes, shows comparisons | No | Yes, cost/temperature plots | Yes - P1 chart + state |
| Mobile support | Yes, 2026 standard | Yes | Mixed | Yes - P1 responsive design |
| Export/download | No | Yes, multiple formats | Mixed | Yes - P2 when validated |
| Parameter tuning | Yes, algorithm-specific | N/A | Yes, SA parameters | Yes - P1 kT/cooling/steps |
| 3D visualization | N/A | Yes, 3D models | N/A | No - anti-feature for 2D graphs |
| Annotation/explanation | Yes, algorithm steps | Yes, properties | Minimal | P3 - complex to do well |

**Key insights:**
- **Algorithm visualization standard:** Play/pause/reset, step-by-step, real-time metrics are expected in all modern algorithm demos
- **Chemistry visualization standard:** Formula input, 2D/3D structures, export are expected in chemistry tools
- **Our unique position:** Intersection of both - chemistry + algorithm visualization has specific requirements from both domains
- **Differentiation opportunity:** Live optimization chart + chemistry structure updates = something neither pure algorithm nor pure chemistry tools do

## Domain-Specific Patterns

### Chemistry Education Apps (2026)

**Must have:**
- Molecular formula as primary input (ChemCollective, MolView pattern)
- 2D structure rendering (ChemDoodle, MolView, RDKit standard)
- Mobile-first or responsive (WebMO BYOD research, Aktiv Chemistry mobile-first)
- Instant feedback (Aktiv Chemistry personalized remediation model)
- Preset examples (reduces barrier, MolPrime+ property presets)

**Valued but not required:**
- Export capabilities (WebMO allows multiple formats)
- Integration with LMS (Aktiv Chemistry syncs gradebooks - NOT relevant for standalone demo)
- Gamification (Interactive Chemistry "Jewels" - NOT appropriate for serious SA demo)

### Algorithm Visualization Education (2026)

**Must have:**
- User control over speed (VisuAlgo research: engagement requires control)
- Ability to use own input (VisuAlgo key feature)
- Visual representation of state (all algorithm viz shows this)
- Step-by-step mode (VisuAlgo research: critical for comprehension)

**Valued but not required:**
- Comparative visualizations (useful but complex)
- Performance metrics (helps understanding but not essential)
- Annotations/explanations (high value when done well)

### Parameter Tuning Education

**Must have:**
- Interactive controls (ICP-MS TuneSim slider-based approach)
- Real-time feedback on adjustments (PhET simulations model)
- Reasonable defaults (students need starting point)
- Visual indication of parameter effects (seeing impact is the point)

**Valued but not required:**
- Presets showing good/bad parameter choices
- Guided exploration (unguided learning is ineffective per research)

## Educational Context Insights

**From research findings:**

1. **Active engagement beats passive viewing** (VisuAlgo research): Students learn more when they control and question the visualization, not just watch it
   - **Implication:** Step-by-step mode and parameter tuning are educationally critical

2. **Own input > provided examples** (VisuAlgo pattern): Students understand better when using their own data
   - **Implication:** Formula input is P1; presets are just for lowering initial barrier

3. **Animation aids understanding but can overwhelm** (chart animation research): Breaking down steps helps, but too much animation = information overload
   - **Implication:** Visual molecule transitions are valuable but should be optional/slow-down is ok

4. **Mobile accessibility is expected** (WebMO BYOD research): Students use personal devices in class
   - **Implication:** Responsive design is not optional in 2026

5. **Parameter tuning requires balance** (ICP-MS TuneSim): Multiple factors need to be adjusted together; UI must show relationships
   - **Implication:** kT, cooling, steps controls should show interdependencies

6. **Export enables assessment** (WebMO research): Instructors want students to save and submit work
   - **Implication:** Download feature enables classroom use; P2 priority

7. **Accessibility is becoming regulatory requirement** (April 2026 WCAG 2.1 AA deadline for higher ed)
   - **Implication:** Color choices, alt text, keyboard navigation should be considered even in MVP

## Sources

**Chemistry Education Web Apps:**
- [ChemCollective](https://chemcollective.org/) - Virtual labs and scenario-based learning
- [Aktiv Chemistry](https://aktiv.com/) - Personalized remediation and mobile-first design
- [Interactive Chemistry](https://interactivechemistry.org/) - Gamified molecular building
- [RSC Education - Chemistry Apps](https://edu.rsc.org/ideas/eight-great-apps-for-chemistry-teachers/2010025.article) - Teacher recommendations

**Molecular Visualization:**
- [Mol* Toolkit](https://molstar.org/) - High-performance WebGL visualization
- [MolView](https://molview.org/) - Web-based structure drawing and 3D visualization
- [Wiley - Molecular Graphics Roadmap](https://onlinelibrary.wiley.com/doi/10.1002/pro.70457) - Interactive molecular graphics best practices
- [RCSB Molecular Graphics Software](https://www.rcsb.org/docs/additional-resources/molecular-graphics-software) - Software overview

**Cheminformatics Education:**
- [Practical Cheminformatics Tutorials](https://github.com/PatWalters/practical_cheminformatics_tutorials) - Jupyter notebooks
- [NFDI4Chem](https://nfdi4chem.de/access-to-open-cheminformatics/) - Open cheminformatics access
- [ChemAxon Platform](https://chemaxon.com/blog/webinar/easily-integrable-instant-cheminformatics-platform) - Web-based tools

**Algorithm Visualization Education:**
- [VisuAlgo](https://visualgo.net/en) - Data structures and algorithms visualization (research shows active engagement critical)
- [PMC - Algorithm Education Practices](https://pmc.ncbi.nlm.nih.gov/articles/PMC6302863/) - Discovery learning research
- [ACM - Algorithm Visualization in CS Education](https://dl.acm.org/doi/10.1145/774833.774846) - Educational impact research
- [Springer - Algorithm Visualization for Data Structures](https://link.springer.com/chapter/10.1007/978-3-031-78561-0_16) - Teaching effectiveness

**Simulated Annealing Demos:**
- [SA Visualization by vgarciasc](https://vgarciasc.github.io/simulated-annealing-viz/) - Clustering application
- [CSE442 - Visualizing SA](https://cse442-17f.github.io/simulated-annealing/) - Interactive controls
- [TSP SA by Inspiaaa](https://github.com/Inspiaaa/TSP-Simulated-Annealing) - Interactive TSP solver
- [Todd Schneider - TSP with SA](https://toddwschneider.com/posts/traveling-salesman-with-simulated-annealing-r-and-shiny/) - R and Shiny implementation
- [CADApps - Floorplanning SA](https://sites.lafayette.edu/cadapps/main-page/floorplanning-simulated-annealing/) - Temperature/cost visualization

**Computational Chemistry Education:**
- [WebMO](https://www.webmo.net/) - Web-based quantum chemistry (BYOD paradigm research)
- [Wiley - WebMO in Education](https://wires.onlinelibrary.wiley.com/doi/epdf/10.1002/wcms.1554) - Barriers minimized through web interface
- [eChem Notebooks](https://pubs.acs.org/doi/10.1021/acs.jchemed.2c01103) - Interactive computational chemistry
- [Wiley - FOSS Computational Chemistry](https://wires.onlinelibrary.wiley.com/doi/10.1002/wcms.1610) - Open source tools

**Molecular Structure Optimization:**
- [PennyLane - Geometry Optimization](https://pennylane.ai/qml/demos/tutorial_mol_geo_opt/) - Quantum approach demo
- [Shodor - Geometry Optimization](http://www.shodor.org/chemviz/optimization/students/background.html) - Educational background

**Chemistry Parameter Tuning:**
- [ACS - ICP-MS TuneSim](https://pubs.acs.org/doi/10.1021/acs.jchemed.8b00343) - Parameter adjustment simulation (slider-based)
- [ACS - Spectroscopic Parameterization](https://pubs.acs.org/doi/10.1021/acs.jchemed.0c00925) - Computer-aided training
- [ResearchGate - PhET Interactive Simulations](https://www.researchgate.net/publication/272272239_PhET_Interactive_Simulations_Transformative_Tools_for_Teaching_Chemistry) - Interactive learning

**Molecular Structure Input:**
- [Journal of Cheminformatics - Structure Input on Web](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-2-1) - Web input patterns
- [ChemDoodle 2D Sketcher](https://web.chemdoodle.com/demos/2d-sketcher) - Drawing tools

**Chart Animation Education:**
- [Flourish - Animated Charts](https://flourish.studio/blog/animated-charts/) - Data storytelling with animation
- [Google Charts - Animation](https://developers.google.com/chart/interactive/docs/animation) - Real-time updates

**Educational Technology 2026:**
- [Cult of Pedagogy - Ed Tech Tools 2026](https://www.cultofpedagogy.com/6-ed-tech-tools-2026/) - Export and share features
- [The 74 - Tech Tools for Educators 2026](https://www.the74million.org/article/10-useful-tech-tools-for-educators-in-2026-a-practical-guide/) - Practical guide

**Accessibility:**
- [Chemistry World - Accessible Science Education](https://www.chemistryworld.com/features/accessible-science-education/3010431.article) - Accessibility features
- [RSC - Inclusive Chemistry Education](https://books.rsc.org/books/edited-volume/2397/Inclusive-and-Accessible-Chemistry-for-Further-and) - 2026 resource
- [Techcolite - Educational Websites 2026](https://www.techcolite.com/building-educational-websites-in-2026-ux-access/) - WCAG 2.1 AA compliance deadline April 24, 2026

---

**Research confidence: MEDIUM**

**High confidence areas:** Algorithm visualization patterns (VisuAlgo research authoritative), chemistry education tool patterns (multiple consistent sources), parameter tuning pedagogy (peer-reviewed research)

**Medium confidence areas:** Specific SA demo feature prevalence (found examples but limited comparative data), mobile usage in chemistry education (WebMO research strong but single source)

**Low confidence areas:** Annotated molecular features implementation complexity (extrapolated from similar tools), comparative runs educational value (no specific research on this feature)

**Validation recommendations:**
- Conduct user testing with chemistry students to validate P1/P2 priority splits
- Survey chemistry instructors about download format preferences (PNG vs SVG, CSV vs JSON)
- Test visual molecule transitions with target users to assess if complexity is justified by educational value

---
*Feature research for: WebFaulon - Interactive SA Molecular Optimization Demo*
*Researched: 2026-02-14*
*Researcher: Claude Opus 4.6*
