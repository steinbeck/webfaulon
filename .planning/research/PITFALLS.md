# Pitfalls Research

**Domain:** Browser-based cheminformatics with Simulated Annealing for constitutional isomer generation
**Researched:** 2026-02-14
**Confidence:** MEDIUM

## Critical Pitfalls

### Pitfall 1: Incorrect Bond Order Redistribution During SA Moves

**What goes wrong:**
The SA displacement equations (eqs 7-11 in Faulon 1996) must correctly redistribute bond orders while preserving atomic valences. If implementation errors exist in these equations, generated structures will have invalid valences (e.g., carbon with 5 bonds, nitrogen with 4 bonds). This invalidates the entire SA search and produces chemically impossible molecules.

**Why it happens:**
The Faulon algorithm involves complex bond order arithmetic where removing a bond order from one location requires adding it elsewhere while respecting valence constraints. Developers often:
- Fail to account for all edge cases in bond redistribution
- Incorrectly implement the multi-step constraint checking (connectivity, then valence)
- Miss the requirement that SA moves must be reversible for detailed balance

**How to avoid:**
- Implement comprehensive unit tests for every SA move type with known valid/invalid transitions
- Verify detailed balance: if move A→B is accepted, reverse move B→A must be possible
- After each move, validate: (1) graph connectivity, (2) all atom valences, (3) total bond order conservation
- Test with edge cases: single atoms, disconnected fragments, saturated vs. unsaturated molecules

**Warning signs:**
- Generated molecules with impossible valences
- SA acceptance ratio > 95% (too permissive) or < 5% (too restrictive)
- Molecules becoming disconnected during SA search
- Total bond order changing between steps (should be invariant)

**Phase to address:**
Phase 1 (Core SA Algorithm Implementation) - Build comprehensive test suite before integration

---

### Pitfall 2: Graph Connectivity Loss After Bond Moves

**What goes wrong:**
When bond orders are redistributed during SA moves, the molecule can become disconnected into separate fragments. For example, removing the only single bond connecting two ring systems creates two separate molecules. Disconnected structures are invalid for constitutional isomer generation but may not be detected if connectivity checking is missing or incorrect.

**Why it happens:**
- Developers assume bond redistribution automatically preserves connectivity
- Connectivity checking (BFS/DFS from arbitrary node) is computationally expensive and skipped for performance
- Edge case: reducing a bond order from 1 to 0 effectively removes the bond
- After complex multi-bond moves, subtle connectivity losses are hard to trace

**How to avoid:**
- ALWAYS perform connectivity check after every SA move that modifies bond orders
- Use efficient BFS from single starting node - O(n+m) for molecular graphs with n atoms, m bonds
- Reject any move that creates disconnected fragments immediately
- For debugging: log which specific bond change caused disconnection

**Warning signs:**
- Generated structures with unexpected "." separators (fragment notation)
- Wiener index becoming undefined (infinite distances between disconnected components)
- Acceptance ratio drops unexpectedly (many moves being rejected late in validation)
- Visual rendering shows separate molecular fragments

**Phase to address:**
Phase 1 (Core SA Algorithm) - Implement as mandatory post-move validation, before valence checking

---

### Pitfall 3: WASM Module Initialization and Memory Management in Workers

**What goes wrong:**
RDKit.js WASM modules cannot be simply imported or passed to Web Workers via postMessage due to structured clone limitations. The WASM binary and JavaScript wrapper must be loaded within the worker context, and developers frequently attempt to share WASM instances across worker boundaries, causing DataCloneError or silent failures.

**Why it happens:**
Based on research, WebAssembly modules are not supported by MessagePort structured clone when crossing process boundaries. Developers expect WASM modules to behave like regular JavaScript objects. The RDKit.js documentation for modern frameworks (React, etc.) notes that "you cannot simply import rdkit from '@rdkit/rdkit'" and requires serving the WASM file from a public location accessible to workers.

**How to avoid:**
- Initialize RDKit.js WASM module INSIDE each Web Worker, not in main thread
- Copy both RDKit_minimal.js and RDKit_minimal.wasm to public assets directory
- Configure worker to load WASM from known public URL path
- Use importScripts() in worker for RDKit_minimal.js, then await initRDKitModule()
- Pass serializable molecular graph representations (adjacency lists, SMILES strings) between worker and main thread, never WASM molecule objects
- Implement WASM initialization error handling with clear messages about missing WASM files

**Warning signs:**
- DataCloneError when posting messages to/from workers
- "RDKit_minimal.wasm not found" errors despite correct npm installation
- WASM module works in main thread but fails in worker
- Browser console shows 404 for .wasm file despite being in node_modules
- Different behavior between dev server (working) and production build (failing)

**Phase to address:**
Phase 2 (Web Worker Integration) - Critical architecture decision, design before implementation

---

### Pitfall 4: Premature Convergence Due to Inappropriate Cooling Schedules

**What goes wrong:**
The Faulon paper describes 33 cooling schedules (f0-f32) with different convergence properties. Using an inappropriate schedule causes either: (a) premature convergence to local minima without exploring diverse structures, or (b) excessive computation time with no convergence. Research shows that if temperature decreases too quickly, the algorithm converges prematurely to suboptimal solutions.

**Why it happens:**
- Developers default to simple linear cooling without understanding the problem landscape
- Cooling schedule parameters from the 1996 paper may not directly transfer to JavaScript/browser context (different performance characteristics)
- Initial temperature not calibrated to achieve ~80% acceptance ratio at start
- Step count and cooling rate not balanced (too few steps with slow cooling, or too many with fast cooling)

**How to avoid:**
- Implement multiple cooling schedules from Faulon paper (geometric, linear, logarithmic variants)
- Measure and log acceptance ratio throughout SA run: should start ~80%, end ~5%
- Adaptive cooling: slow down if acceptance ratio drops too fast, speed up if stuck at plateau
- For initial calibration: run short test with T_initial such that P(accept worse move) ≈ 0.8
- Provide cooling schedule as configurable parameter, not hardcoded
- Benchmark schedules against known test cases (e.g., generate all C6H12 isomers)

**Warning signs:**
- Acceptance ratio drops below 20% in first 10% of steps (too cold too fast)
- Acceptance ratio remains above 60% after 50% of steps (too hot, not converging)
- Same structure appears repeatedly in multiple independent SA runs (stuck in attractor)
- Wiener index shows no improvement over long step sequences
- Generated structures lack diversity even with different random seeds

**Phase to address:**
Phase 3 (SA Parameter Tuning) - After core algorithm works, before production use

---

### Pitfall 5: Inefficient Wiener Index Calculation Creating Performance Bottleneck

**What goes wrong:**
Wiener index (sum of all-pairs shortest paths) is calculated after every SA move to evaluate the objective function. Naive Floyd-Warshall (O(n³)) becomes prohibitively slow even for moderately sized molecules (>30 atoms), causing the browser to freeze or worker to timeout. Research confirms that BFS-based approaches outperform Floyd-Warshall for sparse molecular graphs.

**Why it happens:**
- Floyd-Warshall is conceptually simple and taught in algorithms courses
- Developers don't realize molecular graphs are sparse (average degree ~2-4)
- No performance profiling during development with realistic molecule sizes
- Wiener index calculation happens thousands of times per SA run (once per step)
- Research notes: "Floyd-Warshall is best suited for dense graphs" but molecules are sparse

**How to avoid:**
- Use BFS from each node for unweighted molecular graphs: O(n × (n+m)) = O(n²) for sparse graphs
- Cache results when possible: if move rejected, previous Wiener index still valid
- For molecules >50 atoms, consider incremental updates or approximate Wiener index
- Profile Wiener index calculation separately: if >10ms per call, optimize
- Research shows: for small graphs (<300 vertices), Floyd-Warshall acceptable, but BFS still faster for sparse

**Warning signs:**
- SA step rate < 10 steps/second for molecules with <20 atoms
- Browser DevTools profiler shows >50% time in Wiener index calculation
- Worker timeout errors during SA runs
- UI becomes unresponsive during calculation despite Web Worker usage
- Performance degrades quadratically/cubically with molecule size

**Phase to address:**
Phase 1 (Core SA Algorithm) - Choose correct algorithm from start, changing later is disruptive

---

### Pitfall 6: Invalid Initial Structure Generation from Molecular Formula

**What goes wrong:**
Converting a molecular formula (e.g., C6H12) to a valid initial molecular graph is non-trivial. Naive approaches generate structures with impossible valences, disconnected atoms, or violate chemical constraints. Without a chemically valid starting structure, the SA search cannot produce valid isomers regardless of algorithm correctness.

**Why it happens:**
- Underestimating the complexity of structure generation from formula alone
- Not considering degree of unsaturation (HDI = Hydrogen Deficiency Index)
- Generating random bond assignments without valence checking
- Assuming any connected graph with correct atom counts is valid
- Ignoring that heteroatoms (N, O, S, P) have different valence rules than C and H

**How to avoid:**
- Use deterministic or stochastic structure generators (MAYGEN, OMG-like approaches)
- For simple cases: generate linear chain, then add cycles/branches if HDI indicates unsaturation
- Alternative: start with SMILES string converted via RDKit.js instead of formula
- Validate initial structure: connectivity + all valences correct + matches formula
- Consider using multiple random initial structures to test SA algorithm robustness
- Research shows: deterministic generators (like MOLGEN) vs stochastic (like Faulon's approach) have different trade-offs

**Warning signs:**
- SA runs immediately reject first move (initial structure already invalid)
- Initial structure has zero bonds (just isolated atoms)
- Initial structure violates valence rules but SA proceeds anyway
- Different initial structures for same formula produce vastly different SA results
- Cannot generate initial structure for complex formulas (e.g., with N, O, S)

**Phase to address:**
Phase 1 (Core SA Algorithm) - Initial structure generation is prerequisite for SA testing

---

### Pitfall 7: Implicit vs Explicit Hydrogen Handling Inconsistency

**What goes wrong:**
Cheminformatics tools differ in hydrogen representation: implicit (computed from valence) vs explicit (graph nodes). Mixing approaches or using implicit hydrogens incorrectly leads to valence calculation errors, incorrect Wiener indices (if H atoms included in graph), and RDKit.js conversion failures.

**Why it happens:**
Research shows "a common misconception is that there is a single valence model for chemical informatics; in reality, many vendors and toolkits implement their own interpretations." Developers:
- Switch between implicit/explicit representations without conversion
- Include hydrogen atoms in Wiener index calculation (rarely desired for constitutional isomers)
- Assume RDKit.js auto-handles hydrogens (it doesn't always)
- Don't understand that implicit H count = target_valence - sum(bond_orders)

**How to avoid:**
- Choose ONE representation: implicit hydrogens (recommended for constitutional isomer graphs)
- Document clearly: "Hydrogens are implicit, not graph nodes"
- When calculating Wiener index: use heavy atoms only (C, N, O, etc.), exclude H
- When converting to RDKit.js: explicitly set implicit hydrogen counts on each atom
- Valence checking: bond_order_sum + implicit_H_count must equal standard valence
- Use RDKit.js's addHs()/removeHs() carefully and only when needed for visualization

**Warning signs:**
- Valence errors after converting between internal graph and RDKit molecule
- Wiener index values seem too large (counting H atoms when shouldn't)
- SMILES string from RDKit.js has unexpected [H] explicit notations
- Different results from same structure depending on hydrogen representation
- Aromatic molecules fail kekulization in RDKit.js

**Phase to address:**
Phase 1 (Core Graph Representation) - Architectural decision before algorithm implementation

---

### Pitfall 8: Aromaticity and Kekulization Failures in RDKit.js

**What goes wrong:**
When converting generated molecular graphs to RDKit.js molecules, aromatic systems require kekulization (assigning specific bond orders to aromatic bonds). If the graph representation uses integer bond orders (1, 2, 3) but the molecule is aromatic (benzene, pyridine), RDKit.js kekulization can fail with cryptic errors or produce incorrect structures.

**Why it happens:**
Research shows: "If no perfect matching for the pruned delocalization set exists, the corresponding SMILES is invalid." The Faulon algorithm operates on integer bond orders, but aromatic molecules need special handling. RDKit.js expects either:
- Aromatic SMILES notation (lowercase atoms: c1ccccc1)
- Kekule form with explicit double bonds (C1=CC=CC=C1)

Mixing these or providing impossible configurations causes failures.

**How to avoid:**
- For constitutional isomer generation: use Kekule structures (explicit double bonds), NOT aromatic notation
- Document that generated structures are non-aromatic Kekule forms
- If aromatic detection needed: let RDKit.js sanitize and detect aromaticity AFTER structure generation
- Test with aromatic molecules: benzene (C6H6), pyridine (C5H5N), naphthalene (C10H8)
- Research notes: "kekulization involves finding a perfect matching to satisfy valencies"
- Ensure all atoms have defined formal charges and implicit hydrogen counts before kekulization

**Warning signs:**
- RDKit.js sanitization errors: "Can't kekulize mol"
- Aromatic SMILES input works but generated structures fail
- Molecules with rings all fail conversion while acyclic succeed
- KekulizeException in RDKit.js console logs
- Visual rendering shows incorrect bond orders in aromatic rings

**Phase to address:**
Phase 2 (RDKit.js Integration) - After SA works with graph representation, before visualization

---

## Technical Debt Patterns

Shortcuts that seem reasonable but create long-term problems.

| Shortcut | Immediate Benefit | Long-term Cost | When Acceptable |
|----------|-------------------|----------------|-----------------|
| Skip connectivity check "for performance" | 15-20% faster SA steps | Invalid disconnected structures generated, hard-to-debug failures, wasted computation on invalid paths | Never - O(n+m) BFS is fast enough |
| Use Floyd-Warshall for Wiener index | Simpler to implement (3 nested loops) | O(n³) bottleneck, unusable for >30 atom molecules, browser freezes | Small demos only (<15 atoms), replace before production |
| Hardcode cooling schedule | Faster initial development | No adaptability to different problems, premature convergence, need code changes to experiment | Early prototyping only, must parameterize before phase 3 |
| Pass SMILES strings between worker and main | Avoids graph serialization complexity | RDKit.js WASM must run in both contexts (memory overhead), slower than adjacency list serialization | If RDKit.js needed in both contexts anyway |
| Generate single random initial structure | Simple implementation | SA convergence depends heavily on initial structure quality, low reproducibility | Early testing only, use multiple starts in production |
| Use implicit hydrogens without validation | Smaller graphs, faster computation | Valence errors hidden until RDKit.js conversion, silent failures in complex molecules | Acceptable IF comprehensive valence validation exists |

## Integration Gotchas

Common mistakes when connecting to external services.

| Integration | Common Mistake | Correct Approach |
|-------------|----------------|------------------|
| RDKit.js WASM | Attempting to import module in ES6 style or pass WASM objects via postMessage | Copy .wasm and .js to public assets, use importScripts() in worker, initialize per-worker with initRDKitModule() |
| Web Worker messaging | Sending complex molecular objects assuming they're serializable | Send adjacency list representation {atoms: [...], bonds: [...]} or SMILES strings, reconstruct objects on receiving side |
| RDKit.js molecule creation | Creating molecule without setting implicit hydrogens or formal charges | Explicitly set implicit_H_count and formal_charge for each atom before molecule creation |
| SMILES parsing | Expecting RDKit.js to handle all valid SMILES | Some edge cases fail; validate SMILES independently and handle parse errors gracefully |
| Worker thread count | Spawning one worker per SA run for parallelism | Browser limits workers (typically 4-8); use worker pool pattern, reuse workers across runs |

## Performance Traps

Patterns that work at small scale but fail as usage grows.

| Trap | Symptoms | Prevention | When It Breaks |
|------|----------|------------|----------------|
| Wiener index with Floyd-Warshall | UI freezes, worker timeouts, exponential slowdown | Use BFS-based all-pairs shortest paths O(n²) for sparse graphs | >20 atoms (noticeable), >30 atoms (unusable) |
| Recreating RDKit.js module per SA step | Memory leaks, increasing latency | Initialize WASM once per worker lifetime, reuse module instance | After ~100 SA steps, garbage collection issues |
| Logging every SA step to console | Browser DevTools becomes sluggish | Log summary statistics every N steps (e.g., N=100), full logging only in debug mode | >1000 SA steps, DevTools open |
| Deep copying molecular graphs | Quadratic memory growth, GC pauses | Use structural sharing or immutable data structures, clone only what changes | Molecules >30 atoms, 10K+ SA steps |
| Synchronous RDKit.js operations in worker | Worker becomes unresponsive to termination | Use async/await pattern even in workers, check for cancellation signals periodically | Long-running SA (>30 seconds) |

## Security Mistakes

Domain-specific security issues beyond general web security.

| Mistake | Risk | Prevention |
|---------|------|------------|
| Accepting arbitrary molecular formulas without validation | Resource exhaustion DoS (e.g., C1000H2000 formula crashes browser) | Validate formula: max atoms (e.g., 100), valid elements only, reasonable ratios |
| Loading WASM from CDN without integrity check | WASM binary tampering, arbitrary code execution | Use subresource integrity (SRI) hashes, or bundle WASM in app assets |
| Allowing unlimited SA steps from user input | Worker runs indefinitely, battery drain, thermal throttling | Enforce max steps (e.g., 100K), timeout (e.g., 60s), allow user cancellation |
| Exposing internal graph representation in API | Information leakage, coupling to internal implementation | Export only SMILES/Molfile formats, treat internal graph as private |

## UX Pitfalls

Common user experience mistakes in this domain.

| Pitfall | User Impact | Better Approach |
|---------|-------------|-----------------|
| No progress indication during SA | User doesn't know if browser is frozen or computing | Stream progress updates from worker every N steps (e.g., steps completed, current Wiener index, acceptance ratio) |
| Displaying cryptic RDKit.js error messages | "KekulizeException" means nothing to chemists | Translate technical errors to user-friendly messages: "Could not assign bond orders to aromatic structure" |
| Showing only final structure, not SA trajectory | User can't understand convergence behavior or quality | Visualize Wiener index over time, show intermediate structures, convergence graph |
| Accepting only molecular formulas | Limits users who want to start from specific structure | Also accept SMILES input as initial structure alternative to formula |
| No way to cancel long-running SA | User forced to close browser tab | Implement cancel button that terminates worker gracefully |
| Regenerating same structure without feedback | User clicks "Generate" repeatedly, gets identical results | Show random seed, allow seed input for reproducibility, indicate when result is identical to previous |

## "Looks Done But Isn't" Checklist

Things that appear complete but are missing critical pieces.

- [ ] **SA algorithm:** Often missing post-move connectivity validation - verify every move checks graph remains connected via BFS
- [ ] **Wiener index:** Often missing edge case for disconnected graphs - verify returns infinity or error for disconnected structures
- [ ] **RDKit.js integration:** Often missing explicit hydrogen setup - verify implicit_H_count set on all atoms before molecule creation
- [ ] **Web Worker:** Often missing WASM initialization error handling - verify clear error message when .wasm file not found
- [ ] **Initial structure:** Often missing valence validation - verify starting structure has all valid valences before SA begins
- [ ] **Cooling schedule:** Often missing acceptance ratio monitoring - verify logging acceptance ratio throughout SA run
- [ ] **Bond order moves:** Often missing detailed balance - verify reverse moves possible (A→B implies B→A transition exists)
- [ ] **Aromatic handling:** Often missing kekulization error handling - verify graceful failure when aromatic structure can't be kekulized
- [ ] **Performance testing:** Often missing realistic molecule sizes - verify testing with 20-50 atom molecules, not just C6H6
- [ ] **Worker cleanup:** Often missing proper termination - verify workers are terminated and memory released after SA completion

## Recovery Strategies

When pitfalls occur despite prevention, how to recover.

| Pitfall | Recovery Cost | Recovery Steps |
|---------|---------------|----------------|
| Invalid bond order redistribution | HIGH | Full algorithm rewrite required; comprehensive test suite needed first; expect 2-3 weeks |
| Graph connectivity loss | MEDIUM | Add BFS check after moves; may require move generator redesign; 3-5 days |
| WASM initialization failure | LOW | Refactor asset serving configuration; update worker initialization; 1 day |
| Premature convergence | MEDIUM | Implement multiple cooling schedules; add adaptive cooling; requires re-tuning; 5-7 days |
| Wiener index bottleneck | MEDIUM | Replace Floyd-Warshall with BFS approach; extensive testing needed; 3-5 days |
| Invalid initial structure | MEDIUM | Implement proper structure generator or use RDKit.js; 3-5 days |
| Hydrogen handling errors | LOW | Standardize on implicit hydrogens; update conversion logic; 2-3 days |
| Kekulization failures | MEDIUM | Redesign to avoid aromatic structures or handle properly; test coverage; 3-5 days |

## Pitfall-to-Phase Mapping

How roadmap phases should address these pitfalls.

| Pitfall | Prevention Phase | Verification |
|---------|------------------|--------------|
| Invalid bond redistribution | Phase 1: Core SA Algorithm | Test suite with 100+ move validations, all pass |
| Connectivity loss | Phase 1: Core SA Algorithm | BFS check implemented, no disconnected structures in 1000-step SA runs |
| WASM initialization | Phase 2: Web Worker Integration | WASM loads successfully in worker, error handling tested |
| Premature convergence | Phase 3: SA Parameter Tuning | Acceptance ratio follows expected curve (80% start, 5% end) |
| Wiener index bottleneck | Phase 1: Core SA Algorithm | Wiener calculation <5ms for 50-atom molecules |
| Invalid initial structure | Phase 1: Initial Structure Generation | Generated structures pass RDKit.js sanitization 100% |
| Hydrogen handling | Phase 1: Graph Representation Design | Consistent implicit H throughout, RDKit.js conversion works |
| Kekulization failures | Phase 2: RDKit.js Integration | Aromatic test molecules (benzene, pyridine, naphthalene) convert successfully |

## Sources

### WASM and Web Workers
- [RDKit.js GitHub Repository](https://github.com/rdkit/rdkit-js)
- [RDKit.js npm Package](https://www.npmjs.com/package/@rdkit/rdkit)
- [WebAssembly Structured Clone Issues](https://bugzilla.mozilla.org/show_bug.cgi?id=1412852)
- [WebAssembly MessagePort Limitations](https://bugzilla.mozilla.org/show_bug.cgi?id=1564376)
- [MDN: Structured Clone Algorithm](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Structured_clone_algorithm)
- [MDN: Using Web Workers](https://developer.mozilla.org/en-US/docs/Web/API/Web_Workers_API/Using_web_workers)

### Simulated Annealing
- [Wikipedia: Simulated Annealing](https://en.wikipedia.org/wiki/Simulated_annealing)
- [Simulated Annealing: From Basics to Applications](https://enac.hal.science/hal-01887543v1/document)
- [Computing Initial Temperature of SA](https://www.researchgate.net/publication/227061666_Computing_the_Initial_Temperature_of_Simulated_Annealing)
- [Algorithm Afternoon: Temperature and Cooling Schedules](https://algorithmafternoon.com/books/simulated_annealing/chapter02/)
- [Tuning Parameters in Simulated Annealing](https://aicompetence.org/tuning-parameters-in-simulated-annealing/)

### Molecular Structure Generation
- [MAYGEN: Constitutional Isomer Generator](https://link.springer.com/article/10.1186/s13321-021-00529-9)
- [Faulon's Research: Chemical Structure Generation](https://jfaulon.com/bioretrosynth-background/chemical-structure-and-network-generation/)
- [Fiehn Lab: Molecular Isomer Generators](https://fiehnlab.ucdavis.edu/projects/seven-golden-rules/molecular-isomer-generator)
- [Surge: Constitutional Isomer Generator](https://link.springer.com/article/10.1186/s13321-022-00604-9)

### Molecular Graph Theory
- [Wikipedia: Wiener Index](https://en.wikipedia.org/wiki/Wiener_index)
- [Floyd-Warshall for Molecular Structures](https://towardsdatascience.com/the-floyd-warshall-algorithm-from-graph-theory-applied-to-parsing-molecular-structures-39f8c99c9fe1/)
- [RDKit Cookbook](https://www.rdkit.org/docs/Cookbook.html)
- [Molecular Representations in AI](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00460-5)

### Hydrogen and Aromaticity Handling
- [NextMove: Explicit and Implicit Hydrogens](https://nextmovesoftware.com/blog/2013/02/27/explicit-and-implicit-hydrogens-taking-liberties-with-valence/)
- [Depth-First: Hydrogen Suppression in SMILES](https://depth-first.com/articles/2020/06/08/hydrogen-suppression-in-smiles/)
- [Depth-First: Writing Aromatic SMILES](https://depth-first.com/articles/2021/06/30/writing-aromatic-smiles/)
- [Why Drop Bond Order and Aromaticity](https://blog.omsf.io/why-we-should-drop-graph-based-bond-order-and-no-longer-use-aromaticity/)
- [ChemAxon: Implicit, Explicit and Query Hydrogens](https://docs.chemaxon.com/display/docs/implicit-explicit-and-query-hydrogens.md)

---
*Pitfalls research for: WebFaulon - Browser-based Faulon SA Algorithm Implementation*
*Researched: 2026-02-14*
*Confidence: MEDIUM - Based on web search of current documentation, academic sources, and community discussions. Faulon-specific implementation details are low confidence (original 1996 paper not directly accessed). WASM/Worker integration is high confidence (official MDN docs). SA theory is high confidence (academic consensus).*
