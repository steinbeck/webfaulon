# WebFaulon

Interactive web demonstration of Faulon's simulated annealing algorithm for exploring the space of constitutional isomers.

## Try It Live

**[Launch WebFaulon →](https://steinbeck.github.io/webfaulon/)**

WebFaulon runs entirely in your browser. No installation, no server, no data uploaded. Try it on desktop, tablet, or mobile.

## About

WebFaulon brings to life the simulated annealing algorithm described in Faulon's landmark 1996 paper on computational exploration of constitutional isomer space. Students enter a molecular formula (like C8H18 for octane), configure SA parameters, and watch the algorithm search for molecular structures that optimize the Wiener Index—all in real time.

The tool makes an abstract computational chemistry concept tangible and interactive. Instead of static figures in a textbook, students see the algorithm working: the temperature cooling, the structure evolving, the objective function fluctuating as the Metropolis criterion accepts and rejects proposed changes. It's designed for classroom demonstrations, self-guided learning, and exploration of how stochastic optimization navigates chemical space.

This implementation is based on:

**Faulon, J.-L.** "Stochastic Generator of Chemical Structure. 2. Using Simulated Annealing To Search the Space of Constitutional Isomers." *J. Chem. Inf. Comput. Sci.* **1996**, *36*, 731-740. DOI: [10.1021/ci950093o](https://doi.org/10.1021/ci950093o)

## Features

- **Real-time SA execution with live Wiener Index chart** — Watch the optimization unfold step-by-step, inspired by Figure 4 of the Faulon paper
- **2D molecular structure rendering via RDKit.js** — See the best structure found so far, updated dynamically as the algorithm runs
- **Configurable SA parameters** — Adjust initial temperature, cooling schedule (f0 through f32), steps per cycle, and number of cycles
- **Minimize or maximize Wiener Index** — Explore both optimization directions to understand the structure-property relationship
- **Preset molecules for quick exploration** — Hexane isomers, octane isomers, decane isomers, monoterpene isomers, and more—no chemistry knowledge required to get started
- **Responsive design for classroom projection and mobile devices** — Works on laptops, tablets, and smartphones; optimized for large-screen projection
- **Runs entirely in the browser using Web Workers** — No server required, no data leaves your device, non-blocking execution keeps the UI responsive

## Tech Stack

- **TypeScript + Vite** — Modern build tooling, fast HMR, tree-shaking
- **Alpine.js** — Lightweight reactive UI framework
- **Chart.js** — Live optimization chart rendering
- **RDKit.js (WASM)** — Cheminformatics library for 2D rendering and SMILES generation
- **Web Workers + Comlink** — Non-blocking SA execution in background thread

## Getting Started (Development)

Clone the repository and install dependencies:

```bash
git clone https://github.com/steinbeck/webfaulon.git
cd webfaulon
npm install
npm run dev
```

Then open [http://localhost:5173](http://localhost:5173) in your browser.

## Building for Production

Build the optimized production bundle:

```bash
npm run build
npm run preview
```

The built files will be in `dist/` and can be deployed to any static hosting service (GitHub Pages, Netlify, Vercel, etc.).

## Running Tests

Run the test suite with Vitest:

```bash
npm test
```

For interactive test UI:

```bash
npm run test:ui
```

## How It Works

The implementation follows Faulon's algorithm closely:

1. **Initialization**: The algorithm starts from a random connected molecular graph that satisfies the input molecular formula (e.g., C6H14). Initial bond orders and atom connectivity are generated to respect chemical valence rules.

2. **Displacement Operations**: At each step, the algorithm applies Faulon's displacement operations (equations 7-11 from the paper). Four atoms (x1, y1, x2, y2) are randomly selected, and bond orders among them are redistributed while preserving each atom's valence.

3. **Objective Function**: The Wiener Index—the sum of all pairwise shortest-path distances in the molecular graph—is computed as the objective function. This graph-theoretic measure correlates with molecular properties like boiling point.

4. **Metropolis Criterion**: Proposed moves are accepted or rejected using the Metropolis criterion. Improving moves (lower Wiener Index when minimizing, higher when maximizing) are always accepted. Non-improving moves are accepted probabilistically with probability exp(-ΔE/kT), allowing the algorithm to escape local optima.

5. **Cooling Schedule**: Temperature decreases over time according to a configurable schedule (kT_t = kT_0 - k × kT_0 × t / Δt). The paper found f8 (k=8) with initial kT=100 to be effective.

## Citation

If you use WebFaulon in teaching or research, please cite the original paper:

```bibtex
@article{Faulon1996,
  author = {Faulon, Jean-Loup},
  title = {Stochastic Generator of Chemical Structure. 2. Using Simulated Annealing To Search the Space of Constitutional Isomers},
  journal = {Journal of Chemical Information and Computer Sciences},
  volume = {36},
  number = {4},
  pages = {731--740},
  year = {1996},
  doi = {10.1021/ci950093o}
}
```

## License

MIT License

## Acknowledgments

This project implements the algorithm described by Dr. Jean-Loup Faulon in his 1996 paper. WebFaulon is an independent educational visualization tool and is not affiliated with the original author.

Molecular structure rendering powered by [RDKit.js](https://github.com/rdkit/rdkit-js).
