import { Atom, STANDARD_VALENCE } from './types';

export class MolGraph {
  private atoms: Atom[];
  private bonds: number[][];

  constructor(atoms: Atom[], bonds: number[][]) {
    // Validate bonds matrix is square and matches atoms length
    if (bonds.length !== atoms.length) {
      throw new Error(`Bond matrix size (${bonds.length}) does not match atom count (${atoms.length})`);
    }
    for (let i = 0; i < bonds.length; i++) {
      if (bonds[i]!.length !== atoms.length) {
        throw new Error(`Bond matrix row ${i} has length ${bonds[i]!.length}, expected ${atoms.length}`);
      }
    }

    // Deep copy atoms and bonds to ensure immutability of input
    this.atoms = atoms.map(a => ({ ...a }));
    this.bonds = bonds.map(row => [...row]);

    // Compute implicitH for each atom
    this.recomputeImplicitH();
  }

  private recomputeImplicitH(): void {
    for (let i = 0; i < this.atoms.length; i++) {
      const atom = this.atoms[i]!;
      const bondOrderSum = this.getBondOrderSum(i);
      const standardValence = STANDARD_VALENCE[atom.element];
      if (standardValence === undefined) {
        throw new Error(`Unknown element: ${atom.element}`);
      }
      atom.implicitH = standardValence - bondOrderSum;
    }
  }

  // Getters
  getAtomCount(): number {
    return this.atoms.length;
  }

  getAtom(index: number): Atom {
    const atom = this.atoms[index];
    if (atom === undefined) {
      throw new Error(`Atom index ${index} out of bounds`);
    }
    return { ...atom }; // Return readonly copy
  }

  getAtoms(): readonly Atom[] {
    return this.atoms.map(a => ({ ...a }));
  }

  getBondOrder(i: number, j: number): number {
    if (i < 0 || i >= this.atoms.length || j < 0 || j >= this.atoms.length) {
      throw new Error(`Bond indices (${i}, ${j}) out of bounds`);
    }
    return this.bonds[i]![j]!;
  }

  getBondMatrix(): readonly number[][] {
    return this.bonds.map(row => [...row]);
  }

  getBondOrderSum(atomIndex: number): number {
    if (atomIndex < 0 || atomIndex >= this.atoms.length) {
      throw new Error(`Atom index ${atomIndex} out of bounds`);
    }
    let sum = 0;
    for (let j = 0; j < this.atoms.length; j++) {
      sum += this.bonds[atomIndex]![j]!;
    }
    return sum;
  }

  // Validation methods
  isConnected(): boolean {
    if (this.atoms.length === 0) return true;
    if (this.atoms.length === 1) return true;

    // BFS from atom 0
    const visited = new Array(this.atoms.length).fill(false);
    const queue: number[] = [0];
    visited[0] = true;
    let visitedCount = 1;

    while (queue.length > 0) {
      const current = queue.shift()!;
      for (let neighbor = 0; neighbor < this.atoms.length; neighbor++) {
        if (this.bonds[current]![neighbor]! > 0 && !visited[neighbor]) {
          visited[neighbor] = true;
          queue.push(neighbor);
          visitedCount++;
        }
      }
    }

    return visitedCount === this.atoms.length;
  }

  hasValidValences(): boolean {
    for (let i = 0; i < this.atoms.length; i++) {
      const atom = this.atoms[i]!;
      const standardValence = STANDARD_VALENCE[atom.element];
      if (standardValence === undefined) {
        return false;
      }
      // Valid if bondOrderSum + implicitH == standardValence
      // Since implicitH is computed as standardValence - bondOrderSum,
      // we need to check if implicitH >= 0 (not overvalent)
      if (atom.implicitH < 0) {
        return false;
      }
    }
    return true;
  }

  validate(): { connected: boolean; validValences: boolean } {
    return {
      connected: this.isConnected(),
      validValences: this.hasValidValences(),
    };
  }

  // Mutation
  setBond(i: number, j: number, order: number): void {
    if (i < 0 || i >= this.atoms.length || j < 0 || j >= this.atoms.length) {
      throw new Error(`Bond indices (${i}, ${j}) out of bounds`);
    }
    if (order < 0 || order > 3) {
      throw new Error(`Bond order ${order} out of range [0, 3]`);
    }

    // Set symmetrically
    this.bonds[i]![j] = order;
    this.bonds[j]![i] = order;

    // Update implicitH for both atoms
    const bondOrderSumI = this.getBondOrderSum(i);
    const bondOrderSumJ = this.getBondOrderSum(j);

    const standardValenceI = STANDARD_VALENCE[this.atoms[i]!.element];
    const standardValenceJ = STANDARD_VALENCE[this.atoms[j]!.element];

    if (standardValenceI === undefined || standardValenceJ === undefined) {
      throw new Error(`Unknown element in setBond`);
    }

    this.atoms[i]!.implicitH = standardValenceI - bondOrderSumI;
    this.atoms[j]!.implicitH = standardValenceJ - bondOrderSumJ;
  }

  // Clone
  clone(): MolGraph {
    const atomsCopy = this.atoms.map(a => ({ ...a }));
    const bondsCopy = this.bonds.map(row => [...row]);
    return new MolGraph(atomsCopy, bondsCopy);
  }

  /**
   * Generate a SMILES string representation of this molecular graph
   *
   * Uses DFS-based SMILES generation. This produces valid (but non-canonical) SMILES
   * that can be parsed by RDKit.js. Canonicalization is RDKit's responsibility.
   *
   * @returns SMILES string representation
   */
  toSMILES(): string {
    if (this.atoms.length === 0) {
      return '';
    }

    // Build adjacency list from bond matrix
    const adjacency: number[][] = Array(this.atoms.length).fill(0).map(() => []);
    for (let i = 0; i < this.atoms.length; i++) {
      for (let j = i + 1; j < this.atoms.length; j++) {
        if (this.bonds[i]![j]! > 0) {
          adjacency[i]!.push(j);
          adjacency[j]!.push(i);
        }
      }
    }

    // Track visited atoms and ring closures
    const visited = new Array(this.atoms.length).fill(false);
    const ringClosures: Map<string, number> = new Map(); // key: "i-j" where i < j
    let nextRingNumber = 1;

    const dfs = (atomIdx: number, parentIdx: number = -1): string => {
      visited[atomIdx] = true;
      let smiles = this.atoms[atomIdx]!.element;

      const neighbors = adjacency[atomIdx]!.filter(n => n !== parentIdx);

      // First pass: handle ring closures (back edges to already-visited atoms)
      for (const neighborIdx of neighbors) {
        if (visited[neighborIdx]) {
          const edgeKey = atomIdx < neighborIdx
            ? `${atomIdx}-${neighborIdx}`
            : `${neighborIdx}-${atomIdx}`;

          if (!ringClosures.has(edgeKey)) {
            // First time encountering this ring closure - open it
            const ringNum = nextRingNumber++;
            ringClosures.set(edgeKey, ringNum);
            smiles += ringNum.toString();
          }
          // Note: if the edge is already in ringClosures, it means we're at the
          // closing end, but we already added the digit when opening, so skip
        }
      }

      // Second pass: traverse unvisited neighbors (tree edges)
      let firstBranch = true;
      for (const neighborIdx of neighbors) {
        // Check if visited NOW (not at the start) because recursive calls might have visited it
        if (visited[neighborIdx]) {
          continue;
        }

        const bondOrder = this.bonds[atomIdx]![neighborIdx]!;

        // Add bond symbol for double/triple bonds
        let bondSymbol = '';
        if (bondOrder === 2) {
          bondSymbol = '=';
        } else if (bondOrder === 3) {
          bondSymbol = '#';
        }

        if (firstBranch) {
          // First branch - continue main chain
          smiles += bondSymbol + dfs(neighborIdx, atomIdx);
          firstBranch = false;
        } else {
          // Additional branches - wrap in parentheses
          smiles += '(' + bondSymbol + dfs(neighborIdx, atomIdx) + ')';
        }
      }

      return smiles;
    };

    return dfs(0);
  }

  // Factory methods for testing
  static createLinearAlkane(n: number): MolGraph {
    if (n < 1) {
      throw new Error('Linear alkane must have at least 1 carbon');
    }

    const atoms: Atom[] = [];
    for (let i = 0; i < n; i++) {
      atoms.push({ element: 'C', implicitH: 0 }); // implicitH will be computed
    }

    // Create bond matrix: n x n, all zeros
    const bonds: number[][] = Array(n).fill(0).map(() => Array(n).fill(0));

    // Add single bonds between consecutive carbons
    for (let i = 0; i < n - 1; i++) {
      bonds[i]![i + 1] = 1;
      bonds[i + 1]![i] = 1;
    }

    return new MolGraph(atoms, bonds);
  }

  static createCyclohexane(): MolGraph {
    const atoms: Atom[] = Array(6).fill(0).map(() => ({ element: 'C', implicitH: 0 }));
    const bonds: number[][] = Array(6).fill(0).map(() => Array(6).fill(0));

    // Create ring: 0-1-2-3-4-5-0
    for (let i = 0; i < 6; i++) {
      const next = (i + 1) % 6;
      bonds[i]![next] = 1;
      bonds[next]![i] = 1;
    }

    return new MolGraph(atoms, bonds);
  }

  static createBranched(pattern: 'isobutane' | 'neopentane'): MolGraph {
    if (pattern === 'isobutane') {
      // Central carbon (index 0) bonded to 3 other carbons (indices 1, 2, 3)
      const atoms: Atom[] = Array(4).fill(0).map(() => ({ element: 'C', implicitH: 0 }));
      const bonds: number[][] = Array(4).fill(0).map(() => Array(4).fill(0));

      bonds[0]![1] = 1;
      bonds[1]![0] = 1;
      bonds[0]![2] = 1;
      bonds[2]![0] = 1;
      bonds[0]![3] = 1;
      bonds[3]![0] = 1;

      return new MolGraph(atoms, bonds);
    } else {
      // neopentane: Central carbon (index 0) bonded to 4 other carbons
      const atoms: Atom[] = Array(5).fill(0).map(() => ({ element: 'C', implicitH: 0 }));
      const bonds: number[][] = Array(5).fill(0).map(() => Array(5).fill(0));

      bonds[0]![1] = 1;
      bonds[1]![0] = 1;
      bonds[0]![2] = 1;
      bonds[2]![0] = 1;
      bonds[0]![3] = 1;
      bonds[3]![0] = 1;
      bonds[0]![4] = 1;
      bonds[4]![0] = 1;

      return new MolGraph(atoms, bonds);
    }
  }
}
