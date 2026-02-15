import { describe, it, expect } from 'vitest';
import { MolGraph } from '../MolGraph';
import { Atom } from '../types';

describe('MolGraph', () => {
  describe('Construction', () => {
    it('should create a graph from atoms and bonds', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 1],
        [1, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      expect(graph.getAtomCount()).toBe(2);
      expect(graph.getBondOrder(0, 1)).toBe(1);
    });

    it('should throw if bonds matrix size does not match atoms', () => {
      const atoms: Atom[] = [{ element: 'C', implicitH: 0 }];
      const bonds = [
        [0, 1],
        [1, 0],
      ];
      expect(() => new MolGraph(atoms, bonds)).toThrow('does not match atom count');
    });

    it('should throw if bonds matrix is not square', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 1, 0],
        [1, 0],
      ];
      expect(() => new MolGraph(atoms, bonds)).toThrow('has length');
    });
  });

  describe('Linear Alkane Factory', () => {
    it('should create methane (n=1)', () => {
      const graph = MolGraph.createLinearAlkane(1);
      expect(graph.getAtomCount()).toBe(1);
      expect(graph.getAtom(0).element).toBe('C');
      expect(graph.getAtom(0).implicitH).toBe(4); // No bonds, so 4 hydrogens
    });

    it('should create ethane (n=2)', () => {
      const graph = MolGraph.createLinearAlkane(2);
      expect(graph.getAtomCount()).toBe(2);
      expect(graph.getBondOrder(0, 1)).toBe(1);
      expect(graph.getAtom(0).implicitH).toBe(3); // 1 bond, so 3 hydrogens
      expect(graph.getAtom(1).implicitH).toBe(3);
    });

    it('should create pentane (n=5)', () => {
      const graph = MolGraph.createLinearAlkane(5);
      expect(graph.getAtomCount()).toBe(5);
      // Terminal carbons
      expect(graph.getAtom(0).implicitH).toBe(3);
      expect(graph.getAtom(4).implicitH).toBe(3);
      // Internal carbons
      expect(graph.getAtom(1).implicitH).toBe(2);
      expect(graph.getAtom(2).implicitH).toBe(2);
      expect(graph.getAtom(3).implicitH).toBe(2);
    });

    it('should create decane (n=10)', () => {
      const graph = MolGraph.createLinearAlkane(10);
      expect(graph.getAtomCount()).toBe(10);
      expect(graph.isConnected()).toBe(true);
    });

    it('should throw for n < 1', () => {
      expect(() => MolGraph.createLinearAlkane(0)).toThrow();
    });
  });

  describe('Cyclohexane Factory', () => {
    it('should create cyclohexane with 6 carbons', () => {
      const graph = MolGraph.createCyclohexane();
      expect(graph.getAtomCount()).toBe(6);
    });

    it('should have each carbon bonded to 2 neighbors', () => {
      const graph = MolGraph.createCyclohexane();
      for (let i = 0; i < 6; i++) {
        expect(graph.getBondOrderSum(i)).toBe(2);
        expect(graph.getAtom(i).implicitH).toBe(2); // 4 - 2 = 2
      }
    });

    it('should be connected', () => {
      const graph = MolGraph.createCyclohexane();
      expect(graph.isConnected()).toBe(true);
    });
  });

  describe('Branched Factories', () => {
    it('should create isobutane with central carbon bonded to 3 others', () => {
      const graph = MolGraph.createBranched('isobutane');
      expect(graph.getAtomCount()).toBe(4);
      expect(graph.getBondOrderSum(0)).toBe(3); // Central carbon
      expect(graph.getAtom(0).implicitH).toBe(1); // 4 - 3 = 1
      // Terminal carbons
      expect(graph.getBondOrderSum(1)).toBe(1);
      expect(graph.getAtom(1).implicitH).toBe(3);
    });

    it('should create neopentane with central carbon bonded to 4 others', () => {
      const graph = MolGraph.createBranched('neopentane');
      expect(graph.getAtomCount()).toBe(5);
      expect(graph.getBondOrderSum(0)).toBe(4); // Central carbon
      expect(graph.getAtom(0).implicitH).toBe(0); // 4 - 4 = 0
      // Terminal carbons
      for (let i = 1; i < 5; i++) {
        expect(graph.getBondOrderSum(i)).toBe(1);
        expect(graph.getAtom(i).implicitH).toBe(3);
      }
    });
  });

  describe('Connectivity', () => {
    it('should return true for connected graphs', () => {
      const graph = MolGraph.createLinearAlkane(5);
      expect(graph.isConnected()).toBe(true);
    });

    it('should return false for disconnected graphs', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 0], // Atom 2 is disconnected
      ];
      const graph = new MolGraph(atoms, bonds);
      expect(graph.isConnected()).toBe(false);
    });

    it('should handle single atom as connected', () => {
      const graph = MolGraph.createLinearAlkane(1);
      expect(graph.isConnected()).toBe(true);
    });
  });

  describe('Valence Validation', () => {
    it('should return true for valid alkanes', () => {
      const graph = MolGraph.createLinearAlkane(5);
      expect(graph.hasValidValences()).toBe(true);
    });

    it('should return false for overvalent carbon', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 2, 3], // Carbon 0 has bondOrderSum = 5 (overvalent)
        [2, 0, 0],
        [3, 0, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      expect(graph.hasValidValences()).toBe(false);
      expect(graph.getAtom(0).implicitH).toBe(-1); // 4 - 5 = -1
    });

    it('should validate both connected and validValences together', () => {
      const graph = MolGraph.createLinearAlkane(5);
      const result = graph.validate();
      expect(result.connected).toBe(true);
      expect(result.validValences).toBe(true);
    });
  });

  describe('ImplicitH Computation', () => {
    it('should compute correct implicitH for pentane', () => {
      const graph = MolGraph.createLinearAlkane(5);
      // Terminal carbons: 1 bond -> 3 H
      expect(graph.getAtom(0).implicitH).toBe(3);
      expect(graph.getAtom(4).implicitH).toBe(3);
      // Internal carbons: 2 bonds -> 2 H
      expect(graph.getAtom(1).implicitH).toBe(2);
      expect(graph.getAtom(2).implicitH).toBe(2);
      expect(graph.getAtom(3).implicitH).toBe(2);
    });
  });

  describe('setBond', () => {
    it('should set bond order and update implicitH', () => {
      const graph = MolGraph.createLinearAlkane(3);
      // Initially: C-C-C with single bonds
      expect(graph.getBondOrder(0, 1)).toBe(1);
      expect(graph.getAtom(0).implicitH).toBe(3);

      // Change to double bond
      graph.setBond(0, 1, 2);
      expect(graph.getBondOrder(0, 1)).toBe(2);
      expect(graph.getBondOrder(1, 0)).toBe(2); // Symmetric
      expect(graph.getAtom(0).implicitH).toBe(2); // 4 - 2 = 2
      expect(graph.getAtom(1).implicitH).toBe(1); // 4 - 3 = 1 (2 to atom 0, 1 to atom 2)
    });

    it('should throw for invalid bond order', () => {
      const graph = MolGraph.createLinearAlkane(2);
      expect(() => graph.setBond(0, 1, 4)).toThrow('out of range');
      expect(() => graph.setBond(0, 1, -1)).toThrow('out of range');
    });

    it('should throw for out of bounds indices', () => {
      const graph = MolGraph.createLinearAlkane(2);
      expect(() => graph.setBond(0, 5, 1)).toThrow('out of bounds');
    });
  });

  describe('Clone', () => {
    it('should create independent copy', () => {
      const original = MolGraph.createLinearAlkane(3);
      const clone = original.clone();

      // Modify clone
      clone.setBond(0, 1, 2);

      // Original should be unchanged
      expect(original.getBondOrder(0, 1)).toBe(1);
      expect(clone.getBondOrder(0, 1)).toBe(2);
    });

    it('should clone all atoms independently', () => {
      const original = MolGraph.createLinearAlkane(3);
      const clone = original.clone();

      const originalAtom0 = original.getAtom(0);
      const cloneAtom0 = clone.getAtom(0);

      expect(originalAtom0).toEqual(cloneAtom0);
      expect(originalAtom0).not.toBe(cloneAtom0); // Different objects
    });
  });

  describe('Getters', () => {
    it('should return readonly copies of atoms', () => {
      const graph = MolGraph.createLinearAlkane(2);
      const atoms = graph.getAtoms();

      // Modify returned array should not affect graph
      atoms[0]!.implicitH = 999;
      expect(graph.getAtom(0).implicitH).toBe(3);
    });

    it('should return readonly copy of bond matrix', () => {
      const graph = MolGraph.createLinearAlkane(2);
      const bonds = graph.getBondMatrix();

      // Modify returned matrix should not affect graph
      bonds[0]![1] = 999;
      expect(graph.getBondOrder(0, 1)).toBe(1);
    });

    it('should throw on invalid atom index', () => {
      const graph = MolGraph.createLinearAlkane(2);
      expect(() => graph.getAtom(5)).toThrow('out of bounds');
    });

    it('should throw on invalid bond indices', () => {
      const graph = MolGraph.createLinearAlkane(2);
      expect(() => graph.getBondOrder(0, 5)).toThrow('out of bounds');
    });
  });

  describe('toMolBlock', () => {
    it('should return empty string for empty graph', () => {
      const graph = new MolGraph([], []);
      expect(graph.toMolBlock()).toBe('');
    });

    it('should generate valid MOL block for methane (C1)', () => {
      const graph = MolGraph.createLinearAlkane(1);
      const molBlock = graph.toMolBlock();
      expect(molBlock).toContain('V2000');
      expect(molBlock).toContain('M  END');
      // 1 atom, 0 bonds
      expect(molBlock).toContain('  1  0  0  0  0  0  0  0  0  0999 V2000');
      // Atom line should contain C
      expect(molBlock).toContain(' C ');
    });

    it('should generate valid MOL block for ethane (C2)', () => {
      const graph = MolGraph.createLinearAlkane(2);
      const molBlock = graph.toMolBlock();
      // 2 atoms, 1 bond
      expect(molBlock).toContain('  2  1  0  0  0  0  0  0  0  0999 V2000');
      // Bond line: atom 1 to atom 2, single bond
      expect(molBlock).toContain('  1  2  1  0  0  0  0');
    });

    it('should generate correct atom and bond counts for n-hexane', () => {
      const graph = MolGraph.createLinearAlkane(6);
      const molBlock = graph.toMolBlock();
      // 6 atoms, 5 bonds
      expect(molBlock).toContain('  6  5  0  0  0  0  0  0  0  0999 V2000');
    });

    it('should generate correct bond count for cyclohexane (ring)', () => {
      const graph = MolGraph.createCyclohexane();
      const molBlock = graph.toMolBlock();
      // 6 atoms, 6 bonds (ring has same number of bonds as atoms)
      expect(molBlock).toContain('  6  6  0  0  0  0  0  0  0  0999 V2000');
    });

    it('should encode double bonds correctly (ethene)', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 2],
        [2, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      const molBlock = graph.toMolBlock();
      // Bond type 2 for double bond
      expect(molBlock).toContain('  1  2  2  0  0  0  0');
    });

    it('should encode triple bonds correctly (ethyne)', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 3],
        [3, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      const molBlock = graph.toMolBlock();
      // Bond type 3 for triple bond
      expect(molBlock).toContain('  1  2  3  0  0  0  0');
    });

    it('should include heteroatoms correctly (methanol)', () => {
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'O', implicitH: 0 },
      ];
      const bonds = [
        [0, 1],
        [1, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      const molBlock = graph.toMolBlock();
      expect(molBlock).toContain(' C ');
      expect(molBlock).toContain(' O ');
    });

    it('should use 1-based atom indices in bond block', () => {
      const graph = MolGraph.createBranched('isobutane');
      const molBlock = graph.toMolBlock();
      const lines = molBlock.split('\n');
      // Find bond lines (after atom block, before M  END)
      // Bond lines have exactly the pattern: 3-digit atom1, 3-digit atom2, 3-digit type, then "  0  0  0  0"
      const mEndIdx = lines.findIndex(l => l.trim() === 'M  END');
      const countsIdx = 3; // Line index 3 is the counts line
      // Bond lines are between atom block and M  END
      // Atom block starts at index 4, has 4 atoms, so bonds start at index 8
      const bondLines = lines.slice(countsIdx + 1 + 4, mEndIdx);
      // Isobutane: 3 bonds (0-1, 0-2, 0-3 → 1-2, 1-3, 1-4 in 1-based)
      expect(bondLines.length).toBe(3);
      // All bond indices should be >= 1 (1-based)
      for (const line of bondLines) {
        const parts = line.trim().split(/\s+/).map(Number);
        expect(parts[0]).toBeGreaterThanOrEqual(1);
        expect(parts[1]).toBeGreaterThanOrEqual(1);
      }
    });

    it('should produce consistent output for same graph', () => {
      const graph = MolGraph.createLinearAlkane(4);
      const molBlock1 = graph.toMolBlock();
      const molBlock2 = graph.toMolBlock();
      expect(molBlock1).toBe(molBlock2);
    });
  });
});
