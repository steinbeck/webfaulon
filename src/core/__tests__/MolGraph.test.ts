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

  describe('toSMILES', () => {
    it('should generate SMILES for methane (C1)', () => {
      const graph = MolGraph.createLinearAlkane(1);
      const smiles = graph.toSMILES();
      expect(smiles).toBe('C');
    });

    it('should generate SMILES for ethane (C2)', () => {
      const graph = MolGraph.createLinearAlkane(2);
      const smiles = graph.toSMILES();
      expect(smiles).toBe('CC');
    });

    it('should generate SMILES for n-hexane (C6)', () => {
      const graph = MolGraph.createLinearAlkane(6);
      const smiles = graph.toSMILES();
      expect(smiles).toBe('CCCCCC');
    });

    it('should generate SMILES with branch notation for isobutane', () => {
      const graph = MolGraph.createBranched('isobutane');
      const smiles = graph.toSMILES();
      // Valid SMILES for isobutane: CC(C)C or C(C)CC or C(C)(C)C
      expect(smiles).toBeTruthy();
      expect(smiles.length).toBeGreaterThan(0);
      // Should contain 4 C's
      const cCount = (smiles.match(/C/g) || []).length;
      expect(cCount).toBe(4);
      // Should contain parentheses for branching
      expect(smiles).toMatch(/\(/);
      expect(smiles).toMatch(/\)/);
    });

    it('should generate SMILES with branch notation for neopentane', () => {
      const graph = MolGraph.createBranched('neopentane');
      const smiles = graph.toSMILES();
      expect(smiles).toBeTruthy();
      expect(smiles.length).toBeGreaterThan(0);
      // Should contain 5 C's
      const cCount = (smiles.match(/C/g) || []).length;
      expect(cCount).toBe(5);
      // Should contain parentheses for branching
      expect(smiles).toMatch(/\(/);
      expect(smiles).toMatch(/\)/);
    });

    it('should generate SMILES with ring closure for cyclohexane', () => {
      const graph = MolGraph.createCyclohexane();
      const smiles = graph.toSMILES();
      expect(smiles).toBeTruthy();
      expect(smiles.length).toBeGreaterThan(0);
      // Should contain 6 C's
      const cCount = (smiles.match(/C/g) || []).length;
      expect(cCount).toBe(6);
      // Should contain ring closure digit (1-9)
      expect(smiles).toMatch(/\d/);
    });

    it('should generate SMILES with double bond notation (ethene)', () => {
      // Create ethene: C=C (2 carbons, bond order 2)
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 2],
        [2, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      const smiles = graph.toSMILES();
      expect(smiles).toBe('C=C');
    });

    it('should generate SMILES with triple bond notation (ethyne)', () => {
      // Create ethyne: C#C (2 carbons, bond order 3)
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
      ];
      const bonds = [
        [0, 3],
        [3, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      const smiles = graph.toSMILES();
      expect(smiles).toBe('C#C');
    });

    it('should generate SMILES with heteroatoms (methanol)', () => {
      // Create methanol: C-O
      const atoms: Atom[] = [
        { element: 'C', implicitH: 0 },
        { element: 'O', implicitH: 0 },
      ];
      const bonds = [
        [0, 1],
        [1, 0],
      ];
      const graph = new MolGraph(atoms, bonds);
      const smiles = graph.toSMILES();
      expect(smiles).toBeTruthy();
      expect(smiles.length).toBeGreaterThan(0);
      expect(smiles).toContain('C');
      expect(smiles).toContain('O');
    });

    it('should generate valid SMILES for linear alkanes (round-trip validation)', () => {
      // For linear alkanes n=3 to 8, verify SMILES is non-empty and has correct structure
      for (let n = 3; n <= 8; n++) {
        const graph = MolGraph.createLinearAlkane(n);
        const smiles = graph.toSMILES();
        expect(smiles).toBeTruthy();
        expect(smiles.length).toBe(n); // Linear alkane SMILES should be n C's
        const cCount = (smiles.match(/C/g) || []).length;
        expect(cCount).toBe(n);
      }
    });
  });
});
