import { describe, it, expect } from 'vitest';
import { computeWienerIndex } from '../wiener';
import { MolGraph } from '../MolGraph';
import { Atom } from '../types';

describe('Wiener Index', () => {
  describe('Known values for linear alkanes', () => {
    it('should compute 0 for methane (single atom)', () => {
      const graph = MolGraph.createLinearAlkane(1);
      expect(computeWienerIndex(graph)).toBe(0);
    });

    it('should compute 1 for ethane (C2)', () => {
      const graph = MolGraph.createLinearAlkane(2);
      expect(computeWienerIndex(graph)).toBe(1);
    });

    it('should compute 4 for propane (C3)', () => {
      const graph = MolGraph.createLinearAlkane(3);
      // Formula: n(n^2-1)/6 = 3*8/6 = 4
      expect(computeWienerIndex(graph)).toBe(4);
    });

    it('should compute 10 for n-butane (C4)', () => {
      const graph = MolGraph.createLinearAlkane(4);
      // Formula: 4*15/6 = 10
      expect(computeWienerIndex(graph)).toBe(10);
    });

    it('should compute 20 for n-pentane (C5)', () => {
      const graph = MolGraph.createLinearAlkane(5);
      // Formula: 5*24/6 = 20
      // This is the example from Faulon paper
      expect(computeWienerIndex(graph)).toBe(20);
    });

    it('should compute 35 for n-hexane (C6)', () => {
      const graph = MolGraph.createLinearAlkane(6);
      // Formula: 6*35/6 = 35
      expect(computeWienerIndex(graph)).toBe(35);
    });
  });

  describe('Known values for cyclic structures', () => {
    it('should compute 27 for cyclohexane', () => {
      const graph = MolGraph.createCyclohexane();
      // Pairs at distance 1: 6 (adjacent)
      // Pairs at distance 2: 6 (separated by 1)
      // Pairs at distance 3: 3 (opposite vertices)
      // Total: 6*1 + 6*2 + 3*3 = 6 + 12 + 9 = 27
      expect(computeWienerIndex(graph)).toBe(27);
    });
  });

  describe('Known values for branched structures', () => {
    it('should compute 9 for isobutane', () => {
      const graph = MolGraph.createBranched('isobutane');
      // Central C bonded to 3 others
      // 3 pairs at distance 1 (center to each terminal)
      // 3 pairs at distance 2 (terminal to terminal through center)
      // Total: 3*1 + 3*2 = 3 + 6 = 9
      expect(computeWienerIndex(graph)).toBe(9);
    });

    it('should compute 16 for neopentane', () => {
      const graph = MolGraph.createBranched('neopentane');
      // Central C bonded to 4 others
      // 4 pairs at distance 1 (center to each terminal)
      // 6 pairs at distance 2 (4 choose 2 = 6 terminal pairs through center)
      // Total: 4*1 + 6*2 = 4 + 12 = 16
      expect(computeWienerIndex(graph)).toBe(16);
    });
  });

  describe('Linear alkane formula verification', () => {
    // For linear alkanes, Wiener Index = n(n^2-1)/6
    const testCases = [
      { n: 7, expected: (7 * 48) / 6 },
      { n: 8, expected: (8 * 63) / 6 },
      { n: 10, expected: (10 * 99) / 6 },
    ];

    testCases.forEach(({ n, expected }) => {
      it(`should match formula for n=${n}: ${expected}`, () => {
        const graph = MolGraph.createLinearAlkane(n);
        expect(computeWienerIndex(graph)).toBe(expected);
      });
    });
  });

  describe('Performance', () => {
    it('should compute Wiener Index in under 5ms for 50-atom molecule', () => {
      const graph = MolGraph.createLinearAlkane(50);

      const start = performance.now();
      computeWienerIndex(graph);
      const end = performance.now();

      const duration = end - start;
      expect(duration).toBeLessThan(5);
    });

    it('should handle larger molecules efficiently', () => {
      const graph = MolGraph.createLinearAlkane(100);

      const start = performance.now();
      const result = computeWienerIndex(graph);
      const end = performance.now();

      expect(result).toBeGreaterThan(0);
      expect(end - start).toBeLessThan(20); // Should still be fast
    });
  });

  describe('Edge cases', () => {
    it('should throw error for disconnected graph', () => {
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

      expect(() => computeWienerIndex(graph)).toThrow('disconnected');
    });

    it('should handle two-atom graph correctly', () => {
      const graph = MolGraph.createLinearAlkane(2);
      // Only one pair at distance 1
      expect(computeWienerIndex(graph)).toBe(1);
    });
  });

  describe('Manual distance verification', () => {
    it('should correctly count distances in pentane', () => {
      // Pentane: C0-C1-C2-C3-C4
      const graph = MolGraph.createLinearAlkane(5);

      // Manually count distances:
      // From C0: d(0,1)=1, d(0,2)=2, d(0,3)=3, d(0,4)=4 -> sum = 10
      // From C1: d(1,2)=1, d(1,3)=2, d(1,4)=3 -> sum = 6
      // From C2: d(2,3)=1, d(2,4)=2 -> sum = 3
      // From C3: d(3,4)=1 -> sum = 1
      // Total: 10 + 6 + 3 + 1 = 20

      expect(computeWienerIndex(graph)).toBe(20);
    });

    it('should correctly count distances in isobutane', () => {
      // Isobutane: C1-C0-C2 with C0 also bonded to C3
      //            C3
      const graph = MolGraph.createBranched('isobutane');

      // Manually count distances:
      // From C0: d(0,1)=1, d(0,2)=1, d(0,3)=1 -> sum = 3
      // From C1: d(1,2)=2, d(1,3)=2 -> sum = 4
      // From C2: d(2,3)=2 -> sum = 2
      // Total: 3 + 4 + 2 = 9

      expect(computeWienerIndex(graph)).toBe(9);
    });
  });
});
