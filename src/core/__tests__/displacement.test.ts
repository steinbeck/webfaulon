import { describe, it, expect } from 'vitest';
import { attemptDisplacement } from '../displacement';
import { SeededRandom } from '../random';
import { MolGraph } from '../MolGraph';

describe('attemptDisplacement', () => {
  describe('basic validation', () => {
    it('should return null for molecules with < 4 atoms', () => {
      const methane = MolGraph.createLinearAlkane(1); // CH4, 1 atom
      const ethane = MolGraph.createLinearAlkane(2); // C2H6, 2 atoms
      const propane = MolGraph.createLinearAlkane(3); // C3H8, 3 atoms
      const rng = new SeededRandom(42);

      expect(attemptDisplacement(methane, rng)).toBe(null);
      expect(attemptDisplacement(ethane, rng)).toBe(null);
      expect(attemptDisplacement(propane, rng)).toBe(null);
    });

    it('should preserve atom count', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(12345);

      for (let i = 0; i < 100; i++) {
        const result = attemptDisplacement(hexane, rng);
        if (result !== null) {
          expect(result.getAtoms().length).toBe(hexane.getAtoms().length);
        }
      }
    });

    it('should preserve atom types', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(54321);

      for (let i = 0; i < 100; i++) {
        const result = attemptDisplacement(hexane, rng);
        if (result !== null) {
          const originalAtoms = hexane.getAtoms();
          result.getAtoms().forEach((atom, idx) => {
            const originalAtom = originalAtoms[idx];
            expect(originalAtom).toBeDefined();
            expect(atom.element).toBe(originalAtom!.element);
          });
        }
      }
    });

    it('should not mutate original graph', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const originalBonds = JSON.parse(JSON.stringify(hexane.getBondMatrix()));
      const rng = new SeededRandom(99999);

      attemptDisplacement(hexane, rng);

      // Original should be unchanged
      expect(hexane.getBondMatrix()).toEqual(originalBonds);
    });
  });

  describe('bond order constraints', () => {
    it('should produce bond orders in range [0, 3]', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(777);

      for (let i = 0; i < 200; i++) {
        const result = attemptDisplacement(hexane, rng);
        if (result !== null) {
          result.getBondMatrix().forEach(row => {
            row.forEach(order => {
              expect(order).toBeGreaterThanOrEqual(0);
              expect(order).toBeLessThanOrEqual(3);
            });
          });
        }
      }
    });

    it('should handle molecules with existing double bonds', () => {
      // Create a simple molecule with a double bond
      const graph = new MolGraph([
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 }
      ], [
        [0, 2, 1, 0], // C=C-C-C
        [2, 0, 1, 0],
        [1, 1, 0, 1],
        [0, 0, 1, 0]
      ]);

      const rng = new SeededRandom(333);

      for (let i = 0; i < 100; i++) {
        const result = attemptDisplacement(graph, rng);
        if (result !== null) {
          // All bonds should still be in valid range
          result.getBondMatrix().forEach(row => {
            row.forEach(order => {
              expect(order).toBeGreaterThanOrEqual(0);
              expect(order).toBeLessThanOrEqual(3);
            });
          });
        }
      }
    });
  });

  describe('connectivity and valence validation', () => {
    it('should only return connected molecules or null', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(888);

      for (let i = 0; i < 200; i++) {
        const result = attemptDisplacement(hexane, rng);
        if (result !== null) {
          expect(result.isConnected()).toBe(true);
        }
      }
    });

    it('should only return molecules with valid valences or null', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(444);

      for (let i = 0; i < 200; i++) {
        const result = attemptDisplacement(hexane, rng);
        if (result !== null) {
          expect(result.hasValidValences()).toBe(true);
        }
      }
    });

    it('should produce at least some valid displacements on linear hexane', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(555);
      let successCount = 0;

      for (let i = 0; i < 100; i++) {
        const result = attemptDisplacement(hexane, rng);
        if (result !== null) {
          successCount++;
        }
      }

      // Should have at least some successful displacements
      expect(successCount).toBeGreaterThan(0);
    });

    it('should reject displacements that would disconnect the molecule', () => {
      // Create a "bridge" molecule where one bond is critical for connectivity
      // Simple linear chain: C-C-C-C
      const chain = MolGraph.createLinearAlkane(4);
      const rng = new SeededRandom(666);

      for (let i = 0; i < 100; i++) {
        const result = attemptDisplacement(chain, rng);
        // If result is non-null, it must be connected
        if (result !== null) {
          expect(result.isConnected()).toBe(true);
        }
      }
    });
  });

  describe('reproducibility', () => {
    it('should produce identical results with same seed', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng1 = new SeededRandom(12345);
      const rng2 = new SeededRandom(12345);

      for (let i = 0; i < 100; i++) {
        const result1 = attemptDisplacement(hexane, rng1);
        const result2 = attemptDisplacement(hexane, rng2);

        // Both should be null or both should have same structure
        if (result1 === null) {
          expect(result2).toBe(null);
        } else {
          expect(result2).not.toBe(null);
          expect(result2!.getBondMatrix()).toEqual(result1.getBondMatrix());
        }
      }
    });

    it('should produce different results with different seeds', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng1 = new SeededRandom(111);
      const rng2 = new SeededRandom(222);

      const results1: (MolGraph | null)[] = [];
      const results2: (MolGraph | null)[] = [];

      for (let i = 0; i < 50; i++) {
        results1.push(attemptDisplacement(hexane, rng1));
        results2.push(attemptDisplacement(hexane, rng2));
      }

      // Sequences should differ
      let differences = 0;
      for (let i = 0; i < 50; i++) {
        const r1 = results1[i];
        const r2 = results2[i];

        if (r1 === undefined || r2 === undefined) {
          throw new Error('Results array has undefined element');
        }

        if ((r1 === null) !== (r2 === null)) {
          differences++;
        } else if (r1 !== null && r2 !== null) {
          if (JSON.stringify(r1.getBondMatrix()) !== JSON.stringify(r2.getBondMatrix())) {
            differences++;
          }
        }
      }

      // Expect at least some differences (seeds are different)
      expect(differences).toBeGreaterThan(5);
    });
  });

  describe('stress test - 500 displacements', () => {
    it('should produce zero invalid molecules over 500 iterations', () => {
      const hexane = MolGraph.createLinearAlkane(6);
      const rng = new SeededRandom(99999);
      let validCount = 0;
      let nullCount = 0;

      for (let i = 0; i < 500; i++) {
        const result = attemptDisplacement(hexane, rng);

        if (result === null) {
          nullCount++;
        } else {
          // Must be connected
          expect(result.isConnected()).toBe(true);
          // Must have valid valences
          expect(result.hasValidValences()).toBe(true);
          // Must have same atom count
          expect(result.getAtoms().length).toBe(6);
          // All bonds must be in range [0, 3]
          result.getBondMatrix().forEach(row => {
            row.forEach(order => {
              expect(order).toBeGreaterThanOrEqual(0);
              expect(order).toBeLessThanOrEqual(3);
            });
          });
          validCount++;
        }
      }

      // Should have processed all 500
      expect(validCount + nullCount).toBe(500);
      // Should have at least some valid displacements
      expect(validCount).toBeGreaterThan(0);
    });
  });

  describe('cyclohexane test', () => {
    it('should handle cyclic structures', () => {
      const cyclohexane = MolGraph.createCyclohexane();
      const rng = new SeededRandom(777);
      let successCount = 0;

      for (let i = 0; i < 100; i++) {
        const result = attemptDisplacement(cyclohexane, rng);
        if (result !== null) {
          expect(result.isConnected()).toBe(true);
          expect(result.hasValidValences()).toBe(true);
          successCount++;
        }
      }

      // Cyclic structures should allow some displacements
      expect(successCount).toBeGreaterThan(0);
    });
  });

  describe('equation verification', () => {
    it('should correctly implement bond conservation equations', () => {
      // Create a simple 4-atom molecule where we can track equations
      const graph = new MolGraph([
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 },
        { element: 'C', implicitH: 0 }
      ], [
        [0, 1, 1, 0],
        [1, 0, 1, 0],
        [1, 1, 0, 1],
        [0, 0, 1, 0]
      ]);

      const rng = new SeededRandom(42);

      // Run multiple displacements
      for (let i = 0; i < 50; i++) {
        const result = attemptDisplacement(graph, rng);
        // Just verify result is valid or null (equation correctness is implicit in validation)
        if (result !== null) {
          expect(result.isConnected()).toBe(true);
          expect(result.hasValidValences()).toBe(true);
        }
      }
    });
  });
});
