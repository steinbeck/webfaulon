import { describe, it, expect } from 'vitest';
import { SeededRandom } from '../random';

describe('SeededRandom', () => {
  describe('constructor and determinism', () => {
    it('should produce the same first value for the same seed', () => {
      const rng1 = new SeededRandom(42);
      const rng2 = new SeededRandom(42);

      expect(rng1.next()).toBe(rng2.next());
    });

    it('should produce identical sequences of 1000 values for same seed', () => {
      const rng1 = new SeededRandom(12345);
      const rng2 = new SeededRandom(12345);

      for (let i = 0; i < 1000; i++) {
        expect(rng1.next()).toBe(rng2.next());
      }
    });

    it('should produce different sequences for different seeds', () => {
      const rng1 = new SeededRandom(111);
      const rng2 = new SeededRandom(222);

      const sequence1 = Array.from({ length: 100 }, () => rng1.next());
      const sequence2 = Array.from({ length: 100 }, () => rng2.next());

      // At least some values should be different
      const differences = sequence1.filter((val, idx) => val !== sequence2[idx]);
      expect(differences.length).toBeGreaterThan(90); // Expect most to be different
    });
  });

  describe('next()', () => {
    it('should return values in range [0, 1)', () => {
      const rng = new SeededRandom(999);

      for (let i = 0; i < 10000; i++) {
        const val = rng.next();
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(1);
      }
    });

    it('should produce different values on subsequent calls', () => {
      const rng = new SeededRandom(777);
      const val1 = rng.next();
      const val2 = rng.next();
      const val3 = rng.next();

      expect(val1).not.toBe(val2);
      expect(val2).not.toBe(val3);
      expect(val1).not.toBe(val3);
    });
  });

  describe('nextInt()', () => {
    it('should return integers in range [min, max] inclusive', () => {
      const rng = new SeededRandom(555);

      for (let i = 0; i < 1000; i++) {
        const val = rng.nextInt(0, 5);
        expect(Number.isInteger(val)).toBe(true);
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThanOrEqual(5);
      }
    });

    it('should cover entire range over many calls', () => {
      const rng = new SeededRandom(888);
      const seen = new Set<number>();

      for (let i = 0; i < 1000; i++) {
        seen.add(rng.nextInt(0, 5));
      }

      // Should see all values 0-5
      expect(seen.size).toBe(6);
      expect(seen.has(0)).toBe(true);
      expect(seen.has(5)).toBe(true);
    });

    it('should handle range with same min and max', () => {
      const rng = new SeededRandom(444);

      for (let i = 0; i < 100; i++) {
        expect(rng.nextInt(7, 7)).toBe(7);
      }
    });

    it('should handle negative ranges', () => {
      const rng = new SeededRandom(333);

      for (let i = 0; i < 1000; i++) {
        const val = rng.nextInt(-10, -5);
        expect(val).toBeGreaterThanOrEqual(-10);
        expect(val).toBeLessThanOrEqual(-5);
      }
    });
  });

  describe('selectNDistinct()', () => {
    it('should return exactly n distinct values', () => {
      const rng = new SeededRandom(666);
      const selected = rng.selectNDistinct(4, 10);

      expect(selected.length).toBe(4);
      const unique = new Set(selected);
      expect(unique.size).toBe(4);
    });

    it('should return values in range [0, range)', () => {
      const rng = new SeededRandom(123);
      const selected = rng.selectNDistinct(5, 20);

      selected.forEach(val => {
        expect(val).toBeGreaterThanOrEqual(0);
        expect(val).toBeLessThan(20);
      });
    });

    it('should work when n equals range', () => {
      const rng = new SeededRandom(456);
      const selected = rng.selectNDistinct(10, 10);

      expect(selected.length).toBe(10);
      const unique = new Set(selected);
      expect(unique.size).toBe(10);

      // Should be all values 0-9
      for (let i = 0; i < 10; i++) {
        expect(selected).toContain(i);
      }
    });

    it('should produce different selections for different seeds', () => {
      const rng1 = new SeededRandom(111);
      const rng2 = new SeededRandom(222);

      const sel1 = rng1.selectNDistinct(4, 100);
      const sel2 = rng2.selectNDistinct(4, 100);

      // Arrays should not be identical
      expect(sel1).not.toEqual(sel2);
    });

    it('should produce same selection for same seed', () => {
      const rng1 = new SeededRandom(789);
      const rng2 = new SeededRandom(789);

      const sel1 = rng1.selectNDistinct(6, 50);
      const sel2 = rng2.selectNDistinct(6, 50);

      expect(sel1).toEqual(sel2);
    });
  });
});
