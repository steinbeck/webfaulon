/**
 * Seeded pseudo-random number generator for reproducible SA runs.
 *
 * Uses Mulberry32 algorithm - fast, good distribution, deterministic.
 * Faulon 1996 p.733: "Two runs using the same seed will lead to the same results"
 */
export class SeededRandom {
  private state: number;

  constructor(seed: number) {
    // Initialize state with seed
    this.state = seed | 0;
  }

  /**
   * Generate next random number in range [0, 1).
   * Mulberry32 algorithm by Tommy Ettinger.
   */
  next(): number {
    let t = (this.state += 0x6d2b79f5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  /**
   * Generate random integer in range [min, max] inclusive.
   */
  nextInt(min: number, max: number): number {
    const range = max - min + 1;
    return Math.floor(this.next() * range) + min;
  }

  /**
   * Select n distinct random integers from range [0, range).
   * Used for selecting 4 distinct atoms for displacement.
   */
  selectNDistinct(n: number, range: number): number[] {
    if (n > range) {
      throw new Error(`Cannot select ${n} distinct values from range ${range}`);
    }

    const selected: number[] = [];
    const available = new Set<number>();

    // Initialize available set
    for (let i = 0; i < range; i++) {
      available.add(i);
    }

    // Select n distinct values
    for (let i = 0; i < n; i++) {
      const availableArray = Array.from(available);
      const idx = this.nextInt(0, availableArray.length - 1);
      const value = availableArray[idx];
      if (value === undefined) {
        throw new Error('Internal error: selected undefined value');
      }
      selected.push(value);
      available.delete(value);
    }

    return selected;
  }
}
