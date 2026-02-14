export type FormulaMap = Record<string, number>;

const KNOWN_ELEMENTS = new Set([
  'H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'
]);

const HALOGENS = new Set(['F', 'Cl', 'Br', 'I']);

/**
 * Parse a molecular formula string into a map of element symbols to counts.
 *
 * Examples:
 * - "C6H14" -> {C: 6, H: 14}
 * - "CH4" -> {C: 1, H: 4}
 * - "C2H6O" -> {C: 2, H: 6, O: 1}
 *
 * @param formula - Molecular formula string (e.g., "C6H14")
 * @returns Map of element symbols to counts
 * @throws Error if formula is empty or contains unknown elements
 */
export function parseFormula(formula: string): FormulaMap {
  if (!formula || formula.trim() === '') {
    throw new Error('Empty formula');
  }

  const result: FormulaMap = {};

  // Regex to match element symbol (capital + optional lowercase) followed by optional digits
  // Example: C6H14 -> ["C6", "H14"]
  // Example: CH4 -> ["C", "H4"]
  const pattern = /([A-Z][a-z]?)(\d*)/g;
  let match: RegExpExecArray | null;

  while ((match = pattern.exec(formula)) !== null) {
    const element = match[1];
    const countStr = match[2];

    // Skip empty matches (happens at end of string)
    if (!element) continue;

    // Validate element is known
    if (!KNOWN_ELEMENTS.has(element)) {
      throw new Error(`Unknown element: ${element}`);
    }

    // Parse count (default to 1 if not specified)
    const count = countStr ? parseInt(countStr, 10) : 1;

    // Add to result (accumulate if element appears multiple times)
    result[element] = (result[element] || 0) + count;
  }

  return result;
}

/**
 * Compute the Hydrogen Deficiency Index (HDI) for a molecular formula.
 *
 * HDI = (2C + 2 + N - H - Halogens) / 2
 *
 * HDI indicates the number of rings and/or multiple bonds:
 * - HDI = 0: saturated (no rings or double bonds)
 * - HDI = 1: one ring or one double bond
 * - HDI = 2: two rings, one triple bond, or one ring + one double bond
 * - etc.
 *
 * Examples:
 * - C6H14 (hexane): HDI = 0
 * - C6H12 (cyclohexane or 1-hexene): HDI = 1
 * - C6H6 (benzene): HDI = 4
 *
 * @param formula - Parsed formula map
 * @returns Hydrogen Deficiency Index
 */
export function computeHDI(formula: FormulaMap): number {
  const C = formula.C || 0;
  const H = formula.H || 0;
  const N = formula.N || 0;

  // Count halogens
  let halogens = 0;
  for (const element of HALOGENS) {
    halogens += formula[element] || 0;
  }

  // HDI = (2C + 2 + N - H - Halogens) / 2
  // Oxygen and sulfur don't affect HDI (they're divalent like C in the formula)
  const hdi = (2 * C + 2 + N - H - halogens) / 2;

  return Math.max(0, hdi); // HDI can't be negative
}
