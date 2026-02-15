import { parseFormula } from '../core/formulaParser';

export interface ValidationResult {
  valid: boolean;
  error: string;
  formula: string;
}

/**
 * Two-stage validation:
 * 1. Regex format check (fast, catches obviously wrong input)
 * 2. parseFormula() chemical validity (catches impossible compositions)
 */
export function validateFormula(formula: string): ValidationResult {
  const trimmed = formula.trim();

  if (!trimmed) {
    return { valid: false, error: 'Please enter a molecular formula', formula: trimmed };
  }

  // Stage 1: Format check
  const formatRegex = /^([A-Z][a-z]?\d*)+$/;
  if (!formatRegex.test(trimmed)) {
    return { valid: false, error: 'Invalid format. Use format like C6H14 or C8H10.', formula: trimmed };
  }

  // Stage 2: Chemical validity
  try {
    const parsed = parseFormula(trimmed);

    // Must contain at least carbon
    if (!parsed.C || parsed.C < 1) {
      return { valid: false, error: 'Formula must contain at least one carbon atom.', formula: trimmed };
    }

    // Check HDI is non-negative (already handled by computeHDI clamping, but check raw value)
    const C = parsed.C || 0;
    const H = parsed.H || 0;
    const N = parsed.N || 0;
    let halogens = 0;
    for (const el of ['F', 'Cl', 'Br', 'I']) {
      halogens += parsed[el] || 0;
    }
    const rawHDI = (2 * C + 2 + N - H - halogens) / 2;
    if (rawHDI < 0) {
      return { valid: false, error: `Too many hydrogens for this formula (HDI = ${rawHDI}).`, formula: trimmed };
    }

    // Must have at least 4 heavy atoms for displacement to work
    // (need to pick x1, y1, x2, y2 -- 4 distinct atoms)
    const heavyAtoms = Object.entries(parsed).reduce(
      (sum, [el, count]) => el !== 'H' ? sum + count : sum, 0
    );
    if (heavyAtoms < 4) {
      return { valid: false, error: 'Need at least 4 heavy atoms for SA displacement (e.g., C4H10).', formula: trimmed };
    }

    return { valid: true, error: '', formula: trimmed };
  } catch (e: unknown) {
    const msg = e instanceof Error ? e.message : 'Invalid formula';
    return { valid: false, error: msg, formula: trimmed };
  }
}
