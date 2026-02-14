import { describe, it, expect } from 'vitest';
import { parseFormula, computeHDI } from '../formulaParser';

describe('parseFormula', () => {
  it('parses simple alkane formulas', () => {
    expect(parseFormula('C6H14')).toEqual({ C: 6, H: 14 });
    expect(parseFormula('C8H10')).toEqual({ C: 8, H: 10 });
    expect(parseFormula('CH4')).toEqual({ C: 1, H: 4 });
    expect(parseFormula('C4H10')).toEqual({ C: 4, H: 10 });
  });

  it('handles implicit count of 1', () => {
    expect(parseFormula('CH4')).toEqual({ C: 1, H: 4 });
    expect(parseFormula('CO2')).toEqual({ C: 1, O: 2 });
  });

  it('parses formulas with heteroatoms', () => {
    expect(parseFormula('C2H6O')).toEqual({ C: 2, H: 6, O: 1 });
    expect(parseFormula('C6H5Cl')).toEqual({ C: 6, H: 5, Cl: 1 });
    expect(parseFormula('C3H7N')).toEqual({ C: 3, H: 7, N: 1 });
    expect(parseFormula('C2H6S')).toEqual({ C: 2, H: 6, S: 1 });
  });

  it('parses complex formulas with multiple heteroatoms', () => {
    expect(parseFormula('C2H5NO2')).toEqual({ C: 2, H: 5, N: 1, O: 2 });
    expect(parseFormula('C3H9N3')).toEqual({ C: 3, H: 9, N: 3 });
  });

  it('throws error on empty formula', () => {
    expect(() => parseFormula('')).toThrow('Empty formula');
  });

  it('throws error on unknown elements', () => {
    expect(() => parseFormula('XYZ')).toThrow('Unknown element: X');
    expect(() => parseFormula('C6Q14')).toThrow('Unknown element: Q');
  });

  it('handles formulas with double-digit counts', () => {
    expect(parseFormula('C10H22')).toEqual({ C: 10, H: 22 });
    expect(parseFormula('C12H26')).toEqual({ C: 12, H: 26 });
  });

  it('preserves order independence', () => {
    // Different orderings should produce same result (normalized)
    const result1 = parseFormula('C6H14');
    const result2 = parseFormula('H14C6');
    expect(result1).toEqual({ C: 6, H: 14 });
    expect(result2).toEqual({ C: 6, H: 14 });
  });
});

describe('computeHDI', () => {
  it('computes HDI for saturated alkanes', () => {
    expect(computeHDI({ C: 6, H: 14 })).toBe(0); // hexane
    expect(computeHDI({ C: 4, H: 10 })).toBe(0); // butane
    expect(computeHDI({ C: 1, H: 4 })).toBe(0);  // methane
  });

  it('computes HDI for unsaturated hydrocarbons', () => {
    expect(computeHDI({ C: 6, H: 12 })).toBe(1); // one ring or double bond
    expect(computeHDI({ C: 6, H: 6 })).toBe(4);  // benzene
    expect(computeHDI({ C: 2, H: 4 })).toBe(1);  // ethene
    expect(computeHDI({ C: 2, H: 2 })).toBe(2);  // ethyne
  });

  it('computes HDI for molecules with oxygen', () => {
    expect(computeHDI({ C: 2, H: 6, O: 1 })).toBe(0); // ethanol
    expect(computeHDI({ C: 2, H: 4, O: 1 })).toBe(1); // acetaldehyde
  });

  it('computes HDI for molecules with nitrogen', () => {
    // HDI = (2C + 2 + N - H - Halogens) / 2
    expect(computeHDI({ C: 2, H: 7, N: 1 })).toBe(0); // ethylamine (C2H5NH2)
    expect(computeHDI({ C: 2, H: 5, N: 1 })).toBe(1); // acetonitrile (CH3CN)
  });

  it('computes HDI for molecules with halogens', () => {
    expect(computeHDI({ C: 6, H: 5, Cl: 1 })).toBe(4); // chlorobenzene
    expect(computeHDI({ C: 2, H: 5, Br: 1 })).toBe(0); // bromoethane
  });

  it('handles multiple heteroatoms', () => {
    expect(computeHDI({ C: 2, H: 5, N: 1, O: 2 })).toBe(1); // glycine
  });

  it('returns 0 for edge cases', () => {
    expect(computeHDI({ C: 1, H: 4 })).toBe(0); // methane
    expect(computeHDI({ H: 2 })).toBe(0); // H2
  });
});
