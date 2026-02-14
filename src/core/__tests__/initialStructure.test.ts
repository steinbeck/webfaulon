import { describe, it, expect } from 'vitest';
import { generateInitialStructure } from '../initialStructure';

describe('generateInitialStructure', () => {
  it('generates valid structure for saturated alkanes', () => {
    const graph = generateInitialStructure('C6H14'); // hexane

    // Check atom count
    expect(graph.getAtomCount()).toBe(6);

    // Check all atoms are carbon
    const atoms = graph.getAtoms();
    expect(atoms.every(a => a.element === 'C')).toBe(true);

    // Check connectivity
    expect(graph.isConnected()).toBe(true);

    // Check valid valences
    expect(graph.hasValidValences()).toBe(true);

    // Check total implicit H matches formula
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(14);
  });

  it('generates valid structure for butane', () => {
    const graph = generateInitialStructure('C4H10');

    expect(graph.getAtomCount()).toBe(4);
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(10);
  });

  it('generates valid structure for methane', () => {
    const graph = generateInitialStructure('CH4');

    expect(graph.getAtomCount()).toBe(1);
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    expect(atoms[0]!.implicitH).toBe(4);
  });

  it('generates structure with unsaturation (HDI=1)', () => {
    const graph = generateInitialStructure('C6H12');

    expect(graph.getAtomCount()).toBe(6);
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(12);

    // Should have at least one double bond
    const bonds = graph.getBondMatrix();
    let doubleBondCount = 0;
    for (let i = 0; i < bonds.length; i++) {
      for (let j = i + 1; j < bonds.length; j++) {
        if (bonds[i]![j] === 2) doubleBondCount++;
      }
    }
    expect(doubleBondCount).toBeGreaterThan(0);
  });

  it('generates structure with high unsaturation (HDI=4)', () => {
    const graph = generateInitialStructure('C6H6'); // benzene-like

    expect(graph.getAtomCount()).toBe(6);
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(6);

    // Should have multiple double bonds or triple bonds
    const bonds = graph.getBondMatrix();
    let unsaturationCount = 0;
    for (let i = 0; i < bonds.length; i++) {
      for (let j = i + 1; j < bonds.length; j++) {
        if (bonds[i]![j] === 2) unsaturationCount += 1;
        if (bonds[i]![j] === 3) unsaturationCount += 2;
      }
    }
    expect(unsaturationCount).toBe(4);
  });

  it('generates valid structure with heteroatoms (oxygen)', () => {
    const graph = generateInitialStructure('C2H6O'); // ethanol or dimethyl ether

    expect(graph.getAtomCount()).toBe(3); // 2 carbons + 1 oxygen
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const carbonCount = atoms.filter(a => a.element === 'C').length;
    const oxygenCount = atoms.filter(a => a.element === 'O').length;
    expect(carbonCount).toBe(2);
    expect(oxygenCount).toBe(1);

    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(6);
  });

  it('generates valid structure with nitrogen', () => {
    const graph = generateInitialStructure('C2H7N'); // ethylamine

    expect(graph.getAtomCount()).toBe(3); // 2 C + 1 N
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(7);
  });

  it('handles ethene (C2H4)', () => {
    const graph = generateInitialStructure('C2H4');

    expect(graph.getAtomCount()).toBe(2);
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(4);

    // Should have a double bond
    expect(graph.getBondOrder(0, 1)).toBe(2);
  });

  it('handles ethyne (C2H2)', () => {
    const graph = generateInitialStructure('C2H2');

    expect(graph.getAtomCount()).toBe(2);
    expect(graph.isConnected()).toBe(true);
    expect(graph.hasValidValences()).toBe(true);

    const atoms = graph.getAtoms();
    const totalH = atoms.reduce((sum, a) => sum + a.implicitH, 0);
    expect(totalH).toBe(2);

    // Should have a triple bond
    expect(graph.getBondOrder(0, 1)).toBe(3);
  });

  it('throws error for empty formula', () => {
    expect(() => generateInitialStructure('')).toThrow('Empty formula');
  });

  it('throws error for impossible hydrogen count', () => {
    // C6H99 is impossible with standard valences
    expect(() => generateInitialStructure('C6H99')).toThrow();
  });

  it('produces deterministic output', () => {
    const graph1 = generateInitialStructure('C6H14');
    const graph2 = generateInitialStructure('C6H14');

    // Should produce identical structures
    expect(graph1.getBondMatrix()).toEqual(graph2.getBondMatrix());
    expect(graph1.getAtoms()).toEqual(graph2.getAtoms());
  });

  it('handles formulas with only hydrogen', () => {
    const graph = generateInitialStructure('H2');
    expect(graph.getAtomCount()).toBe(0); // No heavy atoms, H2 is just implicit
  });
});
