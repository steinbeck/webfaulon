import { MolGraph } from './MolGraph';
import { parseFormula, computeHDI } from './formulaParser';
import { Atom, STANDARD_VALENCE } from './types';

/**
 * Generate a deterministic initial molecular structure from a molecular formula.
 *
 * Algorithm:
 * 1. Parse formula to get heavy atom counts (C, N, O, etc.)
 * 2. Create array of heavy atoms (carbons first, then heteroatoms alphabetically)
 * 3. Connect all heavy atoms in a linear chain with single bonds
 * 4. Compute HDI from formula
 * 5. Add unsaturation (double/triple bonds) to satisfy HDI
 * 6. Validate that implicit H matches formula H count
 *
 * This produces a valid, connected starting structure for SA optimization.
 * The specific arrangement doesn't matter - SA will rearrange via displacement.
 *
 * @param formula - Molecular formula string (e.g., "C6H14")
 * @returns Connected MolGraph with valid valences
 * @throws Error if formula is invalid or impossible to satisfy
 */
export function generateInitialStructure(formula: string): MolGraph {
  // Parse formula
  const formulaMap = parseFormula(formula);
  const targetH = formulaMap.H || 0;
  const hdi = computeHDI(formulaMap);

  // Extract heavy atoms (all except H)
  const heavyAtoms: Atom[] = [];

  // Add atoms in deterministic order: C, N, O, S, P, then halogens
  const elementOrder = ['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'];

  for (const element of elementOrder) {
    const count = formulaMap[element] || 0;
    for (let i = 0; i < count; i++) {
      heavyAtoms.push({ element, implicitH: 0 });
    }
  }

  // Handle edge case: no heavy atoms (e.g., H2)
  if (heavyAtoms.length === 0) {
    // Return empty graph (all hydrogens are implicit)
    return new MolGraph([], []);
  }

  // Handle edge case: single atom
  if (heavyAtoms.length === 1) {
    const bonds: number[][] = [[0]];
    const graph = new MolGraph(heavyAtoms, bonds);

    // Verify hydrogen count
    const actualH = graph.getAtom(0).implicitH;
    if (actualH !== targetH) {
      throw new Error(
        `Cannot generate structure for ${formula}: expected ${targetH} H, got ${actualH}`
      );
    }

    return graph;
  }

  // Create bond matrix: n x n, all zeros
  const n = heavyAtoms.length;
  const bonds: number[][] = Array(n)
    .fill(0)
    .map(() => Array(n).fill(0));

  // Step 1: Connect atoms in a linear chain with single bonds
  for (let i = 0; i < n - 1; i++) {
    bonds[i]![i + 1] = 1;
    bonds[i + 1]![i] = 1;
  }

  // Step 2: Add unsaturation to satisfy HDI
  // Each double bond adds 1 HDI, each triple bond adds 2 HDI
  // Strategy: upgrade bonds iteratively until HDI is satisfied
  let remainingHDI = hdi;

  // Keep upgrading bonds until we satisfy HDI or can't upgrade any more
  while (remainingHDI > 0) {
    let upgraded = false;

    // Try to upgrade each bond in the chain
    for (let i = 0; i < n - 1; i++) {
      if (remainingHDI === 0) break;

      const j = i + 1;
      const currentBondOrder = bonds[i]![j]!;

      // Can't upgrade beyond triple bond
      if (currentBondOrder >= 3) continue;

      // Calculate bond order sums for both atoms
      const getBondOrderSum = (atomIndex: number): number => {
        let sum = 0;
        for (let k = 0; k < n; k++) {
          sum += bonds[atomIndex]![k]!;
        }
        return sum;
      };

      const atom1 = heavyAtoms[i]!;
      const atom2 = heavyAtoms[j]!;

      const valence1 = STANDARD_VALENCE[atom1.element];
      const valence2 = STANDARD_VALENCE[atom2.element];

      if (valence1 === undefined || valence2 === undefined) {
        throw new Error(`Unknown element valence`);
      }

      // Try upgrading by 1 (e.g., 1->2 or 2->3)
      const newBondOrder = currentBondOrder + 1;

      // Calculate what the bond sums would be after upgrade
      const currentBondSum1 = getBondOrderSum(i);
      const currentBondSum2 = getBondOrderSum(j);

      const newBondSum1 = currentBondSum1 + 1; // Adding 1 to the bond
      const newBondSum2 = currentBondSum2 + 1;

      const implicitH1 = valence1 - newBondSum1;
      const implicitH2 = valence2 - newBondSum2;

      // Check if both atoms can handle the upgrade
      if (implicitH1 >= 0 && implicitH2 >= 0) {
        bonds[i]![j] = newBondOrder;
        bonds[j]![i] = newBondOrder;
        remainingHDI -= 1;
        upgraded = true;
      }
    }

    // If we couldn't upgrade any bond, we can't satisfy HDI
    if (!upgraded) {
      break;
    }
  }

  // If we couldn't satisfy HDI, the formula might be impossible
  if (remainingHDI > 0) {
    throw new Error(
      `Cannot generate structure for ${formula}: unable to satisfy HDI=${hdi}`
    );
  }

  // Create the graph
  const graph = new MolGraph(heavyAtoms, bonds);

  // Step 3: Verify hydrogen count
  const atoms = graph.getAtoms();
  const actualH = atoms.reduce((sum, a) => sum + a.implicitH, 0);

  if (actualH !== targetH) {
    throw new Error(
      `Cannot generate structure for ${formula}: expected ${targetH} H, got ${actualH}`
    );
  }

  // Step 4: Validate structure
  if (!graph.isConnected()) {
    throw new Error(`Generated structure is not connected`);
  }

  if (!graph.hasValidValences()) {
    throw new Error(`Generated structure has invalid valences`);
  }

  return graph;
}
