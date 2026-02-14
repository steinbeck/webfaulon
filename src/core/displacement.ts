/**
 * Faulon displacement operation (equations 7-11 from 1996 paper).
 *
 * The core SA mutation: select 4 atoms, redistribute bond orders
 * while preserving valences and connectivity.
 *
 * Reference: J. Chem. Inf. Comput. Sci. 1996, 36, 731-740, page 733
 */

import { MolGraph } from './MolGraph';
import { SeededRandom } from './random';

/**
 * Attempt a Faulon displacement on the given graph.
 *
 * @param graph - The molecular graph to displace
 * @param rng - Seeded random number generator for reproducibility
 * @returns A new graph with bonds redistributed, or null if move is invalid
 */
export function attemptDisplacement(
  graph: MolGraph,
  rng: SeededRandom
): MolGraph | null {
  const atomCount = graph.getAtomCount();

  // Need at least 4 atoms for displacement
  if (atomCount < 4) {
    return null;
  }

  // Select 4 distinct atoms (Faulon paper p.733, Table 2, step 2.1)
  const selected = rng.selectNDistinct(4, atomCount);
  const x1 = selected[0]!;
  const y1 = selected[1]!;
  const x2 = selected[2]!;
  const y2 = selected[3]!;

  // Read current bond orders (a notation from paper)
  // Paper p.733: "Let a11, a12, a21, and a22 be the orders of
  // the bonds [x1,y1], [x1,y2], [x2,y1], and [x2,y2]"
  const a11 = graph.getBondOrder(x1, y1);
  const a12 = graph.getBondOrder(y1, y2);
  const a21 = graph.getBondOrder(x1, x2);
  const a22 = graph.getBondOrder(x2, y2);

  // Compute valid range for b11 using equations 10 and 11
  // Equation 10: b11 >= MAX(0, a11-a22, a11-a12, a11+a12-3, a11+a21-3)
  // Equation 11: b11 <= MIN(3, a11+a12, a11+a21, a11-a22+3)
  const b11_min = Math.max(
    0,
    a11 - a22,
    a11 - a12,
    a11 + a12 - 3,
    a11 + a21 - 3
  );

  const b11_max = Math.min(
    3,
    a11 + a12,
    a11 + a21,
    a11 - a22 + 3
  );

  // If no valid range, displacement is impossible for this atom selection
  if (b11_min > b11_max) {
    return null;
  }

  // Choose b11 randomly from valid range
  const b11 = rng.nextInt(b11_min, b11_max);

  // Compute remaining new bond orders using equations 7-9
  // Equation 7: b12 = a11 + a12 - b11
  // Equation 8: b21 = a11 + a21 - b11
  // Equation 9: b22 = a22 - a11 + b11
  const b12 = a11 + a12 - b11;
  const b21 = a11 + a21 - b11;
  const b22 = a22 - a11 + b11;

  // Verify all bond orders are in valid range [0, 3]
  // (Should be guaranteed by constraint equations, but verify)
  if (
    b11 < 0 || b11 > 3 ||
    b12 < 0 || b12 > 3 ||
    b21 < 0 || b21 > 3 ||
    b22 < 0 || b22 > 3
  ) {
    throw new Error(
      `Bond order out of range after displacement: b11=${b11}, b12=${b12}, b21=${b21}, b22=${b22}. ` +
      `This indicates a bug in the displacement equations.`
    );
  }

  // Clone the graph (don't mutate original)
  const newGraph = graph.clone();

  // Apply new bond orders
  newGraph.setBond(x1, y1, b11);
  newGraph.setBond(y1, y2, b12);
  newGraph.setBond(x1, x2, b21);
  newGraph.setBond(x2, y2, b22);

  // Validate result: must be connected and have valid valences
  if (!newGraph.isConnected()) {
    return null;
  }

  if (!newGraph.hasValidValences()) {
    return null;
  }

  return newGraph;
}
