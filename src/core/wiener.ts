import { MolGraph } from './MolGraph';

/**
 * Computes the Wiener Index of a molecular graph using BFS.
 *
 * The Wiener Index is the sum of all pairwise distances between heavy atoms
 * in the molecular graph. It's a topological descriptor used in QSAR studies.
 *
 * Algorithm: O(n^2) where n is the number of atoms.
 * For each atom i, run BFS to compute distances to all other atoms j > i,
 * and sum all distances.
 *
 * @param graph - The molecular graph
 * @returns The Wiener Index (sum of all pairwise distances)
 * @throws Error if graph is disconnected
 */
export function computeWienerIndex(graph: MolGraph): number {
  const n = graph.getAtomCount();

  // Edge case: single atom has Wiener index 0
  if (n <= 1) {
    return 0;
  }

  let totalSum = 0;

  // For each atom i, compute distances to all atoms j > i
  for (let i = 0; i < n; i++) {
    const distances = bfs(graph, i, n);

    // Check for disconnected graph
    for (let j = 0; j < n; j++) {
      if (distances[j] === Infinity) {
        throw new Error('Graph is disconnected');
      }
    }

    // Add distances to atoms j > i (to avoid double-counting)
    for (let j = i + 1; j < n; j++) {
      totalSum += distances[j]!;
    }
  }

  return totalSum;
}

/**
 * BFS to compute shortest path distances from a start atom to all other atoms.
 *
 * @param graph - The molecular graph
 * @param start - Starting atom index
 * @param n - Number of atoms (for array allocation)
 * @returns Array of distances where distances[i] is distance from start to atom i
 */
function bfs(graph: MolGraph, start: number, n: number): number[] {
  const distances = new Array(n).fill(Infinity);
  distances[start] = 0;

  // Array-based queue for O(1) operations
  const queue: number[] = [start];
  let head = 0; // Index of next element to dequeue

  while (head < queue.length) {
    const current = queue[head]!;
    head++;

    const currentDistance = distances[current]!;

    // Check all neighbors
    for (let neighbor = 0; neighbor < n; neighbor++) {
      // Bond order > 0 means connected (treat all bond orders as distance 1)
      if (graph.getBondOrder(current, neighbor) > 0 && distances[neighbor] === Infinity) {
        distances[neighbor] = currentDistance + 1;
        queue.push(neighbor);
      }
    }
  }

  return distances;
}
