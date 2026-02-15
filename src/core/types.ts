export interface Atom {
  element: string;      // 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'
  implicitH: number;    // Computed from valence - (sum of bond orders)
}

export const STANDARD_VALENCE: Record<string, number> = {
  'C': 4, 'N': 3, 'O': 2, 'S': 2, 'P': 3,
  'F': 1, 'Cl': 1, 'Br': 1, 'I': 1
};

export type OptimizationMode = 'MAXIMIZE' | 'MINIMIZE';

export interface SAStepResult {
  step: number;
  currentEnergy: number;  // Current graph's Wiener Index
  bestEnergy: number;     // Best energy found so far
  temperature: number;    // Temperature at this step
  accepted: boolean;      // Whether the move was accepted
}

export interface SAEngineState {
  step: number;
  totalSteps: number;
  cycle: number;
  currentEnergy: number;
  bestEnergy: number;
  temperature: number;
  acceptedMoves: number;
  rejectedMoves: number;
  invalidMoves: number;
  isComplete: boolean;
}

/**
 * Simulated annealing parameters
 */
export interface SAParams {
  formula: string;              // Molecular formula (e.g., "C6H14")
  initialTemp: number;          // kT_0 (default: 100, per paper)
  coolingScheduleK: number;     // k parameter for cooling (default: 8)
  stepsPerCycle: number;        // Steps per temperature cycle (default: 500)
  numCycles: number;            // Number of cooling cycles (default: 4)
  optimizationMode: OptimizationMode; // 'MAXIMIZE' or 'MINIMIZE'
  seed: number;                 // Random seed for reproducibility
}

/**
 * Result of simulated annealing optimization
 */
export interface SAResult {
  bestGraph: import('./MolGraph').MolGraph;          // Graph with best energy found
  bestEnergy: number;           // Best Wiener Index found
  finalGraph: import('./MolGraph').MolGraph;         // Final graph state
  finalEnergy: number;          // Final Wiener Index
  initialEnergy: number;        // Initial Wiener Index
  totalSteps: number;           // Total SA steps executed
  acceptedMoves: number;        // Number of accepted moves
  rejectedMoves: number;        // Number of rejected moves (Metropolis)
  invalidMoves: number;         // Number of invalid displacement attempts
  acceptanceRatio: number;      // acceptedMoves / totalSteps
  history: SAStepResult[];      // Step-by-step results for charting
}
