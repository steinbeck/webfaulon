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
