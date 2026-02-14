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
  energy: number;        // Wiener Index value
  bestEnergy: number;
  temperature: number;
  accepted: boolean;
}
