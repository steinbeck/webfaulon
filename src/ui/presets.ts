export interface PresetMolecule {
  name: string;           // Display name
  formula: string;        // Molecular formula
  description: string;    // Brief chemistry context
}

/**
 * Curated preset molecules for users who don't know chemistry.
 * Selected to show different SA behaviors:
 * - Small saturated: fast convergence
 * - Medium saturated: moderate exploration
 * - Aromatic: interesting isomer space
 * - Larger: demonstrates computation time
 */
export const PRESET_MOLECULES: PresetMolecule[] = [
  {
    name: 'Hexane isomers',
    formula: 'C6H14',
    description: 'Saturated 6-carbon chain. 5 constitutional isomers.'
  },
  {
    name: 'Octane isomers',
    formula: 'C8H18',
    description: 'Saturated 8-carbon chain. 18 constitutional isomers.'
  },
  {
    name: 'Ethylbenzene / Xylenes',
    formula: 'C8H10',
    description: 'Aromatic 8-carbon compound. Includes ring structures.'
  },
  {
    name: 'Decane isomers',
    formula: 'C10H22',
    description: 'Larger saturated chain. 75 constitutional isomers. Longer SA runs.'
  },
  {
    name: 'Naphthalene isomers',
    formula: 'C10H8',
    description: 'Highly unsaturated. Fused ring systems possible.'
  },
];
