/**
 * RDKit.js 2D Molecule Renderer
 *
 * Provides functions for rendering molecular structures to canvas elements.
 * Handles WASM memory management and prevents redundant re-renders.
 */

// Module-level state to track last rendered SMILES
let _lastRenderedSMILES: string | null = null;

/**
 * Render a molecule from SMILES string to canvas using RDKit.js.
 * Skips redundant re-renders if SMILES hasn't changed.
 *
 * @param smiles - The SMILES string to render
 * @param canvas - The HTMLCanvasElement to render on
 * @returns true if render succeeded, false if skipped or failed
 */
export function renderMolecule(smiles: string, canvas: HTMLCanvasElement): boolean {
  // Skip redundant re-render if SMILES hasn't changed
  if (smiles === _lastRenderedSMILES) {
    return false;
  }

  // Get RDKit from global (loaded via script tag)
  const RDKit = (window as any).__rdkit;
  if (!RDKit) {
    console.warn('RDKit.js not loaded - cannot render molecule');
    return false;
  }

  let mol = null;

  try {
    // Parse SMILES to RDKit molecule
    mol = RDKit.get_mol(smiles);

    // Validate molecule
    if (!mol || !mol.is_valid()) {
      throw new Error('Invalid molecule structure');
    }

    // Clear canvas before rendering
    const ctx = canvas.getContext('2d');
    if (ctx) {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    // Render molecule to canvas
    // -1, -1 = use full canvas dimensions
    mol.draw_to_canvas(canvas, -1, -1);

    // Update last rendered SMILES
    _lastRenderedSMILES = smiles;

    return true;
  } catch (error) {
    // Render error fallback on canvas
    const ctx = canvas.getContext('2d');
    if (ctx) {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = 'rgb(239, 68, 68)'; // red-500
      ctx.font = '14px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillText('Could not render molecule', canvas.width / 2, canvas.height / 2);
    }

    console.error('Failed to render molecule:', error);
    return false;
  } finally {
    // CRITICAL: Delete molecule to prevent WASM memory leak
    if (mol) {
      mol.delete();
    }
  }
}

/**
 * Clear the molecule canvas and show placeholder text.
 *
 * @param canvas - The HTMLCanvasElement to clear
 */
export function clearMoleculeCanvas(canvas: HTMLCanvasElement): void {
  // Get 2D context
  const ctx = canvas.getContext('2d');
  if (!ctx) {
    return;
  }

  // Clear entire canvas
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // Reset last rendered SMILES
  _lastRenderedSMILES = null;

  // Draw placeholder text
  ctx.fillStyle = 'rgb(156, 163, 175)'; // gray-400
  ctx.font = '14px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('No structure yet', canvas.width / 2, canvas.height / 2);
}
