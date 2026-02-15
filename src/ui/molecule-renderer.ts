/**
 * RDKit.js 2D Molecule Renderer
 *
 * Renders molecular structures from V2000 MOL blocks to canvas elements.
 * Uses RDKit.js to parse MOL blocks, generate 2D coords, and render.
 * Also derives canonical SMILES from RDKit for display.
 */

// Module-level state to track last rendered MOL block
let _lastRenderedMolBlock: string | null = null;

/**
 * Render a molecule from a V2000 MOL block to canvas using RDKit.js.
 * Also returns canonical SMILES derived by RDKit.
 *
 * @param molBlock - V2000 MOL block string (from MolGraph.toMolBlock())
 * @param canvas - The HTMLCanvasElement to render on
 * @returns canonical SMILES string if render succeeded, null if skipped or failed
 */
export function renderMolecule(molBlock: string, canvas: HTMLCanvasElement): string | null {
  // Skip redundant re-render if MOL block hasn't changed
  if (molBlock === _lastRenderedMolBlock) {
    return null;
  }

  // Get RDKit from global (loaded via script tag)
  const RDKit = (window as any).__rdkit;
  if (!RDKit) {
    console.warn('RDKit.js not loaded - cannot render molecule');
    return null;
  }

  let mol = null;

  try {
    mol = RDKit.get_mol(molBlock);

    if (!mol || !mol.is_valid()) {
      if (mol) { mol.delete(); mol = null; }
      throw new Error('RDKit could not parse MOL block');
    }

    // Get canonical SMILES from RDKit
    const smiles: string = mol.get_smiles();

    // Generate 2D coordinates for rendering (MOL block has all-zero coords
    // since MolGraph stores topology only, not geometry)
    mol.set_new_coords();

    // Clear canvas before rendering
    const ctx = canvas.getContext('2d');
    if (ctx) {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
    }

    // Render molecule to canvas (-1, -1 = use full canvas dimensions)
    mol.draw_to_canvas(canvas, -1, -1);

    _lastRenderedMolBlock = molBlock;
    return smiles;
  } catch (error) {
    const ctx = canvas.getContext('2d');
    if (ctx) {
      ctx.clearRect(0, 0, canvas.width, canvas.height);
      ctx.fillStyle = 'rgb(239, 68, 68)';
      ctx.font = '14px sans-serif';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      const errMsg = error instanceof Error ? error.message : String(error);
      ctx.fillText('Could not render molecule', canvas.width / 2, canvas.height / 2 - 10);
      ctx.font = '11px monospace';
      ctx.fillStyle = 'rgb(180, 80, 80)';
      ctx.fillText(errMsg.slice(0, 60), canvas.width / 2, canvas.height / 2 + 10);
    }
    console.error('[mol-renderer] render failed:', error);
    return null;
  } finally {
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
  const ctx = canvas.getContext('2d');
  if (!ctx) {
    return;
  }

  ctx.clearRect(0, 0, canvas.width, canvas.height);
  _lastRenderedMolBlock = null;

  ctx.fillStyle = 'rgb(156, 163, 175)';
  ctx.font = '14px sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'middle';
  ctx.fillText('No structure yet', canvas.width / 2, canvas.height / 2);
}
