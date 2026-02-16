/**
 * SVG Molecule Renderer
 *
 * Renders molecular structures from backend-generated SVG markup.
 * Replaces RDKit.js WASM canvas rendering with simple innerHTML-based approach.
 * Backend RDKit generates 2D coords and SVG, frontend just displays it.
 */

/**
 * Render a molecule from SVG markup to an HTML element.
 *
 * The backend generates SVG via RDKit, frontend displays via innerHTML.
 * Browser parses SVG natively - no WASM, no canvas operations needed.
 *
 * @param svg - SVG markup string from backend (StatusResponse.best_svg or SSEProgressData.best_svg)
 * @param container - The HTMLElement to render into
 */
export function renderMoleculeSVG(svg: string, container: HTMLElement): void {
  // Guard against missing inputs
  if (!svg || !container) {
    return;
  }

  // Set innerHTML - browser parses SVG markup natively
  // Backend RDKit generates safe SVG, no sanitization needed
  container.innerHTML = svg;
}

/**
 * Clear the molecule display and show placeholder.
 *
 * @param container - The HTMLElement to clear
 */
export function clearMoleculeDisplay(container: HTMLElement): void {
  // Guard against missing container
  if (!container) {
    return;
  }

  // Show centered gray placeholder
  container.innerHTML = `
    <div style="display: flex; align-items: center; justify-content: center; height: 100%;">
      <p style="color: #9ca3af; font-size: 14px; font-family: sans-serif;">No structure yet</p>
    </div>
  `;
}
