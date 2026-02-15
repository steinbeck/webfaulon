import * as Comlink from 'comlink';
import { validateFormula, type ValidationResult } from './validation';
import { PRESET_MOLECULES } from './presets';
import type { ISAWorker, SAProgressData } from '../worker/types';
import type { SAParams, SAResult } from '../core/types';
import { createWienerChart, addChartDataPoint, resetChart } from './chart';
import { renderMolecule, clearMoleculeCanvas } from './molecule-renderer';

// Worker references kept OUTSIDE Alpine reactive state to avoid
// Proxy-of-Proxy issues (Alpine Proxy wrapping Comlink Remote Proxy
// causes infinite recursion / "Maximum call stack size exceeded")
let _rawWorker: Worker | null = null;
let _saWorker: Comlink.Remote<ISAWorker> | null = null;

function createWorker(): Comlink.Remote<ISAWorker> {
  _rawWorker = new Worker(
    new URL('../worker/sa-worker.ts', import.meta.url),
    { type: 'module' }
  );
  _saWorker = Comlink.wrap<ISAWorker>(_rawWorker);
  return _saWorker;
}

function destroyWorker(): void {
  if (_rawWorker) {
    _rawWorker.terminate();
    _rawWorker = null;
    _saWorker = null;
  }
}

export function appComponent() {
  return {
    // Formula input
    formula: '',
    validation: { valid: false, error: '', formula: '' } as ValidationResult,
    selectedPreset: null as number | null,
    presets: PRESET_MOLECULES,

    // SA parameters
    initialTemp: 100,
    coolingScheduleK: 8,
    stepsPerCycle: 500,
    numCycles: 4,
    optimizationMode: 'MINIMIZE' as 'MAXIMIZE' | 'MINIMIZE',

    // Execution state
    state: 'idle' as 'idle' | 'running' | 'paused' | 'complete',
    progress: null as SAProgressData | null,
    result: null as SAResult | null,
    _prevBestEnergy: null as number | null,
    bestSMILES: '' as string,

    // RDKit status
    rdkitReady: false,
    rdkitError: '',

    // Computed properties
    get canStart(): boolean {
      return this.state === 'idle' && this.validation.valid;
    },

    get canPause(): boolean {
      return this.state === 'running';
    },

    get canResume(): boolean {
      return this.state === 'paused';
    },

    get totalSteps(): number {
      return this.stepsPerCycle * this.numCycles;
    },

    get progressPercent(): number {
      if (!this.progress) return 0;
      return (this.progress.step / this.progress.totalSteps) * 100;
    },

    // Lifecycle
    async init() {
      // Load RDKit.js WASM via global initRDKitModule (loaded by script tag in index.html)
      try {
        if (typeof window.initRDKitModule !== 'function') {
          throw new Error('RDKit script not loaded. Check index.html script tag.');
        }
        const RDKit = await window.initRDKitModule();
        this.rdkitReady = true;

        // Smoke test: create a molecule from SMILES
        const mol = RDKit.get_mol('CCCCCC'); // n-hexane
        if (mol) {
          console.log('RDKit smoke test passed: n-hexane loaded');
          mol.delete(); // Clean up WASM memory
        }

        // Store RDKit instance for later use (Phase 3 rendering)
        (window as any).__rdkit = RDKit;

        // Initialize chart on canvas (after Alpine has rendered the DOM)
        // Use queueMicrotask to defer until after Alpine reactive updates
        queueMicrotask(() => {
          const chartCanvas = document.getElementById('wiener-chart') as HTMLCanvasElement;
          if (chartCanvas) {
            createWienerChart(chartCanvas);
          }
          const molCanvas = document.getElementById('molecule-canvas') as HTMLCanvasElement;
          if (molCanvas) {
            clearMoleculeCanvas(molCanvas);
          }
        });
      } catch (e: unknown) {
        const msg = e instanceof Error ? e.message : 'Failed to load RDKit.js';
        this.rdkitError = msg;
        console.error('RDKit initialization failed:', e);
      }
    },

    // Methods
    selectPreset(index: number | string) {
      const idx = typeof index === 'string' ? parseInt(index, 10) : index;
      if (isNaN(idx) || idx < 0 || idx >= this.presets.length) {
        this.selectedPreset = null;
        return;
      }

      const preset = this.presets[idx];
      if (!preset) return;

      this.selectedPreset = idx;
      this.formula = preset.formula;
      this.validateInput();
    },

    validateInput() {
      this.validation = validateFormula(this.formula);
    },

    async start() {
      if (!this.validation.valid) {
        return;
      }

      try {
        // Create worker (references stored outside Alpine reactive scope)
        const worker = createWorker();

        // Create SA parameters
        const params: SAParams = {
          formula: this.formula,
          initialTemp: this.initialTemp,
          coolingScheduleK: this.coolingScheduleK,
          stepsPerCycle: this.stepsPerCycle,
          numCycles: this.numCycles,
          optimizationMode: this.optimizationMode,
          seed: Date.now(),
        };

        // Initialize SA engine in worker
        await worker.initialize(params);

        // Set state to running
        this.state = 'running';
        this.progress = null;
        this.result = null;

        // Capture `this` for use in callback (Alpine proxy context)
        const self = this;

        // Run SA with progress callback
        const result = await worker.run(
          Comlink.proxy((data: SAProgressData) => {
            self.progress = { ...data } as SAProgressData; // Plain copy to avoid proxy issues

            // Update chart with new data point
            addChartDataPoint(data.step, data.bestEnergy);

            // Re-render molecule only when best energy changes (new best structure found)
            if (data.bestMolBlock && data.bestEnergy !== self._prevBestEnergy) {
              const molCanvas = document.getElementById('molecule-canvas') as HTMLCanvasElement;
              if (molCanvas) {
                const smiles = renderMolecule(data.bestMolBlock, molCanvas);
                if (smiles) {
                  self.bestSMILES = smiles;
                }
              }
              self._prevBestEnergy = data.bestEnergy;
            }
          }),
          10 // Report every 10 steps
        );

        // Store result as plain object (strip MolGraph class instances)
        this.result = {
          bestEnergy: result.bestEnergy,
          finalEnergy: result.finalEnergy,
          initialEnergy: result.initialEnergy,
          totalSteps: result.totalSteps,
          acceptedMoves: result.acceptedMoves,
          rejectedMoves: result.rejectedMoves,
          invalidMoves: result.invalidMoves,
          acceptanceRatio: result.acceptanceRatio,
        } as SAResult;
        this.state = 'complete';

        // Render final best molecule
        const progressData = this.progress as SAProgressData | null;
        if (progressData && progressData.bestMolBlock) {
          const molCanvas = document.getElementById('molecule-canvas') as HTMLCanvasElement;
          if (molCanvas) {
            const smiles = renderMolecule(progressData.bestMolBlock, molCanvas);
            if (smiles) {
              this.bestSMILES = smiles;
            }
          }
        }
      } catch (e: unknown) {
        console.error('SA execution failed:', e);
        this.state = 'idle';
        alert('SA execution failed: ' + (e instanceof Error ? e.message : 'Unknown error'));
      }
    },

    pause() {
      if (_saWorker && this.state === 'running') {
        _saWorker.pause();
        this.state = 'paused';
      }
    },

    resume() {
      if (_saWorker && this.state === 'paused') {
        _saWorker.resume();
        this.state = 'running';
      }
    },

    async reset() {
      destroyWorker();
      resetChart();
      this._prevBestEnergy = null;
      this.bestSMILES = '';
      const molCanvas = document.getElementById('molecule-canvas') as HTMLCanvasElement;
      if (molCanvas) {
        clearMoleculeCanvas(molCanvas);
      }
      this.state = 'idle';
      this.progress = null;
      this.result = null;
    },
  };
}
