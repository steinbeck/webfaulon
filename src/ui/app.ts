import * as Comlink from 'comlink';
import { validateFormula, type ValidationResult } from './validation';
import { PRESET_MOLECULES, type PresetMolecule } from './presets';
import type { ISAWorker, SAProgressData } from '../worker/types';
import type { SAParams, SAResult } from '../core/types';

export interface AppState {
  // Formula input
  formula: string;
  validation: ValidationResult;
  selectedPreset: number | null;
  presets: PresetMolecule[];

  // SA parameters
  initialTemp: number;
  coolingScheduleK: number;
  stepsPerCycle: number;
  numCycles: number;
  optimizationMode: 'MAXIMIZE' | 'MINIMIZE';

  // Execution state
  state: 'idle' | 'running' | 'paused' | 'complete';
  progress: SAProgressData | null;
  result: SAResult | null;

  // RDKit status
  rdkitReady: boolean;
  rdkitError: string;

  // Worker reference
  worker: Comlink.Remote<ISAWorker> | null;

  // Computed getters
  canStart: boolean;
  canPause: boolean;
  canResume: boolean;
  totalSteps: number;
  progressPercent: number;

  // Methods
  init(): void;
  selectPreset(index: number | string): void;
  validateInput(): void;
  start(): Promise<void>;
  pause(): void;
  resume(): void;
  reset(): Promise<void>;
  onProgress(data: SAProgressData): void;
}

export function appComponent(): AppState {
  return {
    // Initial state
    formula: '',
    validation: { valid: false, error: '', formula: '' },
    selectedPreset: null,
    presets: PRESET_MOLECULES,

    initialTemp: 100,
    coolingScheduleK: 8,
    stepsPerCycle: 500,
    numCycles: 4,
    optimizationMode: 'MINIMIZE',

    state: 'idle',
    progress: null,
    result: null,

    rdkitReady: false,
    rdkitError: '',

    worker: null,

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
      } catch (e: unknown) {
        const msg = e instanceof Error ? e.message : 'Failed to load RDKit.js';
        this.rdkitError = msg;
        console.error('RDKit initialization failed:', e);
        // RDKit failure is non-blocking -- SA can still run without rendering
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

      // Initialize worker
      const WorkerClass = (await import('../worker/sa-worker?worker')).default;
      const worker = new WorkerClass();
      this.worker = Comlink.wrap<ISAWorker>(worker);

      // Create SA parameters
      const params: SAParams = {
        formula: this.formula,
        initialTemp: this.initialTemp,
        coolingScheduleK: this.coolingScheduleK,
        stepsPerCycle: this.stepsPerCycle,
        numCycles: this.numCycles,
        optimizationMode: this.optimizationMode,
        seed: Date.now(), // Generate seed from timestamp
      };

      // Initialize SA engine
      await this.worker.initialize(params);

      // Set state to running
      this.state = 'running';
      this.progress = null;
      this.result = null;

      try {
        // Run SA with progress callback
        const result = await this.worker.run(
          Comlink.proxy((data: SAProgressData) => this.onProgress(data)),
          10 // Report every 10 steps
        );

        this.result = result;
        this.state = 'complete';
      } catch (e: unknown) {
        console.error('SA execution failed:', e);
        this.state = 'idle';
        alert('SA execution failed: ' + (e instanceof Error ? e.message : 'Unknown error'));
      }
    },

    pause() {
      if (this.worker && this.state === 'running') {
        this.worker.pause();
        this.state = 'paused';
      }
    },

    resume() {
      if (this.worker && this.state === 'paused') {
        this.worker.resume();
        this.state = 'running';
      }
    },

    async reset() {
      if (this.worker) {
        await this.worker.reset();
        // Terminate worker
        const raw = (this.worker as any)[Comlink.proxyMarker];
        if (raw) {
          raw.terminate();
        }
        this.worker = null;
      }

      this.state = 'idle';
      this.progress = null;
      this.result = null;
    },

    onProgress(data: SAProgressData) {
      this.progress = data;

      // If worker reports completion, update state
      if (data.isComplete && this.state === 'running') {
        // State will be set to 'complete' when run() promise resolves
      }
    },
  };
}
