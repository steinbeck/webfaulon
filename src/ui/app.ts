import { validateFormula, type ValidationResult } from './validation';
import { PRESET_MOLECULES } from './presets';
import { createWienerChart, addChartDataPoint, resetChart } from './chart';
import { renderMoleculeSVG, clearMoleculeDisplay } from './molecule-renderer';
import { SAAPIClient } from '../api/client';
import { SSEConnection } from '../api/sse';
import type { SSEProgressData, SSECompleteData } from '../api/types';

// API client and SSE connection kept OUTSIDE Alpine reactive state
// (same pattern as v1.0 worker refs - avoid proxy wrapping issues)
const apiClient = new SAAPIClient();
const sseConnection = new SSEConnection();

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
    progress: null as SSEProgressData | null,
    result: null as SSECompleteData | null,
    sessionId: null as string | null,
    bestSMILES: '' as string,

    // Backend status
    backendError: '',

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
      return (this.progress.step / this.progress.total_steps) * 100;
    },

    // Lifecycle
    async init() {
      // Initialize chart and molecule display placeholder on canvas (after Alpine has rendered the DOM)
      // Use queueMicrotask to defer until after Alpine reactive updates
      queueMicrotask(() => {
        const chartCanvas = document.getElementById('wiener-chart') as HTMLCanvasElement;
        if (chartCanvas) {
          createWienerChart(chartCanvas);
        }
        // Initialize molecule display placeholder
        const molDisplay = document.getElementById('molecule-display');
        if (molDisplay) {
          clearMoleculeDisplay(molDisplay);
        }
      });
    },

    destroy() {
      sseConnection.close();
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
      if (!this.validation.valid) return;

      try {
        this.backendError = '';

        // 1. Configure session via REST API
        const configResponse = await apiClient.configure({
          formula: this.formula,
          initial_temp: this.initialTemp,
          cooling_schedule_k: this.coolingScheduleK,
          steps_per_cycle: this.stepsPerCycle,
          num_cycles: this.numCycles,
          optimization_mode: this.optimizationMode,
          seed: Date.now(),
        });
        this.sessionId = configResponse.session_id;

        // 2. Start SA execution via REST API
        await apiClient.start(this.sessionId);

        // 3. Reset chart for new run
        resetChart();

        // 4. Connect SSE for real-time progress
        const self = this;
        sseConnection.connect(this.sessionId, {
          onProgress(data: SSEProgressData) {
            self.progress = data;
            addChartDataPoint(data.step, data.best_energy);

            // Update molecule SVG
            const molDisplay = document.getElementById('molecule-display');
            if (molDisplay && data.best_svg) {
              renderMoleculeSVG(data.best_svg, molDisplay);
            }

            // Track best SMILES
            if (data.best_smiles) {
              self.bestSMILES = data.best_smiles;
            }
          },
          onComplete(data: SSECompleteData) {
            self.state = 'complete';
            self.result = data;

            // Render final molecule
            const molDisplay = document.getElementById('molecule-display');
            if (molDisplay && data.best_svg) {
              renderMoleculeSVG(data.best_svg, molDisplay);
            }
            if (data.best_smiles) {
              self.bestSMILES = data.best_smiles;
            }
          },
          onError(_event: Event) {
            console.error('SSE connection error');
            // Attempt to re-sync state from backend
            self.resyncState();
          },
        });

        this.state = 'running';
        this.progress = null;
        this.result = null;
      } catch (e: unknown) {
        console.error('Start failed:', e);
        this.backendError = e instanceof Error ? e.message : 'Failed to start SA';
        this.state = 'idle';
      }
    },

    async pause() {
      if (!this.sessionId || this.state !== 'running') return;
      try {
        await apiClient.pause(this.sessionId);
        this.state = 'paused';
      } catch (e: unknown) {
        console.error('Pause failed:', e);
      }
    },

    async resume() {
      if (!this.sessionId || this.state !== 'paused') return;
      try {
        await apiClient.start(this.sessionId);  // Backend treats start on paused session as resume
        this.state = 'running';
      } catch (e: unknown) {
        console.error('Resume failed:', e);
      }
    },

    async reset() {
      sseConnection.close();  // Stop SSE first
      if (this.sessionId) {
        try {
          await apiClient.reset(this.sessionId);
        } catch (e: unknown) {
          console.error('Reset failed:', e);
        }
      }
      resetChart();
      const molDisplay = document.getElementById('molecule-display');
      if (molDisplay) {
        clearMoleculeDisplay(molDisplay);
      }
      this.state = 'idle';
      this.progress = null;
      this.result = null;
      this.sessionId = null;
      this.bestSMILES = '';
      this.backendError = '';
    },

    async resyncState() {
      if (!this.sessionId) return;
      try {
        const status = await apiClient.getStatus(this.sessionId);
        this.state = status.session_state as any;
      } catch (e: unknown) {
        console.error('State resync failed:', e);
      }
    },
  };
}
