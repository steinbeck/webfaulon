import * as Comlink from 'comlink';
import { SAEngine } from '../core/SAEngine';
import type { SAParams, SAResult } from '../core/types';
import type { ISAWorker, SAProgressData } from './types';

class SAWorker implements ISAWorker {
  private engine: SAEngine | null = null;
  private _isPaused = false;

  async initialize(params: SAParams): Promise<void> {
    this.engine = new SAEngine(params);
    this.engine.init();
    this._isPaused = false;
  }

  async run(
    onProgress: (data: SAProgressData) => void,
    reportInterval = 10
  ): Promise<SAResult> {
    if (!this.engine) {
      throw new Error('Engine not initialized. Call initialize() first.');
    }

    while (!this.engine.getState().isComplete) {
      // Check for pause
      if (this._isPaused) {
        // Yield to event loop, then re-check pause flag
        await new Promise<void>(resolve => setTimeout(resolve, 0));
        continue;
      }

      this.engine.step();

      const currentState = this.engine.getState();

      // Report progress at interval
      if (currentState.step % reportInterval === 0 || currentState.isComplete) {
        const progressData: SAProgressData = {
          step: currentState.step,
          totalSteps: currentState.totalSteps,
          cycle: currentState.cycle,
          currentEnergy: currentState.currentEnergy,
          bestEnergy: currentState.bestEnergy,
          temperature: currentState.temperature,
          accepted: currentState.step > 0, // Approximate - use last step's accepted
          acceptedMoves: currentState.acceptedMoves,
          rejectedMoves: currentState.rejectedMoves,
          invalidMoves: currentState.invalidMoves,
          isComplete: currentState.isComplete,
        };
        await onProgress(progressData);
      }

      // Yield to event loop periodically (every 100 steps) to allow message processing
      // This ensures pause/resume messages can be received
      if (currentState.step % 100 === 0) {
        await new Promise<void>(resolve => setTimeout(resolve, 0));
      }
    }

    return this.engine.getResult();
  }

  pause(): void {
    this._isPaused = true;
  }

  resume(): void {
    this._isPaused = false;
  }

  isPaused(): boolean {
    return this._isPaused;
  }

  reset(): void {
    this.engine = null;
    this._isPaused = false;
  }
}

Comlink.expose(new SAWorker());
