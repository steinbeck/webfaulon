import type { SAParams, SAResult } from '../core/types';

/**
 * Progress data sent from worker to main thread during SA execution.
 * Kept small (primitives only) to minimize postMessage overhead.
 */
export interface SAProgressData {
  step: number;
  totalSteps: number;
  cycle: number;
  currentEnergy: number;
  bestEnergy: number;
  bestMolBlock: string;
  temperature: number;
  accepted: boolean;         // Whether the last move was accepted
  acceptedMoves: number;
  rejectedMoves: number;
  invalidMoves: number;
  isComplete: boolean;
}

/**
 * Worker API interface exposed via Comlink.
 * All methods are async from the main thread's perspective.
 */
export interface ISAWorker {
  /**
   * Initialize SA engine with parameters and begin step-by-step execution.
   * Must be called before run().
   */
  initialize(params: SAParams): Promise<void>;

  /**
   * Run SA to completion with progress callbacks.
   * Calls onProgress every `reportInterval` steps.
   * Respects pause/resume signals.
   *
   * @param onProgress - Callback receiving progress data (use Comlink.proxy)
   * @param reportInterval - Report progress every N steps (default: 10)
   */
  run(
    onProgress: (data: SAProgressData) => void,
    reportInterval?: number
  ): Promise<SAResult>;

  /**
   * Pause execution (takes effect after current step completes).
   */
  pause(): void;

  /**
   * Resume paused execution.
   */
  resume(): void;

  /**
   * Check if currently paused.
   */
  isPaused(): boolean;

  /**
   * Reset engine to uninitialized state.
   * Allows reuse of the worker for a new run.
   */
  reset(): void;
}
