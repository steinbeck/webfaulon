/**
 * Simulated Annealing Engine for molecular graph optimization
 *
 * Implements the SA algorithm from Faulon's 1996 paper (Table 2).
 * Optimizes Wiener Index by iteratively applying displacement operations
 * with Metropolis acceptance criterion.
 */

import { MolGraph } from './MolGraph';
import { generateInitialStructure } from './initialStructure';
import { computeWienerIndex } from './wiener';
import { attemptDisplacement } from './displacement';
import { SeededRandom } from './random';
import { computeTemperature } from './cooling';
import type { SAStepResult, SAEngineState, SAParams, SAResult } from './types';

/**
 * Simulated Annealing Engine
 *
 * Usage:
 * ```typescript
 * const engine = new SAEngine({
 *   formula: 'C6H14',
 *   initialTemp: 100,
 *   coolingScheduleK: 8,
 *   stepsPerCycle: 500,
 *   numCycles: 4,
 *   optimizationMode: 'MINIMIZE',
 *   seed: 42
 * });
 * const result = engine.run();
 * ```
 */
export class SAEngine {
  private params: SAParams;
  private rng: SeededRandom;
  private currentGraph!: MolGraph;
  private currentEnergy!: number;
  private bestGraph!: MolGraph;
  private bestEnergy!: number;
  private initialEnergy!: number;
  private acceptedMoves = 0;
  private rejectedMoves = 0;
  private invalidMoves = 0;
  private history: SAStepResult[] = [];
  private totalSteps = 0;
  private globalStep = 0;
  private currentTemperature = 0;
  private initialized = false;
  private completed = false;

  constructor(params: SAParams) {
    this.params = params;
    this.rng = new SeededRandom(params.seed);
  }

  /**
   * Initialize the SA engine for step-by-step execution
   *
   * Sets up initial structure, computes initial energy, and prepares state.
   * Must be called before step().
   */
  init(): void {
    // Generate initial structure
    this.currentGraph = generateInitialStructure(this.params.formula);
    this.currentEnergy = computeWienerIndex(this.currentGraph);
    this.initialEnergy = this.currentEnergy;
    this.bestGraph = this.currentGraph.clone();
    this.bestEnergy = this.currentEnergy;

    // Reset state
    this.totalSteps = this.params.stepsPerCycle * this.params.numCycles;
    this.globalStep = 0;
    this.acceptedMoves = 0;
    this.rejectedMoves = 0;
    this.invalidMoves = 0;
    this.history = [];
    this.currentTemperature = this.params.initialTemp;
    this.initialized = true;
    this.completed = false;
  }

  /**
   * Execute a single SA iteration
   *
   * Advances the algorithm by one step. Can be called repeatedly with
   * arbitrary delays between calls (enabling pause/resume).
   *
   * @throws Error if init() hasn't been called or execution is already complete
   */
  step(): void {
    if (!this.initialized) {
      throw new Error('SAEngine.step() called before init()');
    }

    if (this.completed) {
      throw new Error('SAEngine.step() called after completion');
    }

    // Increment step counter
    this.globalStep++;

    // Compute current temperature
    this.currentTemperature = computeTemperature(
      this.globalStep - 1, // 0-indexed for temperature calculation
      this.totalSteps,
      this.params.initialTemp,
      this.params.coolingScheduleK
    );

    // Execute one SA iteration
    this.iterate(this.currentTemperature, this.globalStep);

    // Check if complete
    if (this.globalStep === this.totalSteps) {
      this.completed = true;
    }
  }

  /**
   * Get current execution state
   *
   * Returns a snapshot of the current SA execution state.
   *
   * @throws Error if init() hasn't been called
   * @returns Current state including step number, energy values, and completion status
   */
  getState(): SAEngineState {
    if (!this.initialized) {
      throw new Error('SAEngine.getState() called before init()');
    }

    const cycle = this.globalStep === 0
      ? 0
      : Math.floor((this.globalStep - 1) / this.params.stepsPerCycle) + 1;

    return {
      step: this.globalStep,
      totalSteps: this.totalSteps,
      cycle,
      currentEnergy: this.currentEnergy,
      bestEnergy: this.bestEnergy,
      bestSMILES: this.bestGraph.toSMILES(),
      temperature: this.currentTemperature,
      acceptedMoves: this.acceptedMoves,
      rejectedMoves: this.rejectedMoves,
      invalidMoves: this.invalidMoves,
      isComplete: this.completed,
    };
  }

  /**
   * Get final optimization result
   *
   * Returns the complete SAResult with best/final graphs and statistics.
   * Can only be called after all steps are complete.
   *
   * @throws Error if execution is not complete
   * @returns Final optimization result
   */
  getResult(): SAResult {
    if (!this.completed) {
      throw new Error('SAEngine.getResult() called before completion');
    }

    return {
      bestGraph: this.bestGraph,
      bestEnergy: this.bestEnergy,
      finalGraph: this.currentGraph,
      finalEnergy: this.currentEnergy,
      initialEnergy: this.initialEnergy,
      totalSteps: this.totalSteps,
      acceptedMoves: this.acceptedMoves,
      rejectedMoves: this.rejectedMoves,
      invalidMoves: this.invalidMoves,
      acceptanceRatio: this.acceptedMoves / this.totalSteps,
      history: this.history,
    };
  }

  /**
   * Run the full SA optimization to completion
   *
   * Convenience method that delegates to init(), step(), and getResult().
   * Maintains backward compatibility with existing code.
   *
   * @returns Complete optimization result
   */
  run(): SAResult {
    this.init();
    while (!this.completed) {
      this.step();
    }
    return this.getResult();
  }

  /**
   * Execute a single SA iteration
   *
   * @param temperature Current temperature (kT)
   * @param stepNumber Current step number (1-indexed for display)
   */
  private iterate(temperature: number, stepNumber: number): void {
    // Attempt displacement
    const proposedGraph = attemptDisplacement(this.currentGraph, this.rng);

    // If displacement returned null, it's invalid
    if (proposedGraph === null) {
      this.invalidMoves++;
      this.recordStep(stepNumber, temperature, false);
      return;
    }

    // Compute energy of proposed graph
    const proposedEnergy = computeWienerIndex(proposedGraph);

    // Compute energy delta based on optimization mode
    // For MINIMIZE: positive delta = worsening move
    // For MAXIMIZE: positive delta = worsening move
    const deltaE =
      this.params.optimizationMode === 'MINIMIZE'
        ? proposedEnergy - this.currentEnergy
        : this.currentEnergy - proposedEnergy;

    // Metropolis acceptance criterion
    const accepted = this.metropolisAccept(deltaE, temperature);

    if (accepted) {
      this.acceptedMoves++;
      this.currentGraph = proposedGraph;
      this.currentEnergy = proposedEnergy;

      // Update best if this is better
      if (this.isBetter(this.currentEnergy, this.bestEnergy)) {
        this.bestGraph = this.currentGraph.clone();
        this.bestEnergy = this.currentEnergy;
      }
    } else {
      this.rejectedMoves++;
    }

    this.recordStep(stepNumber, temperature, accepted);
  }

  /**
   * Check if proposed energy is better than current best based on optimization mode
   *
   * @param proposedEnergy Energy to test
   * @param currentBest Current best energy
   * @returns true if proposed is better than current best
   */
  private isBetter(proposedEnergy: number, currentBest: number): boolean {
    return this.params.optimizationMode === 'MINIMIZE'
      ? proposedEnergy < currentBest
      : proposedEnergy > currentBest;
  }

  /**
   * Metropolis acceptance criterion (ALG-05)
   *
   * Always accepts improving moves (delta_e <= 0)
   * Probabilistically accepts worsening moves: P = exp(-delta_e / temperature)
   *
   * @param deltaE Energy change (positive = worsening for both modes)
   * @param temperature Current temperature (kT)
   * @returns true if move should be accepted
   */
  private metropolisAccept(deltaE: number, temperature: number): boolean {
    // Always accept improving moves
    if (deltaE <= 0) {
      return true;
    }

    // Probabilistically accept worsening moves
    const probability = Math.exp(-deltaE / temperature);
    return this.rng.next() < probability;
  }

  /**
   * Record step result in history array
   *
   * @param stepNumber Step number (1-indexed)
   * @param temperature Current temperature
   * @param accepted Whether move was accepted
   */
  private recordStep(
    stepNumber: number,
    temperature: number,
    accepted: boolean
  ): void {
    this.history.push({
      step: stepNumber,
      currentEnergy: this.currentEnergy,
      bestEnergy: this.bestEnergy,
      temperature,
      accepted,
    });
  }
}
