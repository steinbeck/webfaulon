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
import type { OptimizationMode, SAStepResult } from './types';

/**
 * Parameters for simulated annealing optimization
 */
export interface SAParams {
  formula: string;              // Molecular formula (e.g., "C6H14")
  initialTemp: number;          // kT_0 (default: 100, per paper)
  coolingScheduleK: number;     // k parameter for cooling (default: 8)
  stepsPerCycle: number;        // Steps per temperature cycle (default: 500)
  numCycles: number;            // Number of cooling cycles (default: 4)
  optimizationMode: OptimizationMode; // 'MAXIMIZE' or 'MINIMIZE'
  seed: number;                 // Random seed for reproducibility
}

/**
 * Result of simulated annealing optimization
 */
export interface SAResult {
  bestGraph: MolGraph;          // Graph with best energy found
  bestEnergy: number;           // Best Wiener Index found
  finalGraph: MolGraph;         // Final graph state
  finalEnergy: number;          // Final Wiener Index
  initialEnergy: number;        // Initial Wiener Index
  totalSteps: number;           // Total SA steps executed
  acceptedMoves: number;        // Number of accepted moves
  rejectedMoves: number;        // Number of rejected moves (Metropolis)
  invalidMoves: number;         // Number of invalid displacement attempts
  acceptanceRatio: number;      // acceptedMoves / totalSteps
  history: SAStepResult[];      // Step-by-step results for charting
}

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

  constructor(params: SAParams) {
    this.params = params;
    this.rng = new SeededRandom(params.seed);
  }

  /**
   * Run the full SA optimization to completion
   *
   * Algorithm (Faulon paper Table 2):
   * 1. Parse formula and generate initial structure
   * 2. Compute initial Wiener Index
   * 3. For each cycle:
   *    a. Reset temperature to initialTemp
   *    b. For each step:
   *       - Attempt displacement
   *       - Compute energy delta
   *       - Metropolis acceptance
   *       - Update current/best graphs
   *       - Record step result
   *       - Update temperature
   * 4. Return results
   */
  run(): SAResult {
    // Step 1: Initialize
    this.currentGraph = generateInitialStructure(this.params.formula);
    this.currentEnergy = computeWienerIndex(this.currentGraph);
    this.initialEnergy = this.currentEnergy;
    this.bestGraph = this.currentGraph.clone();
    this.bestEnergy = this.currentEnergy;

    // Step 2: Run SA cycles
    const totalSteps = this.params.stepsPerCycle * this.params.numCycles;
    let globalStep = 0;

    for (let cycle = 0; cycle < this.params.numCycles; cycle++) {
      for (let step = 0; step < this.params.stepsPerCycle; step++) {
        globalStep++;
        this.totalSteps = globalStep;

        // Compute current temperature
        const temperature = computeTemperature(
          globalStep - 1, // 0-indexed for temperature calculation
          totalSteps,
          this.params.initialTemp,
          this.params.coolingScheduleK
        );

        // Attempt displacement
        this.iterate(temperature, globalStep);
      }
    }

    // Step 3: Return results
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
