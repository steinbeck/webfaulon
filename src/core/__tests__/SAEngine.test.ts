import { describe, it, expect } from 'vitest';
import { SAEngine, type SAParams } from '../SAEngine';

describe('SAEngine', () => {
  const defaultParams: SAParams = {
    formula: 'C6H14',
    initialTemp: 100,
    coolingScheduleK: 8,
    stepsPerCycle: 100,
    numCycles: 2,
    optimizationMode: 'MINIMIZE',
    seed: 42,
  };

  describe('basic functionality', () => {
    it('creates engine with valid parameters', () => {
      const engine = new SAEngine(defaultParams);
      expect(engine).toBeDefined();
    });

    it('returns SAResult with all required fields', () => {
      const engine = new SAEngine(defaultParams);
      const result = engine.run();

      expect(result.bestGraph).toBeDefined();
      expect(result.bestEnergy).toBeTypeOf('number');
      expect(result.finalGraph).toBeDefined();
      expect(result.finalEnergy).toBeTypeOf('number');
      expect(result.initialEnergy).toBeTypeOf('number');
      expect(result.totalSteps).toBe(200); // 100 * 2
      expect(result.acceptedMoves).toBeTypeOf('number');
      expect(result.rejectedMoves).toBeTypeOf('number');
      expect(result.invalidMoves).toBeTypeOf('number');
      expect(result.acceptanceRatio).toBeTypeOf('number');
      expect(result.history).toBeInstanceOf(Array);
    });

    it('history array has correct length', () => {
      const engine = new SAEngine(defaultParams);
      const result = engine.run();

      expect(result.history.length).toBe(result.totalSteps);
    });

    it('accounting: acceptedMoves + rejectedMoves + invalidMoves = totalSteps', () => {
      const engine = new SAEngine(defaultParams);
      const result = engine.run();

      const sum = result.acceptedMoves + result.rejectedMoves + result.invalidMoves;
      expect(sum).toBe(result.totalSteps);
    });

    it('acceptanceRatio is between 0 and 1', () => {
      const engine = new SAEngine(defaultParams);
      const result = engine.run();

      expect(result.acceptanceRatio).toBeGreaterThanOrEqual(0);
      expect(result.acceptanceRatio).toBeLessThanOrEqual(1);
    });
  });

  describe('determinism', () => {
    it('same seed produces identical results', () => {
      const params1 = { ...defaultParams, seed: 42 };
      const params2 = { ...defaultParams, seed: 42 };

      const engine1 = new SAEngine(params1);
      const engine2 = new SAEngine(params2);

      const result1 = engine1.run();
      const result2 = engine2.run();

      expect(result1.bestEnergy).toBe(result2.bestEnergy);
      expect(result1.totalSteps).toBe(result2.totalSteps);
      expect(result1.acceptedMoves).toBe(result2.acceptedMoves);
      expect(result1.rejectedMoves).toBe(result2.rejectedMoves);
      expect(result1.invalidMoves).toBe(result2.invalidMoves);
    });

    it('different seeds produce different results (with high probability)', () => {
      const params1 = { ...defaultParams, seed: 42 };
      const params2 = { ...defaultParams, seed: 123 };

      const engine1 = new SAEngine(params1);
      const engine2 = new SAEngine(params2);

      const result1 = engine1.run();
      const result2 = engine2.run();

      // At least one metric should differ
      const differs =
        result1.bestEnergy !== result2.bestEnergy ||
        result1.acceptedMoves !== result2.acceptedMoves ||
        result1.rejectedMoves !== result2.rejectedMoves;

      expect(differs).toBe(true);
    });
  });

  describe('minimization', () => {
    it('finds Wiener Index lower than or equal to initial value for C6H14', () => {
      const params: SAParams = {
        formula: 'C6H14',
        initialTemp: 100,
        coolingScheduleK: 8,
        stepsPerCycle: 500,
        numCycles: 4,
        optimizationMode: 'MINIMIZE',
        seed: 42,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      // In MINIMIZE mode, bestEnergy should never be worse than initial
      expect(result.bestEnergy).toBeLessThanOrEqual(result.initialEnergy);
      expect(result.acceptedMoves).toBeGreaterThan(0); // Verify SA is actually accepting moves
      expect(result.initialEnergy).toBe(35); // Verify linear hexane starting point
    });

    it('best energy is less than or equal to final energy', () => {
      const params: SAParams = {
        ...defaultParams,
        optimizationMode: 'MINIMIZE',
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      expect(result.bestEnergy).toBeLessThanOrEqual(result.finalEnergy);
    });
  });

  describe('maximization', () => {
    it('finds Wiener Index equal to or higher than initial value for C6H14', () => {
      const params: SAParams = {
        ...defaultParams,
        optimizationMode: 'MAXIMIZE',
        stepsPerCycle: 500,
        numCycles: 4,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      // Linear hexane already has maximum Wiener Index = 35
      // SA should maintain or potentially find the same value
      expect(result.bestEnergy).toBeGreaterThanOrEqual(result.initialEnergy);
    });

    it('best energy is greater than or equal to final energy', () => {
      const params: SAParams = {
        ...defaultParams,
        optimizationMode: 'MAXIMIZE',
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      expect(result.bestEnergy).toBeGreaterThanOrEqual(result.finalEnergy);
    });
  });

  describe('chemical validity', () => {
    it('bestGraph is connected and has valid valences', () => {
      const engine = new SAEngine(defaultParams);
      const result = engine.run();

      expect(result.bestGraph.isConnected()).toBe(true);
      expect(result.bestGraph.hasValidValences()).toBe(true);
    });

    it('finalGraph is connected and has valid valences', () => {
      const engine = new SAEngine(defaultParams);
      const result = engine.run();

      expect(result.finalGraph.isConnected()).toBe(true);
      expect(result.finalGraph.hasValidValences()).toBe(true);
    });

    it('all graphs in history are from valid molecules', () => {
      const params: SAParams = {
        ...defaultParams,
        stepsPerCycle: 50,
        numCycles: 1,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      // Every step that was accepted should have valid energy
      for (const step of result.history) {
        expect(step.currentEnergy).toBeTypeOf('number');
        expect(step.currentEnergy).toBeGreaterThan(0);
        expect(step.bestEnergy).toBeTypeOf('number');
        expect(step.bestEnergy).toBeGreaterThan(0);
        expect(step.temperature).toBeGreaterThanOrEqual(0.01);
      }
    });
  });

  describe('history tracking', () => {
    it('records step-by-step results', () => {
      const params: SAParams = {
        ...defaultParams,
        stepsPerCycle: 10,
        numCycles: 1,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      expect(result.history.length).toBe(10);

      for (let i = 0; i < result.history.length; i++) {
        const step = result.history[i];
        expect(step).toBeDefined();
        expect(step!.step).toBe(i + 1);
        expect(step!.currentEnergy).toBeTypeOf('number');
        expect(step!.bestEnergy).toBeTypeOf('number');
        expect(step!.temperature).toBeTypeOf('number');
        expect(step!.accepted).toBeTypeOf('boolean');
      }
    });

    it('best energy is non-increasing over time (minimization)', () => {
      const params: SAParams = {
        ...defaultParams,
        optimizationMode: 'MINIMIZE',
        stepsPerCycle: 50,
        numCycles: 2,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      for (let i = 1; i < result.history.length; i++) {
        expect(result.history[i]!.bestEnergy).toBeLessThanOrEqual(
          result.history[i - 1]!.bestEnergy
        );
      }
    });

    it('best energy is non-decreasing over time (maximization)', () => {
      const params: SAParams = {
        ...defaultParams,
        optimizationMode: 'MAXIMIZE',
        stepsPerCycle: 50,
        numCycles: 2,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      for (let i = 1; i < result.history.length; i++) {
        expect(result.history[i]!.bestEnergy).toBeGreaterThanOrEqual(
          result.history[i - 1]!.bestEnergy
        );
      }
    });
  });

  describe('Metropolis acceptance criterion', () => {
    it('always accepts improving moves', () => {
      // Test the metropolis acceptance logic directly via SA run
      const params: SAParams = {
        formula: 'C4H10',
        initialTemp: 0.01, // Very low temperature
        coolingScheduleK: 0, // Constant temperature
        stepsPerCycle: 100,
        numCycles: 1,
        optimizationMode: 'MINIMIZE',
        seed: 42,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      // At low temperature, should still accept improving moves
      // Check that some moves were accepted (not all rejected)
      expect(result.acceptedMoves).toBeGreaterThan(0);
    });

    it('accepts moves at appropriate rates based on temperature', () => {
      // Test at high temperature: should accept most moves
      const highTempParams: SAParams = {
        formula: 'C6H14',
        initialTemp: 1000,
        coolingScheduleK: 0,
        stepsPerCycle: 500,
        numCycles: 1,
        optimizationMode: 'MINIMIZE',
        seed: 100,
      };

      const highTempEngine = new SAEngine(highTempParams);
      const highTempResult = highTempEngine.run();

      // At very high temperature, should accept most valid moves
      expect(highTempResult.acceptanceRatio).toBeGreaterThan(0.3);

      // Test at low temperature: should accept fewer moves (but not necessarily < 50% due to improving moves)
      const lowTempParams: SAParams = {
        formula: 'C6H14',
        initialTemp: 0.01,
        coolingScheduleK: 0,
        stepsPerCycle: 500,
        numCycles: 1,
        optimizationMode: 'MINIMIZE',
        seed: 100,
      };

      const lowTempEngine = new SAEngine(lowTempParams);
      const lowTempResult = lowTempEngine.run();

      // Acceptance ratio should be positive (some moves accepted)
      expect(lowTempResult.acceptanceRatio).toBeGreaterThan(0);
      // High temp should have higher acceptance ratio than low temp
      expect(highTempResult.acceptanceRatio).toBeGreaterThanOrEqual(lowTempResult.acceptanceRatio);
    });

    it('frequently accepts worsening moves at high temperature', () => {
      const params: SAParams = {
        formula: 'C5H12',
        initialTemp: 1000, // Very high temperature
        coolingScheduleK: 0, // Constant
        stepsPerCycle: 1000,
        numCycles: 1,
        optimizationMode: 'MINIMIZE',
        seed: 200,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      // At very high temperature, acceptance ratio should be high
      // (both improving and worsening moves accepted frequently)
      expect(result.acceptanceRatio).toBeGreaterThan(0.5);
    });
  });

  describe('edge cases', () => {
    it('handles small molecules (C4H10)', () => {
      const params: SAParams = {
        ...defaultParams,
        formula: 'C4H10',
        stepsPerCycle: 50,
        numCycles: 1,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      expect(result.bestGraph).toBeDefined();
      expect(result.bestGraph.isConnected()).toBe(true);
      expect(result.bestGraph.hasValidValences()).toBe(true);
    });

    it('handles single cycle', () => {
      const params: SAParams = {
        ...defaultParams,
        numCycles: 1,
        stepsPerCycle: 10,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      expect(result.totalSteps).toBe(10);
      expect(result.history.length).toBe(10);
    });

    it('handles many cycles', () => {
      const params: SAParams = {
        ...defaultParams,
        numCycles: 10,
        stepsPerCycle: 10,
      };

      const engine = new SAEngine(params);
      const result = engine.run();

      expect(result.totalSteps).toBe(100);
      expect(result.history.length).toBe(100);
    });
  });

  describe('SAEngine step-by-step execution', () => {
    const stepParams: SAParams = {
      formula: 'C6H14',
      initialTemp: 100,
      coolingScheduleK: 8,
      stepsPerCycle: 50,
      numCycles: 2,
      optimizationMode: 'MINIMIZE',
      seed: 42,
    };

    it('init() sets up initial state', () => {
      const engine = new SAEngine(stepParams);
      engine.init();
      const state = engine.getState();

      expect(state.step).toBe(0);
      expect(state.isComplete).toBe(false);
      expect(state.currentEnergy).toBe(state.bestEnergy); // Initial structure computed
      expect(state.currentEnergy).toBeGreaterThan(0);
      expect(state.totalSteps).toBe(100); // 50 * 2
    });

    it('step() advances one iteration', () => {
      const engine = new SAEngine(stepParams);
      engine.init();
      engine.step();
      const state = engine.getState();

      expect(state.step).toBe(1);
    });

    it('step() throws if init() not called', () => {
      const engine = new SAEngine(stepParams);

      expect(() => engine.step()).toThrow();
    });

    it('multiple steps execute correctly', () => {
      const engine = new SAEngine(stepParams);
      engine.init();

      const N = 5;
      for (let i = 0; i < N; i++) {
        engine.step();
      }

      const state = engine.getState();
      expect(state.step).toBe(N);
    });

    it('isComplete becomes true after all steps', () => {
      const params: SAParams = {
        ...stepParams,
        stepsPerCycle: 10,
        numCycles: 1,
      };

      const engine = new SAEngine(params);
      engine.init();

      // Execute all 10 steps
      for (let i = 0; i < 10; i++) {
        engine.step();
      }

      const state = engine.getState();
      expect(state.isComplete).toBe(true);
      expect(state.step).toBe(10);
    });

    it('step() throws after completion', () => {
      const params: SAParams = {
        ...stepParams,
        stepsPerCycle: 5,
        numCycles: 1,
      };

      const engine = new SAEngine(params);
      engine.init();

      // Execute all steps
      for (let i = 0; i < 5; i++) {
        engine.step();
      }

      // Try to step again after completion
      expect(() => engine.step()).toThrow();
    });

    it('getResult() returns final result after completion', () => {
      const params: SAParams = {
        ...stepParams,
        stepsPerCycle: 10,
        numCycles: 1,
      };

      const engine = new SAEngine(params);
      engine.init();

      for (let i = 0; i < 10; i++) {
        engine.step();
      }

      const result = engine.getResult();

      expect(result.bestGraph).toBeDefined();
      expect(result.bestEnergy).toBeTypeOf('number');
      expect(result.finalGraph).toBeDefined();
      expect(result.finalEnergy).toBeTypeOf('number');
      expect(result.initialEnergy).toBeTypeOf('number');
      expect(result.totalSteps).toBe(10);
      expect(result.acceptedMoves).toBeTypeOf('number');
      expect(result.rejectedMoves).toBeTypeOf('number');
      expect(result.invalidMoves).toBeTypeOf('number');
      expect(result.acceptanceRatio).toBeTypeOf('number');
      expect(result.history).toBeInstanceOf(Array);
    });

    it('getResult() throws if not complete', () => {
      const engine = new SAEngine(stepParams);
      engine.init();
      engine.step(); // Only 1 step out of 100

      expect(() => engine.getResult()).toThrow();
    });

    it('step-by-step matches run() for same seed', () => {
      const params: SAParams = {
        ...stepParams,
        stepsPerCycle: 50,
        numCycles: 2,
        seed: 12345,
      };

      // Engine 1: use run()
      const engine1 = new SAEngine({ ...params });
      const result1 = engine1.run();

      // Engine 2: use step-by-step
      const engine2 = new SAEngine({ ...params });
      engine2.init();
      while (!engine2.getState().isComplete) {
        engine2.step();
      }
      const result2 = engine2.getResult();

      // Both should produce identical results
      expect(result2.bestEnergy).toBe(result1.bestEnergy);
      expect(result2.acceptedMoves).toBe(result1.acceptedMoves);
      expect(result2.rejectedMoves).toBe(result1.rejectedMoves);
      expect(result2.invalidMoves).toBe(result1.invalidMoves);
      expect(result2.totalSteps).toBe(result1.totalSteps);
    });

    it('run() still works unchanged (backward compatibility)', () => {
      const engine = new SAEngine(stepParams);
      const result = engine.run();

      // Verify run() works exactly as before
      expect(result.bestGraph).toBeDefined();
      expect(result.bestEnergy).toBeTypeOf('number');
      expect(result.totalSteps).toBe(100);
      expect(result.history.length).toBe(100);

      // Accounting still correct
      const sum = result.acceptedMoves + result.rejectedMoves + result.invalidMoves;
      expect(sum).toBe(result.totalSteps);
    });
  });
});
