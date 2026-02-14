import { describe, it, expect } from 'vitest';
import { computeTemperature, type CoolingScheduleType } from '../cooling';

describe('cooling schedules', () => {
  describe('computeTemperature', () => {
    const initialTemp = 100;
    const totalSteps = 1000;

    it('k=0 (constant temperature) maintains initial temperature at all steps', () => {
      const k: CoolingScheduleType = 0;
      expect(computeTemperature(0, totalSteps, initialTemp, k)).toBe(100);
      expect(computeTemperature(500, totalSteps, initialTemp, k)).toBe(100);
      expect(computeTemperature(1000, totalSteps, initialTemp, k)).toBe(100);
    });

    it('k=1 (linear cooling) at step 0 returns initial temperature', () => {
      const k: CoolingScheduleType = 1;
      expect(computeTemperature(0, totalSteps, initialTemp, k)).toBe(100);
    });

    it('k=1 (linear cooling) at step totalSteps reaches minimum temperature', () => {
      const k: CoolingScheduleType = 1;
      // T = T0 - k * T0 * step / totalSteps = 100 - 1 * 100 * 1000 / 1000 = 0
      // But clamped to 0.01
      expect(computeTemperature(totalSteps, totalSteps, initialTemp, k)).toBe(0.01);
    });

    it('k=1 (linear cooling) at halfway returns half initial temperature', () => {
      const k: CoolingScheduleType = 1;
      // T = 100 - 1 * 100 * 500 / 1000 = 50
      expect(computeTemperature(500, totalSteps, initialTemp, k)).toBe(50);
    });

    it('k=1 (linear cooling) at 25% returns 75% of initial temperature', () => {
      const k: CoolingScheduleType = 1;
      // T = 100 - 1 * 100 * 250 / 1000 = 75
      expect(computeTemperature(250, totalSteps, initialTemp, k)).toBe(75);
    });

    it('k=8 (fast decay) at step 0 returns initial temperature', () => {
      const k: CoolingScheduleType = 8;
      expect(computeTemperature(0, totalSteps, initialTemp, k)).toBe(100);
    });

    it('k=8 (fast decay) reaches minimum temperature early', () => {
      const k: CoolingScheduleType = 8;
      // At step = totalSteps/8 = 125:
      // T = 100 - 8 * 100 * 125 / 1000 = 100 - 100 = 0, clamped to 0.01
      expect(computeTemperature(125, totalSteps, initialTemp, k)).toBe(0.01);
    });

    it('k=8 (fast decay) stays at minimum after reaching it', () => {
      const k: CoolingScheduleType = 8;
      expect(computeTemperature(500, totalSteps, initialTemp, k)).toBe(0.01);
      expect(computeTemperature(1000, totalSteps, initialTemp, k)).toBe(0.01);
    });

    it('k=32 (very fast decay) reaches minimum very early', () => {
      const k: CoolingScheduleType = 32;
      // At step = totalSteps/32 = 31.25:
      // T = 100 - 32 * 100 * 32 / 1000 = 100 - 102.4 = negative, clamped to 0.01
      expect(computeTemperature(32, totalSteps, initialTemp, k)).toBe(0.01);
    });

    it('clamps to minimum temperature of 0.01', () => {
      const k: CoolingScheduleType = 1;
      // Should never go below 0.01
      expect(computeTemperature(totalSteps, totalSteps, initialTemp, k)).toBeGreaterThanOrEqual(0.01);

      const k2: CoolingScheduleType = 100;
      expect(computeTemperature(1, totalSteps, initialTemp, k2)).toBeGreaterThanOrEqual(0.01);
    });

    it('works with different initial temperatures', () => {
      const k: CoolingScheduleType = 1;
      const temp200 = 200;
      expect(computeTemperature(0, totalSteps, temp200, k)).toBe(200);
      expect(computeTemperature(500, totalSteps, temp200, k)).toBe(100);
      expect(computeTemperature(1000, totalSteps, temp200, k)).toBe(0.01);
    });

    it('handles edge case of totalSteps = 1', () => {
      const k: CoolingScheduleType = 1;
      expect(computeTemperature(0, 1, initialTemp, k)).toBe(100);
      expect(computeTemperature(1, 1, initialTemp, k)).toBe(0.01);
    });
  });
});
