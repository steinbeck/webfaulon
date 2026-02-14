/**
 * Cooling schedule implementations from Faulon paper Tables 3-4
 *
 * The paper uses: kT_t = kT_0 - k * kT_0 * t / delta_t
 * where:
 * - kT_t: temperature at step t
 * - kT_0: initial temperature
 * - k: cooling rate parameter (0, 1, 8, 32, etc.)
 * - t: current step
 * - delta_t: total steps
 *
 * The parameter k controls cooling rate:
 * - k=0 (f0): constant temperature (kT stays at kT_0)
 * - k=1 (f1): linear cooling to 0
 * - k=8 (f8): fast decay (paper's best schedule)
 * - k=32 (f32): very fast decay
 */

/**
 * Cooling schedule type (k parameter)
 *
 * Common values:
 * - 0: constant temperature
 * - 1: linear cooling
 * - 8: fast decay (recommended in paper)
 * - 32: very fast decay
 */
export type CoolingScheduleType = number;

const MIN_TEMPERATURE = 0.01; // Minimum to avoid division by zero in Metropolis

/**
 * Compute temperature at a given step using Faulon cooling schedule
 *
 * Formula: T = max(MIN_TEMP, T0 - k * T0 * step / totalSteps)
 *
 * @param step Current step number (0 to totalSteps)
 * @param totalSteps Total number of steps
 * @param initialTemp Initial temperature (kT_0)
 * @param scheduleK Cooling rate parameter (k)
 * @returns Temperature at current step (clamped to MIN_TEMPERATURE)
 *
 * @example
 * // Constant temperature (k=0)
 * computeTemperature(500, 1000, 100, 0) // => 100
 *
 * @example
 * // Linear cooling (k=1)
 * computeTemperature(0, 1000, 100, 1)    // => 100 (start)
 * computeTemperature(500, 1000, 100, 1)  // => 50 (halfway)
 * computeTemperature(1000, 1000, 100, 1) // => 0.01 (end, clamped)
 *
 * @example
 * // Fast decay (k=8)
 * computeTemperature(0, 1000, 100, 8)    // => 100 (start)
 * computeTemperature(125, 1000, 100, 8)  // => 0.01 (reaches min at 1/8)
 */
export function computeTemperature(
  step: number,
  totalSteps: number,
  initialTemp: number,
  scheduleK: CoolingScheduleType
): number {
  // Formula from Faulon paper: T = T0 - k * T0 * t / delta_t
  const temperature = initialTemp - (scheduleK * initialTemp * step) / totalSteps;

  // Clamp to minimum to avoid division by zero in Metropolis
  return Math.max(MIN_TEMPERATURE, temperature);
}
