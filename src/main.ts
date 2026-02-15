import * as Comlink from 'comlink';
import type { ISAWorker } from './worker/types';

let worker: Worker | null = null;
let saWorker: Comlink.Remote<ISAWorker> | null = null;

export async function initWorker(): Promise<Comlink.Remote<ISAWorker>> {
  if (saWorker) return saWorker;

  worker = new Worker(
    new URL('./worker/sa-worker.ts', import.meta.url),
    { type: 'module' }
  );
  saWorker = Comlink.wrap<ISAWorker>(worker);
  return saWorker;
}

export function terminateWorker(): void {
  if (worker) {
    worker.terminate();
    worker = null;
    saWorker = null;
  }
}

// Cleanup on page exit
window.addEventListener('beforeunload', terminateWorker);

// Log initialization for verification
console.log('WebFaulon main.ts loaded');

// Temporary smoke test (will be replaced by Alpine.js UI in Plan 03)
async function smokeTest() {
  const proxy = await initWorker();
  await proxy.initialize({
    formula: 'C6H14',
    initialTemp: 100,
    coolingScheduleK: 8,
    stepsPerCycle: 50,
    numCycles: 2,
    optimizationMode: 'MINIMIZE',
    seed: 42
  });

  const result = await proxy.run(
    Comlink.proxy((progress) => {
      console.log(`Step ${progress.step}/${progress.totalSteps} | Energy: ${progress.currentEnergy} | Best: ${progress.bestEnergy}`);
    }),
    10
  );

  console.log('SA complete:', result.bestEnergy, 'accepted:', result.acceptedMoves);
  terminateWorker();
}

// Run smoke test on load
smokeTest().catch(console.error);
