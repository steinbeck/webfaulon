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
