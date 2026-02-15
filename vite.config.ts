import { defineConfig } from 'vite';
import path from 'path';

export default defineConfig({
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
  worker: {
    format: 'es',
  },
  build: {
    target: 'es2020',
  },
  optimizeDeps: {
    exclude: ['@rdkit/rdkit'], // Don't pre-bundle WASM
  },
});
