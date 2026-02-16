import { defineConfig } from 'vite';
import path from 'path';

export default defineConfig({
  base: '/webfaulon/',
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
  server: {
    proxy: {
      '/api': {
        target: 'http://localhost:8000',
        changeOrigin: true,
      },
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
