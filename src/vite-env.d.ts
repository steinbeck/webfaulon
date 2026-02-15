/// <reference types="vite/client" />

// Vite worker module declarations
declare module '*?worker' {
  const WorkerFactory: new () => Worker;
  export default WorkerFactory;
}
