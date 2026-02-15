/// <reference types="vite/client" />
/// <reference types="@rdkit/rdkit" />

// Vite worker module declarations
declare module '*?worker' {
  const WorkerFactory: new () => Worker;
  export default WorkerFactory;
}
