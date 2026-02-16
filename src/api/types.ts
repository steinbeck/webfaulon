/**
 * TypeScript interfaces mirroring backend Pydantic models.
 *
 * IMPORTANT: All fields use snake_case to match backend JSON exactly.
 * Do NOT convert to camelCase - FastAPI serializes Pydantic with snake_case.
 */

/**
 * SA configuration parameters.
 * Mirrors: backend/app/models/sa_params.py::SAParams
 */
export interface SAConfigParams {
  formula: string;
  initial_temp: number;
  cooling_schedule_k: number;
  steps_per_cycle: number;
  num_cycles: number;
  optimization_mode: 'MINIMIZE' | 'MAXIMIZE';
  seed: number;
}

/**
 * Response from configure endpoint.
 * Mirrors: backend/app/api/sa_configure.py::ConfigureResponse
 */
export interface ConfigureResponse {
  session_id: string;
  status: string;
}

/**
 * Response from control endpoints (start/pause/reset).
 * Mirrors: backend/app/api/sa_control.py::ControlResponse
 */
export interface ControlResponse {
  session_id: string;
  status: string;
}

/**
 * Response from status endpoint.
 * Mirrors: backend/app/api/sa_status.py::StatusResponse
 */
export interface StatusResponse {
  session_id: string;
  session_state: string;
  step: number;
  total_steps: number;
  cycle: number;
  current_energy: number;
  best_energy: number;
  best_smiles: string;
  best_svg: string;
  temperature: number;
  accepted_moves: number;
  rejected_moves: number;
  invalid_moves: number;
  is_complete: boolean;
}

/**
 * SSE "progress" event data payload.
 * Mirrors: backend/app/api/sa_stream.py lines 84-97 (progress event data field)
 */
export interface SSEProgressData {
  step: number;
  total_steps: number;
  cycle: number;
  temperature: number;
  current_energy: number;
  best_energy: number;
  best_smiles: string;
  best_svg: string;
  accepted_moves: number;
  rejected_moves: number;
  invalid_moves: number;
  is_complete: boolean;
}

/**
 * SSE "complete" event data payload.
 * Mirrors: backend/app/api/sa_stream.py lines 60-70 (complete event data field)
 */
export interface SSECompleteData {
  step: number;
  total_steps: number;
  best_energy: number;
  best_smiles: string;
  best_svg: string;
  accepted_moves: number;
  rejected_moves: number;
  invalid_moves: number;
}

/**
 * SSE "waiting" event data payload.
 * Mirrors: backend/app/api/sa_stream.py lines 141-146 (waiting event data field)
 */
export interface SSEWaitingData {
  session_state: string;
  message: string;
}

/**
 * SSE "error" event data payload.
 * Mirrors: backend/app/api/sa_stream.py lines 151-154 (error event data field)
 */
export interface SSEErrorData {
  error: string;
}
