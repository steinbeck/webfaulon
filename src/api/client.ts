/**
 * Typed REST API client for backend SA endpoints.
 *
 * Provides methods for all 5 SA REST endpoints with typed requests/responses.
 */

import type {
  SAConfigParams,
  ConfigureResponse,
  ControlResponse,
  StatusResponse,
} from './types';

/**
 * API client for Simulated Annealing backend endpoints.
 */
export class SAAPIClient {
  private baseURL = import.meta.env.VITE_API_URL
    ? `${import.meta.env.VITE_API_URL}/api/sa`
    : '/api/sa';

  /**
   * Create a new SA session with provided configuration.
   *
   * @param params - SA configuration parameters
   * @returns Session ID and status
   * @throws Error if request fails
   */
  async configure(params: SAConfigParams): Promise<ConfigureResponse> {
    const response = await fetch(`${this.baseURL}/configure`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(params),
    });

    if (!response.ok) {
      await this.handleError(response, 'configure');
    }

    return response.json();
  }

  /**
   * Start an SA session.
   *
   * @param sessionId - Session ID to start
   * @returns Session ID and status
   * @throws Error if request fails
   */
  async start(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'start');
  }

  /**
   * Pause a running SA session.
   *
   * @param sessionId - Session ID to pause
   * @returns Session ID and status
   * @throws Error if request fails
   */
  async pause(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'pause');
  }

  /**
   * Reset an SA session to initial state.
   *
   * @param sessionId - Session ID to reset
   * @returns Session ID and status
   * @throws Error if request fails
   */
  async reset(sessionId: string): Promise<ControlResponse> {
    return this.controlRequest(sessionId, 'reset');
  }

  /**
   * Get current status of an SA session.
   *
   * @param sessionId - Session ID to query
   * @returns Complete session state snapshot
   * @throws Error if request fails
   */
  async getStatus(sessionId: string): Promise<StatusResponse> {
    const response = await fetch(`${this.baseURL}/${sessionId}/status`);

    if (!response.ok) {
      await this.handleError(response, 'getStatus');
    }

    return response.json();
  }

  /**
   * Execute a control action (start/pause/reset).
   *
   * @param sessionId - Session ID
   * @param action - Control action name
   * @returns Session ID and status
   * @throws Error if request fails
   */
  private async controlRequest(
    sessionId: string,
    action: string
  ): Promise<ControlResponse> {
    const response = await fetch(`${this.baseURL}/${sessionId}/${action}`, {
      method: 'POST',
    });

    if (!response.ok) {
      await this.handleError(response, action);
    }

    return response.json();
  }

  /**
   * Handle HTTP error responses.
   *
   * Attempts to parse JSON error body and extract meaningful error message.
   * Falls back to generic error if parsing fails.
   *
   * @param response - Failed fetch response
   * @param action - Action name for error message
   * @throws Error with extracted or generic message
   */
  private async handleError(response: Response, action: string): Promise<never> {
    let errorMessage = `HTTP ${response.status}: ${action} failed`;

    try {
      const errorBody = await response.json();
      // FastAPI returns errors in "detail" field, custom errors might use "error"
      if (errorBody.detail) {
        errorMessage = errorBody.detail;
      } else if (errorBody.error) {
        errorMessage = errorBody.error;
      }
    } catch {
      // JSON parsing failed, use default message
    }

    throw new Error(errorMessage);
  }
}
