/**
 * SSE (Server-Sent Events) wrapper for SA stream endpoint.
 *
 * Provides typed EventSource connection with automatic cleanup and typed handlers.
 */

import type {
  SSEProgressData,
  SSECompleteData,
  SSEWaitingData,
} from './types';

/**
 * Typed handlers for SSE events.
 */
export interface SSEHandlers {
  onProgress?: (data: SSEProgressData) => void;
  onComplete?: (data: SSECompleteData) => void;
  onWaiting?: (data: SSEWaitingData) => void;
  onError?: (event: Event) => void;
}

/**
 * SSE connection wrapper for SA streaming endpoint.
 *
 * Handles EventSource lifecycle, typed event parsing, and automatic cleanup.
 */
export class SSEConnection {
  private eventSource: EventSource | null = null;

  /**
   * Connect to SSE stream endpoint for a session.
   *
   * Automatically closes any existing connection before creating new one.
   * Auto-closes on "complete" event.
   *
   * @param sessionId - Session ID to stream
   * @param handlers - Typed event handlers
   */
  connect(sessionId: string, handlers: SSEHandlers): void {
    // Clean up any existing connection first
    this.close();

    // Create new EventSource connection
    this.eventSource = new EventSource(`/api/sa/${sessionId}/stream`);

    // Register typed event listeners
    if (handlers.onProgress) {
      this.eventSource.addEventListener('progress', (e: MessageEvent) => {
        try {
          const data: SSEProgressData = JSON.parse(e.data);
          handlers.onProgress?.(data);
        } catch (error) {
          console.error('[SSE] Failed to parse progress event:', error);
        }
      });
    }

    if (handlers.onComplete) {
      this.eventSource.addEventListener('complete', (e: MessageEvent) => {
        try {
          const data: SSECompleteData = JSON.parse(e.data);
          handlers.onComplete?.(data);
        } catch (error) {
          console.error('[SSE] Failed to parse complete event:', error);
        } finally {
          // Auto-close on completion
          this.close();
        }
      });
    }

    if (handlers.onWaiting) {
      this.eventSource.addEventListener('waiting', (e: MessageEvent) => {
        try {
          const data: SSEWaitingData = JSON.parse(e.data);
          handlers.onWaiting?.(data);
        } catch (error) {
          console.error('[SSE] Failed to parse waiting event:', error);
        }
      });
    }

    // Register error handler
    if (handlers.onError) {
      this.eventSource.onerror = (event: Event) => {
        handlers.onError?.(event);
      };
    }
  }

  /**
   * Close the SSE connection if active.
   *
   * Safe to call multiple times. Sets eventSource to null after closing.
   */
  close(): void {
    if (this.eventSource) {
      this.eventSource.close();
      this.eventSource = null;
    }
  }

  /**
   * Check if connection is currently open.
   *
   * @returns true if EventSource is in OPEN state, false otherwise
   */
  isConnected(): boolean {
    return this.eventSource?.readyState === EventSource.OPEN;
  }
}
