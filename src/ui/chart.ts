import {
  Chart,
  LineController,
  LineElement,
  PointElement,
  LinearScale,
  Filler,
  Tooltip,
} from 'chart.js';

// Register only needed components (tree-shake Chart.js)
// Note: Decimation is built-in, doesn't need registration
Chart.register(
  LineController,
  LineElement,
  PointElement,
  LinearScale,
  Filler,
  Tooltip
);

// Module-level chart instance (CRITICAL: NOT inside Alpine reactive data)
let _chartInstance: Chart | null = null;

/**
 * Create a new Wiener Index line chart on the provided canvas.
 * Destroys any existing chart first.
 *
 * @param canvas - The HTMLCanvasElement to render the chart on
 */
export function createWienerChart(canvas: HTMLCanvasElement): void {
  // Clean up old chart if exists
  if (_chartInstance) {
    _chartInstance.destroy();
    _chartInstance = null;
  }

  // Create new chart instance
  _chartInstance = new Chart(canvas, {
    type: 'line',
    data: {
      datasets: [
        {
          label: 'Wiener Index',
          data: [],
          borderColor: 'rgb(59, 130, 246)', // blue-500
          backgroundColor: 'rgba(59, 130, 246, 0.1)',
          fill: true,
          borderWidth: 2,
          pointRadius: 0, // disable point rendering for performance
          tension: 0, // straight lines (required for decimation)
        },
      ],
    },
    options: {
      animation: false, // CRITICAL for real-time performance
      parsing: false, // CRITICAL: provide pre-parsed {x, y} data
      responsive: true,
      maintainAspectRatio: true,
      aspectRatio: 2.5, // wide chart for step timeline
      plugins: {
        decimation: {
          enabled: true,
          algorithm: 'min-max', // preserves SA peaks/valleys
        },
        legend: {
          display: false, // single dataset, label not needed
        },
        tooltip: {
          enabled: true,
        },
      },
      scales: {
        x: {
          type: 'linear',
          title: {
            display: true,
            text: 'Step',
            font: {
              size: 14,
            },
          },
        },
        y: {
          type: 'linear',
          title: {
            display: true,
            text: 'Wiener Index',
            font: {
              size: 14,
            },
          },
        },
      },
    },
  });
}

/**
 * Add a data point to the chart and update it.
 *
 * @param step - The step number (x-axis)
 * @param wienerIndex - The Wiener Index value (y-axis)
 */
export function addChartDataPoint(step: number, wienerIndex: number): void {
  if (!_chartInstance || !_chartInstance.data.datasets[0]) {
    return; // silently return if no chart instance or dataset
  }

  // Push pre-parsed {x, y} data point
  _chartInstance.data.datasets[0].data.push({ x: step, y: wienerIndex });

  // Update with 'none' mode = no animation, fastest
  _chartInstance.update('none');
}

/**
 * Clear all data from the chart and update it.
 */
export function resetChart(): void {
  if (!_chartInstance || !_chartInstance.data.datasets[0]) {
    return; // silently return if no chart instance or dataset
  }

  // Clear dataset
  _chartInstance.data.datasets[0].data = [];

  // Update with 'none' mode = no animation
  _chartInstance.update('none');
}

/**
 * Destroy the chart instance and free resources.
 */
export function destroyChart(): void {
  if (_chartInstance) {
    _chartInstance.destroy();
    _chartInstance = null;
  }
}
