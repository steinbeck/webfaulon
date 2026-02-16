"""
Cooling schedule implementations from Faulon paper Tables 3-4.

The paper uses: kT_t = kT_0 - k * kT_0 * t / delta_t
where:
- kT_t: temperature at step t
- kT_0: initial temperature
- k: cooling rate parameter (0, 1, 8, 32, etc.)
- t: current step
- delta_t: total steps

The parameter k controls cooling rate:
- k=0 (f0): constant temperature (kT stays at kT_0)
- k=1 (f1): linear cooling to 0
- k=8 (f8): fast decay (paper's best schedule)
- k=32 (f32): very fast decay
"""

from typing import Union

# Cooling schedule type (k parameter)
CoolingScheduleType = Union[int, float]

# Minimum temperature to avoid division by zero in Metropolis
MIN_TEMPERATURE = 0.01


def compute_temperature(
    step: int,
    total_steps: int,
    initial_temp: float,
    schedule_k: CoolingScheduleType
) -> float:
    """
    Compute temperature at a given step using Faulon cooling schedule.

    Formula: T = max(MIN_TEMP, T0 - k * T0 * step / totalSteps)

    Args:
        step: Current step number (0 to total_steps)
        total_steps: Total number of steps
        initial_temp: Initial temperature (kT_0)
        schedule_k: Cooling rate parameter (k)

    Returns:
        Temperature at current step (clamped to MIN_TEMPERATURE)

    Examples:
        Constant temperature (k=0):
        >>> compute_temperature(500, 1000, 100, 0)
        100

        Linear cooling (k=1):
        >>> compute_temperature(0, 1000, 100, 1)
        100
        >>> compute_temperature(500, 1000, 100, 1)
        50
        >>> compute_temperature(1000, 1000, 100, 1)
        0.01

        Fast decay (k=8):
        >>> compute_temperature(0, 1000, 100, 8)
        100
        >>> compute_temperature(125, 1000, 100, 8)  # reaches min at 1/8
        0.01
    """
    # Formula from Faulon paper: T = T0 - k * T0 * t / delta_t
    temperature = initial_temp - (schedule_k * initial_temp * step) / total_steps

    # Clamp to minimum to avoid division by zero in Metropolis
    return max(MIN_TEMPERATURE, temperature)
