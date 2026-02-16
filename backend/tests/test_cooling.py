"""Tests for cooling schedule - ported from src/core/__tests__/cooling.test.ts"""

import pytest
import sys
from pathlib import Path

# Add backend directory to path to import app modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.core.cooling import compute_temperature


class TestComputeTemperature:
    """Test compute_temperature function"""

    initial_temp = 100
    total_steps = 1000

    def test_k0_constant_temperature(self):
        """k=0 maintains initial temperature at all steps"""
        k = 0
        assert compute_temperature(0, self.total_steps, self.initial_temp, k) == 100
        assert compute_temperature(500, self.total_steps, self.initial_temp, k) == 100
        assert compute_temperature(1000, self.total_steps, self.initial_temp, k) == 100

    def test_k1_initial_temperature(self):
        """k=1 at step 0 returns initial temperature"""
        k = 1
        assert compute_temperature(0, self.total_steps, self.initial_temp, k) == 100

    def test_k1_final_temperature(self):
        """k=1 at step totalSteps reaches minimum temperature"""
        k = 1
        # T = T0 - k * T0 * step / totalSteps = 100 - 1 * 100 * 1000 / 1000 = 0
        # But clamped to 0.01
        assert compute_temperature(self.total_steps, self.total_steps, self.initial_temp, k) == 0.01

    def test_k1_halfway(self):
        """k=1 at halfway returns half initial temperature"""
        k = 1
        # T = 100 - 1 * 100 * 500 / 1000 = 50
        assert compute_temperature(500, self.total_steps, self.initial_temp, k) == 50

    def test_k1_quarter(self):
        """k=1 at 25% returns 75% of initial temperature"""
        k = 1
        # T = 100 - 1 * 100 * 250 / 1000 = 75
        assert compute_temperature(250, self.total_steps, self.initial_temp, k) == 75

    def test_k8_initial(self):
        """k=8 at step 0 returns initial temperature"""
        k = 8
        assert compute_temperature(0, self.total_steps, self.initial_temp, k) == 100

    def test_k8_reaches_min_early(self):
        """k=8 reaches minimum temperature early"""
        k = 8
        # At step = totalSteps/8 = 125:
        # T = 100 - 8 * 100 * 125 / 1000 = 100 - 100 = 0, clamped to 0.01
        assert compute_temperature(125, self.total_steps, self.initial_temp, k) == 0.01

    def test_k8_stays_at_min(self):
        """k=8 stays at minimum after reaching it"""
        k = 8
        assert compute_temperature(500, self.total_steps, self.initial_temp, k) == 0.01
        assert compute_temperature(1000, self.total_steps, self.initial_temp, k) == 0.01

    def test_k32_reaches_min_very_early(self):
        """k=32 reaches minimum very early"""
        k = 32
        # At step = totalSteps/32 = 31.25:
        # T = 100 - 32 * 100 * 32 / 1000 = 100 - 102.4 = negative, clamped to 0.01
        assert compute_temperature(32, self.total_steps, self.initial_temp, k) == 0.01

    def test_clamps_to_minimum(self):
        """Never goes below 0.01"""
        k = 1
        # Should never go below 0.01
        assert compute_temperature(self.total_steps, self.total_steps, self.initial_temp, k) >= 0.01

        k2 = 100
        assert compute_temperature(1, self.total_steps, self.initial_temp, k2) >= 0.01

    def test_different_initial_temperatures(self):
        """Works with different initial temperatures"""
        k = 1
        temp200 = 200
        assert compute_temperature(0, self.total_steps, temp200, k) == 200
        assert compute_temperature(500, self.total_steps, temp200, k) == 100
        assert compute_temperature(1000, self.total_steps, temp200, k) == 0.01

    def test_total_steps_1(self):
        """Handle edge case of totalSteps = 1"""
        k = 1
        assert compute_temperature(0, 1, self.initial_temp, k) == 100
        assert compute_temperature(1, 1, self.initial_temp, k) == 0.01
