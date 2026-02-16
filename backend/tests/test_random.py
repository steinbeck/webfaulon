"""Tests for SeededRandom - ported from src/core/__tests__/random.test.ts"""

import pytest
import sys
from pathlib import Path

# Add backend directory to path to import app modules
sys.path.insert(0, str(Path(__file__).parent.parent))

from app.core.random import SeededRandom


class TestSeededRandomDeterminism:
    """Test constructor and determinism"""

    def test_same_seed_produces_same_first_value(self):
        """Two SeededRandom(42) instances produce identical first value"""
        rng1 = SeededRandom(42)
        rng2 = SeededRandom(42)

        assert rng1.next() == rng2.next()

    def test_identical_sequences_for_same_seed(self):
        """1000 values from seed 12345 match"""
        rng1 = SeededRandom(12345)
        rng2 = SeededRandom(12345)

        for _ in range(1000):
            assert rng1.next() == rng2.next()

    def test_different_sequences_for_different_seeds(self):
        """Seeds 111 vs 222, >90 of 100 values differ"""
        rng1 = SeededRandom(111)
        rng2 = SeededRandom(222)

        sequence1 = [rng1.next() for _ in range(100)]
        sequence2 = [rng2.next() for _ in range(100)]

        differences = sum(1 for v1, v2 in zip(sequence1, sequence2) if v1 != v2)
        assert differences > 90  # Expect most to be different


class TestNext:
    """Test next() method"""

    def test_next_range_0_to_1(self):
        """10000 values all in [0, 1)"""
        rng = SeededRandom(999)

        for _ in range(10000):
            val = rng.next()
            assert 0 <= val < 1

    def test_next_produces_different_values(self):
        """3 consecutive calls are distinct"""
        rng = SeededRandom(777)
        val1 = rng.next()
        val2 = rng.next()
        val3 = rng.next()

        assert val1 != val2
        assert val2 != val3
        assert val1 != val3


class TestNextInt:
    """Test next_int() method"""

    def test_next_int_range_inclusive(self):
        """1000 calls of next_int(0,5) all in [0,5]"""
        rng = SeededRandom(555)

        for _ in range(1000):
            val = rng.next_int(0, 5)
            assert isinstance(val, int)
            assert 0 <= val <= 5

    def test_next_int_covers_range(self):
        """Over 1000 calls, all values 0-5 appear"""
        rng = SeededRandom(888)
        seen = set()

        for _ in range(1000):
            seen.add(rng.next_int(0, 5))

        # Should see all values 0-5
        assert len(seen) == 6
        assert 0 in seen
        assert 5 in seen

    def test_next_int_same_min_max(self):
        """next_int(7,7) always returns 7"""
        rng = SeededRandom(444)

        for _ in range(100):
            assert rng.next_int(7, 7) == 7

    def test_next_int_negative_range(self):
        """next_int(-10,-5) stays in range"""
        rng = SeededRandom(333)

        for _ in range(1000):
            val = rng.next_int(-10, -5)
            assert -10 <= val <= -5


class TestSelectNDistinct:
    """Test select_n_distinct() method"""

    def test_select_n_distinct_count(self):
        """select_n_distinct(4,10) returns exactly 4 unique values"""
        rng = SeededRandom(666)
        selected = rng.select_n_distinct(4, 10)

        assert len(selected) == 4
        unique = set(selected)
        assert len(unique) == 4

    def test_select_n_distinct_range(self):
        """Values in [0, range)"""
        rng = SeededRandom(123)
        selected = rng.select_n_distinct(5, 20)

        for val in selected:
            assert 0 <= val < 20

    def test_select_n_distinct_n_equals_range(self):
        """select_n_distinct(10,10) returns all 0-9"""
        rng = SeededRandom(456)
        selected = rng.select_n_distinct(10, 10)

        assert len(selected) == 10
        unique = set(selected)
        assert len(unique) == 10

        # Should be all values 0-9
        for i in range(10):
            assert i in selected

    def test_select_n_distinct_different_seeds(self):
        """Different seeds produce different selections"""
        rng1 = SeededRandom(111)
        rng2 = SeededRandom(222)

        sel1 = rng1.select_n_distinct(4, 100)
        sel2 = rng2.select_n_distinct(4, 100)

        # Arrays should not be identical
        assert sel1 != sel2

    def test_select_n_distinct_same_seed(self):
        """Same seed produces same selection"""
        rng1 = SeededRandom(789)
        rng2 = SeededRandom(789)

        sel1 = rng1.select_n_distinct(6, 50)
        sel2 = rng2.select_n_distinct(6, 50)

        assert sel1 == sel2
