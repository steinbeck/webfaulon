"""
Seeded pseudo-random number generator for reproducible SA runs.

Uses Mulberry32 algorithm - fast, good distribution, deterministic.
Faulon 1996 p.733: "Two runs using the same seed will lead to the same results"

This implementation must produce identical output to the TypeScript version
in src/core/random.ts for cross-language determinism.
"""


class SeededRandom:
    """Seeded PRNG using Mulberry32 algorithm by Tommy Ettinger"""

    def __init__(self, seed: int):
        """
        Initialize state with seed.

        Args:
            seed: Initial seed value (will be converted to 32-bit signed int)
        """
        # JavaScript: this.state = seed | 0
        # The | 0 converts to 32-bit signed integer
        # In Python, we simulate this by masking to 32 bits
        self.state = seed & 0xFFFFFFFF

    def next(self) -> float:
        """
        Generate next random number in range [0, 1).

        Uses Mulberry32 algorithm with careful 32-bit unsigned arithmetic
        to match JavaScript's Math.imul behavior.

        Returns:
            Random float in [0, 1)
        """
        # Mulberry32 algorithm from TypeScript:
        # let t = (this.state += 0x6d2b79f5);
        # t = Math.imul(t ^ (t >>> 15), t | 1);
        # t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        # return ((t ^ (t >>> 14)) >>> 0) / 4294967296;

        # Step 1: Add magic constant and update state
        self.state = (self.state + 0x6d2b79f5) & 0xFFFFFFFF
        t = self.state

        # Step 2: First mixing step
        # Math.imul(a, b) does 32-bit signed multiplication
        # In Python: multiply and mask to 32 bits, handling sign
        a = (t ^ (t >> 15)) & 0xFFFFFFFF
        b = (t | 1) & 0xFFFFFFFF
        t = self._imul32(a, b)

        # Step 3: Second mixing step
        a = (t ^ (t >> 7)) & 0xFFFFFFFF
        b = (t | 61) & 0xFFFFFFFF
        mixed = self._imul32(a, b)
        t = (t ^ (t + mixed)) & 0xFFFFFFFF

        # Step 4: Final mixing and normalization
        t = (t ^ (t >> 14)) & 0xFFFFFFFF

        # Convert to [0, 1) - divide by 2^32
        return t / 4294967296

    def _imul32(self, a: int, b: int) -> int:
        """
        Emulate JavaScript's Math.imul - 32-bit signed integer multiplication.

        JavaScript's Math.imul returns the low 32 bits of the product,
        treating inputs and output as signed 32-bit integers.

        Args:
            a: First operand
            b: Second operand

        Returns:
            Low 32 bits of product as unsigned integer
        """
        # Convert to signed 32-bit range if needed
        if a & 0x80000000:
            a = a - 0x100000000
        if b & 0x80000000:
            b = b - 0x100000000

        # Multiply and take low 32 bits
        product = (a * b) & 0xFFFFFFFF
        return product

    def next_int(self, min_val: int, max_val: int) -> int:
        """
        Generate random integer in range [min, max] inclusive.

        Args:
            min_val: Minimum value (inclusive)
            max_val: Maximum value (inclusive)

        Returns:
            Random integer in [min_val, max_val]
        """
        range_size = max_val - min_val + 1
        return int(self.next() * range_size) + min_val

    def select_n_distinct(self, n: int, range_size: int) -> list[int]:
        """
        Select n distinct random integers from range [0, range_size).

        Used for selecting 4 distinct atoms for displacement in SA algorithm.

        Args:
            n: Number of distinct values to select
            range_size: Upper bound of range (exclusive)

        Returns:
            List of n distinct integers from [0, range_size)

        Raises:
            ValueError: If n > range_size
        """
        if n > range_size:
            raise ValueError(f"Cannot select {n} distinct values from range {range_size}")

        selected: list[int] = []
        available: set[int] = set(range(range_size))

        # Select n distinct values
        for _ in range(n):
            available_list = list(available)
            idx = self.next_int(0, len(available_list) - 1)
            value = available_list[idx]
            selected.append(value)
            available.remove(value)

        return selected
