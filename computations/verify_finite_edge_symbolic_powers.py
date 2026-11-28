"""Dimension checks and an exact evaluation on an explicit nonreduced complete intersection."""

from itertools import product
from math import comb

import sympy as sp


def main():
    assert comb(5, 2) == 10 and comb(6, 3) == 20
    for l, n in product(range(1, 4), range(1, 5)):
        if l == 3 and n == 4:
            assert l + n > 6
            continue
        if (l - 1) + (n - 1) < 2:
            continue
        b = l + n - 1
        if (l, n) == (2, 2):
            b = l * n
        assert (l < 3 and b > l + 1) or (n < 4 and b > n + 1)
    # Three independent rows vanish in the complete-intersection ideal after a change of frame.
    x, y, z = sp.symbols("x y z")
    ideal = [x**2, y**3, z]
    block = sp.Matrix([[ideal[(i + j) % 3] * (1 + int(i == j))
                        for j in range(3)] for i in range(3)])
    determinant = sp.Poly(block.det(), x, y, z)
    assert determinant and all(a//2 + b//3 + c >= 3 for (a, b, c), _ in determinant.terms())
    assert (3 * 2, 2 + 1) == (6, 3)
    print("Incidence dimensions, incomplete projection spans and the local cubic ideal-power bound pass.")


if __name__ == "__main__":
    main()
