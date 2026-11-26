"""Quadratic projection and degree obstructions for paired binary centers."""
from itertools import product

import sympy as sp


def main():
    a, b = sp.symbols("a b")
    matrix = sp.Matrix([[a, b, 0], [0, a, b]])
    normal = matrix.row(0).cross(matrix.row(1))
    assert normal == sp.Matrix([[b*b, -a*b, a*a]])
    assert normal[0]*normal[2] == normal[1]**2
    survivors = []
    for k, m, b2 in product(range(1, 4), range(2), range(6)):
        d3 = 2*k + m
        d4 = 3*d3 - 8 + b2
        if 4 <= d3 <= 6 and d4 <= 5:
            survivors.append((k, m, d3, d4))
    assert set(survivors) == {(2, 0, 4, 4), (2, 0, 4, 5)}
    assert all(3*d4-3*d3+3+1 > 3 for _, _, d3, d4 in survivors)
    conic_cases = [(2-k, 2*k+m) for k, m in product((1, 2), range(2)) if 2*k+m <= 3]
    assert conic_cases == [(1, 2), (1, 3)]
    print("PASS: binary cone, cubic contradiction and remaining conic degrees")


if __name__ == "__main__":
    main()
