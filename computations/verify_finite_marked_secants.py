"""Check reduced and doubled rank-one secants scheme-theoretically."""
from itertools import combinations

import sympy as sp


def unit(i, j):
    m = sp.zeros(3)
    m[i, j] = 1
    return m


a, b, c, d = variables = sp.symbols("a b c d")
first = sp.Matrix([[1, 0, 2], [3, 1, -1], [1, 4, 2]])
second = sp.Matrix([[0, 3, 1], [1, -2, 4], [3, 1, 1]])
pairs = tuple(combinations(range(3), 2))
for direction, divisor in [(unit(1, 1), a*b), (unit(0, 1) + unit(1, 0), b**2)]:
    matrix = a*unit(0, 0) + b*direction + c*first + d*second
    minors = [matrix.extract(rows, cols).det() for rows in pairs for cols in pairs]
    for coord in variables:
        others = [v for v in variables if v != coord]
        actual = sp.groebner([m.subs(coord, 1) for m in minors], *others)
        expected = sp.groebner([m.subs(coord, 1) for m in (c, d, divisor)], *others)
        assert actual.polys == expected.polys
print("PASS: both rank-one secant schemes, including the tangent double point")
