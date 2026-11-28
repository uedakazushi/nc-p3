"""Intersection numbers and the two branch conditions in the cubic cone proof."""
from itertools import combinations_with_replacement
from math import comb

from sympy import Matrix, factor, symbols


def blowup_intersection(first, second, third):
    # Coordinates are the coefficients of H, E_1, ..., E_20.
    return sum(a * b * c for a, b, c in zip(first, second, third))


polarization = [3] + [-1] * 20
degree, marks = symbols("d m")
assert blowup_intersection(polarization, polarization, [degree, -marks] + [0] * 19) == 9 * degree - marks
assert blowup_intersection(polarization, polarization, [6] + [-3] * 20) == -6
assert [comb(6, 3) - comb(6 - e, 3) for e in (1, 2, 3)] == [10, 16, 19]

# On the normalized nodal base choose the two branch lifts (1:0),(0:1).
# Differences in each positive degree are surjective using the first pure power.
w, z = symbols("w z")
jumps = Matrix([[w**2, w, 1], [z**2, z, 1]])
assert factor(jumps[:, 1:].det()) == w - z
assert sum(3 * j + 1 for j in range(4)) - 19 - 2 == 1

bounds = {1: 7, 2: 14, 3: 22}
partitions = [p for k in range(3, 9) for p in combinations_with_replacement(bounds, k) if sum(p) == 8]
maximum = max(sum(bounds[e] for e in p) for p in partitions)
assert maximum == 58 < 60
print(f"PASS: cubic dimensions 10,16,19; two independent conductor values; maximum multiplicity {maximum}/60")
