"""Degree, intersection and residue budgets for a nodal marked cubic cone."""
from itertools import product

from sympy import expand, symbols


def intersect(first, second, k):
    a, b = first
    c, d = second
    return -k * a * c + a * d + b * c


assert [u for u in range(16, 21) if (50 - 3 * u) % 4 in (0, 3)] == [17, 18]
t = symbols("t")
assert expand((1 + 2*t + 3*t**2 + t**3) * (1-t)**2) == 1 - 3*t**3 + t**4 + t**5
for l, u in product((1, 2), (17, 18)):
    k = 3 - l
    a, j = (1, k), (2, 6 - l)
    length = sum(intersect(x, y, k) for x, y in ((a, a), (a, j), (j, j)))
    assert length == 21 - 2*l
    assert 0 <= length - (u-1) <= 3
    assert (-k - 2) + k + 6 - l == 4 - l
print("PASS: length-seven vertex, u=17 or 18, at most three residual values separated on the base")
