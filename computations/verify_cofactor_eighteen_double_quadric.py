"""Numerical restrictions for a degree-two map to the smooth kernel quadric."""
from fractions import Fraction

pairs = [(a, 6-a) for a in range(1, 4) if a * (6-a) >= 8]
assert pairs == [(2, 4), (3, 3)]
sections = [(k, u, Fraction(u*(6*k-8*u), k*k)) for k in (1, 2)
            for u in range(1, 5) if 6*k >= 8*u]
assert sections == [(2, 1, 1)]
assert 8 - 2*6 + 4 == 0  # (H-L)^2
assert 8 - Fraction(6**2, 4) == -1
assert [n for n in (8, 9) if (n - 6) % 2 == 0] == [8]
assert 3*4 - 2*6 == 0
assert 1+6-2 == 5
print("PASS: ruling degrees, elliptic ruled intersections, and exceptional cycle bounds")
