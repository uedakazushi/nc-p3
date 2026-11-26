"""Artin directions of a length-three point and the length-four trace model."""
from itertools import combinations

import sympy as sp

a, b, c, d = coeffs = sp.symbols("a b c d")
# At the quadruple point, trace extracts the cubic coefficient in k[x]/(x^4).
# In the trace hyperplane with basis 1,x,x^2, the trace matrices are Hankel.
matrix = sp.Matrix([[d, c, b], [c, b, a], [b, a, 0]])
pairs = list(combinations(range(3), 2))
minors = [matrix.extract(i, j).det() for i in pairs for j in pairs]
ideal = sp.groebner([m.subs(d, 1) for m in minors], a, b, c)
assert ideal.polys == sp.groebner([a-c**3, b-c**2, c**4], a, b, c).polys
for coord in (a, b, c):
    others = [v for v in (a, b, c) if v != coord]
    at_infinity = [m.subs({d: 0, coord: 1}) for m in minors]
    assert sp.groebner(at_infinity, *others).contains(sp.Integer(1))

x, y, z, t = sp.symbols("x y z t")
for quadric in (x*z-y*y, y*z, z*z):
    value = quadric.subs({x: 1, y: t, z: t*t})
    assert sp.rem(value, t**3, t) == 0
print("PASS: quadruple trace model on all charts and triple-point direction quadrics")
