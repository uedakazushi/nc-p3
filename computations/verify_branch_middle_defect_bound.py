"""Check the primitive quartic line and its simple and double local defects."""
from itertools import combinations

import sympy as sp

s, t = sp.symbols("s t")
linear = sp.Matrix([
    [s, 0, 0, 0],
    [t, s, 0, 0],
    [0, t, s, 0],
    [0, 0, t, s],
    [0, 0, 0, t],
])
for loss in range(3):
    cols = [linear[:, i] if i >= loss else sp.eye(5)[:, i]*t for i in range(4)]
    matrix = sp.Matrix.hstack(*cols)
    minors = [matrix[list(rows), :].det() for rows in combinations(range(5), 4)]
    gcd = sp.gcd_list(minors)
    assert gcd == t**loss
    degrees = {sp.Poly(sp.cancel(m/gcd), s, t).total_degree() for m in minors if m != 0}
    assert degrees == {4-loss}
print("PASS: primitive degrees at zero, one and two elementary losses")
