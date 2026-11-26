"""Check the ruling section and the first-order common-kernel equations."""
from itertools import combinations

import sympy as sp

a, b, s, t = variables = sp.symbols("a b s t")
matrix = sp.Matrix([[a, 0, b], [s, t, 0], [0, s, t]])
pairs = list(combinations(range(3), 2))
minors = [matrix.extract(i, j).det() for i in pairs for j in pairs]
for x in variables:
    others = [v for v in variables if v != x]
    actual = sp.groebner([m.subs(x, 1) for m in minors], *others)
    expected = sp.groebner([v.subs(x, 1) for v in (s, t)], *others)
    assert actual.polys == expected.polys

u, v = sp.symbols("u v")
left = sp.Matrix([[u, v, 0], [0, u, v]])
assert [left[:, j].det() for j in pairs] == [u*u, u*v, v*v]
right = sp.Matrix([[0, 1], [1, 0]])
assert right.det() != 0
print("PASS: reduced ruling line on all charts and full binary contraction pencil")
