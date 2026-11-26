"""The Hankel evaluation at every direction and the rank-eleven splitting list."""
from itertools import combinations, combinations_with_replacement

import sympy as sp

s, t = sp.symbols("s t")
matrix = sp.Matrix([[s, t, 0, 0], [0, s, t, 0], [0, 0, s, t]])
assert [matrix[:, j].det() for j in combinations(range(4), 3)] == [
    s**3, s*s*t, s*t*t, t**3
]
assert [v for v in combinations_with_replacement(range(4), 3) if sum(v) == 4] == [
    (0, 1, 3), (0, 2, 2), (1, 1, 2)
]
block = sp.diag(sp.Matrix([[s, 0], [t, s], [0, t]]), sp.Matrix([s, t]), sp.Matrix([s, t]))
annihilator = sp.diag(sp.Matrix([[t*t, -s*t, s*s]]), sp.Matrix([[t, -s]]), sp.Matrix([[t, -s]]))
assert annihilator * block == sp.zeros(3, 4)
print("PASS: Hankel contractions and the degree (2,1,1) quotient at rank eleven")
