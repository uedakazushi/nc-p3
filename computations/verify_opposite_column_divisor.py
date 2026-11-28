"""Check the residual Hilbert-function budget, including nonsemisimple layers."""
from itertools import combinations_with_replacement

from sympy import Matrix, expand, symbols

survivors = []
for count in range(1, 11):
    for lengths in combinations_with_replacement(range(1, 8), count):
        if 13 <= sum(lengths) <= 20 and sum(min(5, a) for a in lengths) <= 10:
            survivors.append(lengths)
assert survivors == [(6, 7), (7, 7)]
x, y, r0, r1, r2, r3 = symbols("x y r0 r1 r2 r3")
matrix = Matrix([[x*x, x*y, x, 2*x], [x*y, r0, r1, 2], [x, r1, r2, 3], [2*x, 2, 3, r3]])
assert expand(matrix.det()).subs(x, 0) == 0
assert expand(matrix.det().diff(x)).subs(x, 0) == 0
print("PASS: only layers (6,7),(7,7) fit rank ten; proportional contraction has double factor")
