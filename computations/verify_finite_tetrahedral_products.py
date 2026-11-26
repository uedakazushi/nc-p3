"""Reproduce the dimension of the simultaneous off-diagonal cubic identities."""
import sympy as sp
from itertools import combinations

x = sp.symbols("x0:4")
pairs = list(combinations(range(4), 2))
quadrics = [x[i]*x[j] for i in range(4) for j in range(i, 4)]
coeffs = sp.symbols("c0:60")
p = {pair: sum(c*m for c, m in zip(coeffs[10*n:10*n+10], quadrics))
     for n, pair in enumerate(pairs)}
equations = []
for k in range(4):
    for i, j in combinations([n for n in range(4) if n != k], 2):
        form = x[i]*p[tuple(sorted((j, k)))] - x[j]*p[tuple(sorted((i, k)))]
        equations.extend(sp.Poly(form, *x).coeffs())
matrix, _ = sp.linear_eq_to_matrix(equations, coeffs)
expected = sp.Matrix([int(m == x[i]*x[j]) for i, j in pairs for m in quadrics])
assert matrix.nullspace() == [expected]
print("PASS: all six off-diagonal quadrics have one common scalar")
