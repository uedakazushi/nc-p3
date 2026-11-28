"""Check component intersection budgets and ruled cubic normalization invariants."""
from itertools import product

from sympy import symbols, expand

e, m = symbols("e m")
assert expand(4 * (3 * e) - 2 * (9 * e - m) - m) == m - 6 * e
assert 4 * 2**3 - 20 == 60 - 6 * 8 == 12
isolated = [(n, b) for n, b in product(range(16), range(5)) if n + b >= 18]
assert isolated == [(14, 4), (15, 3), (15, 4)]
ruled = [(j, k) for k in range(4) for j in range(k + 1) if 2 * k - j == 3]
assert ruled == [(1, 2), (3, 3)]
assert all(2 * k - j - 2 == 1 and k - j <= 1 for j, k in ruled)
assert 2 * (15 + 3) - 6 * 3 == 18
print("PASS: 32 points with residual 12, three cubic budgets, both ruled normalizations")
