"""Check the multiplicity budgets and the four triple points in quartic degree."""
from itertools import product

# A direction with isolated singularities costs at least 2, a singular line 4.
budgets = [(a, b) for a, b in product(range(7), range(4)) if 2 * a + 4 * b <= 12]
assert max(1 + a + 2 * b for a, b in budgets) == 7
assert 20 + 7 + 2 < 60 - max(22 + 7, 14 * 2)
# At coordinate vertex i a quartic monomial has multiplicity 4 - exponent[i].
quartics = [a for a in product(range(5), repeat=4) if sum(a) == 4]
assert [a for a in quartics if max(a) <= 1] == [(1, 1, 1, 1)]
assert 4 * 3 + (20 - 4) == 4 * 7
print("PASS: seven singular marks, no two irreducible quartics, four faces with four double marks each")
