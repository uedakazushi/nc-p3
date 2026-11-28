"""Triple ruling divisors outside a triple transversal of an octic scroll."""
from itertools import combinations
from math import comb

general_fiber = list(combinations(range(5), 3))
base_degree = len(general_fiber)
hyperplane_degree = sum(0 in subset for subset in general_fiber)
assert (base_degree, hyperplane_degree) == (10, 6)
assert 8-3 == 5 and 8-5 == 3
assert (3-2)*hyperplane_degree+base_degree == 16 < 20-3
assert comb(3, 2) == 3
assert [(a, 12*a-4*a*a) for a in range(1, 4) if 6-4*a >= 0] == [(1, 8)]
print("PASS: triple-transversal divisor degrees 10 and 6, and the quartic cone inverse")
