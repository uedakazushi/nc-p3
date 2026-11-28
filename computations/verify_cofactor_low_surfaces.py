"""Weighted cone adjunction and degree alternatives for an integral intersection."""
from math import comb
from sympy import Matrix, eye, symbols
from itertools import combinations

def weighted_forms(n):
    return sum(n-2*j+1 for j in range(n//2+1)) if n >= 0 else 0

for degree in range(1, 13):
    for twist in range(10, 15):
        dimension = weighted_forms(2*twist)-weighted_forms(2*twist-degree)
        if degree == 10:
            assert dimension == 10*twist-15
    assert (weighted_forms(degree-10) != 0) == (degree >= 10)
assert 6*10-2*15 > 2*14
pairs = [(a,b) for a in range(3,7) for b in range(a,7) if 2*a <= 8 and a+b >= 7]
assert pairs == [(3,4),(3,5),(3,6),(4,4),(4,5),(4,6)]
for a,b in pairs:
    if b >= 5:
        assert 2*b > 8 and 2*(6-b) < a
    elif a+b == 7:
        assert comb(a+b-7+3,3) == 1
    else:
        assert a == b == 4
# Verify that contractions of the exterior square recover its supporting space.
for rank in range(2,5):
    columns = eye(4)[:, :rank]
    wedges = [columns[:,i]*columns[:,j].T-columns[:,j]*columns[:,i].T
              for i,j in combinations(range(rank),2)]
    assert Matrix.hstack(*wedges).rank() == rank
print("PASS: cone Hilbert polynomial, mark bound, complete-intersection cases and contractions")
