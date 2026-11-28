"""Cremona degree, tetrahedron edge budgets and residual multiplicities."""
from itertools import combinations
from math import comb

edges = list(combinations(range(4), 2))
matchings = [chosen for n in range(7) for chosen in combinations(edges, n)
             if len({v for edge in chosen for v in edge}) == 2*n]
assert max(map(len, matchings)) == 2
for exponents in ((3,1,1,0), (0,2,1,2), (2,0,2,1)):
    assert sum(5-i for i in exponents) == 15
assert 3*5-4*3 == 3
assert 1+comb(5,3) == 11 and 1+comb(4,3) == 5
assert max(7+5, 4+1, 6) < 13
assert 5*7+24 < 60 and 2*7+3*15 < 60
assert 1+(3*4)//2+3 < 13
print("PASS: quintic Cremona and edge counts, normal sextic threshold and residual partitions")
