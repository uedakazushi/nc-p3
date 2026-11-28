"""Check the edge saturation, disjoint plane/line budget, and regulus interpolation."""
from itertools import combinations
from sympy import factor, symbols

edges = list(combinations(range(4), 2))
for chosen in combinations(edges, 5):
    assert min(sum(v in e for e in chosen) for v in range(4)) >= 2
assert 20-9 > 6+4
assert 6+3+3*3 < 20 and 14-1+3+3 < 20
s, t = symbols("s t")
basis = [s*(s-t), t*(t-s), s*t]
for i, (sp, tp) in enumerate([(1,0), (0,1), (1,1)]):
    assert [q.subs({s:sp,t:tp}) for q in basis] == [int(i==j) for j in range(3)]
assert factor(basis[0]*basis[1]*basis[2]) == -s**2*t**2*(s-t)**2
print("PASS: five saturated edges, both marked bounds and the three-line interpolation coefficients")
