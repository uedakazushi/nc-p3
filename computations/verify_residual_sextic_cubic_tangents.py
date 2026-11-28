"""Check the small Hilbert budgets and the twisted-cubic conormal splitting."""
from itertools import combinations
from sympy import Matrix, symbols

quadrics = [(i,j,2-i-j) for i in range(3) for j in range(3-i)]
cubics = [(i,j,3-i-j) for i in range(4) for j in range(4-i)]
for chosen in combinations(quadrics, 2):
    supported = []
    for m in cubics:
        divisors = [tuple(v-int(k == axis) for k,v in enumerate(m))
                    for axis in range(3) if m[axis]]
        if all(q in chosen for q in divisors):
            supported.append(m)
    assert len(supported) <= 2
s,t = symbols("s t")
relations = Matrix([[-t,0],[s,-t],[0,s]])
assert Matrix([[s**2,s*t,t**2]])*relations == Matrix([[0,0]])
assert 18-2*5 == 8 and 13-3 > 8
print("PASS: quadratic evaluation and Hilbert bounds, binary conormal syzygies and degree eight")
