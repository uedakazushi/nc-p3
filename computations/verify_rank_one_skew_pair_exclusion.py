"""Algebraic identities and intersection budgets for two skew rank-one lines."""
from sympy import Matrix, expand, factor, symbols

u, w, a, b = symbols("u w a b")
first = Matrix([u**2, u*w, w**2])
second = Matrix([a**2, a*b, b**2])
delta = u*b-w*a
assert first.cross(second).applyfunc(expand) == (
    delta*Matrix([w*b, -(u*b+w*a), u*a])
).applyfunc(expand)
t, alpha, c3, c4 = symbols("t alpha c3 c4", nonzero=True)
assert factor(c3*((t+alpha)*c4-1).subs(t, -1/c3)) == alpha*c3*c4-c3-c4

def intersection(left, right):
    return left[0]*right[1]+left[1]*right[0]

for d in range(1, 5):
    for e1 in range(2*d+1):
        e2 = 2*d-e1
        assert d+max(e1, e2)+min(e1, e2) == 3*d
        assert 2*d+min(e1, e2) <= 3*d

assert intersection((2, 1), (0, 1)) == 2
assert intersection((3, 1), (1, 3)) == 10
assert intersection((1, 1), (1, 1)) == 2
for i,j in ((1, 0), (0, 1), (1, 1)):
    residual = intersection((3-i, 1-j), (1-i, 3-j))
    assert residual == (4 if i == j else 6)
assert 3+6+4 < 14 and 4+4+4 < 14
assert 11+1 < 14 and 4+4+3 < 14
print("PASS: binary square minors, scalar compatibility, and direction/divisor counts")
