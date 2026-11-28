"""Sectional genera and intersection lattices of the quartic normalization."""
from fractions import Fraction

ruled = []
for q in range(1, 5):
    for k in range(1, 5):
        for l in range(1, 5):
            if (2*q-2)*k-Fraction(8, k) > -4 or (2*q-2)*l-Fraction(4, l) > -2:
                continue
            square = 12*Fraction(l, k)-8*Fraction(l*l, k*k)
            if square >= 4:
                ruled.append((q, k, l, square))
assert ruled == [(1, 1, 1, 4), (1, 2, 1, 4), (1, 2, 2, 4), (2, 1, 1, 4)]
assert 2*6/2-4/4 == 5
assert 8+4-2*6-2*2+2*2 == 0
for invariant in (0, 2):
    def meet(x, y):
        return -invariant*x[0]*y[0]+x[0]*y[1]+x[1]*y[0]
    L = (1, (4+invariant)//2)
    K = (-2, -2-invariant)
    for a in (1, 2):
        H = (a, 6-(2-invariant//2)*a)
        assert (meet(L,L), meet(H,H), meet(H,L)) == (4, 8, 6)
        trace = tuple(2*l-k-2*h for l,k,h in zip(L,K,H))
        assert trace == ((2, invariant-2) if a == 1 else (0, 2))
print("PASS: all irregular ruled, adjoint pencil and minimal-degree quartic cases")
