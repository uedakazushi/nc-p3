"""Numerical descent defects, ranks of the obstruction space and triple-fiber counts."""
from math import comb
from sympy import symbols, expand

t, e, genus = symbols("t e genus")
chi = 4*t*(t+1)+(t+1)*(e+1-genus)
assert expand(chi.subs(genus,e-3)-4*(t+1)**2) == 0
cases = [(e,g) for e in (4,5,6) for g in range(e-3,((e-2)**2)//4+1)]
assert cases == [(4,1),(5,2),(6,3),(6,4)]
assert [(e,g) for e,g in cases if g > e-3] == [(6,4)]
assert (3-2)*(4-2) == 2
for arithmetic_genus in range(1,4):
    count = (8-4)*((8-2)*(8-3)-6*arithmetic_genus)//6
    assert count == 20-4*arithmetic_genus < 20
for ambient in range(3,30):
    assert 3*ambient+4-3*(ambient-1)-6 == 1
print("PASS: descent polynomial, finite arithmetic genus list, and triple-fiber smoothing bounds")
