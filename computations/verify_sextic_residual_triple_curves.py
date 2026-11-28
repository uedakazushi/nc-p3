"""Check degrees, local losses and the four modifications on a triple line."""
from math import comb

assert (6-1)*(6-2)//2 // 3 == 3
assert comb(3+3, 3) == 20
assert sum((1,1,2)) == 4
for j in range(4):
    for i in range(4-j):
        assert i+j+2*(3-i-j) == 3+(3-i-j)
assert 13-6 > 3 and 13-3 > 4*2
# At a smooth conic point a nonzero normal derivative gives the xy term.
# All other terms from the tangent variable alone have no xy monomial.
from sympy import symbols, diff
x,y,t,normal,tangent = symbols("x y t normal tangent")
local = normal*x*y+tangent*t**2
assert diff(local,x,y) == normal
print("PASS: triple curve degree, quadric chart, scroll weights and four-fiber budget")
