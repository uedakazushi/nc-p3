"""Euler polynomial, degree-eight syzygies and curve-length budgets."""
from math import comb
from sympy import expand, symbols

t = symbols("t")
def hilbert(s):
    return (s+1)*(s+2)*(s+3)/6

assert expand(hilbert(t)-10*hilbert(t-6)+15*hilbert(t-8)-6*hilbert(t-10)) == 80
assert 10*comb(5,3)-15 == 85
d, chi, b, a = symbols("d chi b a")
assert expand((d*t+chi+a)+(d*(6-t)+chi+b)) == 6*d+2*chi+b+a
for components in range(1, 9):
    for small in range(components+1):
        for connected in range(1, components+1):
            assert 2*(components-connected)-2*small < 2*(components-small)
for curve_degree in range(1, 30):
    for plucker_degree in range(curve_degree, 2*curve_degree+1):
        upper = 2*(4*curve_degree-2*plucker_degree)-6*curve_degree+3*plucker_degree-plucker_degree
        assert upper == 2*(curve_degree-plucker_degree) <= 0
print("PASS: cofactor Euler polynomial, syzygy counts and component-degree obstruction")
