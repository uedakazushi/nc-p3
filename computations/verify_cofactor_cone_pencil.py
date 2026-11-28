"""Test first nonzero arc coefficients for the two sheets of a cone ruling."""
from sympy import Symbol, series, sqrt

s = Symbol("s")
for ramified in (False, True):
    orders = range(3) if ramified else range(1, 3)
    for m in orders:
        base = s**2 if ramified else s
        unit = s if ramified else 1
        z = 2
        w = unit*sqrt(1+base+z)
        z2 = z + 2*base**m*w
        w2 = -unit*sqrt(1+base+z2)
        order = 2*m+1 if ramified else m
        def coeff(expr):
            return series(expr, s, 0, order+1).removeO().expand()
        # B=(1,s_base) and C=s_base^m(1+s_base,3s_base) have parallel leading terms.
        assert coeff(z-z2 + (w-w2)*base**m*(1+base)) == 0
        assert coeff((z-z2)*base + (w-w2)*base**m*3*base) == 0
        assert coeff(z-z2) != 0

assert all(2+2*(3-b) < 8 for b in (1, 2))
print("PASS: unramified and vertically ramified arc collisions, and plane-pencil degrees")
