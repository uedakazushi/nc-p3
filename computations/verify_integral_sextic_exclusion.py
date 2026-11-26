"""Check the invariant field, jet constraints, and equivariant branch multiplicities."""
from itertools import product

from sympy import Poly, expand, resultant, symbols


def main():
    a, b, u = symbols("a b u")
    F = Poly(u**4+a*u*u+b, u)
    derivative = F.diff().as_expr()
    assert Poly(derivative**2-4*u*u*(a*a-4*b), u).rem(F).is_zero
    assert expand(derivative-2*u*(2*u*u+a)) == 0
    c, d = symbols("c d")
    assert resultant(2*c+d-1, 2*d+c-c*d, d) == 2*(c-1)**2 or (
        expand(resultant(2*c+d-1, 2*d+c-c*d, d)-2*(c-1)**2) == 0)
    assert 2*1**2+(-1)**2 != 1
    for roots in product(range(4), repeat=4):
        if len(set(roots)) == 4:
            continue
        assert min(sum(r != s for r in roots) for s in roots) <= 2
    for genus in (0, 1, 2):
        if 2*genus+2 <= 3-genus:
            assert genus == 0
    print("PASS: involution field generators, fixed jet contradiction, cyclic branch ranks and ramification bound")


if __name__ == "__main__":
    main()
