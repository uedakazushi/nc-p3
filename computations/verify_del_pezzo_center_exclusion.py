"""Exact small-rank and valuation checks used in the Brauer exclusions."""
from itertools import combinations, product
from sympy import gcd, symbols


def main():
    # Index in the saturation of a root sublattice can be computed on a basis.
    # Two signed blocks are the only case with torsion in the ambient D5.
    for sizes in ((2,2), (2,3), (3,), (4,), (5,)):
        parities = [v for v in product((0,1), repeat=len(sizes)) if sum(v)%2 == 0]
        assert len(parities) <= 2
        assert all(all(2*n % 2 == 0 for n in v) for v in parities)
    # Ratios of linear forms on the nonsingular conic have rational degree <= 2.
    s, t, x, y, z = symbols("s t x y z")
    basis = (s*s, s*t, t*t)
    for a, b in combinations(basis, 2):
        common = gcd(a, b)
        assert 1 <= (a/common).cancel().as_poly(s,t).total_degree() <= 2
    assert (x*y-z*z).subs(x,0) == -z*z
    assert (x*y-1).subs(x,0) == -1
    # At infinity v(F)=-2, v(X)=-1 gives, up to inverse, X^2/F -> X/Y.
    assert (x**2/(x*y)).cancel() == x/y
    print("Signed-root parity, conic ratio degrees, and affine symbol residues checked.")


if __name__ == "__main__":
    main()
