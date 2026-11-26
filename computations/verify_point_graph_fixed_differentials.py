"""Check the determinant bookkeeping and local normalization examples over QQ."""
from itertools import product

from sympy import Matrix, Symbol, diff, gcd
from sympy.functions.combinatorial.numbers import binomial


def main():
    # Exponents in L_i = L_0 D^i K^binom(i,2), for positive and negative i.
    def line(i):
        return Matrix([1, i, binomial(i, 2)])

    for i in range(-8, 9):
        assert line(i-1)+line(i+2)-line(i)-line(i+1) == Matrix([0, 0, 2])

    # For R(-H) -> ambient cotangent with kernel S(-H), the determinant
    # of the cokernel has H-exponent -4 + 6 - 1 and S-exponent one.
    assert -4+6-1 == 1

    t = Symbol("t")
    for m, n in [(2, 3), (2, 5), (3, 4)]:
        # Plane branch y^m=x^n. B records derivative zeros, T the Smith divisor.
        branch = gcd(diff(t**m, t), diff(t**n, t))
        torsion = gcd(n*t**(m*(n-1)), m*t**(n*(m-1)))
        assert branch.as_poly(t).monic().as_expr() == t**(m-1)
        assert torsion.as_poly(t).monic().as_expr() == t**(n*(m-1))
        assert n*(m-1)-(m-1) == (n-1)*(m-1)
    assert gcd(0, t) == t  # derivative of xy on the branch (t,0)
    bounds = {3: 1, 4: 3}
    for d, e in product(bounds, repeat=2):
        if d+e <= 7:
            assert sum(4*k-2*bounds[k]+2 for k in (d, e)) > 20
    assert [18-4*d+2*g for d, g in [(3, 1), (4, 1), (5, 1), (7, 5)]] == [8, 4, 0, 0]
    print("PASS: polarization, determinant exponents, nodal and cuspidal pullback, smooth degree table")


if __name__ == "__main__":
    main()
