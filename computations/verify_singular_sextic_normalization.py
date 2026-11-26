"""Local fixed ideals and the two-fixed-point genus-two trace obstruction."""
from sympy import I, Matrix, series, symbols


def main():
    t, b = symbols("t b", nonzero=True)
    transported = t/(1+b*t)
    for k in (2, 3):
        assert series(transported**k-t**k, t, 0, k+2).removeO().expand() == -k*b*t**(k+1)
    # A cusp with multiplier -1 has fixed ideal (t^3), of colength three.
    semigroup = {0, *range(2, 13)}
    assert {n for n in semigroup if n-3 not in semigroup} == {0, 2, 4}
    for xi in (I, -I):
        s = Matrix.diag(1, xi)
        tau = Matrix.diag(1, -1)
        g = -tau
        assert (Matrix.eye(2)-g)*s == Matrix.eye(2)+tau
        assert s**2 == tau
        assert all(2-m.trace()-m.conjugate().trace() == 0 for m in (s, tau*s))
    assert (2+4*7) % 4 == 2   # one Kummer branch exponent two
    print("PASS: cusp fixed-ideal orders, genus-two eigenvalues and cyclic-cover obstruction")


if __name__ == "__main__":
    main()
