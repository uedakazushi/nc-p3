"""Rank and divisor calculations of the two incidence scroll charts."""
import sympy as sp


def main():
    s, t, m, n, b = sp.symbols("s t m n b")
    # Coefficients of (s*y1-t*y0)(m*y0+n*y1) and (s*y1-t*y0)b.
    z = sp.Matrix([-t*m, s*m-t*n, s*n, -t*b, s*b])
    monomials = (s*m, s*n, s*b, t*m, t*n, t*b)
    matrix = sp.Matrix([[sp.expand(poly).coeff(term) for term in monomials] for poly in z])
    assert matrix.rank() == 5
    incidence = z[0] + z[4]
    assert sp.expand(incidence+t*m-s*b) == 0
    # Its exceptional section is m=b=0; its special fiber is t=b=0.
    assert z.subs({m: 0, b: 0, n: 1}) == sp.Matrix([0, -t, s, 0, 0])
    assert z.subs({t: 0, b: 0, s: 1}) == sp.Matrix([0, m, n, 0, 0])
    for lost in (0, 1):
        residual = 2+lost
        assert 1-residual < 0
    l1, l2, l3 = sp.symbols("l1 l2 l3")
    left = sp.Matrix([l1, l2, l3])
    v = sp.Matrix(sp.symbols("v1:4"))
    w = sp.Matrix(sp.symbols("w1:4"))
    pencil = sp.Matrix.hstack(s*left, t*left, v, w)
    delta = sp.Matrix.hstack(left, v, w).det()
    assert sp.expand(pencil[:, (0, 2, 3)].det()-s*delta) == 0
    assert sp.expand(pencil[:, (1, 2, 3)].det()-t*delta) == 0
    print("PASS: scroll parametrization, exceptional lines, residual and quintic factors")


if __name__ == "__main__":
    main()
