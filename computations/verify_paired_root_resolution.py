"""Check the eight local coordinate systems and compact exceptional companion curves."""

from itertools import combinations

import sympy as sp

from verify_all_zero_pairing_defects import parameters
from verify_klein_elliptic_completion import coeffs
from verify_kummer_polar import polarity


def check():
    b, c, d, z, l, x = sp.symbols("b c d z l x")
    coordinates = (b, c, d)
    _, p, q = coeffs(sp.S.One, b, c, d, z, l)
    excluded = sp.Poly(
        z * l * (z * z - 1) * (l * l * z - 1) * (l * l - 1) * (l * l * z * z - 1), z, l
    )

    def unit(value):
        top, bottom = sp.factor(value).as_numer_denom()
        return value != 0 and all(
            excluded.rem(sp.Poly(term, z, l)).is_zero
            for poly in (top, bottom)
            for term, _ in sp.factor_list(poly)[1]
        )

    jacobians = [sp.Matrix(triple).jacobian(coordinates) for triple in (p, q)]
    inverse = polarity(z, l).inv()
    for i in range(8):
        first, companion, residual = parameters(i, z, l, x)
        base = dict(zip(coordinates, list(first)[1:]))
        assert unit(residual)
        for jac in jacobians:
            assert unit(jac.subs(base).applyfunc(sp.factor).det())
        normal = (inverse * sp.Matrix(list(first))).applyfunc(sp.factor)
        assert sp.factor(normal.dot(sp.Matrix(list(companion)))) == 0
        _, other, _ = coeffs(*companion, z, l)
        assert sp.factor(other[1] ** 2 - 4 * other[0] * other[2]) == 0
        # A base-point-free basis of binary quadratics gives a conic embedding.
        den = sp.lcm_list([sp.denom(sp.cancel(value)) for value in companion])
        polys = [sp.Poly(sp.cancel(den * value), x) for value in companion]
        matrix = sp.Matrix([[poly.nth(k) for k in range(3)] for poly in polys])
        assert max(poly.degree() for poly in polys) == 2
        assert any(
            unit(matrix[list(rows), :].det()) for rows in combinations(range(4), 3)
        )
    print(
        "Both triples are local coordinates and every exceptional companion map embeds a conic at all eight nodes."
    )


if __name__ == "__main__":
    check()
