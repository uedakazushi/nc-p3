"""Saturated coefficient loci and bilinear fibers when both determinant quadratics vanish."""

from itertools import product

import sympy as sp

from verify_klein_elliptic_completion import coeffs, numerator


def audit():
    a, b, c, d, z, l, t, x, y = sp.symbols("a b c d z l t x y")
    _, P, Q = coeffs(a, b, c, d, z, l)
    assert sp.cancel(P[2] - (a - z * d) * (a * z - d) / (z + 1) ** 2) == 0
    assert sp.cancel(Q[2] - (b - c) * (z * z * b - c) / (z + 1) ** 2) == 0
    for variable, value, poly in [
        (d, z * a, P),
        (d, a / z, P),
        (c, b, Q),
        (c, z * z * b, Q),
    ]:
        sample = {a: 1, b: 3, c: 4, d: 5, z: 2, l: 3}
        sample[variable] = value.subs(sample)
        assert poly[2].subs(sample) == 0 and poly[0].subs(sample) != 0

    cases = {
        (0, 0, 0, 0): 1 / l,
        (0, 0, 1, 1): l * z,
        (0, 1, 0, 1): l,
        (0, 1, 1, 0): 1 / (l * z),
        (1, 0, 0, 1): 1 / (l * z),
        (1, 0, 1, 0): l,
        (1, 1, 0, 0): l / z,
        (1, 1, 1, 1): 1 / (l * z * z),
    }
    one = sp.Integer(1)
    for flags in product(range(2), repeat=4):
        dd, cc, ff, ee = flags
        d = (z, 1 / z)[dd]
        c = (b, z * z * b)[cc]
        f0 = sp.cancel((1 + d) / (z + 1))
        g0 = z * f0
        e0 = sp.cancel((z * b + c) / (z + 1))
        target = [
            1 / f0 if not ff else 1 / g0,
            1 / e0 if not ee else z / e0,
            z / e0 if not ee else 1 / e0,
            1 / f0 if not ff else z * z / g0,
        ]
        det = d - b * c
        equations = [
            sp.factor(n - det * v).as_numer_denom()[0]
            for n, v in zip(numerator(one, b, c, d, z, l), target)
        ]
        unit = sp.together(
            z * b * l * (z * z - 1) * (1 - l * l * z) * det
        ).as_numer_denom()[0]
        ideal = sp.groebner(equations + [t * unit - 1], t, l, b, z, order="grevlex")
        if flags not in cases:
            assert list(ideal) == [1]
            continue
        value = cases[flags]
        assert ideal.reduce(sp.together(b - value).as_numer_denom()[0])[1] == 0
        c, det = c.subs(b, value), det.subs(b, value)
        coeff, P, Q = coeffs(one, value, c, d, z, l)
        assert all(sp.cancel(v) == 0 for v in (*P, *Q))
        f0, f1, g0, g1, e0, e1, h0, h1 = coeff
        e, f, g, h = [
            sp.cancel(v) for v in (e0 + e1 * y, f0 + f1 * x, g0 + g1 * x, h0 + h1 * y)
        ]
        assert sp.cancel((f1, g1)[ff]) == sp.cancel((e1, h1)[ee]) == 0
        A = numerator(e, f, g, h, z, l)
        C = numerator(one, f, g, d, z, one)
        D = numerator(e, value, c, h, z, one)
        det1 = e * h - f * g
        residual = [
            sp.cancel(
                A[i] * (y if i in (0, 3) else x)
                - det1 * (D[i] if i in (0, 3) else C[i])
            )
            for i in range(4)
        ]
        active = [sp.Poly(v, x, y) for v in residual if v != 0]
        pivot = active[0].coeff_monomial(x * y)
        assert pivot != 0

        # Removing precisely the inverted factors leaves a nonzero constant.
        def is_unit(v):
            num, den = sp.factor(v).as_numer_denom()
            return all(
                sp.rem(unit**8, p, z, l) == 0
                for q in (num, den)
                for p, _ in sp.factor_list(q)[1]
            )

        # 'unit' contains b before specialization.
        unit = sp.cancel(unit.subs(b, value)).as_numer_denom()[0] * l * z
        assert is_unit(pivot)
        for poly in active:
            ratio = sp.cancel(poly.coeff_monomial(x * y) / pivot)
            assert (
                is_unit(ratio)
                and sp.cancel(poly.as_expr() - ratio * active[0].as_expr()) == 0
            )
        assert all(poly.degree(x) <= 1 and poly.degree(y) <= 1 for poly in active)
    print("Sixteen saturated systems: eight empty, eight unit-leading bilinear fibers.")


if __name__ == "__main__":
    audit()
