"""Branch quartics and length-four defects of the two flip eigenspaces."""

import sympy as sp

from flint_rational_function import rational_function_field
from two_window_ext import PathRow
from verify_all_zero_pairing_defects import relation_bases
from verify_klein_elliptic_completion import coeffs


def ordinary_unit(value, z, l):
    factors = (l, z, z - 1, z + 1, l - 1, l + 1, l * z - 1, l * z + 1, l * l * z - 1)
    value = sp.factor(value)
    return value != 0 and all(
        f in factors for side in sp.fraction(value) for f, _ in sp.factor_list(side)[1]
    )


def branch_polynomial(sign, a, z, l):
    _, (p, m, n), _ = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    scale = l * (z * z - 1) * (a * a - z)
    return sp.Poly(sp.cancel(scale**2 * (m * m - 4 * p * n)), a)


def check_quartic(sign, a, z, l):
    quartic = branch_polynomial(sign, a, z, l)
    A, B, C, D, E = quartic.all_coeffs()
    scale = l**2 * (z * z - 1) ** 2
    assert sp.expand(A - scale) == sp.expand(E - z * z * scale) == 0
    resultant = sp.resultant(quartic.as_expr(), a * a - z, a)
    assert (
        sp.factor(resultant - 16 * z**4 * (z - sign) ** 4 * (l * l * z - 1) ** 4) == 0
    )
    disc = sp.factor(quartic.discriminant())
    common = (
        4096
        * l**4
        * z**6
        * (l * l - 1) ** 2
        * (l * l * z * z - 1) ** 2
        * (l * l * z - 1) ** 4
    )
    assert sp.factor(disc / common - (z - 1) ** 4 * (z - sign) ** 12) == 0
    invariant = sp.factor(12 * A * E - 3 * B * D + C * C)
    u = (l * l - 1) * (l * l * z * z - 1)
    v = (l * l * z - 1) ** 2
    assert (
        sp.expand(invariant - 16 * z * z * (z - sign) ** 4 * (u * u - u * v + v * v))
        == 0
    )
    assert (
        sp.factor(disc - 16**3 * z**6 * (z - sign) ** 12 * u * u * (u - v) ** 2 * v * v)
        == 0
    )

    co, (p, m, _), (q, k, _) = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    x, y = -m / (2 * p), -k / (2 * q)
    f, e = co[0] + co[1] * x, co[4] + co[5] * y
    factors = set()
    for value in (p, q, x, y, e, f, e * e - z * f * f):
        value = sp.factor(value)
        assert value != 0
        for side in sp.fraction(value):
            factors.update(term for term, _ in sp.factor_list(side)[1])
    for term in factors:
        assert ordinary_unit(sp.resultant(quartic.as_expr(), term, a), z, l), term


def bad_sections(z, l):
    yield from [
        (1, l * z, 1, l),
        (1, l * z, z, l * z),
        (1, 1 / l, z, 1 / l),
        (1, 1 / l, 1, 1 / (l * z)),
    ]
    yield from [
        (-1, l * z, -1, l),
        (-1, l * z, z, -l * z),
        (-1, -l * z, 1, l),
        (-1, -l * z, -z, -l * z),
        (-1, 1 / l, -z, 1 / l),
        (-1, 1 / l, 1, -1 / (l * z)),
        (-1, -1 / l, z, 1 / l),
        (-1, -1 / l, -1, -1 / (l * z)),
    ]


def check_sections(z, l):
    field = rational_function_field(("z", "l"))
    for sign, a, e, f in bad_sections(z, l):
        A = sp.Matrix([[a, 1], [sign * z, sign * a]])
        B = sp.Matrix([[e, f], [sign * z * f, sign * e]])
        co, (p, m, n), (q, k, o) = coeffs(*A, z, l)
        x = sign * (a * a - z * f * f)
        y = sign * (e * e - z)
        assert all(
            ordinary_unit(v, z, l)
            for v in (A.det(), B.det(), p, q, x, y, 2 * p * x + m)
        )
        assert all(
            sp.cancel(v) == 0
            for v in (
                p * x * x + m * x + n,
                q * y * y + k * y + o,
                2 * p * x + m + 2 * q * y + k,
            )
        )
        f0, f1, g0, g1, e0, e1, h0, h1 = co
        assert all(
            sp.cancel(v) == 0
            for v in (
                e0 + e1 * y - e,
                f0 + f1 * x - f,
                g0 + g1 * x - sign * z * f,
                h0 + h1 * y - sign * e,
            )
        )
        windows = relation_bases(A, B, z, l, {z: 3, l: 2}, field)
        profiles = [PathRow(windows, i, 5, field.one).dims[-2:] for i in range(4)]
        dims = [quartic for quartic, _ in profiles]
        assert sorted(dims) == ([36] * 4 if sign == 1 else [35, 35, 36, 36])
        assert sorted(quintic for _, quintic in profiles) == (
            [68] * 4 if sign == 1 else [56, 60, 60, 68]
        )
        print(sign, a, e, f, profiles, flush=True)


if __name__ == "__main__":
    a, z, l = sp.symbols("a z l")
    for sign in (1, -1):
        check_quartic(sign, a, z, l)
    check_sections(z, l)
