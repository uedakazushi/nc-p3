"""Eliminate the single-zero locus and certify a seventh-degree defect over each rational chart."""

from itertools import product
import sys

import sympy as sp

from flint_rational_function import rational_function_field
from two_window_ext import PathRow
from verify_all_zero_pairing_defects import relation_bases
from verify_klein_elliptic_completion import coeffs, numerator


def boundary_equations(which, second_ratio, second_slope, z, l):
    a, b, c, d = sp.symbols("a b c d")
    if which == "P":
        a = sp.S.One
        d = 1 / z if second_ratio else z
        f0 = sp.cancel((1 + d) / (1 + z))
        g0 = z * f0
        values = (1 / g0, z * z / g0) if second_slope else (1 / f0, 1 / f0)
        unknowns = (b, c)
        indices = (0, 3)
    else:
        b = sp.S.One
        c = z * z if second_ratio else sp.S.One
        e0 = sp.cancel((z + c) / (1 + z))
        values = (z / e0, 1 / e0) if second_slope else (1 / e0, z / e0)
        unknowns = (a, d)
        indices = (1, 2)
    n = numerator(a, b, c, d, z, l)
    determinant = a * d - b * c
    e, f = [sp.factor(n[i] - v * determinant) for v, i in zip(values, indices)]
    return (a, b, c, d), unknowns, e, f, sp.factor(values[1] * e - values[0] * f)


def chart_unit(value, z, l):
    inverted = z * l * (z * z - 1) * (l * l * z - 1) * (l * l - 1) * (l * l * z * z - 1)
    for poly in sp.factor(value).as_numer_denom():
        if not poly:
            return False
        if any(sp.rem(inverted, p, z, l) != 0 for p, _ in sp.factor_list(poly)[1]):
            return False
    return True


def ordinary_slopes():
    z, l = sp.symbols("z l")
    for which, t, k in product(("P", "Q"), range(2), range(2)):
        entries, (b, c), e, f, linear = boundary_equations(which, t, k, z, l)
        coefficient = sp.diff(linear, c)
        assert chart_unit(coefficient, z, l)
        root = sp.factor(-linear.subs(c, 0) / coefficient)
        remaining = sp.factor(e.subs(c, root))
        polynomial = sp.Poly(remaining, b)
        assert polynomial.degree() == 2 and chart_unit(polynomial.LC(), z, l)
        factors = sp.factor_list(sp.factor(remaining).as_numer_denom()[0], b)[1]
        roots = []
        for fac, power in factors:
            poly = sp.Poly(fac, b)
            assert poly.degree() == 1
            roots.extend([sp.cancel(-poly.TC() / poly.LC())] * power)
        assert len(roots) == 2
        for value in roots:
            sub = {b: value, c: root.subs(b, value)}
            args = [v.subs(sub) for v in entries]
            _, p, q = coeffs(*args, z, l)
            assert all(sp.cancel(eq) == 0 for eq in (*p, *q))
    print(
        "Outside the four exceptional slopes either zero polynomial forces both to vanish.",
        flush=True,
    )


def exceptional(which, case):
    z, u, v = sp.symbols("z u v")
    t, k, reciprocal, negative = [case >> j & 1 for j in (3, 2, 1, 0)]
    sign = (-1) ** negative
    l = sign / z if reciprocal else sp.Integer(sign)
    entries, unknowns, e, f, linear = boundary_equations(which, t, k, z, l)
    (fixed,) = [q for q in unknowns if linear.has(q)]
    (free,) = [q for q in unknowns if q != fixed]
    multiplier = sp.diff(linear, fixed)
    for poly in sp.factor(multiplier).as_numer_denom():
        assert all(
            sp.rem(z * (z * z - 1), fac, z) == 0 for fac, _ in sp.factor_list(poly)[1]
        )
    substitution = {fixed: sp.factor(-linear.subs(fixed, 0) / multiplier), free: u}
    assert sp.cancel(e.subs(substitution)) == sp.cancel(f.subs(substitution)) == 0
    a, b, c, d = [sp.factor(q.subs(substitution)) for q in entries]
    cs, p, q = coeffs(a, b, c, d, z, l)
    zero, other = (p, q) if which == "P" else (q, p)
    assert all(sp.cancel(value) == 0 for value in zero)
    assert (
        sp.cancel(other[0]) != 0
        and sp.cancel(other[1] ** 2 - 4 * other[0] * other[2]) == 0
    )
    double = sp.factor(-other[1] / (2 * other[0]))
    x, y = (v, double) if which == "P" else (double, v)
    f0, f1, g0, g1, e0, e1, h0, h1 = cs
    E, F, G, H = [
        sp.factor(value)
        for value in (e0 + e1 * y, f0 + f1 * x, g0 + g1 * x, h0 + h1 * y)
    ]
    a0 = sp.Matrix([[a, b], [c, d]])
    a1 = sp.Matrix([[E, F], [G, H]])
    n = numerator(E, F, G, H, z, l)
    diagonal = numerator(E, b, c, H, z, sp.S.One)
    off = numerator(a, F, G, d, z, sp.S.One)
    for i in range(4):
        coord, component = (y, diagonal[i]) if i in (0, 3) else (x, off[i])
        assert sp.cancel(n[i] * coord - a1.det() * component) == 0
    field = rational_function_field(("z", "u", "v"))
    windows = relation_bases(a0, a1, z, l, {z: 3, u: 5, v: 7}, field)
    dims = PathRow(windows, 0, 7, field.one).dims
    assert dims[:7] == [1, 4, 10, 20, 35, 56, 84] and dims[7] > 120
    print(which, case, "start 0, degree seven", dims[7], flush=True)


if __name__ == "__main__":
    if not sys.argv[1:] or sys.argv[1:] == ["eliminate"]:
        ordinary_slopes()
    jobs = sys.argv[1:] or [
        which + str(i) for which, i in product(("P", "Q"), range(16))
    ]
    for job in jobs:
        if job != "eliminate":
            exceptional(job[0], int(job[1:]))
