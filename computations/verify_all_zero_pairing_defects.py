"""Hilbert defects on the eight irreducible total spaces of the simultaneous-zero chart."""

import sys

import sympy as sp

from flint_rational_function import rational_function_field
from two_window_ext import PathRow
from verify_klein_elliptic_completion import gram_twist, numerator


def parameters(case, z, l, x):
    # Table of the simultaneous-zero chart in row-major order, left then right.
    reciprocal, scaled, second = (case // 4, case // 2 % 2, case % 2)
    b = (1 / l, l * z, l, 1 / (l * z), 1 / (l * z), l, l / z, 1 / (l * z * z))[case]
    c = z * z * b if scaled else b
    d = 1 / z if reciprocal else z
    f0 = sp.cancel((1 + d) / (z + 1))
    g0 = z * f0
    e0 = sp.cancel((z * b + c) / (z + 1))
    y = sp.Symbol("y")
    f, g = (f0 - x / g0, g0) if second else (f0, g0 - x / f0)
    e, h = (e0 + y / e0, e0) if reciprocal ^ scaled ^ second else (e0, e0 + y / e0)
    a0 = sp.Matrix([[1, b], [c, d]])
    a1 = sp.Matrix([[e, f], [g, h]])
    n = numerator(e, f, g, h, z, l)
    diagonal = numerator(e, b, c, h, z, sp.S.One)
    off = numerator(sp.S.One, f, g, d, z, sp.S.One)
    determinant = a1.det()
    equations = [
        sp.factor(
            n[i] * (y if i in (0, 3) else x)
            - determinant * (diagonal[i] if i in (0, 3) else off[i])
        )
        for i in range(4)
    ]
    form = sp.Poly(next(eq for eq in equations if eq != 0), x, y)
    leading = form.coeff_monomial(x * y)
    b1, c1, d1 = [sp.factor(form.coeff_monomial(m) / leading) for m in (x, y, 1)]
    residual = sp.factor(b1 * c1 - d1)
    assert residual != 0
    root = sp.factor(residual / (x + c1) - b1)
    assert all(sp.cancel(eq.subs(y, root)) == 0 for eq in equations)
    return a0, a1.subs(y, root).applyfunc(sp.factor), residual


def relation_bases(a0, a1, z, l, sample, field):
    w, _, rev = gram_twist(a0, a1, z, l)
    windows = []
    for i, cut in enumerate((w, rev, w.T, rev.T)):
        fixed = cut.subs(sample)
        cols = list(fixed.to_DM().rref()[1])
        assert len(cols) == 6
        rules = []
        for j in cols:
            rule = {}
            for k, coeff in enumerate(cut[:, j]):
                if coeff:
                    pair = divmod(k, 4)
                    value = field.from_sympy(coeff)
                    if value:
                        rule[pair[::-1] if i % 2 else pair] = value
            rules.append(rule)
        windows.append(rules)
    return windows


def check(case):
    z, l, x = sp.symbols("z l x")
    a, b, residual = parameters(case, z, l, x)
    field = rational_function_field(("z", "l", "x"))
    sample = {z: 3, l: 2, x: 4}
    assert residual.subs(sample) != 0
    windows = relation_bases(a, b, z, l, sample, field)
    start = 2 if case in (2, 3, 4, 5) else 0
    dims = PathRow(windows, start, 6, field.one).dims
    assert dims == [1, 4, 10, 20, 35, 56, 92], (case, dims)
    print(
        f"Row {case // 2 + 1}, choice {case % 2 + 1}, start {start}: {dims}", flush=True
    )


if __name__ == "__main__":
    for i in map(int, sys.argv[1:]) if sys.argv[1:] else range(8):
        check(i)
