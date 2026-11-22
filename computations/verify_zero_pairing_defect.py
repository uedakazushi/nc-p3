"""Exact function-field Hilbert obstructions on the first zero-quadratics boundary."""

import sympy as sp
from flint_rational_function import rational_function_field
from two_window_ext import PathRow
from verify_klein_elliptic_completion import gram_twist, four_relations


def check(names, parameters, sample, length, expected):
    field = rational_function_field(names)
    z, l, x = parameters
    y = sp.cancel(
        (z - 1) ** 2
        * (l**2 * x - (l**2 * z - 1))
        / (l**2 * ((l**2 * z - 1) * x - (z - 1) ** 2))
    )
    a = sp.Matrix([[1, 1 / l], [1 / l, z]])
    # In this row f1=e1=0, g1=-1/f0 and h1=1/e0.
    f = (1 + z) / (z + 1)
    e = 1 / l
    b = sp.Matrix([[e, f], [z * f - x / f, e + y / e]])
    w, _, rev = gram_twist(a, b, z, l)
    windows = []
    for i, m in enumerate((w, rev, w.T, rev.T)):
        selected = list(m.subs(sample).to_DM().rref()[1])
        assert len(selected) == 6
        rules = []
        for j in selected:
            rule = {}
            for k, c in enumerate(m[:, j]):
                if c:
                    pair = divmod(k, 4)
                    val = field.from_sympy(c)
                    if val:
                        rule[pair[::-1] if i % 2 else pair] = val
            rules.append(rule)
        windows.append(rules)
    dims = [PathRow(windows, i, length, field.one).dims[-1] for i in range(4)]
    assert dims == expected, (names, dims)
    print("Function field", names, "degree", length, "dimensions", dims, flush=True)


if __name__ == "__main__":
    z, l, x = sp.symbols("z l x")
    check(("x",), (sp.Integer(3), sp.Integer(2), x), {x: 4}, 6, [92, 84, 85, 84])
    check(("z", "l"), (z, l, z - 1), {z: 3, l: 2}, 4, [35, 36, 36, 35])
    a = sp.Matrix([[1, sp.Rational(1, 2)], [sp.Rational(1, 2), 3]])
    b = sp.Matrix([[sp.Rational(1, 2), 1], [-1, sp.Rational(3, 4)]])
    w, _, rev = gram_twist(a, b, sp.Integer(3), sp.Integer(2))
    windows = [
        [{divmod(i, 4): c for i, c in enumerate(row) if c} for _, row in rules]
        for rules in four_relations(w, rev)
    ]
    assert PathRow(windows, 0, 6).dims == [1, 4, 10, 20, 35, 56, 92]
