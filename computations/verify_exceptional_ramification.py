"""Exact change of exceptional slope, opposite pencils and evaluation blocks."""

from itertools import combinations, combinations_with_replacement, product

import sympy as sp

from verify_klein_elliptic_completion import gram_twist
from verify_klein_two_doublets import opposite
from verify_kummer_discriminant import discriminant_form


def check():
    variables, form = discriminant_form()
    a, b, c, d, z, l = variables
    A = sp.Matrix([[a, b], [c, d]])
    D = A.det()
    quadric = (a + d) * (z * b + c) - (z + 1) * (a * d + b * c)
    change = {a: z * b, b: a, c: z * d, d: c}
    assert (
        sp.factor(sp.hessian(quadric, variables[:4]).det())
        == (z - 1) ** 2 * (z + 1) ** 2
    )
    for slope, sign, scale in ((1, 1, 1), (-1, -1, 1), (1 / z, 1, z), (-1 / z, -1, z)):
        target = quadric.xreplace({b: sign * b, c: sign * c})
        assert sp.cancel(form.subs(l, slope) - (z - 1) ** 2 * target**2 / scale**2) == 0
    assert sp.expand(quadric.xreplace(change) - z * quadric) == 0
    no = (a - z * d) * (z * a - d) * (b - c) * (z * z * b - c)
    assert sp.cancel(no.xreplace(change) / no) == z * z
    X = sp.Matrix([[1, 1], [1, z]])
    Y = sp.Matrix([[1, 1], [z, 1]])
    L = -A.inv().T * Y.T
    K = -A.inv() * X
    assert sp.cancel(L.det() + K.det()) == 0 and sp.cancel(K.det() - (z - 1) / D) == 0
    factors = (a - b, a - c, a - z * b, z * a - c, b - d, c - d, z * b - d, c - z * d)
    for factor in factors:
        coord = next(v for v in variables[:4] if sp.diff(factor, v) != 0)
        sub = {coord: sp.solve(factor, coord)[0]}
        _, residual = sp.div(
            sp.expand(no.subs(sub)),
            sp.expand(quadric.subs(sub)),
            *[v for v in variables[:4] if v != coord],
            domain=sp.QQ.frac_field(z),
        )
        assert residual == 0
    w, _, _ = gram_twist(A, A, z, sp.S.One)
    opp = opposite(w, True)
    pairs = list(combinations_with_replacement(range(4), 2))
    blocks = []
    for h in range(4):
        row = [i for i, (j, k) in enumerate(pairs) if j ^ k == h]
        col = [i for i, (j, k) in enumerate(pairs) if j ^ k == (h ^ 1)]
        blocks.append(opp[row, col])
    assert [sp.factor(m.det()) for m in blocks[2:]] == [z - 1, 1 - z]
    differences = [
        (
            a - d,
            (a - b) * (a - c) / D,
            (a - z * b) * (z * a - c) / D,
            (b - d) * (c - d) / D,
            (c - z * d) * (z * b - d) / D,
            -(z - 1) * (z * b - c) / D,
        ),
        (
            z * b - c,
            (a - b) * (z * b - d) / D,
            (a - z * b) * (b - d) / D,
            (c - d) * (z * a - c) / D,
            -(a - c) * (z * d - c) / D,
            -(a - d) * (z - 1) / D,
        ),
    ]
    for block, known in zip([blocks[0], blocks[1].T], differences):
        for rows, target in zip(combinations(range(4), 2), known):
            minor = sp.factor(block[list(rows), :].det())
            assert sp.cancel(minor / target) in (1, -1)
    new = A.xreplace(change)
    positive, _, _ = gram_twist(A, sp.Matrix([[d, c / z], [z * b, a]]), z, 1 / z)
    standard, _, _ = gram_twist(new, new, z, sp.S.One)
    for i, j, k, n in product(range(4), repeat=4):
        original = positive[4 * i + j, 4 * (k ^ 1) + (n ^ 1)]
        assert (
            sp.cancel(
                standard[4 * i + j, 4 * k + n]
                - (1, 1, 1 / z, 1)[k] * (z, 1, 1, 1)[n] * original
            )
            == 0
        )
    print(
        "Four double quadrics, both diagonal pencils, opposite determinants, and all 256 coefficients of the positive slope change verified."
    )


if __name__ == "__main__":
    check()
