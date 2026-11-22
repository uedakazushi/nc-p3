"""Two-position exchange on the fixed elliptic curves and their degree-two quotients."""

import sympy as sp

from verify_fixed_line_elliptic_curves import (
    branch_polynomial,
    bad_sections,
    ordinary_unit,
)
from verify_klein_elliptic_completion import coeffs


def exchange(sign, z, l):
    return (
        sp.Matrix([[2 * l * z, -z * (l * l * z + 1)], [l * l * z + 1, -2 * l * z]])
        if sign == 1
        else sp.Matrix([[0, -z], [1, 0]])
    )


def check(sign, a, z, l):
    u, f = sp.symbols("u f")
    m = exchange(sign, z, l)
    A, B, C = m[0, 0], m[0, 1], m[1, 0]
    N, D = m * sp.Matrix([a, 1])
    b = N / D
    quartic = branch_polynomial(sign, a, z, l)
    co, _, _ = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    f0, f1 = co[:2]
    assert sp.cancel(b + z * f1 / co[5]) == 0
    assert (m * m + m.det() * sp.eye(2)).applyfunc(sp.expand) == sp.zeros(2)
    transform = sum(c * N**i * D ** (4 - i) for (i,), c in quartic.terms())
    assert sp.expand(transform - quartic.as_expr() * m.det() ** 2) == 0
    assert ordinary_unit(sp.resultant(quartic.as_expr(), a * D - N, a), z, l)
    # Dividing the exchanged pair by f normalizes the first off-diagonal to 1.
    # The second off-diagonal becomes 1/f. Compare the two signed coordinates.
    scale_ratio = sp.factor((b * b - z) * D**2 / (a * a - z))
    difference = sp.factor(
        scale_ratio * (f + 2 * sign * z * f1.subs(a, b))
        - m.det() * f * (1 + 2 * sign * z * f1 * f)
    )
    relation = sign * z * f1 * f**2 + f - f0 - sign * f1 * a * a
    assert sp.cancel(sp.rem(difference, relation, f)) == 0
    # On C != 0 recover the quadratic H directly by comparing its top three terms.
    U_num = C * a * a + B
    R = u * u - 4 * (A * u + B) / C
    residual = sp.Poly(C * C * quartic.as_expr(), a)
    H = 0
    for power in (2, 1, 0):
        term = sp.Poly(U_num**power * D ** (2 - power), a)
        coeff = sp.factor(residual.nth(term.degree()) / term.LC())
        H += coeff * u**power
        residual = sp.Poly(sp.expand(residual.as_expr() - coeff * term.as_expr()), a)
    assert residual.is_zero
    assert sp.factor(R.subs(u, U_num / D) - (2 * a - U_num / D) ** 2) == 0
    rows = [(a0, e, h) for s, a0, e, h in bad_sections(z, l) if s == sign]
    for a0, e, h in rows:
        assert any(
            all(sp.cancel(d) == 0 for d in (aa - e / h, ee - a0 / h, hh - 1 / h))
            for aa, ee, hh in rows
        )


if __name__ == "__main__":
    a, z, l = sp.symbols("a z l")
    for sign in (1, -1):
        check(sign, a, z, l)
    print(
        "Both signed lifts preserve the differential; all twelve sections are paired."
    )
