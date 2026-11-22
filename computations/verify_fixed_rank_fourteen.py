"""Classify every admissible rank-fourteen point by two binary quadratics on a slope divisor."""

import sympy as sp

from verify_klein_elliptic_completion import coeffs


def verify():
    X, b, r, w = sp.symbols("X b r w")
    I = sp.I
    F = (1 + I) * X - I * b
    E = (1 - I) * b + I * X
    # a=rX, e=rE, f=F. Changing the root r or conjugating I covers all four lines.
    assert sp.expand(X + E - b - F) == 0 and sp.expand(X - E - I * (b - F)) == 0
    K = sp.expand(E * ((w * w + 1) * X - 2 * w * b) - F * (2 * w * X - (w * w + 1) * b))
    T = r * w * w + w * r * r - 4 * r * w + w + r
    affine_K = K.subs(b, 1)
    assert sp.factor(sp.discriminant(affine_K, X)) == 16 * w * (w - 1) ** 2
    for C, target in (
        (X * X - b * b, -4 * (w * w - 1) ** 2),
        (E * E - F * F, -4 * (w * w - 1) ** 2),
        (X * X - F * F, -4 * (w - 1) ** 2 * (w * w - 6 * w + 1)),
        (E * E - b * b, -4 * (w - 1) ** 2 * (w * w - 6 * w + 1)),
        (X - E, -4 * I * (w - 1) ** 2),
        (b - F, 4 * I * (w - 1) ** 2),
        (X + E, -8 * w),
        (b + F, -8 * w),
    ):
        assert sp.expand(sp.resultant(affine_K, C.subs(b, 1), X) - target) == 0
    assert sp.expand(T - r * (w * w - 6 * w + 1) - w * (r + 1) ** 2) == 0
    co, _, _ = coeffs(r * X, b, r * r * b, r * X, r * r, w / r)
    conditions = (
        (
            F - co[0] - co[1] * r * r * (X * X - F * F),
            -2 * I,
            w * w + (-2 + 2 * I) * w + 1,
        ),
        (
            r * E - co[4] - co[5] * r * r * (E * E - b * b),
            2 * I * r * r,
            w * w + (-2 - 2 * I) * w + 1,
        ),
    )
    base = sp.Poly(K, X, b)
    for equation, unit, extra in conditions:
        num, den = sp.fraction(sp.cancel(equation))
        assert sp.expand(den - w * (r * r + 1) * (X + b)) == 0
        form = sp.Poly(num, X, b)
        assert base.monoms() == form.monoms() == [(2, 0), (1, 1), (0, 2)]
        res = sp.resultant(affine_K, num.subs(b, 1), X)
        assert sp.expand(res - unit * (w + 1) ** 2 * extra * T * T) == 0
        A = base.coeffs()
        B = form.coeffs()
        for i in range(3):
            for j in range(i):
                assert sp.rem(A[i] * B[j] - A[j] * B[i], T, r) == 0
    t = sp.Symbol("t")
    root = -2 * (t + 3) / (t * t - 1)
    lam = -((t - 1) ** 2) / (t + 3) ** 2
    ww = sp.factor(lam * root)
    assert sp.cancel(T.subs({r: root, w: ww})) == 0
    assert (
        sp.cancel(
            (2 * root**2 * lam + root**2 - 5 * root + 2) / (root * (root - 1)) - t
        )
        == 0
    )
    assert (
        sp.cancel(ww * (t + 1) ** 2 * (t + 3) ** 2 - 2 * (t - 1) * (t + 1) * (t + 3))
        == 0
    )
    assert (root.subs(t, -sp.Rational(1, 3)), lam.subs(t, -sp.Rational(1, 3))) == (
        6,
        -sp.Rational(1, 4),
    )
    print(
        "All ordinary rank-14 fibers have four admissible points exactly on the stated divisor.",
        flush=True,
    )


if __name__ == "__main__":
    verify()
