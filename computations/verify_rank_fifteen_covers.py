"""Biquadratic rank-drop covers, their elliptic ramification curves and split fibers."""

from itertools import product

import sympy as sp
from flint import nmod
from verify_fixed_branch_regular import check_window
from verify_fixed_exchange_isogeny import exchange
from verify_fixed_line_elliptic_curves import ordinary_unit
from verify_klein_elliptic_completion import coeffs


def projective_equations():
    a, b, z, lam, u = sp.symbols("a b z l u")
    N, D = exchange(1, z, lam) * sp.Matrix([a, b])
    e, f = u * N, u * D
    det = a * a - z * b * b
    G = lam * (z + 1) * u * det - a * b - e * f
    co = coeffs(a, b, z * b, a, z, lam)[0]
    assert (
        sp.cancel(
            f - co[0] - co[1] * (a * a - z * f * f) - D * G / (lam * det * (z + 1))
        )
        == 0
    )
    assert (
        sp.cancel(
            e - co[4] - co[5] * (e * e - z * b * b) - N * G / (lam * det * (z + 1))
        )
        == 0
    )
    assert sp.cancel(z * co[1] * f + co[5] * e) == 0
    assert sp.expand(a * f + b * e - (lam * lam * z + 1) * u * det) == 0
    assert sp.expand(a * e + z * b * f - 2 * lam * z * u * det) == 0


def verify_covers():
    n, z, l = sp.symbols("n z l")
    D = 4 * l * z
    K = (l + 1) * (l * z + 1)
    M = -(l - 1) * (l * z - 1)
    poly = D**2 * n**4 + (D**2 - z * K**2) * n * n + z * M * M
    a, e = D * n * (1 + n), D * n * (1 - n)
    b, f = K * n + M, K * n - M
    for x, y, s, t in ((a, b, e, f), (z * b, a, z * f, e)):
        for eq in (
            2 * l * z * (x * t + y * s) - (l * l * z + 1) * (x * s + z * y * t),
            2 * l * z * (x * y + s * t) - l * (z + 1) * (x * s + z * y * t),
        ):
            assert sp.rem(eq, poly, n).expand() == 0
        for numerator in (
            x * x - z * y * y,
            s * s - z * t * t,
            x * x - z * t * t,
            s * s - z * y * y,
        ):
            assert ordinary_unit(sp.resultant(poly, numerator, n), z, l)
    assert sp.rem(a * a + e * e - 2 * z * b * f, poly, n).expand() == 0
    assert (
        sp.rem(
            z * (b * b + f * f) - 2 * a * e - 4 * n * n * (z * K * K - D * D), poly, n
        ).expand()
        == 0
    )
    A = 1 - z * (K / D) ** 2
    B = z * (M / D) ** 2
    disc = sp.factor(A * A - 4 * B)
    numerator = sp.fraction(sp.factor(disc.subs(l, 2)))[0]
    assert sp.gcd(numerator, sp.diff(numerator, z)) == 1 and sp.degree(numerator, z) > 0
    assert numerator.subs(z, 0) != 0
    assert sp.factor(sp.discriminant(n**4 + A * n * n + B, n) - 16 * B * disc**2) == 0


def ramification_parameters(t, y):
    u = (t - 1 / t) / 2
    v = (t + 1 / t) / 4
    k = 4 / (1 - u * u) - 1
    n = (k + 1) * v
    r = (2 - (k - 1) ** 2 + (k - 1) * (k + 1) * u) / 2
    w = ((k + 1) ** 2 - 2 + (k + 1) * n * y / (2 * t)) / 2
    return r, w, k, n


def frame_on_branch(t, y, h):
    r, w, k, n = ramification_parameters(t, y)
    z = r * r
    a, b, e, f = 1 + n, (k - n) / r, 1 - n, (k + n) / r
    if h:
        a, b, e, f = z * b, a, z * f, e
    return z, w / r, a, b, e, f


def check_quartic_model():
    t, y = sp.symbols("t y")
    quartic = -(t**4) + 14 * t * t - 1
    r, w, k, n = ramification_parameters(t, y)

    def zero(expr):
        num, den = sp.fraction(sp.cancel(expr))
        assert (
            not sp.rem(num, y * y - quartic, y) and sp.rem(den, y * y - quartic, y) != 0
        )

    zero(r * r - (2 - (k - 1) ** 2) * r + 1)
    zero(w * w - ((k + 1) ** 2 - 2) * w + 1)
    zero(2 * n * n - k * k + 1)
    zero(t - (2 * n + (2 * r - 2 + (k - 1) ** 2) / (k - 1)) / (k + 1))
    zero(y - 2 * t * (2 * w + 2 - (k + 1) ** 2) / ((k + 1) * n))
    assert 256 * (14**2 + 12) ** 3 / sp.discriminant(quartic, t) == sp.Rational(
        35152, 9
    )


def check_sample():
    prime = 139
    k = nmod(2, prime)
    R, S = 2 - (k - 1) ** 2, (k + 1) ** 2 - 2
    signs = [
        sp.sqrt_mod(int(value), prime, all_roots=True)
        for value in (R * R - 4, S * S - 4, (k * k - 1) / 2)
    ]
    assert all(len(roots) == 2 for roots in signs)
    for rr, ww, nn in product(*signs):
        r, w, n = (R + nmod(rr, prime)) / 2, (S + nmod(ww, prime)) / 2, nmod(nn, prime)
        z, lam = r * r, w / r
        assert (
            z
            * lam
            * (z * z - 1)
            * (lam * lam - 1)
            * (lam * lam * z - 1)
            * (lam * lam * z * z - 1)
        )
        for h in (0, 1):
            a, b, e, f = 1 + n, (k - n) / r, 1 - n, (k + n) / r
            if h:
                a, b, e, f = z * b, a, z * f, e
            aa, bb, ee, ff, zz, ll = (sp.Integer(int(v)) for v in (a, b, e, f, z, lam))
            A = sp.Matrix([[aa, bb], [zz * bb, aa]])
            B = sp.Matrix([[ee, ff], [zz * ff, ee]])
            check_window(
                A, B, zz, ll, (2, int(r), int(w), int(n), h), prime, opposite_rank=15
            )


if __name__ == "__main__":
    projective_equations()
    verify_covers()
    check_quartic_model()
    check_sample()
    print("Biquadratic covers, their branch models and all 16 sample windows verified.")
