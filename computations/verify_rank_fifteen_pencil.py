"""The rank-fifteen pencil, its branch divisors and split fibers on both genus-three levels."""

from itertools import product

import sympy as sp
from flint import nmod
from verify_fixed_branch_regular import check_window


def identities():
    c, K = sp.symbols("c K")
    v = (K * K - 1) / (1 + c * c)
    s = c + 1 / c
    for sign in (1, -1):
        midpoint = K + sign * c * v
        assert (
            sp.cancel(
                midpoint**2 - 1 - v * c * c / (1 + c * c) * ((K + sign * s) ** 2 - 1)
            )
            == 0
        )
    assert (
        sp.cancel((1 - K * K) ** 2 - 4 * c * c * v * v - (1 - c * c) ** 2 * v * v) == 0
    )
    x = sp.Symbol("x")
    pair = (x * x - 1, (x + s) ** 2 - 1, (x - s) ** 2 - 1)
    for a, b, expected in (
        (0, 1, s * s * (s * s - 4)),
        (0, 2, s * s * (s * s - 4)),
        (1, 2, 16 * s * s * (s * s - 1)),
    ):
        assert sp.cancel(sp.resultant(pair[a], pair[b], x) - expected) == 0
    for value, genus in ((3, 5), (2, 1), (1, 3), (-1, 3), (-2, 1)):
        roots = [{-1, 1}, {-value - 1, -value + 1}, {value - 1, value + 1}]
        for mask in range(1, 8):
            support = set()
            for j, row in enumerate(roots):
                if mask >> j & 1:
                    support.symmetric_difference_update(row)
            assert support
        pts = set.union(*roots)
        assert 2 * genus - 2 == -16 + 4 * len(pts)
        assert not set.intersection(*roots)


def ordinary_boundary():
    K, s = sp.symbols("K s")
    det = K * K + s * K - 1
    assert sp.discriminant(det, K) == s * s + 4
    for eq, expected in (
        (K * K - 1, -s * s),
        ((K + s) ** 2 - 1, -s * s),
        ((K - s) ** 2 - 1, s * s * (4 * s * s - 9)),
    ):
        assert sp.expand(sp.resultant(eq, det, K) - expected) == 0
    for val, count in (
        (1, 44),
        (2, 48),
        (sp.I * 2, 48),
        (sp.Rational(3, 2), 48),
        (3, 56),
    ):
        branched = {sp.Integer(-1), sp.Integer(1), -val - 1, -val + 1, val - 1, val + 1}
        extra = {sp.oo, sp.Integer(0), *sp.solve(det.subs(s, val), K)}
        assert 4 * len(branched) + 8 * len(extra - branched) == count


def genus_three_model():
    t, x = sp.symbols("t x")
    denominator = 2 * t * (t * t - 1)
    K = -((t * t + 1) ** 2) / denominator
    U = (t * t - 2 * t - 1) / (t * t + 1)
    V = (1 - 2 * t - t * t) / (t * t + 1)
    F = t**8 + 14 * t**4 + 1
    A, B = K * U, K * V
    assert sp.cancel(K * K - 1 - F / denominator**2) == 0
    assert sp.cancel(A * A - K * (K + 2)) == 0
    assert sp.cancel(B * B - K * (K - 2)) == 0
    assert sp.cancel((V - 1) / (U - 1) - t) == 0
    assert sp.cancel((V + 1) / (U + 1) + 1 / t) == 0
    assert sp.cancel(denominator.subs(t, -1 / t) - denominator / t**4) == 0
    quotient = (x * x + 4) * (x**4 + 4 * x * x + 16)
    assert sp.cancel(quotient.subs(x, t - 1 / t) - (t + 1 / t) ** 2 * F / t**4) == 0
    assert sp.gcd(F, sp.diff(F, t)) == sp.gcd(quotient, sp.diff(quotient, x)) == 1


def sample(prime):
    cs = sp.sqrt_mod(-3, prime, all_roots=True)
    assert len(cs) == 2
    for cc in cs:
        c = (1 + nmod(cc, prime)) / 2
        assert c * c - c + 1 == 0
        K = nmod(3, prime)
        n2 = (K * K - 1) / (1 + c * c)
        R, W = K + c * n2, K - c * n2
        choices = [
            sp.sqrt_mod(int(v), prime, all_roots=True)
            for v in (n2, R * R - 1, W * W - 1)
        ]
        assert all(len(row) == 2 for row in choices)
        for nn, rr, ww in product(*choices):
            n = nmod(nn, prime)
            r = R + nmod(rr, prime)
            w = W + nmod(ww, prime)
            z = r * r
            l = w / r
            assert (
                z
                * l
                * (z * z - 1)
                * (l * l - 1)
                * (l * l * z - 1)
                * (l * l * z * z - 1)
            )
            k = (l + 1) * (l * z + 1) / (4 * l * z)
            m = -(l - 1) * (l * z - 1) / (4 * l * z)
            assert (
                K == r * k
                and c * n * n == r * m
                and n**4 + (1 - z * k * k) * n * n + z * m * m == 0
            )
            a, b, e, f = 1 + n, k + m / n, 1 - n, k - m / n
            for h, frame in enumerate(((a, b, e, f), (z * b, a, z * f, e))):
                aa, bb, ee, ff = map(lambda x: sp.Integer(int(x)), frame)
                zz, ll = sp.Integer(int(z)), sp.Integer(int(l))
                A = sp.Matrix([[aa, bb], [zz * bb, aa]])
                B = sp.Matrix([[ee, ff], [zz * ff, ee]])
                check_window(
                    A,
                    B,
                    zz,
                    ll,
                    ("genus-three", int(c), int(r), int(w), int(n), h),
                    prime,
                    opposite_rank=15,
                )


if __name__ == "__main__":
    identities()
    ordinary_boundary()
    genus_three_model()
    sample(139)
    print("Pencil, genus-three octic and all 32 split windows verified.", flush=True)
