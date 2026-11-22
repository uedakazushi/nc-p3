"""A single homogeneous frame on the rank-14 double cover and all ordinary F101 residue points."""

from collections import Counter
import sympy as sp
from flint import nmod

from verify_fixed_branch_regular import check_window


def frame(s, i):
    X = i * (s * s - (1 + i) * s - 1)
    b = s * s + (1 - i) * s - 1
    E, F = (1 - i) * b + i * X, (1 + i) * X - i * b
    return X, b, E, F


def slopes(t):
    r = -2 * (t + 3) / (t * t - 1)
    lam = -((t - 1) ** 2) / (t + 3) ** 2
    return r, lam


def verify_frame():
    s = sp.Symbol("s")
    i = sp.I
    X, b, E, F = frame(s, i)
    w = s * s
    assert (
        sp.expand(E * ((w * w + 1) * X - 2 * w * b) - F * (2 * w * X - (w * w + 1) * b))
        == 0
    )
    assert sp.expand(X / i - b + 2 * s) == 0 and X.subs(s, 0) != 0
    assert sp.expand(E + X.subs(s, -s)) == sp.expand(F + b.subs(s, -s)) == 0
    for value, expected in (
        (X * X - b * b, -2 * (s * s - 1) * (s - i) ** 2),
        (E * E - F * F, -2 * (s * s - 1) * (s + i) ** 2),
        (X * X - F * F, -2 * (s * s - 1) * (s * s - 2 * s - 1)),
        (E * E - b * b, -2 * (s * s - 1) * (s * s + 2 * s - 1)),
    ):
        assert sp.expand(value - expected) == 0


def check_residue_points():
    p = 101
    fourth_roots = [v for v in range(p) if v * v % p == p - 1]
    assert fourth_roots == [10, 91]
    # The first choice of two degree-zero characters fails the Artinian test
    # at these four residues; the second pair certifies the same center.
    alternatives = {(10, 58, 0, 1), (91, 24, 20, 1), (91, 43, 63, 1), (91, 97, 74, 1)}
    used = set()
    counts = Counter()
    for t0 in range(p):
        t = nmod(t0, p)
        if not (t * t - 1) * (t + 3):
            continue
        r, lam = slopes(t)
        z = r * r
        if (
            not z
            * lam
            * (z * z - 1)
            * (lam * lam * z - 1)
            * (lam * lam - 1)
            * (lam * lam * z * z - 1)
        ):
            continue
        for s0 in range(p):
            s = nmod(s0, p)
            if s * s != lam * r:
                continue
            for i0 in fourth_roots:
                X, b, E, F = frame(s, nmod(i0, p))
                scale = b if b else X
                X, b, E, F = [v / scale for v in (X, b, E, F)]
                key = (i0, t0, int(X), int(b))
                if key in alternatives:
                    options = dict(degree_zero_choices=((5, 3, 1), (7, 2, 13)))
                    used.add(key)
                else:
                    options = {}
                # Lift the coefficients to ordinary integers before creating
                # the four cuts; the only target of this test is the residue field.
                aa, bb, ee, ff = map(int, (r * X, b, r * E, F))
                zz, ll = map(sp.Integer, (int(z), int(lam)))
                A = sp.Matrix([[aa, bb], [zz * bb, aa]])
                B = sp.Matrix([[ee, ff], [zz * ff, ee]])
                check_window(A, B, zz, ll, key, p, opposite_rank=14, **options)
                counts[i0] += 1
    assert used == alternatives and counts == {10: 96, 91: 96}
    print(dict(counts), "ordinary rational residue windows certified", flush=True)


if __name__ == "__main__":
    verify_frame()
    check_residue_points()
