"""Outer-pair quotient test on the finite invertible charts at (3,2) mod 101."""

from itertools import product
from collections import Counter

import sympy as sp
from flint import nmod, nmod_mat

from verify_klein_elliptic_completion import coeffs, gram_twist
from verify_fixed_line_elliptic_curves import bad_sections


def relation_blocks(A, B, z, l, p):
    w, _, rev = gram_twist(A, B, z, l)
    blocks = []
    for side, cut in enumerate((w, rev, w.T, rev.T)):
        groups = []
        for color in range(4):
            indices = [4 * i + (i ^ color) for i in range(4)]
            pairs = [4 * (i ^ color) + i for i in range(4)] if side % 2 else indices
            matrix = cut.extract(pairs, indices).T
            rows = nmod_mat(
                4, 4, [int(v.p) * pow(int(v.q), -1, p) % p for v in matrix], p
            )
            reduced, rank = rows.rref()
            assert rank == (2, 2, 1, 1)[color]
            kernel, _ = rows.nullspace()
            groups.append((reduced.tolist()[:rank], kernel))
        blocks.append(groups)
    return blocks


def path_dimensions(blocks, start, p):
    left, middle, right = [blocks[(start + j) % 4] for j in range(3)]
    dims = []
    for h in range(4):
        columns = [
            (u, h ^ u, j, k)
            for u in range(4)
            for j in range(4 - len(left[u][0]))
            for k in range(4 - len(right[h ^ u][0]))
        ]
        rows = []
        for head, tail in product(range(4), repeat=2):
            for rule in middle[head ^ tail ^ h][0]:
                rows.append(
                    [
                        int(
                            rule[head ^ u]
                            * left[u][1][head, j]
                            * right[v][1][tail ^ v, k]
                        )
                        for u, v, j, k in columns
                    ]
                )
        dims.append(len(columns) - nmod_mat(rows, p).rank())
    return tuple(dims)


def audit(sign):
    p = 101
    z, l = sp.Integer(3), sp.Integer(2)
    zz, ll = nmod(3, p), nmod(2, p)
    defects = set()
    counts = Counter()
    for aa in range(p):
        a = nmod(aa, p)
        if a * a == zz:
            continue
        co, (P, m, n), (Q, k, _) = coeffs(a, nmod(1, p), sign * zz, sign * a, zz, ll)
        if not P or not Q:
            continue
        for s in range(p):
            if s * s % p != int(m * m - 4 * P * n):
                continue
            x, y = (s - m) / (2 * P), (-s - k) / (2 * Q)
            e, f = co[4] + co[5] * y, co[0] + co[1] * x
            if not x or not y or e * e == zz * f * f:
                continue
            A = sp.Matrix([[aa, 1], [sign * z, sign * aa]])
            B = sp.Matrix([[int(e), int(f)], [sign * z * int(f), sign * int(e)]])
            blocks = relation_blocks(A, B, z, l, p)
            dims = [path_dimensions(blocks, i, p) for i in range(4)]
            assert all(h in ((9, 10, 8, 8), (10, 10, 8, 8)) for h in dims)
            jumps = sum(h[0] == 10 for h in dims)
            counts[jumps] += 1
            if jumps:
                defects.add((aa, int(e), int(f)))
    expected = {
        (int(a) % p, int(e) % p, int(f) % p)
        for sig, a, e, f in bad_sections(zz, ll)
        if sig == sign
    }
    assert defects == expected
    assert counts == ({0: 97, 4: 4} if sign == 1 else {0: 96, 2: 8})
    print(sign, dict(counts), "exceptional signed triples", sorted(defects), flush=True)


if __name__ == "__main__":
    for sign in (1, -1):
        audit(sign)
