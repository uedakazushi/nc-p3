"""Exact set-theoretic length-four degeneracy on both fixed curves at slopes (3,2)."""

import random
from fractions import Fraction

import sympy as sp
from flint import nmod_mat

from verify_fixed_exchange_isogeny import exchange
from verify_fixed_quartic_ranks import projected_entries
from verify_klein_elliptic_completion import coeffs, gram_twist, four_relations
from two_window_ext import PathRow

a, e, f = sp.symbols("a e f")
z, l = sp.Integer(3), sp.Integer(2)


def symbolic_blocks(sign, z=z, l=l, b=1):
    A = sp.Matrix([[a, b], [sign * z * b, sign * a]])
    B = sp.Matrix([[e, f], [sign * z * f, sign * e]])
    w, _, rev = gram_twist(A, B, z, l)
    blocks = []
    for side, cut in enumerate((w, rev, w.T, rev.T)):
        groups = []
        for color, r in enumerate((2, 2, 1, 1)):
            indices = [4 * i + (i ^ color) for i in range(4)]
            pairs = [4 * (i ^ color) + i for i in range(4)] if side % 2 else indices
            rows = cut.extract(pairs, indices[:r]).T
            high = rows[:, 4 - r :]
            assert not (high.free_symbols & {a, e, f}) and high.det() != 0
            rows = high.inv() * rows
            kernel = sp.eye(4 - r).col_join(-rows[:, : 4 - r])
            groups.append((rows.tolist(), kernel))
        blocks.append(groups)
    return blocks


def reduced_snapshot(matrix, sign):
    # The split branch fiber supplies a nonzero pivot; all minors below are over Q.
    prime = 307

    def red(v):
        return int(v.p) * pow(int(v.q), -1, prime) % prime

    aa = sp.Integer(102 if sign == 1 else 80)
    co, P, Q = coeffs(aa, sp.S.One, sign * z, sign * aa, z, l)
    values = {
        a: aa,
        f: red(co[0] - co[1] * P[1] / (2 * P[0])),
        e: red(co[4] - co[5] * Q[1] / (2 * Q[0])),
    }
    return [[red(v) for v in row] for row in matrix.subs(values).tolist()]


def pivot_columns(rows):
    mat, rank = nmod_mat(rows, 307).rref()
    return [next(j for j in range(mat.ncols()) if mat[i, j]) for i in range(rank)]


def minor_indices(matrix, snapshot, seed):
    rng = random.Random(seed)
    rows = list(range(matrix.rows))
    cols = list(range(matrix.cols))
    rng.shuffle(rows)
    rng.shuffle(cols)
    if seed % 3:
        rows.sort(key=lambda i: sum(v != 0 for v in matrix.row(i)))
        cols.sort(key=lambda j: sum(v != 0 for v in matrix.col(j)))
    J = [cols[j] for j in pivot_columns([[snapshot[i][j] for j in cols] for i in rows])]
    I = [rows[i] for i in pivot_columns([[snapshot[i][j] for i in rows] for j in J])]
    return I, J


def verify_rank(matrix, target, sign, N, D, curve, excluded, boundary):
    snapshot = reduced_snapshot(matrix, sign)
    common = None
    minors = []
    for seed in range(8):
        sub = matrix.extract(*minor_indices(matrix, snapshot, seed))
        assert sub.shape == (target, target)
        dm = sub.to_DM()
        poly = sp.Poly(dm.domain.to_sympy(dm.det()), a, e, f)
        degree = poly.degree(e)
        W = sp.Poly(
            sum(
                c * a**i * f ** (j + k) * N**j * D ** (degree - j)
                for (i, j, k), c in poly.terms()
            ),
            f,
            a,
        )
        norm = sp.Poly(curve.resultant(W).as_expr(), a, domain=sp.QQ)
        assert not norm.is_zero
        common = norm if common is None else common.gcd(norm)
        minors.append(W)
        if seed == 0:
            continue
        support = common.sqf_part()
        residual = support.exquo(support.gcd(excluded))
        if residual.degree() == 0:
            return len(minors)
        base = sp.groebner([residual.as_expr(), curve.as_expr()], f, a, domain=sp.QQ)
        ideal = sp.groebner(
            [p.as_expr() for p in base.polys]
            + [base.reduce(w.as_expr())[1] for w in minors],
            f,
            a,
            domain=sp.QQ,
        )
        rem = ideal.reduce(boundary)[1]
        if ideal.reduce(rem**2)[1] == 0:
            return len(minors)
    raise AssertionError((sign, target, residual, rem))


def finite_chart(sign):
    co, P, Q = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    N, D = exchange(sign, z, l) * sp.Matrix([a, 1])
    curve = sp.Poly(
        sp.factor(
            sign * z * co[1] * f * f + f - co[0] - sign * co[1] * a * a
        ).as_numer_denom()[0],
        f,
        a,
    )
    excluded = (a * a - z) * (2 * a - 1) * (a - 6) * D
    if sign == -1:
        excluded *= (2 * a + 1) * (a + 6)
    excluded *= (
        sp.factor(P[0]).as_numer_denom()[0] * sp.factor(Q[0]).as_numer_denom()[0]
    )
    excluded = sp.Poly(excluded, a, domain=sp.QQ)
    # xy det(B), with denominators cleared and nonzero scalar factors removed.
    boundary = (
        (a * a - z * f * f) * (N * N * f * f - z * D * D) * (N * N - z * D * D) * f * f
    )
    if sign == 1:
        boundary *= (a * a - z) * D * D + (N * N - z * D * D) * f * f
    blocks = symbolic_blocks(sign)
    for i in range(4):
        counts = [
            verify_rank(
                sp.Matrix(projected_entries(blocks, i, h)),
                target,
                sign,
                N,
                D,
                curve,
                excluded,
                boundary,
            )
            for h, target in enumerate((17, 16, 16, 16))
        ]
        print(sign, i, "maximal minors per character", counts, flush=True)


def boundary_fibers(sign):
    points = (
        [(sp.Rational(13, 4), 1), (sp.Rational(12, 13), 1)] if sign == 1 else [(0, 1)]
    ) + [(1, 0)]
    admissible = []
    for aa, bb in points:
        A = sp.Matrix([[aa, bb], [sign * z * bb, sign * aa]])
        co, P, Q = coeffs(*A, z, l)

        def roots(coeff):
            poly = sp.Poly.from_list(list(coeff), a)
            assert not poly.is_zero
            result = sp.roots(poly)
            assert sum(result.values()) == poly.degree()
            return result

        for x in roots(P):
            for y in roots(Q):
                if 2 * P[0] * x + P[1] + 2 * Q[0] * y + Q[1] != 0:
                    continue
                e0, f0 = co[4] + co[5] * y, co[0] + co[1] * x
                B = sp.Matrix([[e0, f0], [sign * z * f0, sign * e0]])
                if x * y * B.det() == 0:
                    continue
                admissible.append((aa, bb))
                w, _, rev = gram_twist(A, B, z, l)
                windows = [
                    [
                        {divmod(j, 4): Fraction(v) for j, v in enumerate(row) if v}
                        for _, row in group
                    ]
                    for group in four_relations(w, rev)
                ]
                assert all(
                    PathRow(windows, i, 5).dims[-2:] == [35, 56] for i in range(4)
                )
    assert admissible == (points if sign == 1 else [])


if __name__ == "__main__":
    for sign in (1, -1):
        finite_chart(sign)
        boundary_fibers(sign)
