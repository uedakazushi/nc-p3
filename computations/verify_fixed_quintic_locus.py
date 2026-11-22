"""All length-five ranks off the quartic defect in the (3,2) fiber."""

from itertools import product

import sympy as sp
from flint import nmod_mat

from verify_fixed_exchange_isogeny import exchange
from verify_fixed_geometric_quartic_locus import (
    a,
    e,
    f,
    symbolic_blocks,
    minor_indices,
    reduced_snapshot,
    boundary_fibers,
)
from verify_klein_elliptic_completion import coeffs
from verify_uniform_fixed_quartic_support import (
    u,
    polynomial,
    determinant,
    cleared_rows,
)


def middle_relations(blocks, start, h):
    first, left, right, last = [blocks[(start + j) % 4] for j in range(4)]
    columns = [
        (u, v, j, k)
        for u, v in product(range(4), repeat=2)
        for j in range(4 - len(first[u][0]))
        for k in range(4 - len(last[v][0]))
    ]
    rows = []
    for arrow, degree in product(range(4), repeat=2):
        outer = h ^ arrow ^ degree
        for rule in left[degree][0]:
            for t in range(4 - len(last[outer][0])):
                rows.append(
                    [
                        rule[arrow ^ u] * first[u][1][arrow, j]
                        if v == outer and k == t
                        else 0
                        for u, v, j, k in columns
                    ]
                )
        for rule in right[degree][0]:
            for t in range(4 - len(first[outer][0])):
                rows.append(
                    [
                        rule[arrow ^ v ^ degree] * last[v][1][arrow ^ v, k]
                        if u == outer and j == t
                        else 0
                        for u, v, j, k in columns
                    ]
                )
    return sp.Matrix(rows)


def expression(poly):
    return sum(int(c) * a**i * u**j for (i, j, _, _), c in poly.to_dict().items())


def projection_table(det, curve, excluded):
    table = {}
    for g, _ in det.factor()[1]:
        norm = curve.resultant(g, 1) if g.degrees()[1] else g
        assert norm
        for factor, _ in norm.factor()[1]:
            key = expression(factor)
            if key not in excluded:
                table.setdefault(key, []).append(expression(g))
    return table


def check_residual(point, irreducible, matrix):
    # The residual scheme projects onto a single irreducible number-field orbit.
    # One nonzero modular minor suffices to certify the rank on that entire orbit.
    curve = sp.Poly(irreducible, a, domain=sp.QQ)
    assert curve.degree() == 8 and curve.is_irreducible
    linear, base = point.polys
    assert sp.Poly(base.as_expr(), a, domain=sp.QQ).monic() == curve.monic()
    U = sp.expand(u - linear.as_expr())
    assert not U.has(u)
    p = 1009
    aa = next(i for i in range(p) if curve.eval(i) % p == 0)
    uu = U.subs(a, aa)

    def red(v):
        return int(v.p) * pow(int(v.q), -1, p) % p

    at = matrix.subs({a: aa, u: red(uu)})
    assert nmod_mat(at.rows, at.cols, [red(v) for v in at], p).rank() == 86
    return curve.monic().as_expr(), U


def verify_position(matrix, curve, N, D, excluded, boundary, sign):
    assert matrix.shape == (120, 100)
    snapshot = reduced_snapshot(matrix, sign)
    cache = {v: sp.factor(v.subs({e: N * u, f: D * u})) for v in set(matrix)}
    restricted = matrix.applyfunc(cache.get)
    tables = []
    for seed in (1, 2):
        sub = restricted.extract(*minor_indices(matrix, snapshot, seed))
        assert sub.shape == (86, 86)
        det = determinant(cleared_rows(sub, lambda den: not den.free_symbols))
        assert det
        tables.append(projection_table(det, curve, excluded))
    extra = set()
    for key in tables[0].keys() & tables[1].keys():
        base = sp.groebner([key, expression(curve)], u, a, domain=sp.QQ)
        conditions = []
        for table in tables:
            value = sp.S.One
            for factor in table[key]:
                value = base.reduce(value * base.reduce(factor)[1])[1]
            conditions.append(value)
        ideal = sp.groebner(
            [g.as_expr() for g in base.polys] + conditions, u, a, domain=sp.QQ
        )
        residue = ideal.reduce(boundary)[1]
        if ideal.reduce(residue**2)[1] != 0:
            extra.add(check_residual(ideal, key, restricted))
    return extra


def check(sign):
    z, l = sp.Integer(3), sp.Integer(2)
    N, D = exchange(sign, z, l) * sp.Matrix([a, 1])
    co, _, _ = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    relation = sp.factor(
        sign * z * co[1] * (D * u) ** 2 + D * u - co[0] - sign * co[1] * a * a
    )
    (curve,) = [
        polynomial(g)
        for g, _ in sp.factor_list(relation.as_numer_denom()[0])[1]
        if g.has(u)
    ]
    exceptions = [
        a,
        a * a - z,
        a - 6,
        2 * a - 1,
        D,
        sp.factor(co[1]).as_numer_denom()[0],
        sp.factor(co[5]).as_numer_denom()[0],
    ]
    if sign == -1:
        exceptions.extend((a + 6, 2 * a + 1))
    excluded = {
        expression(g) for poly in exceptions for g, _ in polynomial(poly).factor()[1]
    }
    H = (a * a - z * (D * u) ** 2) * ((N * u) ** 2 - z) * (N * N - z * D * D) * u
    if sign == 1:
        H *= a * a - z + (N * N - z * D * D) * u * u
    blocks = symbolic_blocks(sign)
    residual = set()
    for start in (0, 1):
        cases = [
            verify_position(
                middle_relations(blocks, start, h), curve, N, D, excluded, H, sign
            )
            for h in range(4)
        ]
        residual.update(*cases)
        print(
            sign,
            start,
            "residual number fields checked per character",
            [len(points) for points in cases],
            flush=True,
        )
    # The p or q endpoints and infinity were excluded in projection. The positive
    # a=0 point exchanges with the p=0 endpoint and the negative one is inadmissible.
    boundary_fibers(sign)
    return residual


if __name__ == "__main__":
    for sign in (1, -1):
        check(sign)
