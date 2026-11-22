"""Eliminate the signed coordinate from 32 generic minors over the integer polynomial ring."""

import sympy as sp
from flint import fmpz_mpoly_ctx

from verify_fixed_exchange_isogeny import exchange
from verify_fixed_geometric_quartic_locus import (
    a,
    e,
    f,
    minor_indices,
    reduced_snapshot,
    symbolic_blocks,
)
from verify_fixed_line_elliptic_curves import ordinary_unit
from verify_fixed_quartic_ranks import projected_entries
from verify_klein_elliptic_completion import coeffs
from verify_opposite_fixed_determinants import section

u, z, l = sp.symbols("u z l")
ring = fmpz_mpoly_ctx.get(("a", "u", "z", "l"))


def polynomial(value):
    return ring.from_dict(
        {exp: int(c) for exp, c in sp.Poly(value, a, u, z, l, domain=sp.ZZ).terms()}
    )


def determinant(rows):
    size = len(rows)
    previous = ring.constant(1)
    sign = 1
    for k in range(size - 1):
        choices = [(i, j) for i in range(k, size) for j in range(k, size) if rows[i][j]]
        if not choices:
            return ring.constant(0)

        def cost(cell):
            i, j = cell
            return (
                len(rows[i][j]),
                sum(bool(x) for x in rows[i][k:])
                * sum(bool(rows[t][j]) for t in range(k, size)),
            )

        i, j = min(choices, key=cost)
        if i != k:
            rows[i], rows[k] = rows[k], rows[i]
            sign = -sign
        if j != k:
            for row in rows:
                row[j], row[k] = row[k], row[j]
            sign = -sign
        pivot = rows[k][k]
        for i in range(k + 1, size):
            for j in range(k + 1, size):
                value = rows[i][j] * pivot - rows[i][k] * rows[k][j]
                rows[i][j] = value / previous if value else value
        previous = pivot
    return sign * rows[-1][-1]


def cleared_rows(mat, allowed_denominator):
    rows = []
    for row in mat.tolist():
        den = sp.lcm([sp.denom(v) for v in row])
        assert allowed_denominator(den)
        rows.append([polynomial(sp.factor(den * v)) for v in row])
    return rows


def candidate_factors(sign):
    factors = [
        z,
        l,
        z - 1,
        z + 1,
        l - 1,
        l + 1,
        z * l - 1,
        z * l + 1,
        z * l * l - 1,
        a,
        a * a - z,
        a - z * l,
        a * l - 1,
    ]
    if sign == -1:
        factors.extend((a + z * l, a * l + 1))
    else:
        factors.extend((2 * a * l - z * l * l - 1, a * z * l * l + a - 2 * z * l))
        r = sp.Symbol("r")
        for epsilon in (1, -1):
            A, _ = section(r, l, epsilon)
            factors.append(sp.resultant(A[0, 1] * a - A[0, 0], r * r - z, r))
    return {str(g) for factor in factors for g, _ in polynomial(factor).factor()[1]}


def check(sign):
    N, D = exchange(sign, z, l) * sp.Matrix([a, 1])
    co, _, _ = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    relation = sp.factor(
        sign * z * co[1] * (D * u) ** 2 + D * u - co[0] - sign * co[1] * a * a
    )
    numer, denom = relation.as_numer_denom()
    (signed,) = [g for g, _ in sp.factor_list(numer)[1] if g.has(u)]
    curve = polynomial(signed)
    candidates = candidate_factors(sign)
    # Removed scalar factors are invertible on the open used for elimination.
    assert all(
        str(g) in candidates
        for factor in (denom, sp.cancel(numer / signed))
        for g, _ in polynomial(factor).factor()[1]
    )
    blocks = symbolic_blocks(sign, z, l)
    for i in range(4):
        for h, target in enumerate((17, 16, 16, 16)):
            matrix = sp.Matrix(projected_entries(blocks, i, h))
            special = matrix.subs({z: 3, l: 2})
            selected = minor_indices(special, reduced_snapshot(special, sign), 1)
            sub = matrix.extract(*selected)
            assert sub.shape == (target, target)
            cache = {v: sp.factor(v.subs({e: N * u, f: D * u})) for v in set(sub)}
            det = determinant(
                cleared_rows(sub.applyfunc(cache.get), lambda v: ordinary_unit(v, z, l))
            )
            assert det
            for factor, _ in det.factor()[1]:
                norm = curve.resultant(factor, 1) if factor.degrees()[1] else factor
                assert norm and all(str(g) in candidates for g, _ in norm.factor()[1])
        print(
            sign, i, "all four projected rank certificates over Z[a,u,z,l]", flush=True
        )


if __name__ == "__main__":
    for sign in (1, -1):
        check(sign)
