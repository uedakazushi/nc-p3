"""Quintic rank certificates on the entire negative line and on a common slope open."""

from collections import Counter

import sympy as sp

from verify_fixed_exchange_isogeny import exchange
from verify_fixed_geometric_quartic_locus import (
    a,
    e,
    f,
    symbolic_blocks,
    minor_indices,
    reduced_snapshot,
)
from verify_fixed_line_elliptic_curves import ordinary_unit
from verify_fixed_quintic_locus import middle_relations
from verify_klein_elliptic_completion import coeffs
from verify_uniform_fixed_quartic_support import (
    u,
    z,
    l,
    polynomial,
    cleared_rows,
    candidate_factors,
)


def determinant_factors(rows):
    """Return nonconstant irreducible factors, with multiplicities, up to a rational unit."""
    powers = Counter()
    polynomials = {}

    def record(value, multiplicity):
        if not value.is_constant():
            label = str(value)
            powers[label] += multiplicity
            polynomials[label] = value

    def primitive(row):
        content = None
        for value in row:
            if value:
                content = value if content is None else content.gcd(value)
                if content.is_one():
                    return row
        if content is None:
            return row
        record(content, 1)
        return [value / content if value else value for value in row]

    rows = [primitive(row) for row in rows]
    while rows:
        width = len(rows)
        column_sizes = [sum(bool(row[j]) for row in rows) for j in range(width)]
        choices = []
        for i, row in enumerate(rows):
            count = sum(bool(v) for v in row)
            for j, v in enumerate(row):
                if v:
                    choices.append(((count - 1) * (column_sizes[j] - 1), len(v), i, j))
        assert choices, "the selected minor must be nonzero"
        _, _, i, j = min(choices)
        rows[0], rows[i] = rows[i], rows[0]
        for row in rows:
            row[0], row[j] = row[j], row[0]
        top = rows[0]
        pivot = top[0]
        record(pivot, 1)
        reduced = []
        for row in rows[1:]:
            if not row[0]:
                reduced.append(row[1:])
                continue
            common = row[0].gcd(pivot)
            alpha, beta = pivot / common, row[0] / common
            # Scaling the row by alpha accounts for its inverse in the original determinant.
            record(alpha, -1)
            reduced.append(
                primitive([alpha * v - beta * w for v, w in zip(row[1:], top[1:])])
            )
        rows = reduced
    factors = Counter()
    irreducibles = {}
    for label, n in powers.items():
        if n:
            for g, k in polynomials[label].factor()[1]:
                key = str(g)
                irreducibles[key] = g
                factors[key] += k * n
    assert all(n >= 0 for n in factors.values())
    return [(irreducibles[key], n) for key, n in factors.items() if n]


def specialized(poly):
    coeff = {}
    for (i, j, k, n), c in poly.to_dict().items():
        assert j == 0
        coeff[i] = coeff.get(i, 0) + int(c) * 3**k * 2**n
    return sp.Poly.from_dict({(i,): v for i, v in coeff.items()}, (a,), domain=sp.QQ)


def forbidden_slice(sign, N, D, curve):
    subs = {z: 3, l: 2}
    ns, ds = N.subs(subs), D.subs(subs)
    expr = sum(
        int(c) * a**i * u**j * 3**k * 2**n
        for (i, j, k, n), c in curve.to_dict().items()
    )
    bad = a * (a * a - 3) * ns * ds * (a - 6) * (2 * a - 1)
    if sign == 1:
        bad *= (a * a - 78 * a + 69) * (23 * a * a - 78 * a + 3)
    else:
        bad *= (a + 6) * (2 * a + 1)
    invalid = (
        (a * a - 3 * (ds * u) ** 2) * ((ns * u) ** 2 - 3) * (ns * ns - 3 * ds * ds) * u
    )
    bad *= sp.Poly(expr, u).LC() * sp.resultant(expr, invalid, u)
    return sp.Poly(bad, a, domain=sp.QQ)


def check(sign):
    N, D = exchange(sign, z, l) * sp.Matrix([a, 1])
    co, _, _ = coeffs(a, sp.S.One, sign * z, sign * a, z, l)
    equation = sp.factor(
        sign * z * co[1] * D**2 * u**2 + D * u - co[0] - sign * co[1] * a * a
    )
    (curve,) = [
        polynomial(g)
        for g, _ in sp.factor_list(equation.as_numer_denom()[0])[1]
        if g.has(u)
    ]
    candidates = candidate_factors(sign)
    forbidden = forbidden_slice(sign, N, D, curve)
    blocks = symbolic_blocks(sign, z, l)
    for i in (0, 1):
        for h in (0, 2):
            matrix = middle_relations(blocks, i, h)
            special = matrix.subs({z: 3, l: 2})
            I, J = minor_indices(special, reduced_snapshot(special, sign), 1)
            assert len(I) == len(J) == 86
            sub = matrix.extract(I, J)
            cache = {v: sp.factor(v.subs({e: N * u, f: D * u})) for v in set(sub)}
            factors = determinant_factors(
                cleared_rows(sub.applyfunc(cache.get), lambda d: ordinary_unit(d, z, l))
            )
            extra = {}
            for g, _ in factors:
                norm = curve.resultant(g, 1) if g.degrees()[1] else g
                assert norm
                for fac, _ in norm.factor()[1]:
                    if str(fac) in candidates:
                        continue
                    assert sign == 1 and fac.degrees()[0] > 0
                    p = specialized(fac)
                    assert p.degree() == p.sqf_part().degree() == fac.degrees()[0]
                    assert p.gcd(forbidden).degree() == 0
                    extra[str(fac)] = p.degree()
            print(
                sign,
                i,
                h,
                "additional finite factors by degree",
                sorted(extra.values()),
                flush=True,
            )


if __name__ == "__main__":
    for sign in (-1, 1):
        check(sign)
