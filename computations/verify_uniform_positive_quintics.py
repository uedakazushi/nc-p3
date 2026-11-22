"""Verify four positive quintic rank minors inside the localized quadratic curve algebra."""

import json
from pathlib import Path

import sympy as sp
from flint import fmpz_mpoly_ctx

from verify_fixed_exchange_isogeny import exchange
from verify_fixed_geometric_quartic_locus import a, e, f, symbolic_blocks
from verify_fixed_line_elliptic_curves import ordinary_unit
from verify_fixed_quintic_locus import middle_relations
from verify_klein_elliptic_completion import coeffs
from verify_uniform_fixed_quartic_support import (
    candidate_factors,
    cleared_rows,
    polynomial,
    u,
    z,
    l,
)

base = fmpz_mpoly_ctx.get(("a", "z", "l"))


def coefficient(poly, j):
    return base.from_dict(
        {(i, k, n): v for (i, power, k, n), v in poly.to_dict().items() if power == j}
    )


def verify(certificate, blocks, N, D, q0, q1, q2, allowed):
    def is_unit(value):
        if not value:
            return False
        for fac in allowed:
            while True:
                quotient, remainder = divmod(value, fac)
                if remainder:
                    break
                value = quotient
        return value.is_constant()

    c = q0 * q2

    def mul(x, y):
        v = x[1] * y[1]
        return x[0] * y[0] - c * v, x[0] * y[1] + x[1] * y[0] - q1 * v

    def strip(row):
        content = base.constant(0)
        for x, y in row:
            content = content.gcd(x).gcd(y)
            if content.is_one():
                return row
        assert is_unit(content)
        return [(x / content, y / content) for x, y in row]

    def lift(poly):
        assert poly.degrees()[1] <= 2
        x, y, t = [coefficient(poly, j) for j in range(3)]
        return q2 * q2 * x - c * t, q2 * y - q1 * t

    i, h = certificate["start"], certificate["character"]
    I, J = certificate["rows"], certificate["columns"]
    assert len(I) == len(set(I)) == len(J) == len(set(J)) == 86
    assert set(I) <= set(range(120)) and set(J) <= set(range(100))
    sub = middle_relations(blocks, i, h).extract(I, J)
    cache = {v: sp.factor(v.subs({e: N * u, f: D * u})) for v in set(sub)}
    rows = cleared_rows(
        sub.applyfunc(cache.get), lambda value: ordinary_unit(value, z, l)
    )
    assert is_unit(q2)
    rows = [strip([lift(value) for value in row]) for row in rows]
    while rows:
        top = rows[0]
        x, y = top[0]
        conjugate = x - q1 * y, -y
        norm = x * x - q1 * x * y + c * y * y
        assert is_unit(norm)
        reduced = []
        for row in rows[1:]:
            if not any(row[0]):
                reduced.append(row[1:])
                continue
            product = mul(row[0], conjugate)
            common = norm.gcd(product[0]).gcd(product[1])
            alpha = norm / common
            beta = product[0] / common, product[1] / common
            new = []
            for value, above in zip(row[1:], top[1:]):
                p, q = mul(beta, above)
                new.append((alpha * value[0] - p, alpha * value[1] - q))
            reduced.append(strip(new))
        rows = reduced
    print(i, h, "86 invertible pivots on the positive curve", flush=True)


if __name__ == "__main__":
    N, D = exchange(1, z, l) * sp.Matrix([a, 1])
    co, _, _ = coeffs(a, sp.S.One, z, a, z, l)
    eq = sp.factor(z * co[1] * D**2 * u**2 + D * u - co[0] - co[1] * a * a)
    (curve,) = [
        polynomial(g) for g, _ in sp.factor_list(eq.as_numer_denom()[0])[1] if g.has(u)
    ]
    q0, q1, q2 = [coefficient(curve, j) for j in range(3)]
    allowed = sorted(
        [
            coefficient(
                polynomial(
                    sp.sympify(s.replace("^", "**"), locals={"a": a, "z": z, "l": l})
                ),
                0,
            )
            for s in candidate_factors(1)
        ],
        key=len,
    )
    path = (
        Path(__file__).resolve().parents[1]
        / "research"
        / "fixed_quintic_pivot_certificate.json"
    )
    cases = json.loads(path.read_text())["minors"]
    assert [(v["start"], v["character"]) for v in cases] == [
        (0, 0),
        (0, 2),
        (1, 0),
        (1, 2),
    ]
    blocks = symbolic_blocks(1, z, l)
    for case in cases:
        verify(case, blocks, N, D, q0, q1, q2, allowed)
