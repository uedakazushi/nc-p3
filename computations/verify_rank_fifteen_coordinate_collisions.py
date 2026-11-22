"""Exact necessary conditions for exceptional indexed equivalences among coordinate ends."""

from functools import cache
from itertools import combinations, product

import sympy as sp

from verify_rank_fifteen_coordinate_ends import (
    EllipticChart,
    GaussianChart,
    ONE,
    ZERO,
    VARIABLE,
    elliptic_frame,
    gaussian_frame,
)
from verify_rank_fifteen_endpoint_orbits import multiply, trace, trace_blocks

t, s, u, U, v = sp.symbols("t s u U v")
I = sp.I


def oriented_cubic(S):
    def triangle(h, g, e, x, y, z):
        return sum(S[h][a][x] * S[g][a ^ h][y] * S[e][a ^ h ^ g][z] for a in range(4))

    total = 0
    for i, j, k in product(range(4), repeat=3):
        h, g, e = i ^ j, j ^ k, k ^ i
        total += triangle(h, g, e, j, k, i) + triangle(g, h, e, k, j, i)
    return 2 * total


def fifth_contraction(T, field):
    """Contract the two cycles of Q after clearing a common denominator."""
    denominator = ONE
    for block in T:
        for row in block:
            for x in row:
                denominator = (
                    denominator // denominator.gcd(x.denominator) * x.denominator
                )

    def numerator(x):
        scale = denominator // x.denominator
        return x.coordinates[0] * scale, x.coordinates[2] * scale

    entries = [[[numerator(x) for x in row] for row in block] for block in T]
    P = (4, 0, 3, 1, 2)
    Q = (2, 4, 3, 0, 1)
    real = radical = ZERO
    for a in product(range(4), repeat=5):
        shifts = [a[j] ^ a[P[j]] for j in range(5)]
        if shifts[1] != shifts[4]:
            continue
        for b0, b1 in product(range(4), repeat=2):
            b = (b0, b1, b0 ^ shifts[0], b0 ^ shifts[0] ^ shifts[2], b1 ^ shifts[1])
            x, y = ONE, ZERO
            for j in range(5):
                c, d = entries[shifts[j]][a[j]][b[Q[j]]]
                x, y = x * c + field.radicand * y * d, x * d + y * c
            real += x
            radical += y
    return field((real, ZERO, radical, ZERO), denominator**5)


def quotient_rational(expression, sign):
    field = sp.QQ_I.frac_field(s)

    def poly(a):
        return sp.Poly(a, t, domain=field)

    equation = poly(t * t - s * t + sign)
    a, b = sp.fraction(sp.cancel(expression, extension=I))
    result = ((poly(a) % equation) * (poly(b) % equation).invert(equation)) % equation
    assert result.degree() <= 0
    return result.as_expr()


@cache
def invariants(branch):
    if branch == -1:
        return [f.xreplace({I: -I}) for f in invariants(1)]
    field, frame = (
        (EllipticChart, elliptic_frame)
        if branch == 0
        else (GaussianChart, gaussian_frame)
    )
    actual = frame(field.lift(VARIABLE), field((ZERO, ZERO, ONE, ZERO)))
    T, S = trace_blocks(actual)
    den = sum(map(trace, S))
    square = [multiply(A, A) for A in S]
    cube = [multiply(A, B) for A, B in zip(square, S)]
    values = (
        sum(map(trace, square)) / den**2,
        sum(map(trace, cube)) / den**3,
        oriented_cubic(S) / den**3,
        fifth_contraction(T, field) ** 2 / den**5,
        den,
    )
    out = []
    for entry in values:
        a, b, d = [
            sum(int(c) * t**i for i, c in enumerate(g))
            for g in (entry.coordinates[0], entry.coordinates[2], entry.denominator)
        ]
        if branch:
            out.append(quotient_rational((a + I * b) / d, -1))
        else:
            re = quotient_rational(a / d, 1)
            if b:
                re += v * quotient_rational(b * t * t / ((t * t - 1) * d), 1)
            out.append(re)
    return out


def weight_equation(branch):
    if branch:
        return branch * I * (s * s + 4) - 8 * s * U
    return (s + 10) ** 2 - 16 * (s + 2) * U


def projected_bound(left, right):
    a, b = invariants(left), invariants(right)
    indices = (0, 1, 2, 3) if left and right else (0, 1, 3)
    equations = []
    for index in indices:
        an, ad = sp.fraction(sp.cancel(a[index], extension=I))
        bn, bd = sp.fraction(sp.cancel(b[index], extension=I))
        row = sp.Poly(an * bd.subs(s, u) - ad * bn.subs(s, u), u, s, domain=sp.QQ_I)
        if left == right:
            row = row.exquo(sp.Poly(s - u, u, s, domain=sp.QQ_I))
        equations.append(row)
    support = None
    for A, B in combinations(equations, 2):
        resultant = A.resultant(B)
        if resultant:
            support = resultant if support is None else support.gcd(resultant)
    support = support.sqf_part()
    weight = weight_equation(left)
    boundary = (
        s * (s * s - 4) * (s * s + 4)
        if left
        else (s * s - 4) * (s - 14) * (s - 6) * (s + 10)
    )
    excluded = sp.Poly(
        boundary * weight.subs(U, 1) * weight.subs(U, -3), s, domain=sp.QQ_I
    )
    support = support.exquo(support.gcd(excluded))
    return (
        sp.Poly(sp.resultant(support.as_expr(), weight, s), U, domain=sp.QQ_I)
        .sqf_part()
        .monic()
    )


def verify():
    for left, right, expected in (
        (0, 0, (U - 9) * (U * U - 74 * U + 169)),
        (1, 1, (U + 1) * (U * U - 42 * U + 9)),
        (1, -1, U - 9),
        (0, 1, U**3 + 13 * U * U + 259 * U - 529),
        (1, 0, U**3 - 147 * U * U - 61 * U - 49),
    ):
        assert projected_bound(left, right) == sp.Poly(expected, U, domain=sp.QQ_I)
    E = invariants(0)
    amplitude = sp.cancel(E[2].diff(v))
    assert (
        sp.expand(sp.fraction(amplitude)[0] + 48 * (s + 2) * (s * s - 28 * s + 4)) == 0
    )
    assert sp.expand(weight_equation(0).subs(U, 3) - (s * s - 28 * s + 4)) == 0

    for branch, expected in ((0, sp.Rational(64, 375)), (1, -sp.Rational(64, 81))):
        J, K, _, _, den = invariants(branch)
        pole = sp.Poly(weight_equation(branch).subs(U, -3), s, domain=sp.QQ_I).monic()
        dn, dd = sp.fraction(sp.cancel(den, extension=I))
        assert sp.Poly(dn, s, domain=sp.QQ_I).monic() == pole
        assert pole.gcd(sp.Poly(dd, s, domain=sp.QQ_I)).degree() == 0
        # tr(S^2) stays nonzero at these poles, so the ratio is indexed-invariant there.
        a, b = sp.fraction(sp.cancel(J * den**2, extension=I))
        assert all(
            pole.gcd(sp.Poly(c, s, domain=sp.QQ_I)).degree() == 0 for c in (a, b)
        )
        a, b = sp.fraction(sp.cancel(K * K / J**3, extension=I))
        assert not (sp.Poly(a - expected * b, s, domain=sp.QQ_I) % pole)
        assert pole.gcd(sp.Poly(b, s, domain=sp.QQ_I)).degree() == 0
    print(
        "Five collision projections verified, including sign separation and the trace-zero comparison."
    )


if __name__ == "__main__":
    verify()
