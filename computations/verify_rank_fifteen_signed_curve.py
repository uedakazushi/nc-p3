"""The entire signed residual curve, including every common exceptional horizontal fiber."""

from functools import cache
from math import lcm

from flint import fq_default_ctx, fq_default_poly_ctx, nmod_poly
from rank_fifteen_signed_curve_field import C, ONE, ZERO, ZERO_TRIPLE, Curve
from verify_klein_center_lift import add_multiple, echelon
from verify_rank_fifteen_bivariate_sixth_cut import quotient_profile
from verify_rank_fifteen_bivariate_word_paths import joint_square, tensor
from verify_rank_fifteen_fixed_center_obstruction import (
    specialize_sections,
    universal_sections,
)

PROFILE = (4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0)


def curve_equation(c, v):
    u = c * v
    return u**6 - (c * c + 1) * u**4 + (c**4 - c**3 + c) * u**2 + c**3 - c * c - c + 1


def sparse_echelon(rows):
    pending = [r for r in rows if r]
    pivots = []
    while pending:
        width = min(map(len, pending))
        _, column, index = min(
            (z.pivot_degree(), j, i)
            for i, r in enumerate(pending)
            if len(r) == width
            for j, z in r.items()
        )
        row = pending.pop(index)
        inverse = row.pop(column).inverse()
        row = {j: z * inverse for j, z in row.items()}
        pivots.append((column, row))
        for other in pending:
            if column in other:
                add_multiple(other, row, -other.pop(column))
        pending = [r for r in pending if r]
    result = {}
    for column, row in reversed(pivots):
        for j in list(row):
            if j in result:
                add_multiple(row, result[j], -row[j])
        row[column] = Curve.constant(1)
        result[column] = row
    return result


def divisor(centers, omitted, shift):
    Curve.divisor = ONE
    c = Curve.constant(C)
    u = Curve(ZERO_TRIPLE, (ONE, ZERO, ZERO))
    v = u / c
    assert not curve_equation(c, v) and not (u * (1 / u) - 1)
    v += shift
    k = 1 / v
    h = c * (k - v)

    @cache
    def monomial(i, j):
        return k**i * h**j

    def evaluate(poly):
        return sum(
            (int(a) * monomial(int(i), int(j)) for (i, j), a in poly.to_dict().items()),
            Curve.constant(0),
        )

    def specialize(value):
        return evaluate(value.n) / evaluate(value.d)

    entries = {w: (e, specialize(z)) for w, (e, z) in tensor.items()}
    joint = specialize(joint_square)

    @cache
    def power(e):
        assert e % 2 == 0
        return joint ** (e // 2)

    cuts = [
        [
            {w: specialize(z) for w, z in row.items()}
            for j, row in enumerate(side)
            if j != omitted
        ]
        for side in centers
    ]
    rows = quotient_profile(
        cuts, 0, 13, Curve.constant(1), entries, power, sparse_echelon
    )
    assert tuple(sum(row) for _, row in rows) == PROFILE
    return Curve.divisor


def good_presentation(field, c, v, centers):
    for shift in (0, 1, -1):
        for sign in (1, -1):
            sections, entries, power = specialize_sections(
                field, c, sign * v + shift, centers
            )
            for omitted in (0, 1, 3):
                cuts = [
                    [row for j, row in enumerate(side) if j != omitted]
                    for side in sections
                ]
                rows = quotient_profile(cuts, 0, 13, power(0), entries, power, echelon)
                if all(sum(row) == PROFILE[n - 1] for n, row in rows):
                    return True
    return False


def exceptional_fiber(modulus, centers):
    d = modulus.degree()
    base = fq_default_ctx(3, d, "a")
    ring = fq_default_poly_ctx(base)
    coeff = list(map(int, modulus))
    c = ring(coeff).roots()[0][0]
    polynomial = curve_equation(c, ring.gen())
    relative = lcm(*(f.degree() for f, _ in polynomial.factor()[1]))
    field = fq_default_ctx(3, d * relative, "a")
    ring = fq_default_poly_ctx(field)
    c = ring(coeff).roots()[0][0]
    x = c * c
    roots = curve_equation(c, ring.gen()).roots()
    assert sum(e for _, e in roots) == 6
    remaining = {(v**3 - v) ** 2 for v, _ in roots}
    new_orbits = 0
    while remaining:
        b = next(iter(remaining))
        w = b
        while w in remaining:
            remaining.remove(w)
            w = w ** (3**d)
        if x**9 == x or b**9 == b:
            continue
        ratio = b * x**3 / (1 - x) ** 2
        if ratio**9 == ratio:
            continue
        v = next(v for v, _ in roots if (v**3 - v) ** 2 == b)
        assert good_presentation(field, c, v, centers)
        new_orbits += 1
    return new_orbits


def verify(centers):
    first = divisor(centers, 0, 0)
    second = divisor(centers, 3, 1)
    assert (first.degree(), second.degree()) == (6986, 17641)
    p = first.gcd(second)
    t = nmod_poly([0, 1], 3)
    while (g := p.gcd(t**18 - t * t)).degree() > 0:
        p = p // g
    factors = p.factor()[1]
    assert p.degree() == 237 and len(factors) == 41
    assert sum(exceptional_fiber(m, centers) for m, _ in factors) == 100
    print("PASS: signed curve, both generic presentations and 100 exceptional orbits")


if __name__ == "__main__":
    verify(universal_sections(range(5)))
