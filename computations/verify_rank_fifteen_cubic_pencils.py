"""All geometric horizontal coordinates for the eight residual cubic levels."""

from functools import cache
from math import lcm

import sympy as sp
from flint import fq_default_ctx, fq_default_poly_ctx, nmod_poly

from verify_klein_center_lift import add_multiple
from verify_rank_fifteen_bivariate_sixth_cut import quotient_profile
from verify_rank_fifteen_bivariate_word_paths import joint_square, tensor
from verify_rank_fifteen_fixed_center_obstruction import universal_sections
from verify_rank_fifteen_quadratic_verticals import (
    coefficient_field,
    split_unit_coordinates,
)
from verify_rank_fifteen_unit_residues import verify_unit_point

HILBERT = (4, 10, 20, 31, 40, 44, 40, 31, 20, 10, 4, 1, 0)
CENTRAL_COLUMNS = (0, 1, 2, 4)


def level(field, constant, imaginary):
    b = field(constant)
    return b + (-field.one()).sqrt() if imaginary else b


def splitting_field(constant, imaginary):
    degrees = (2, 4, 6, 12) if imaginary else ((3,) if constant == 1 else (2,))
    for degree in degrees:
        field = fq_default_ctx(3, degree, "a")
        ring = fq_default_poly_ctx(field)
        t = ring.gen()
        b = level(field, constant, imaginary)
        if not imaginary:
            v = (t**3 - t - 1).roots()[0][0] if constant == 1 else (-field.one()).sqrt()
        else:
            roots = ((t**3 - t) ** 2 - ring([b])).roots()
            if not roots:
                continue
            assert len(roots) == 6
            v = roots[0][0]
        assert (v**3 - v) ** 2 == b
        return field, ring, b, v
    raise AssertionError((constant, imaginary))


def sparse_echelon(rows):
    """Order pivots by support, rational degree, column and then original row."""
    pending = [row for row in rows if row]
    pivots = []
    while pending:
        width = min(map(len, pending))
        _, column, index = min(
            (z.coordinates[0].degree() + z.denominator.degree(), j, i)
            for i, row in enumerate(pending)
            if len(row) == width
            for j, z in row.items()
        )
        row = pending.pop(index)
        inverse = row.pop(column).inverse()
        row = {j: z * inverse for j, z in row.items()}
        pivots.append((column, row, inverse.constant(1)))
        for other in pending:
            if column in other:
                add_multiple(other, row, -other.pop(column))
        pending = [other for other in pending if other]
    reduced = {}
    for column, row, one in reversed(pivots):
        for j in list(row):
            if j in reduced:
                add_multiple(row, reduced[j], -row[j])
        row[column] = one
        reduced[column] = row
    return reduced


def chart_divisor(centers, field, ring, v):
    curve = coefficient_field(ring)
    k = 1 / v

    def polynomial(poly):
        terms = {}
        for (i, j), a in poly.to_dict().items():
            j = int(j)
            terms[j] = terms.get(j, field.zero()) + int(a) * k ** int(i) * (k - v) ** j
        return ring(
            [terms.get(i, field.zero()) for i in range(max(terms, default=-1) + 1)]
        )

    def specialize(z):
        return curve.constant(polynomial(z.n)) / curve.constant(polynomial(z.d))

    relations = {w: (e, specialize(z)) for w, (e, z) in tensor.items()}
    joint = specialize(joint_square)
    joint.inverse()  # All normalizations of word parity require this unit.

    @cache
    def power(exponent):
        assert exponent % 2 == 0
        return joint ** (exponent // 2)

    cuts = [
        [{w: specialize(z) for w, z in form.items()} for form in side]
        for side in centers
    ]
    rows = quotient_profile(
        cuts, 0, 13, curve.constant(1), relations, power, sparse_echelon
    )
    assert tuple(sum(dims) for _, dims in rows) == HILBERT
    return curve.divisor.monic()


def exceptional_bound(divisors, ring, degree, b, relative_degree):
    core = ring.zero()
    for chart in divisors:
        for _ in range(degree // relative_degree):
            core = core.gcd(chart)
            chart = ring([a ** (3**relative_degree) for a in chart])
    c = ring.gen()
    numerator = ring([b]) * c**6
    denominator = (1 - c * c) ** 2
    for known in (
        c**18 - c * c,
        numerator**9 * denominator - numerator * denominator**9,
    ):
        while (common := core.gcd(known)).degree() > 0:
            core = core // common
    conjugate = ring([a**3 for a in core])
    if relative_degree == 2:
        core *= conjugate // core.gcd(conjugate)
    assert all(not any(map(int, a.to_list()[1:])) for a in core)
    return nmod_poly([int(a.to_list()[0]) for a in core], 3)


def certify_factor(f, constant, imaginary):
    q = 9 if imaginary else 3
    degree = lcm(2 if imaginary else 1, f.degree())
    for multiplier in (1, 2, 3, 4, 6, 12):
        field = fq_default_ctx(3, degree * multiplier, "a")
        ring = fq_default_poly_ctx(field)
        b = level(field, constant, imaginary)
        roots = [r for r, _ in ring(list(map(int, f))).roots()]
        assert len(roots) == f.degree()
        points = []
        while roots:
            c = roots[0]
            w = c
            while w in roots:
                roots.remove(w)
                w = w**q
            x = c * c
            if x**9 == x:
                continue
            y = x**3 * b
            ratio = y / (1 - x) ** 2
            if ratio**9 == ratio:
                continue
            coordinates = split_unit_coordinates(ring.gen(), x, y, c)
            if coordinates is None:
                break
            points.append((x, y, coordinates))
        else:
            for x, y, coordinates in points:
                verify_unit_point(field, x, y, coordinates, ((0, 1, 0), (0, 0, 1)))
            return
    raise AssertionError((f, constant, imaginary))


def verify_level(centers, constant, imaginary):
    field, ring, b, v = splitting_field(constant, imaginary)
    second = v + 1 if (constant, imaginary) == (-1, False) else -v
    divisors = [chart_divisor(centers, field, ring, w) for w in (v, second)]
    bound = exceptional_bound(divisors, ring, field.degree(), b, 2 if imaginary else 1)
    factors = bound.factor()[1]
    print(
        (constant, imaginary),
        [f.degree() for f in divisors],
        bound.degree(),
        [f.degree() for f, _ in factors],
        flush=True,
    )
    for f, _ in factors:
        certify_factor(f, constant, imaginary)
    print("PASS: all geometric fibers at", (constant, imaginary), flush=True)


def verify_quotient_curve():
    x, multiplier = sp.symbols("x multiplier")
    cubic = (1 - x) ** 2 - 4 * multiplier * x**3
    assert (
        sp.factor(sp.discriminant(cubic, x) - 16 * multiplier * (1 - 27 * multiplier))
        == 0
    )
    special = cubic.subs(multiplier, sp.Rational(1, 27))
    assert sp.expand(27 * special + (x - 3) ** 2 * (4 * x - 3)) == 0


if __name__ == "__main__":
    from concurrent.futures import ProcessPoolExecutor
    from multiprocessing import get_context

    verify_quotient_curve()
    centers = universal_sections(CENTRAL_COLUMNS)

    def run(pair):
        verify_level(centers, *pair)

    levels = [(1, False), (-1, False), (-1, True), (0, True), (1, True)]
    with ProcessPoolExecutor(max_workers=3, mp_context=get_context("fork")) as pool:
        list(pool.map(run, levels))
