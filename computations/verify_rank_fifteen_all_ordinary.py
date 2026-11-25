"""Maximal-minor descent for every ordinary positive fixed rank-fifteen window."""

from functools import reduce

from central_border_paths import central_border
from verify_rank_fifteen_bivariate_border import (
    EXPECTED,
    independent_sources,
    reduce_words,
)
from verify_rank_fifteen_bivariate_seventh_cut import (
    fiber,
    finite_certificate,
    horizontal_bound,
)
from verify_rank_fifteen_bivariate_sixth_cut import (
    base,
    c,
    exceptional_product,
    reduced_power,
    v,
)
from verify_rank_fifteen_bivariate_word_paths import (
    adjust,
    factors,
    field,
    open_factors,
    tensor,
)
from verify_rank_fifteen_center_border import border_data
from verify_rank_fifteen_fixed_center_obstruction import universal_sections
from verify_rank_fifteen_signed_curve import curve_equation


def minor_divisors(centers):
    relations, weighted = border_data(centers, tensor, adjust, 0)
    parts = {}
    block = 0

    def choose(rows):
        nonlocal block
        degree, h = divmod(block, 4)
        degree += 1
        block += 1
        if degree < 9:
            assert not rows
            return []
        determinant = field.one

        def capture(pivot):
            nonlocal determinant
            determinant *= pivot

        selected = independent_sources(rows, capture=capture)
        if degree in (10, 12) or (degree in (11, 13) and h == 0):
            part = parts.setdefault(degree, set())
            part.update(str(f) for f, _ in determinant.n.factor()[1])
            part.update(factors)
        return selected

    direction = [field.one] + [field.zero] * 4
    profile = central_border(
        relations,
        weighted,
        direction,
        field.one,
        reduce_words=reduce_words,
        select=choose,
    )
    assert tuple(sum(dims) for _, dims in profile) == EXPECTED
    return parts


def six_coefficients(p):
    by_power = {}
    for (j, i), value in p.to_dict().items():
        by_power.setdefault(int(j), {})[int(i), 0] = int(value)
    out = [base.constant(0) for _ in range(6)]
    for j, coefficients in by_power.items():
        a = base.from_dict(coefficients)
        for k, r in enumerate(reduced_power(j)):
            out[k] += a * r
    return out


def signed_curve():
    curve = curve_equation(c, v)
    assert curve.factor()[1] == [(curve, 1)]
    product = curve * curve_equation(c, v + 1) * curve_equation(c, v - 1)
    polynomial, *rest = six_coefficients(product)
    assert not any(rest) and tuple(polynomial.degrees()) == (18, 3)
    return polynomial


def exceptional_pairs(divisors, shared_curve, expected):
    polynomial = exceptional_product(divisors)
    assert tuple(polynomial.degrees()) == expected[0]
    coefficients = six_coefficients(polynomial)
    common = reduce(lambda a, b: a.gcd(b), coefficients)
    assert common == shared_curve
    coefficients = [a // common for a in coefficients]
    bound = horizontal_bound(coefficients, signed_root=True)
    assert bound.degree() == expected[1]
    pairs = set()
    for f, _ in bound.factor()[1]:
        if f.degree() > 2:
            modulus = tuple(map(int, f))
            pairs.update((modulus, b) for b in fiber(coefficients, modulus))
    assert len(pairs) == expected[2]
    return pairs


def verify():
    factors.clear()
    centers = universal_sections(range(5))
    factors.update(open_factors)
    parts = minor_divisors(centers)
    common = signed_curve()
    table = {
        10: ((120, 114), 6, 2),
        11: ((84, 80), 0, 0),
        12: ((171, 170), 9, 3),
        13: ((88, 86), 3, 1),
    }
    points = set()
    for n, expected in table.items():
        points.update(exceptional_pairs(parts[n], common, expected))
        print("Certified descent at length", n, flush=True)
    assert points == {
        ((1, 1, 2, 1), (1, 1, 2)),
        ((2, 2, 2, 1), (1, 0, 0)),
        ((2, 2, 2, 1), (2, 2, 1)),
    }
    for point in sorted(points):
        finite_certificate(point, signed_root=True)
    print("PASS: off the signed curve, all unit pairs through thirteen")


if __name__ == "__main__":
    verify()
