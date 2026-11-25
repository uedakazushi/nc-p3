"""Certify the ninth character that suffices by translation of the four-cut."""

from collections import Counter
from functools import reduce
from itertools import product

from verify_rank_fifteen_bivariate_seventh_cut import (
    fiber,
    finite_certificate,
    horizontal_bound,
    primitive_echelon,
)
from verify_rank_fifteen_bivariate_sixth_cut import (
    X,
    base,
    exceptional_product,
    quotient_profile,
    reduced_power,
    sections,
)
from verify_rank_fifteen_bivariate_word_paths import (
    adjust,
    centers0,
    factors,
    field,
    paths,
    tensor,
)


def walsh_symmetries():
    permutation = (0, 3, 2, 1)
    for i, j, k, t in product(range(4), repeat=4):
        h = i ^ k
        if h != j ^ t:
            continue
        for s, u in product(range(2), range(h & 1, 4, 2)):

            def phase(a, b, d, e):
                return (
                    ((a ^ b) & 1) * s
                    + (u >> 1) * (d >> 1)
                    + (permutation[u] >> 1) * (e >> 1)
                ) % 2

            initial = phase(i, j, k, t)
            assert phase(i ^ 1, j ^ 1, k ^ 1, t ^ 1) == initial
            assert phase(i ^ 2, j ^ 2, k ^ 2, t ^ 2) == (initial ^ (h & 1))


def ninth_divisors():
    basis, image = paths()
    cuts = sections(basis, centers0(basis, image), (0, 1, 3, 4))
    image.cache_clear()
    block = 0

    def reducer(rows):
        nonlocal block
        block += 1
        # At the last length only character three is queried. No later
        # calculation reads the three skipped quotient spaces.
        if 33 <= block <= 35:
            return {}
        return primitive_echelon(rows, rational_order=True)

    for n, dims in quotient_profile(cuts, 0, 9, field.one, tensor, adjust, reducer):
        if n < 9:
            assert sum(dims) == (0, 4, 10, 20, 31, 40, 44, 40, 31)[n]
        else:
            assert dims[3] == 5
    return set(factors)


def signed_exceptions(divisors):
    p = exceptional_product(divisors)
    assert tuple(p.degrees()) == (181, 160)
    coefficients = [base.constant(0) for _ in range(6)]
    for (j, i), a in p.to_dict().items():
        for k, r in enumerate(reduced_power(int(j))):
            coefficients[k] += int(a) * X ** int(i) * r
    assert reduce(lambda a, b: a.gcd(b), coefficients) == 1
    bound = horizontal_bound(coefficients, signed_root=True)
    assert bound.degree() == 67
    pairs = set()
    for f, _ in bound.factor()[1]:
        if f.degree() > 2:
            m = tuple(map(int, f))
            pairs.update((m, b) for b in fiber(coefficients, m))
    assert Counter(len(m) - 1 for m, b in pairs) == {3: 8, 4: 3}
    return sorted(pairs)


if __name__ == "__main__":
    walsh_symmetries()
    for pair in signed_exceptions(ninth_divisors()):
        finite_certificate(pair, signed_root=True)
    print("PASS: signed ninth-character bound and all eleven residual pairs")
