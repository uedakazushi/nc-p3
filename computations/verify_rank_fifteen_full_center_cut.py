"""Generic quotient profiles and pivot loci on the characteristic-three cover."""

from functools import partial, reduce

from flint import nmod_poly

from verify_rank_fifteen_bivariate_seventh_cut import (
    fiber,
    finite_certificate,
    horizontal_bound,
    primitive_echelon,
)
from verify_rank_fifteen_bivariate_sixth_cut import (
    exceptional_product,
    orbit_coefficients,
    quotient_profile,
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


def first_start(cuts, expected):
    profile = quotient_profile(
        cuts,
        0,
        len(expected),
        field.one,
        tensor,
        adjust,
        partial(primitive_echelon, rational_order=True),
    )
    for (n, dims), total in zip(profile, expected, strict=True):
        assert sum(dims) == total
        if n == 9:
            assert dims == [6, 6, 6, 5]
        print("Length", n, "dimensions", dims, flush=True)
    return set(factors)


def descend(divisors):
    polynomial = exceptional_product(divisors)
    assert tuple(polynomial.degrees()) == (177, 140)
    coefficients = orbit_coefficients(polynomial)
    assert len(coefficients) == 12
    assert reduce(lambda x, y: x.gcd(y), coefficients) == 1
    # Both six-row blocks use X=c^2; neither treats c as a signed invariant.
    bound = horizontal_bound(coefficients[:6]).gcd(horizontal_bound(coefficients[6:]))
    t = nmod_poly([0, 1], 3)
    assert bound == (t + 1) ** 10 * (t * t + 1) ** 2 * (t**3 + 2 * t + 2) * (
        t**5 + 2 * t**4 + 2 * t + 2
    )
    residual = set()
    for f, _ in bound.factor()[1]:
        if f.degree() <= 2:  # Previously certified quadratic verticals.
            continue
        modulus = tuple(map(int, f))
        residual.update((modulus, b) for b in fiber(coefficients, modulus))
    assert residual == {
        ((2, 2, 0, 1), (1, 2, 1)),
        ((2, 2, 0, 0, 2, 1), (0, 0, 1, 2, 0)),
    }
    for pair in sorted(residual):
        finite_certificate(pair)


def verify():
    basis, image = paths()
    pivots = centers0(basis, image)
    saved = set(factors)
    image.cache_clear()
    cuts = sections(basis, pivots, range(5))
    divisors = first_start(cuts, [4, 10, 20, 30, 36, 34, 20, 0])
    descend(divisors)
    factors.clear()
    factors.update(saved)
    fixed = sections(basis, pivots, (0, 1, 3, 4))
    first_start(fixed, [4, 10, 20, 31, 40, 44, 40, 31, 23])


if __name__ == "__main__":
    verify()
