"""Check length eight at the first start on the cover with fixed signed root c."""

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from math import comb

import sympy as sp
from verify_rank_fifteen_bivariate_seventh_cut import (
    fiber,
    finite_certificate,
    horizontal_bound,
    primitive_echelon,
)
from verify_rank_fifteen_bivariate_sixth_cut import (
    X as C,
)
from verify_rank_fifteen_bivariate_sixth_cut import (
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


def lift_identities():
    assert comb(15, 3) - 4 * 31 - 10 * 31 - 20 == 1
    c, v, n, h, u = sp.symbols("c v n h u")
    k = 1 / v
    residue_h = c * (k - v)
    assert sp.cancel(residue_h**3 * k**3 / (k * k - 1) ** 3 - c**3) == 0
    sign_sum = n * n * ((2 * n * n + 1 - u) + (2 * (h / n) ** 2 + 1 - u))
    assert sp.expand(sign_sum - 2 * (n**4 + (1 - u) * n * n + h * h)) == 0


def first_start():
    basis, image = paths()
    cuts = sections(basis, centers0(basis, image), (0, 1, 3, 4))
    image.cache_clear()
    reducer = partial(primitive_echelon, rational_order=True)
    for n, dims in quotient_profile(cuts, 0, 8, field.one, tensor, adjust, reducer):
        assert sum(dims) == (0, 4, 10, 20, 31, 40, 44, 40, 31)[n]
        print("Length", n, "first-start dimensions", dims, flush=True)
    assert dims == [7, 8, 8, 8]
    return set(factors)


def signed_fibers(divisors):
    p = exceptional_product(divisors)
    assert tuple(p.degrees()) == (167, 146)
    coefficients = [base.constant(0) for i in range(6)]
    for (j, i), a in p.to_dict().items():
        for k, r in enumerate(reduced_power(j)):
            coefficients[k] += int(a) * C**i * r
    common = coefficients[0]
    for r in coefficients[1:]:
        common = common.gcd(r)
    assert common == 1
    bound = horizontal_bound(coefficients, signed_root=True)
    assert bound.degree() == 33
    pairs = set()
    for modulus, _ in bound.factor()[1]:
        if modulus.degree() > 2:
            m = tuple(map(int, modulus))
            pairs.update((m, b) for b in fiber(coefficients, m))
    assert len(pairs) == 8
    return sorted(pairs)


if __name__ == "__main__":
    lift_identities()
    pairs = signed_fibers(first_start())
    with ProcessPoolExecutor(max_workers=3) as pool:
        list(pool.map(partial(finite_certificate, signed_root=True), pairs))
