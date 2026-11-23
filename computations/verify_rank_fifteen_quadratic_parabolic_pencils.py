"""Certify every quadratic residual ratio Y/(1-X)^2 at arbitrary X."""

import gc
from math import gcd, lcm
from flint import fq_default_ctx, fq_default_poly_ctx
from verify_rank_fifteen_genus_two_unit_line import central_divisor, unit_model
from verify_rank_fifteen_quadratic_verticals import (
    coefficient_field,
    descend,
    split_unit_coordinates,
)
from verify_rank_fifteen_unit_quartic_towers import encode
from verify_rank_fifteen_unit_residues import verify_unit_point


def verify_line(s):
    field = fq_default_ctx(3, 4, "a")
    ring = fq_default_poly_ctx(field)
    verify_pencil(s, field, ring)
    # Collect the dynamic coefficient class before python-flint clears its native field.
    gc.collect()


def verify_pencil(s, field, ring):
    functions = coefficient_field(ring)
    q = functions.constant(ring.gen())
    i = (-field.one()).sqrt()
    ratio = s + i
    r = ratio.sqrt()
    a = r**27
    assert a**3 == r and a != r
    b = (q - q**3) / functions.constant(r - a)
    v = q + functions.constant(a) * b
    d = b**3 - b
    assert v**3 - v == functions.constant(r) * d
    windows = unit_model(
        v, functions.constant(i), 1 / b, functions.constant((ratio - 1).sqrt()) * d
    )
    degrees, totals, exception_degrees = {
        0: ((1284, 766, 971, 1142), (3831, 66, 84), [4, 4, 6]),
        1: ((1476, 743, 1755, 1723), (5416, 60, 96), [4, 4, 8]),
        -1: ((1426, 1136, 1640, 1844), (5768, 60, 72), [4, 4, 4]),
    }[s]
    divisor = central_divisor(
        windows, functions.constant(1), functions, degrees, selected=(0, 1, 3, 4)
    )
    stable, prime_stable, levels = descend(ring, divisor)
    assert (divisor.degree(), stable, prime_stable) == totals
    assert sorted(len(g) - 1 for g in levels) == exception_degrees
    for g in levels:
        check(g, s)
    print(
        "Verified all geometric horizontal residues for ratio",
        s,
        "+ i and its conjugate.",
        flush=True,
    )


def check(g, s):
    d = len(g) - 1
    deg = lcm(4, d)
    for mult in (1, 2, 3, 4, 6, 12):
        F = fq_default_ctx(3, deg * mult, "a")
        ring = fq_default_poly_ctx(F)
        t = ring.gen()
        i = (-F.one()).sqrt()
        R = s + i
        r = R.sqrt()
        # Inverse Frobenius on the F81 constants.
        a = r**27
        kappa = -((r - a) ** 2)
        assert kappa**9 == kappa
        roots = sorted([v for v, e in ring(list(g)).roots()], key=encode)
        seen = set()
        points = []
        for B in roots:
            if encode(B) in seen:
                continue
            t0 = B
            while encode(t0) not in seen:
                seen.add(encode(t0))
                t0 = t0**9
            x = kappa / B
            if not x.is_square():
                continue
            y = R * (1 - x) ** 2
            coords = split_unit_coordinates(t, x, y, x.sqrt())
            if coords is not None:
                points.append((x, y, coords))
        if len(points) != gcd(d, 2):
            continue
        for x, y, coords in points:
            verify_unit_point(F, x, y, coords, ((0, 1, 0), (0, 0, 1)))
        gc.collect()
        print("Verified", s, g, deg * mult, len(points), flush=True)
        return
    raise AssertionError((g, s))


def check_descent_parity():
    from itertools import combinations, product
    from verify_rank_fifteen_genus_two_unit_line import normalize_support

    def beta(g, h):
        return (g & 1) * (h >> 1)

    def parity(p, word):
        local = sum((p + j) % 4 == 2 and a == 3 for j, a in enumerate(word))
        return (local + sum(beta(a, b) for a, b in combinations(word, 2))) % 2

    for word, (_, powers) in normalize_support().items():
        for p in range(4):
            shifted = word[p:] + word[:p]
            assert parity(p, shifted) == powers[3] % 2
    for g, h, k in product(range(4), repeat=3):
        difference = beta(h, k) + beta(g, h ^ k) - beta(g, h) - beta(g ^ h, k)
        assert difference % 2 == 0 and (0 not in (g, h, k) or difference == 0)
    print(
        "Verified the descent character for all 36 coefficients and all cyclic positions."
    )


if __name__ == "__main__":
    check_descent_parity()
    for s in (0, 1, -1):
        verify_line(s)
