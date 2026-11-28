"""Finite arithmetic and polynomial identities used in the opposite-contact bound."""

from itertools import combinations_with_replacement
from math import isqrt

import sympy as sp


def planar_defect(n):
    return (n - isqrt(n - 1) - 1) // 2


def partitions(n, lower=1):
    if n == 0:
        yield ()
    for first in range(lower, n + 1):
        for tail in partitions(n - first, first):
            yield (first, *tail)


def check_lengths():
    for lengths in partitions(20):
        if len(lengths) < 10:
            continue
        assert sum(map(planar_defect, lengths)) <= 3
        for i, n in enumerate(lengths):
            if n < 8:
                continue
            others = lengths[:i] + lengths[i + 1:]
            assert max(others) < 8
            if n <= 9:
                assert n - 6 + sum(map(planar_defect, others)) <= 3
            else:
                assert n <= 11 and max(others) <= 2 and len(others) >= 9


def check_derivatives():
    coords = sp.symbols("w:3")
    conic = sum(w**2 for w in coords)
    monomials = [a*b for a, b in combinations_with_replacement(coords, 2)]
    rows = []
    for w in coords:
        for deriv in coords:
            form = sp.Poly(sp.diff(conic*w, deriv), *coords)
            rows.append([form.coeff_monomial(monomial) for monomial in monomials])
    assert sp.Matrix(rows).rank() == 6


if __name__ == "__main__":
    check_lengths()
    check_derivatives()
    print("Verified all length-20 partitions on at least ten sites and the quadratic derivative span.")
