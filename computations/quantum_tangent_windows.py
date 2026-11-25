"""Finite-edge normal forms over a ternary quantum base with a maximal tangent point."""

from itertools import permutations

import sympy as sp

from curvilinear_window_orbits import torus_invariants


def finite_edges(c):
    if any(not any(row) for row in c):
        return False
    # At vertex j the two incident binary polynomials must not both vanish.
    return all(c[(j - 1) % 3][2] or c[(j - 2) % 3][0] for j in range(3))


def local_lengths(zeta, c):
    """Lengths at the unique distinguished point and at the other geometric points."""
    if not finite_edges(c):
        raise ValueError("edge is positive-dimensional")
    others = []
    for a, b, d in c:
        exponents = [s for s, v in enumerate((a, b, d)) if v]
        span = max(exponents) - min(exponents)
        if span == 2 and sp.expand(zeta**2 * b**2 - 4 * zeta * a * d) == 0:
            others.append(2)
        else:
            others.extend([1] * span)
    return 20 - sum(others), tuple(sorted(others))


def shifted(zeta, c, step=1):
    return tuple(
        tuple(sp.expand(zeta ** (step * (s + 1)) * v) for s, v in enumerate(row))
        for row in c
    )


def permuted(zeta, c, order):
    direct = order[1] == (order[0] + 1) % 3
    target = zeta if direct else 1 / zeta
    out = [[0] * 3 for _ in range(3)]
    for a, row in enumerate(c):
        for s, v in enumerate(row):
            j = s if direct else 2 - s
            factor = 1 if direct else target ** (-s * (2 - s))
            out[order[a]][j] = sp.expand(v * factor)
    return target, out


def labelled_coordinates(c):
    support, weights, values = [], [], []
    for a, row in enumerate(c):
        for s, v in enumerate(row):
            if not v:
                continue
            w = [0] * 3
            w[a], w[(a + 1) % 3], w[(a + 2) % 3] = -1, s, 2 - s
            support.append((a, s))
            weights.append(w)
            values.append(sp.sympify(v))
    return tuple(support), tuple(
        sp.cancel(v) for v in torus_invariants(weights, values)
    )


def equivalent(zeta, c, other, d):
    zeta, other = sp.sympify(zeta), sp.sympify(other)
    key = labelled_coordinates(d)
    for order in permutations(range(3)):
        target, cc = permuted(zeta, c, order)
        if target == other and labelled_coordinates(cc) == key:
            return True
    return False


def period(zeta, c):
    if not finite_edges(c):
        raise ValueError("this period test assumes finite edges")
    return next(n for n in (1, 2, 4) if equivalent(zeta, c, zeta, shifted(zeta, c, n)))
