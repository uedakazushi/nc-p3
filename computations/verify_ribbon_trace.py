#!/usr/bin/env python3
"""Exact products on the genus-one split ribbon on a conic."""
from itertools import combinations, combinations_with_replacement

import sympy as sp


def sections(n, parameter, t):
    # A section on the affine line is f + eta*g, with eta^2=0 and eta transforming as dt.
    return ([(t**i, 0) for i in range(n)]
            + [(t**n, parameter*t**(n-1))]
            + [(0, t**i) for i in range(n-1)])


def multiply(x, y):
    a, b = x
    c, d = y
    return a*c, a*d+b*c


def coordinates(section, n, parameter, t):
    f, g = (sp.Poly(a, t) for a in section)
    assert f.degree() <= n and g.degree() <= n-1
    assert g.nth(n-1) == parameter*f.nth(n)
    return sp.Matrix([f.nth(i) for i in range(n+1)]+[g.nth(i) for i in range(n-1)])


def check(l, m):
    t = sp.symbols("t")
    n = l+m
    x, y, middle = sections(2, l, t), sections(2, m, t), sections(4, n, t)
    w = [middle[0], middle[4]]
    def column(a, b):
        return coordinates(multiply(a, b), 8, 2*n, t)
    multiples = sp.Matrix.hstack(*(column(a, b) for a in w for b in middle))
    trace, = multiples.T.nullspace()

    def pairing(left, right):
        return sp.Matrix([[trace.dot(column(a, b)) for b in right] for a in left])
    mixed = [multiply(a, b) for a in x for b in y]
    left = [multiply(a, b) for a, b in combinations_with_replacement(x, 2)]
    right = [multiply(a, b) for a, b in combinations_with_replacement(y, 2)]
    assert pairing(mixed, mixed).rank() == 6
    assert pairing(left, right).rank() == (8 if l != m else 6)
    # Nonzero values of the trace on all nilpotent products are not required,
    # but the restriction to their span must be nonzero.
    assert any(trace[-7:, 0])
    gaussian = [sp.Poly((j-i)*t**(i+j-1), t) for i, j in combinations(range(5), 2)]
    assert sp.Matrix([[g.nth(k) for g in gaussian] for k in range(7)]).rank() == 7


if __name__ == "__main__":
    check(1, 3)
    check(1, 1)
    print("PASS two ribbon pencils, adjacent/opposite ranks, socle on the nilradical and full Gaussian map")
