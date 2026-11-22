"""Exact nodal and cuspidal controls for the noninvertible degree-four limit."""
from fractions import Fraction
from itertools import product

import sympy as sp

from two_window_ext import PathRow

t = sp.symbols("t")


def coefficients(expressions, degree):
    polys = [sp.Poly(p, t) for p in expressions]
    return sp.Matrix(degree+1, len(polys), lambda i, j: polys[j].nth(i))


def section(polynomial, degree, nodal):
    p = sp.Poly(polynomial, t)
    return polynomial, p.nth(0)+t*p.nth(degree if nodal else 1)


def multiply(pair, other):
    return pair[0]*other[0], pair[1]*other[1]


def fourth_degree(sections):
    # Core degree 14 and exceptional degree 2. Its other two coordinates
    # are fixed by the conductor data, but keeping them makes both cases uniform.
    return coefficients([p for p, _ in sections], 14).col_join(
        coefficients([q for _, q in sections], 2))


def verify(nodal):
    a = [1+t**4, t, t**2, t**3] if nodal else [1, t**2, t**3, t**4]
    b = [1, t, t**2, t**3]
    w = [1+2*t+4*t**2+7*t**6+3*t**7, t+4*t**4+t**7]
    wa = [section(p, 7, nodal) for p in w]
    assert sp.resultant(*w, t)
    assert coefficients([q for _, q in wa], 1).rank() == 2
    left = [(p, sp.Poly(p, t).nth(0)) for p in a]
    right = [section(p, 3, nodal) for p in b]
    middle = [section(t**i, 7, nodal) for i in range(8)]
    products = [multiply(p, q) for p, q in product(wa, middle)]
    fourth = [multiply(p, q) for p, q in product(middle, repeat=2)]
    columns = fourth_degree(products+fourth)
    assert columns[:, :16].rank() == 15 and columns.rank() == 16
    trace = next(v for v in columns[:, :16].T.nullspace() if (v.T*columns).rank())

    def pairing(first, second):
        all_products = [multiply(p, q) for p, q in product(first, second)]
        return (trace.T*fourth_degree(all_products)).reshape(len(first), len(second))

    def all_pairs(first, second):
        return [multiply(p, q) for p, q in product(first, second)]

    mixed = all_pairs(left, right)
    assert pairing(mixed, mixed).rank() == 6
    assert pairing(all_pairs(left, left), all_pairs(right, right)).rank() == 8

    pencil_quotient = sp.Matrix.hstack(*coefficients(w, 7).T.nullspace()).T
    windows = []
    for first, second, expected in [(a, b, 12), (b, a, 13)]:
        forward = coefficients([p*q for p, q in product(first, second)], 7)
        back = coefficients([p*q for p, q in product(second, first)], 7)
        kernel = sp.Matrix.hstack(*forward.nullspace())
        reversed_kernel = sp.Matrix.hstack(*back.nullspace())
        quad_relations = sp.kronecker_product(kernel, sp.eye(4)).row_join(
            sp.kronecker_product(sp.eye(4), reversed_kernel))
        assert 64-quad_relations.rank() == expected
        r = (pencil_quotient*forward).T
        windows.append([{divmod(j, 4): Fraction(r[j, i]) for j in range(16) if r[j, i]}
                        for i in range(6)])
        r_back = (pencil_quotient*back).T
        cubic = sp.kronecker_product(r, sp.eye(4)).row_join(sp.kronecker_product(sp.eye(4), r_back))
        assert 48-cubic.rank() == 4
    for start, seventh in [(0, 120), (1, 121)]:
        assert PathRow(windows*2, start, 7).dims == [1, 4, 10, 20, 35, 56, 84, seventh]
    print("nodal" if nodal else "cuspidal", "flattenings 6,8; cubic hulls 12,13; window overlaps 4,4; degree seven 120,121")


if __name__ == "__main__":
    verify(True)
    verify(False)
