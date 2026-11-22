#!/usr/bin/env python3
"""Exact trace ranks for the degree-four pencils glued at one elliptic multiple point."""
from itertools import combinations_with_replacement, product

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def sections(degrees, parameter):
    # Polynomial coordinates on each branch at zero, concatenated in branch order.
    offsets = [sum(degrees[:i])+i for i in range(len(degrees))]
    width = sum(degrees)+len(degrees)
    constant = sp.zeros(width, 1)
    for start in offsets:
        constant[start] = 1
    constant[offsets[0]+1] = parameter
    result = [constant]
    for start in offsets[1:]:
        v = sp.zeros(width, 1)
        v[start+1], v[offsets[0]+1] = 1, -1
        result.append(v)
    for start, degree in zip(offsets, degrees):
        for j in range(2, degree+1):
            v = sp.zeros(width, 1)
            v[start+j] = 1
            result.append(v)
    return sp.Matrix.hstack(*result)


def multiplication(d, e, left, right, symmetric=False):
    pairs = list(combinations_with_replacement(range(left.cols), 2)) if symmetric else list(product(range(left.cols), range(right.cols)))
    output = []
    for col, other in pairs:
        i = j = 0
        coefficients = []
        for a, b in zip(d, e):
            values = [0]*(a+b+1)
            for s in range(a+1):
                for t in range(b+1):
                    values[s+t] += left[i+s, col]*right[j+t, other]
            coefficients.extend(values)
            i += a+1
            j += b+1
        output.append(sp.Matrix(coefficients))
    return sp.Matrix.hstack(*output)


def check(d, e):
    n = tuple(a+b for a, b in zip(d, e))
    dd, ee = [2*a for a in d], [2*b for b in e]
    x, y, h = sections(d, 1), sections(e, 3), sections(n, 4)
    xy = multiplication(d, e, x, y)
    xx = multiplication(d, d, x, x, True)
    yy = multiplication(e, e, y, y, True)
    assert xy.rank() == xx.rank() == yy.rank() == 8
    # b has only its highest coefficient on every branch; a is nonzero where b vanishes.
    b = sp.zeros(h.rows, 1)
    offset = 0
    for deg in n:
        b[offset+deg] = 1
        offset += deg+1
    w = h[:, 0].row_join(b)
    multiples = multiplication(n, n, w, h)
    assert multiples.rank() == 15
    products = multiplication(dd, ee, xx, yy)
    trace = next(v for v in multiples.T.nullspace() if not (v.T*products).is_zero_matrix)
    opposite = (trace.T*products).reshape(xx.cols, yy.cols)
    expected = 6 if sorted(zip(d, e)) == [(1, 3), (3, 1)] else 8
    assert DomainMatrix.from_Matrix(opposite).rank() == expected
    mixed = (trace.T*multiplication(n, n, xy, xy)).reshape(xy.cols, xy.cols)
    assert mixed.rank() == 6


if __name__ == "__main__":
    for d, e in [((3, 1), (3, 1)), ((2, 2), (2, 2)), ((3, 1), (1, 3)),
                 ((2, 1, 1), (2, 1, 1)), ((2, 1, 1), (1, 2, 1)), ((1,)*4, (1,)*4)]:
        check(d, e)
    print("PASS six multiple-point models: quadratic ranks, pencil codimension and both trace flattenings")
