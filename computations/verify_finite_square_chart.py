#!/usr/bin/env python3
"""Length-twenty symmetric CI example: traces, square spaces and their Grassmann tangent."""
from itertools import combinations, combinations_with_replacement

import sympy as sp


def products(left, right, symmetric=False):
    pairs = (combinations_with_replacement(range(left.cols), 2) if symmetric else
             ((i, j) for i in range(left.cols) for j in range(right.cols)))
    return sp.Matrix.hstack(*(left[:, i].multiply_elementwise(right[:, j]) for i, j in pairs))


def frame():
    forms = [sp.eye(4)[:, i] for i in range(4)] + [sp.ones(4, 1)]
    pairs = []
    for a, b in combinations(forms, 2):
        pairs.extend([(a+b, a-b), (a-b, a+b)])
    return tuple(sp.Matrix.vstack(*(point[i].T for point in pairs)) for i in (0, 1))


def tangent_dimension(points, squares):
    # Linearize U*U in the fixed ten-space of squares, allowing arbitrary lifted columns.
    rows = []
    for functional in squares.T.nullspace():
        for i, j in combinations_with_replacement(range(4), 2):
            row = [0]*80
            for k in range(20):
                row[4*k+i] += functional[k]*points[k, j]
                row[4*k+j] += functional[k]*points[k, i]
            rows.append(row)
    return 80 - sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(rows)).rank() - 16


def verify_ratio_chart(points, square_basis, genuine):
    square_basis = sp.Matrix.hstack(*square_basis.columnspace())
    unit = points*sp.Matrix([1, 3, 9, 27])
    other = points*sp.Matrix([2, -1, 4, 7])
    assert all(unit)
    ratio = sp.diag(*(1/v for v in unit))*other
    multiplier = sp.diag(*ratio)
    annihilator = sp.Matrix.hstack(*square_basis.T.nullspace()).T
    matrix = annihilator*multiplier*square_basis
    kernel = square_basis*sp.Matrix.hstack(*matrix.nullspace())
    assert kernel.cols == 4
    kernel_squared = products(kernel, kernel, True)
    assert sp.polys.matrices.DomainMatrix.from_Matrix(kernel_squared).rank() == 10
    normals = sp.Matrix.hstack(*square_basis.row_join(multiplier*square_basis).T.nullspace())
    tangent_equations = products(kernel, normals)
    assert sp.polys.matrices.DomainMatrix.from_Matrix(tangent_equations).rank() == 13
    annihilator2 = sp.Matrix.hstack(*kernel_squared.T.nullspace()).T
    equations = sp.Matrix.vstack(*(annihilator2*sp.diag(*q) for q in square_basis.columnspace()))
    assert sp.polys.matrices.DomainMatrix.from_Matrix(equations).rank() == (19 if genuine else 20)
    if genuine:
        assert kernel.row_join(sp.diag(*unit)*points).rank() == 4
        assert equations*unit.multiply_elementwise(unit) == sp.zeros(equations.rows, 1)


def verify_support_component(squares):
    q = sp.Matrix.hstack(*squares.columnspace())
    positions = [0, 2, 3, 7, 8, 16]
    diagonal = sp.diag(*(int(i in positions) for i in range(20)))
    kernel = q*sp.Matrix.hstack(*q[positions, :].nullspace())
    normals = sp.Matrix.hstack(*q.row_join(diagonal*q).T.nullspace())
    assert kernel.cols == normals.cols == 4
    assert kernel[positions, :] == sp.zeros(6, 4)
    assert sp.polys.matrices.DomainMatrix.from_Matrix(products(kernel, normals)).rank() == 13


def verify_nonreduced_ratio():
    # A separate square-space control in Q[t]/t^20, in divided-power coordinates.
    # exp(a*t) exp(b*t) = exp((a+b)*t) for the truncated exponentials.
    rates = (0, 1, 4, 16)
    def columns(values):
        return sp.Matrix(20, len(values), lambda k, i: values[i]**k)
    quadratic = [a+b for a, b in combinations_with_replacement(rates, 2)]
    cubic = [a+b+c for a, b, c in combinations_with_replacement(rates, 3)]
    assert len(set(cubic)) == 20  # Vandermonde gives the full cubic evaluation.
    square = columns(quadratic)
    # s=exp(4*t), t'=exp(16*t), hence their ratio has exponent 12.
    ratio_map = sp.Matrix(20, 20, lambda k, j:
                          sp.binomial(k, j)*12**(k-j) if j <= k else 0)
    assert ratio_map*square == columns([12+a for a in quadratic])
    annihilator = sp.Matrix.hstack(*square.T.nullspace()).T
    matrix = annihilator*ratio_map*square
    kernel = square*sp.Matrix.hstack(*matrix.nullspace())
    assert kernel.cols == 4 and kernel.row_join(columns([4+a for a in rates])).rank() == 4


def main():
    x, y = frame()
    xx, yy = products(x, x, True), products(y, y, True)
    xy = products(x, y)
    assert xx.rank() == yy.rank() == xy.rank() == 10
    for linear, square in ((x, xx), (y, yy)):
        assert products(linear, square).rank() == 20
    trace, = products(xx, yy).T.nullspace()
    assert all(trace)
    assert xx.T*sp.diag(*trace)*yy == sp.zeros(10)
    assert xy.T*sp.diag(*trace)*xy == sp.zeros(16)
    verify_ratio_chart(x, xx, True)
    verify_ratio_chart(y, xx, False)
    verify_support_component(xx)
    verify_nonreduced_ratio()
    tangent = tangent_dimension(x, xx)
    assert tangent == 0
    print("PASS 20 points, cubic ranks 20, Frobenius trace and orthogonality;",
          "square-space Grassmann tangent dimension", tangent)
    return tangent


if __name__ == "__main__":
    main()
