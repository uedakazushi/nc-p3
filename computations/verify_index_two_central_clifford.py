"""Adjugate polarization and unique simultaneous symmetric identifications in rank two."""
from itertools import combinations_with_replacement

import sympy as sp


def balanced_symmetric_space(matrix):
    basis = []
    for i, j in combinations_with_replacement(range(4), 2):
        form = sp.zeros(4)
        form[i, j] = form[j, i] = 1
        basis.append((matrix*form - form*matrix.T).reshape(16, 1))
    return 10 - sp.Matrix.hstack(*basis).rank()


def main():
    a = sp.Matrix(2, 2, sp.symbols("a:4"))
    b = sp.Matrix(2, 2, sp.symbols("b:4"))
    scalar = sp.expand((a+b).det() - a.det() - b.det())
    target = sp.eye(2) * scalar
    for polarized in [a.adjugate()*b + b.adjugate()*a, a*b.adjugate() + b*a.adjugate()]:
        assert (polarized - target).applyfunc(sp.expand) == sp.zeros(2)
    rank_one = sp.zeros(4)
    rank_one[0,2] = 1
    rank_two = rank_one.copy()
    rank_two[1,3] = 1
    for matrix, expected in [(sp.diag(0,0,0,1), 7), (sp.diag(0,0,1,1), 6),
                             (rank_one, 7), (rank_two, 6)]:
        assert balanced_symmetric_space(matrix) == expected
    t = sp.symbols("t")
    assert sp.cancel((1+t)**4/(1-t**2)**4 - 1/(1-t)**4) == 0
    print("Both adjugate products, all critical block dimensions, and the Clifford Hilbert series.")


if __name__ == "__main__":
    main()
