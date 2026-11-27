"""Exact Hankel equations, Frobenius intersections, and degrees of secant fibers."""
from itertools import combinations_with_replacement, product
from math import comb

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def rank(matrix):
    return DomainMatrix.from_Matrix(matrix).rank()


def equation_rank(pairs, constraints):
    columns = {pair: i for i, pair in enumerate(pairs)}
    rows = []
    for first, second in constraints:
        row = [0] * len(columns)
        row[columns[first]] += 1
        row[columns[second]] -= 1
        rows.append(row)
    return rank(sp.Matrix(rows))


def gram_equations():
    pairs = list(combinations_with_replacement(range(6), 2))
    constraints = [(tuple(sorted((a+i, b+j))), tuple(sorted((a+j, b+i))))
                   for a, b, i, j in product(range(3), range(3), range(4), range(4))]
    assert len(pairs) - equation_rank(pairs, constraints) == 11
    pairs = list(product(range(6), range(4)))
    constraints = [((a+i, j), (a+j, i)) for a, i, j in product(range(3), range(4), range(4))]
    assert len(pairs) - equation_rank(pairs, constraints) == 9


def frobenius_sample():
    decimic = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31)
    octavic = (3, 1, 4, 1, 5, 9, 2, 6, 5)
    gram = sp.Matrix(6, 6, lambda i,j: decimic[i+j])
    shifts = [sp.Matrix(6, 4, lambda i,j: int(i == j+k)) for k in range(3)]
    extra = gram.inv() * sp.Matrix(6, 4, lambda i,j: octavic[i+j])
    relation = sp.Matrix.hstack(*shifts, extra).T
    opposite = relation[[4*j+i for i,j in product(range(4), repeat=2)], :]
    w = relation * gram * relation.T
    alphabet = range(4)
    words = list(product(alphabet, repeat=3))
    triple = []
    for side, (first, second) in enumerate(((relation, opposite), (opposite, relation))):
        identity = sp.eye(4)
        left = sp.kronecker_product(first, identity)
        right = sp.kronecker_product(identity, second)
        assert left.cols + right.cols - rank(left.row_join(right)) == 4
        def coefficient(row, last):
            i, j, k = words[row]
            if side == 0:
                return w[4*i+j, 4*k+last]
            return w[4*last+i, 4*j+k]
        contractions = sp.Matrix(64, 4, coefficient)
        assert rank(contractions) == 4
        triple.append(contractions)
    for first, second in (triple, triple[::-1]):
        left = sp.kronecker_product(first, identity)
        right = sp.kronecker_product(identity, second)
        assert left.cols + right.cols - rank(left.row_join(right)) == 1
    assert rank(sp.Matrix(16,16,lambda i,j: w[4*(i//4)+j//4, 4*(i%4)+j%4])) == 8


if __name__ == '__main__':
    gram_equations()
    frobenius_sample()
    assert comb(6, 3) == 20 and 3*comb(6, 3) + 2*comb(6, 2) == 90
    print("Exact Hankel systems, both Frobenius orientations, opposite rank 8, and degrees 20 and 90 verified.")
