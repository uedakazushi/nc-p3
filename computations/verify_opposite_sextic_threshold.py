"""Check reciprocal linear systems symbolically and double sextic ranks on both Clifford frames."""

from math import prod
from random import Random

import numpy as np

import sympy as sp

from verify_finite_square_chart import frame
from verify_opposite_root_triple_octics import monomials, quotient_projection
from verify_semisimple_period_obstruction import P, reduced
from verify_single_edge_instanton_extension import coefficients_to_field


def reciprocal_form():
    h, a, b, c, coefficient = sp.symbols("h a b c t")
    v = sp.Matrix([a, b, c])
    linear = sp.diag(coefficient*h, h, h, h)
    linear[1:, 0] = v
    linear[0, 1:] = v.T
    u = sp.Matrix([h, -a, -b, -c])
    q = coefficient*h**2 - v.dot(v)
    quadratic = u*u.T + q*sp.diag(0, 1, 1, 1)
    assert (linear.adjugate() - h*quadratic).applyfunc(sp.expand).is_zero_matrix


def double_sextic_nullity(points):
    words = monomials(6)
    jets = []
    for point in points.tolist():
        for i in range(4):
            row = []
            for word in words:
                letters = list(word)
                if i in letters:
                    letters.remove(i)
                    row.append(word.count(i)*prod(point[j] for j in letters))
                else:
                    row.append(0)
            jets.append(row)
    return len(words) - len(reduced(coefficients_to_field(sp.Matrix(jets)))[1])


def generic_edge():
    rng = Random(92871)
    relations = sp.Matrix(16, 6, lambda i, j: rng.randrange(-3, 4))
    quintics, projection = quotient_projection(relations, 5)
    words = monomials(6)
    blocks = []
    for i in range(4):
        partial = np.zeros((len(quintics), len(words)), dtype=np.int64)
        for j, word in enumerate(words):
            if i in word:
                letters = list(word)
                letters.remove(i)
                partial[quintics[tuple(letters)], j] = word.count(i)
        blocks.append(projection @ partial % P)
    assert len(reduced(np.vstack(blocks))[1]) == 80

    pairs = monomials(2)
    wedge = [pair for pair in pairs if pair[0] < pair[1]]
    matrix = sp.zeros(24)
    for k in range(6):
        for i in range(4):
            for j in range(4):
                for row, (a, b) in enumerate(wedge):
                    if i == b:
                        matrix[4*row+j, 4*k+i] += relations[4*a+j, k]
                    if i == a:
                        matrix[4*row+j, 4*k+i] -= relations[4*b+j, k]
    assert len(reduced(coefficients_to_field(matrix))[1]) == 24


def main():
    reciprocal_form()
    generic_edge()
    for points in frame():
        assert double_sextic_nullity(points) == 10
    print("Reciprocal identity and generic edge ranks 80 and 24 verified; Clifford double sextic nullities 10.")


if __name__ == "__main__":
    main()
