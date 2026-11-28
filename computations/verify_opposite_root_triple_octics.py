"""Degree-eight triple-point equations from a finite edge's maximal minors."""

from itertools import combinations, combinations_with_replacement, product
from random import Random

import numpy as np
import sympy as sp

from verify_finite_square_chart import frame
from verify_semisimple_period_obstruction import P, reduced
from verify_single_edge_instanton_extension import coefficients_to_field


def monomials(degree):
    return list(combinations_with_replacement(range(4), degree))


def quotient_projection(rel, degree):
    x = sp.symbols("x:4")
    matrix = sp.Matrix(6, 4, lambda k, j: sum(x[i]*rel[4*i+j, k] for i in range(4)))
    quartics = []
    for rows in combinations(range(6), 4):
        polynomial = sp.Poly(matrix[list(rows), :].det(method="domain-ge"), *x)
        quartics.append({tuple(i for i, n in enumerate(exp) for _ in range(n)): coefficient
                         for exp, coefficient in polynomial.terms()})
    basis = {word: i for i, word in enumerate(monomials(degree))}
    multipliers = monomials(degree-4)
    multiples = sp.zeros(len(quartics)*len(multipliers), len(basis))
    for r, (quartic, pair) in enumerate(product(quartics, multipliers)):
        for word, coefficient in quartic.items():
            multiples[r, basis[tuple(sorted(word+pair))]] = coefficient
    rows, pivots = reduced(coefficients_to_field(multiples))
    free = [i for i in range(len(basis)) if i not in pivots]
    assert len(free) == 20
    projection = np.zeros((20, len(basis)), dtype=np.int64)
    projection[:, free] = np.eye(20, dtype=np.int64)
    projection[:, pivots] = -rows[:len(pivots), free].T
    return basis, projection


def octic_nullity(rel):
    six, projection = quotient_projection(rel, 6)
    eight = monomials(8)
    multipliers = monomials(2)
    blocks = []
    for i, j in multipliers:
        derivative = np.zeros((len(six), len(eight)), dtype=np.int64)
        for k, word in enumerate(eight):
            count = word.count(i)*(word.count(j)-int(i == j))
            if count:
                remaining = list(word)
                remaining.remove(i)
                remaining.remove(j)
                derivative[six[tuple(remaining)], k] = count
        blocks.append(projection@derivative % P)
    return len(eight)-len(reduced(np.vstack(blocks))[1])


def main():
    rng = Random(92871)
    rel = sp.Matrix(16, 6, lambda i, j: rng.randrange(-3, 4))
    assert octic_nullity(rel) == 0
    for points in frame():
        rows = []
        for point in points.tolist():
            for i, j in monomials(2):
                row = []
                for word in monomials(8):
                    value = word.count(i)*(word.count(j)-int(i == j))
                    if value:
                        remaining = list(word)
                        remaining.remove(i)
                        remaining.remove(j)
                        for letter in remaining:
                            value *= point[letter]
                    row.append(value)
                rows.append(row)
        assert len(reduced(coefficients_to_field(sp.Matrix(rows)))[1]) == 165
    print("Full degree-eight triple evaluation ranks for the general-edge test and both Clifford projections.")


if __name__ == "__main__":
    main()
