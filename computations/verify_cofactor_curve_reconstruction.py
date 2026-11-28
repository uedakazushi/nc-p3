"""Postulation of the general trisecant source and square-root embeddings of an elliptic octic."""

from itertools import combinations_with_replacement
from random import Random

import numpy as np
import sympy as sp
from flint import nmod_mat

from verify_opposite_root_triple_octics import monomials, octic_nullity, quotient_projection
from verify_single_edge_instanton_extension import P
from verify_trisecant_middle_simple_obstruction import middle_relations


def source_postulation():
    transform = sp.Matrix([[-1, 5, -5, 1, -4], [1, 2, 1, 5, 4], [5, 0, -4, -1, 2],
                           [-3, -1, -1, 4, -5], [1, -2, -5, -2, 5]])
    rel = middle_relations(transform)[[4*(i % 4)+i//4 for i in range(16)], :]
    degree5, projection = quotient_projection(rel, 5)
    six = monomials(6)
    rows = []
    for i in range(4):
        derivative = np.zeros((len(degree5), len(six)), dtype=np.int64)
        for j, word in enumerate(six):
            if i in word:
                rest = list(word)
                rest.remove(i)
                derivative[degree5[tuple(rest)], j] = word.count(i)
        rows.extend((projection@derivative % P).tolist())
    assert len(six)-nmod_mat(rows, P).rank() == 7
    assert octic_nullity(rel) == 1


def elliptic_matrices():
    x, y = sp.symbols("x y")
    cubic = sp.prod(x-r for r in (0, 1, 3))
    octic_basis = [1, x, x**2, x**3, x**4, y, x*y, x**2*y]
    pairs = list(combinations_with_replacement(range(4), 2))
    matrices = []
    for branch in [None, 0, 1, 3]:
        if branch is None:
            sections, multiplier = [1, x, y, x**2], 1
        else:
            sections, multiplier = [1, x, y, y/(x-branch)], x-branch
        coefficients = []
        for i, j in pairs:
            element = sp.Poly(multiplier*sections[i]*sections[j], y).rem(sp.Poly(y**2-cubic, y))
            polynomial = sp.Poly(sp.cancel(element.as_expr()), x, y)
            row = [polynomial.coeff_monomial(monomial) for monomial in octic_basis]
            assert sp.expand(sum(c*m for c, m in zip(row, octic_basis))-polynomial.as_expr()) == 0
            coefficients.append(row)
        matrix = sp.Matrix(coefficients)
        assert matrix.rank() == 8
        matrices.append(matrix)
    # Exact evaluations distinguish each of the four determinantal quartics in their pencil.
    rng = Random(579)
    values = []
    for _ in range(12):
        point = sp.Matrix([rng.randrange(-5, 6) for _ in octic_basis])
        row = []
        for matrix in matrices:
            forms = matrix*point
            symmetric = sp.Matrix(4, 4, lambda i, j: forms[pairs.index(tuple(sorted((i, j))))])
            row.append(symmetric.det())
        values.append(row)
    values = sp.Matrix(values)
    assert values.rank() == 2
    assert all(values[:, [i, j]].rank() == 2 for i in range(4) for j in range(i))


if __name__ == "__main__":
    source_postulation()
    elliptic_matrices()
    print("Verified the sextic and octic source dimensions, and the four elliptic square-root embeddings.")
