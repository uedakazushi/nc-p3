"""Differentiate the forgetful maps to the first edge and then to its projected scheme."""

from itertools import combinations_with_replacement

from flint import nmod_mat
import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix

from verify_klein_elliptic_completion import rational_example
from verify_klein_full_rank_configurations import partner
from verify_klein_two_doublets import modular_points
from verify_unrestricted_doublet_deformations import (
    K,
    exchange,
    kernel,
    normal_data,
    rank,
    symmetric_example,
    tensor_basis,
)


def fixed_edge(w):
    second_normal, _ = normal_data(exchange(w))
    first_left = sp.Matrix.hstack(*kernel(w.T)).T
    dimensions = []
    for h in K:
        columns = [
            (first_left * t).reshape(160, 1).col_join(second_normal(exchange(t)))
            for t in tensor_basis(h, 0)
        ]
        dimensions.append(len(columns) - rank(sp.Matrix.hstack(*columns)))
    assert dimensions == [7, 8, 8, 8]


def fixed_projection(w, prime):
    first = modular_points(w, prime)
    second = partner(first, w, False, prime)
    relations = sp.Matrix.hstack(*w.columnspace()).applyfunc(
        lambda x: int(x.p) * pow(int(x.q), -1, prime) % prime
    )
    equations = []
    for a, b in zip(first, second):
        evaluation = sp.Matrix(
            [list(sp.Matrix(1, 4, a) * relations[:, i].reshape(4, 4)) for i in range(6)]
        )
        left_kernel = (
            DomainMatrix.from_Matrix(evaluation.T).convert_to(GF(prime)).nullspace()
        )
        assert left_kernel.shape == (3, 6)
        for row in left_kernel.to_Matrix().tolist():
            equations.append(
                [int(c) * x * y % prime for c in row for x in a for y in b]
            )
    assert nmod_mat(equations, prime).rank() == 45
    print(
        f"Fixed first edge: 31 tangent directions; fixed first projection: codimension 45 at p={prime}."
    )


def repeated_first_space(w):
    flip = [4 * j + i for i in K for j in K]
    n = exchange(w)
    overlaps = []
    for edge in (w, n[flip, :], w.T, n.T[flip, :]):
        relations = sp.Matrix.hstack(*edge.columnspace())
        left = sp.kronecker_product(relations, sp.eye(4))
        right = sp.kronecker_product(sp.eye(4), relations[flip, :])
        overlaps.append(48 - rank(left.row_join(right)))
    return overlaps


def square_roots(points, prime, doubled):
    domain = GF(prime)

    def matrix(data):
        return DomainMatrix.from_Matrix(sp.Matrix(data)).convert_to(domain)

    pairs = list(combinations_with_replacement(K, 2))
    square = matrix([[q[i] * q[j] for i, j in pairs] for q in points])
    normal = square.transpose().nullspace().to_Matrix().tolist()
    pivots = matrix(points).transpose().rref()[1]
    free = [i for i in range(20) if i not in pivots]
    jac = matrix(
        [
            [
                int(row[n]) * ((k == i) * points[n][j] + (k == j) * points[n][i])
                for n in free
                for k in K
            ]
            for i, j in pairs
            for row in normal
        ]
    )
    assert jac.rank() == (63 if doubled else 64)
    if doubled:
        (tangent,) = jac.nullspace().to_Matrix().tolist()
        rhs = [
            sum(
                int(row[n]) * int(tangent[4 * m + i]) * int(tangent[4 * m + j])
                for m, n in enumerate(free)
            )
            for i, j in pairs
            for row in normal
        ]
        assert jac.hstack(matrix(rhs)).rank() == 64


if __name__ == "__main__":
    for w, prime in [(symmetric_example(), 12619), (rational_example()[0], 2657)]:
        assert repeated_first_space(w) == ([4, 0, 4, 0] if prime == 12619 else [0] * 4)
        fixed_edge(w)
        fixed_projection(w, prime)
        points = modular_points(w, prime)
        first = points
        n = exchange(w)
        for k, (mat, flip) in enumerate(
            [(w, False), (n, True), (w.T, False), (n.T, True)]
        ):
            if prime == 12619 and k == 2:
                assert points == first
            square_roots(points, prime, prime == 12619 and k % 2 == 0)
            points = partner(points, mat, flip, prime)
        assert points == first
        print(
            "Checked all four square-root tangents and the collision where applicable."
        )
