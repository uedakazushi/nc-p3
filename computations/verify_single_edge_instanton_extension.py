"""Recover a full regular four-window from a single rational symplectic edge."""

from itertools import combinations, combinations_with_replacement, product

import numpy as np
import sympy as sp

from verify_finite_square_tangent_instantons import edge, finite_projection, rank
from verify_semisimple_period_obstruction import P, reduced


def coefficients_to_field(matrix):
    def residue(x):
        numerator, denominator = x.as_numer_denom()
        return int(numerator)*pow(int(denominator), -1, P) % P

    return np.asarray(matrix.applyfunc(residue).tolist(), dtype=np.int64)


def companion(rel):
    pairs = list(combinations(range(4), 2))
    flip = sp.Matrix([[int(i == t)*rel[4*j+k, c]-int(j == t)*rel[4*i+k, c]
                       for t, c in product(range(4), range(6))]
                      for i, j in pairs for k in range(4)])
    overlaps = flip.to_DM().to_field().nullspace().to_Matrix().T
    assert overlaps.shape == (24, 4)
    tensor = sp.Matrix(16, 16, lambda iv, jz:
                       sum(rel[iv, c]*overlaps[6*(jz//4)+c, jz%4] for c in range(6)))
    pivots = tensor.T.to_DM().rref()[1]
    other = tensor.T[:, list(pivots)]
    assert other.shape == (16, 6)
    opposite = sp.Matrix(16, 16, lambda ij, vz: tensor[4*(ij//4)+vz//4, 4*(ij%4)+vz%4])
    assert rank(opposite) == 10
    return tensor, other


def symplectic_tangent(rel):
    basis = rel.T.to_DM().to_field().nullspace().to_Matrix().T
    coordinates = list(combinations(range(10), 2))
    equations = sp.Matrix([[
        sum(basis[4*u+i, k]*basis[4*v+j, l]-basis[4*v+j, k]*basis[4*u+i, l]
            for u, v in ((a, b), (b, a))) for k, l in coordinates]
        for a, b in combinations_with_replacement(range(4), 2) for i, j in combinations(range(4), 2)])
    solution, = equations.to_DM().to_field().nullspace().to_Matrix().tolist()
    form = sp.zeros(10)
    for (k, l), value in zip(coordinates, solution):
        form[k, l], form[l, k] = value, -value
    assert rank(form) == 10


def standard_frobenius(tensor):
    data = coefficients_to_field(tensor).reshape(4, 4, 4, 4)
    for n in range(1, 5):
        spaces = []
        expected = [4, 6, 4, 1][n-1]
        for i in range(4):
            entries = data.transpose(tuple((i+j)%4 for j in range(4))).reshape(4**n, -1)
            pivots = reduced(entries)[1]
            assert len(pivots) == expected
            spaces.append(entries[:, pivots])
        if n == 1:
            continue
        for i in range(4):
            intersection = np.column_stack((np.kron(spaces[i], np.eye(4, dtype=np.int64)),
                                            np.kron(np.eye(4, dtype=np.int64), spaces[(i+1)%4])))
            assert 8*expected-len(reduced(intersection)[1]) == {2: 4, 3: 1, 4: 0}[n]


def main():
    rel = edge()
    tensor, other = companion(rel)
    for rels in (rel, other):
        finite_projection(rels)
        finite_projection(rels[[4*(i%4)+i//4 for i in range(16)], :])
        symplectic_tangent(rels)
    standard_frobenius(tensor)
    print("Full Frobenius window, four finite projections and both symplectic tangent lines verified.")


if __name__ == "__main__":
    main()
