"""Smooth K3 trilinear slice and the differential of its completion map."""
from itertools import combinations_with_replacement, permutations

import numpy as np
from flint import nmod_mat
import sympy as sp

from verify_klein_elliptic_completion import rational_example

PRIME = 1009


def matrix(rows):
    return nmod_mat(np.asarray(rows, dtype=np.int64).tolist(), PRIME)


def null_columns(m):
    basis, dimension = m.nullspace()
    return np.asarray(basis.tolist(), dtype=np.int64)[:, :dimension]


def exponents(degree):
    return [tuple(indices.count(i) for i in range(4))
            for indices in combinations_with_replacement(range(4), degree)]


def smooth_quartic(tensor):
    vars_ = sp.symbols('x:4')
    m = sp.Matrix(4, 4, lambda i, j: sum(int(tensor[k, i, j])*x
                                          for k, x in enumerate(vars_)))
    polynomial = sp.Poly(m.det(method='domain-ge'), *vars_, modulus=PRIME)
    target = {h: i for i, h in enumerate(exponents(9))}
    rows = []
    for x in vars_:
        derivative = polynomial.diff(x).terms()
        for h in exponents(6):
            row = [0]*len(target)
            for q, value in derivative:
                row[target[tuple(a+b for a, b in zip(h, q))]] = int(value)
            rows.append(row)
    return matrix(rows).rank()


def main():
    w, _, _ = rational_example()
    data = np.array([int(x.p)*pow(int(x.q), -1, PRIME) % PRIME for x in w],
                    dtype=np.int64).reshape((4,)*4)
    covector = np.array([1, 2, 3, 5], dtype=np.int64)
    t = np.tensordot(covector, data, axes=1) % PRIME
    smooth = [smooth_quartic(np.moveaxis(t, axis, 0)) for axis in range(3)]
    equations = []
    for axes in [(0, 1, 2, 3), (1, 2, 3, 0)]:
        m = matrix(data.transpose(axes).reshape(16, 16))
        assert m.rank() == 6
        left, right = null_columns(m.transpose()), null_columns(m)
        for u in left.T:
            for v in right.T:
                row = np.outer(u, v).reshape((4,)*4).transpose(np.argsort(axes))
                equations.append(row.reshape(256) % PRIME)
    relations = matrix(equations)
    tangent = null_columns(relations)
    image = np.tensordot(covector, tangent.reshape(4, 64, -1), axes=1)
    differential = matrix(image % PRIME).rank()
    print('Jacobian quartic ranks', smooth, 'cut tangent dimension', tangent.shape[1],
          'slice differential rank', differential)
    assert smooth == [220, 220, 220]
    assert tangent.shape[1] == 76 and differential == 64


if __name__ == '__main__':
    main()
