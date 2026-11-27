"""Construct the charge-four linear monad of the finite companion window."""
from itertools import combinations, combinations_with_replacement, permutations, product

import numpy as np
import sympy as sp

from verify_semisimple_period_obstruction import P, reduced
from verify_twenty_point_companion import gram_half


def kernel(matrix):
    return matrix.to_DM().to_field().nullspace().to_Matrix().T


def monad(relations):
    generators = kernel(relations)
    columns = [generators[:, j].reshape(4, 4) for j in range(10)]
    pairs = list(combinations_with_replacement(range(4), 2))
    quad = sp.Matrix(40, 40, lambda r, c: sum(columns[c//4][i, r % 4] for i in range(4)
                     if tuple(sorted((i, c % 4))) == pairs[r//4]))
    linear = kernel(quad)
    assert linear.cols == 4
    x = sp.symbols('x:4')
    b = sp.Matrix(4, 10, lambda i, j: sum(x[k]*columns[j][k, i] for k in range(4)))
    a = sp.Matrix(10, 4, lambda i, j: sum(x[k]*linear[4*i+k, j] for k in range(4)))
    assert b*a == sp.zeros(4)
    return x, a


def minor_span(variables, matrix):
    def residue(c):
        return int(c.p)*pow(int(c.q), -1, P) % P

    coeffs = np.array([matrix.applyfunc(lambda f: residue(f.coeff(x))).tolist()
                       for x in variables], dtype=np.int64)
    # Evaluation on these 35 points tests homogeneous quartics; rank 35 alone
    # suffices for the conclusion even without using unisolvence of the grid.
    grid = [(1, a, b, c) for a, b, c in product(range(5), repeat=3) if a+b+c <= 4]
    values = (np.array(grid)@coeffs.reshape(4, -1) % P).reshape(35, 10, 4)
    submatrices = values[:, list(combinations(range(10), 4)), :]
    det = np.zeros(submatrices.shape[:2], dtype=np.int64)
    for perm in permutations(range(4)):
        sign = (-1)**sum(a > b for a, b in combinations(perm, 2))
        det = (det+sign*np.prod(submatrices[:, :, range(4), perm], axis=2)) % P
    return len(reduced(det)[1])


def main():
    z = sp.Symbol('z')
    assert ((4+6*z+4*z**2)*(1-z)**4).expand() == 4-10*z+4*z**2+4*z**3+4*z**4-10*z**5+4*z**6
    variables, a = monad(gram_half(1))
    assert minor_span(variables, a) == 35
    print('PASS: four linear syzygies, monad composition zero and full quartic minor span')


if __name__ == '__main__':
    main()
