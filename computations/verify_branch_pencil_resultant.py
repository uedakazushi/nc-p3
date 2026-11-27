"""Exact complementary cut identities, on all pairs of biquadratic basis forms."""
from itertools import product

import sympy as sp

from verify_branch_cubic_contact import flattenings


def evaluate(form, x, t):
    return sum(form[i, j]*x[0]**(form.rows-1-i)*x[1]**i
               *t[0]**(form.cols-1-j)*t[1]**j
               for i, j in product(range(form.rows), range(form.cols)))


def verify(B):
    # A grid of five affine values interpolates each binary quartic.
    v = [sp.Matrix([1, a]) for a in range(5)]
    J = sp.Matrix([[0, -1], [1, 0]])
    T, U = J*B.T, J*B
    assert U*T == -B.det()*sp.eye(2)
    rows = ((2, 4, 6, 8, 10, 12), (1, 2, 3, 8, 9, 10))
    extra = ((3, 7, 11, 15), (12, 13, 14, 15))
    core = flattenings(B, sp.zeros(3), sp.zeros(3), sp.zeros(4))
    inv = [c.extract(r, r).inv() for c, r in zip(core, rows)]
    basis = [sp.eye(9)[:, i].reshape(3, 3) for i in range(9)]
    blocks = [flattenings(B, f, f, sp.zeros(4)) for f in basis]
    grid = sp.Matrix(5, 4, lambda i, j: i**j)
    bgrid = grid[:, :2]*B*grid[:, :2].T
    values = [grid[:, :3]*f*grid[:, :3].T for f in basis]
    left = [sp.Matrix([evaluate(f, x, T*x) for x in v]) for f in basis]
    right = [sp.Matrix([evaluate(f, U*t, t) for t in v]) for f in basis]
    for i, j in product(range(9), repeat=2):
        for k, (r, z) in enumerate(zip(rows, extra)):
            # The even variables change sides in the second cut.
            d = blocks[i][k].extract(z, r)*inv[k]*blocks[j][k].extract(r, z)
            if k:
                d = d.T
            discrepancy = values[i].multiply_elementwise(values[j]) - bgrid.multiply_elementwise(grid*d*grid.T)
            a, b = ((i, j) if k == 0 else (j, i))
            assert discrepancy*B.det()**2 == left[a]*right[b].T


if __name__ == '__main__':
    for matrix in (sp.Matrix([[1, 2], [3, 4]]), sp.eye(2), sp.Matrix([[0, 1], [-1, 0]])):
        verify(matrix)
    print("PASS: 486 Schur-product identities, including both marked-point collisions")
