"""An elliptic Bezout matrix recovered from its multiplication pairings at a good prime."""

from itertools import combinations
from random import Random

import numpy as np
import sympy as sp
from flint import nmod_mat

from verify_configuration_intersections import kernel
from verify_opposite_root_triple_octics import monomials
from verify_single_edge_instanton_extension import P


def multiplication(first, second, multiplier, x, y, cubic=None):
    basis = [1, x, x**2, x**3, x**4, y, x*y, x**2*y]
    equation = sp.Poly(y**2-(x**3-x+1 if cubic is None else cubic), y)
    rows = []
    for f in first:
        block = []
        for g in second:
            normal = sp.Poly(multiplier*f*g, y).rem(equation)
            polynomial = sp.Poly(sp.cancel(normal.as_expr()), x, y)
            row = [int(polynomial.coeff_monomial(w)) for w in basis]
            assert sp.expand(sum(c*w for c, w in zip(row, basis))-polynomial.as_expr()) == 0
            block.append(row)
        rows.append(block)
    return np.array(rows, dtype=np.int64) % P


def inverse_pairing(c, transform, mixed, symmetric):
    # The two sections annihilate the same octic moment. Compute it up to scale.
    second = transform@c
    constraints = np.vstack((np.einsum("i,ijk->jk", c, mixed),
                             np.einsum("j,ijk->ik", second, mixed))) % P
    return inverse_moment(constraints, symmetric)


def inverse_moment(constraints, symmetric):
    moment = kernel(nmod_mat(constraints.tolist(), P), P)
    if moment.ncols() != 1:
        return None
    coeff = np.asarray(moment.tolist(), dtype=np.int64)[:, 0]
    matrix = nmod_mat((np.einsum("ijk,k->ij", symmetric, coeff) % P).tolist(), P)
    if not matrix.det():
        return None
    return np.asarray(matrix.inv().tolist(), dtype=np.int64)


def main():
    x, y = sp.symbols("x y")
    base = [1, x, y, x**2]
    mixed = multiplication([1, x, y, (y+1)/x], [1, x, y, (y-1)/x], x, x, y)
    symmetric = multiplication(base, base, 1, x, y)
    def evaluate(c, transform):
        return inverse_pairing(c, transform, mixed, symmetric)

    verify_quadratic_presentation(evaluate, 580)
    print("Verified ten independent quadratic entries, nineteen quartic minors and eight sextic cofactors.")


def verify_quadratic_presentation(evaluate, seed):
    rng = Random(seed)
    transform = np.array([[rng.randrange(-4, 5) for _ in range(4)] for _ in range(4)], dtype=np.int64)
    pairs = monomials(2)
    rows = []
    for _ in range(20):
        c = np.array([rng.randrange(P) for _ in range(4)], dtype=np.int64)
        value = evaluate(c, transform)
        if value is None:
            continue
        entries = np.array([value[i, j] for i, j in pairs])
        pivot = next(i for i, v in enumerate(entries) if v)
        quad = np.array([c[i]*c[j] % P for i, j in pairs])
        equations = np.kron(np.eye(10, dtype=np.int64), quad)*entries[pivot]
        equations[:, 10*pivot:10*(pivot+1)] -= entries[:, None]*quad[None, :]
        rows.extend((equations % P).tolist())
    solution = kernel(nmod_mat(rows, P), P)
    assert solution.ncols() == 1
    change = np.asarray(solution.tolist(), dtype=np.int64).reshape(10, 10)
    assert nmod_mat(change.tolist(), P).rank() == 10
    coords = sp.symbols("z:4")
    quadrics = [sum(int(c)*coords[i]*coords[j] for c, (i, j) in zip(row, pairs)) for row in change]
    matrix = sp.Matrix(4, 4, lambda i, j: quadrics[pairs.index(tuple(sorted((i, j))))])

    def coefficients(expression, degree):
        polynomial = sp.Poly(expression, *coords, modulus=P)
        return [int(polynomial.coeff_monomial(tuple(w.count(i) for i in range(4)))) % P
                for w in monomials(degree)]

    minors = [coefficients(matrix.extract(i, j).det(), 4)
              for i in combinations(range(4), 2) for j in combinations(range(4), 2)]
    cofactors = [coefficients(matrix.minor_submatrix(i, j).det(method="domain-ge")*(-1)**(i+j), 6)
                 for i, j in pairs]
    assert nmod_mat(minors, P).rank() == 19 and nmod_mat(cofactors, P).rank() == 8


if __name__ == "__main__":
    main()
