"""Exact finite edge from a charge-four symplectic monad, and its two tangent ranks."""

from itertools import combinations, combinations_with_replacement

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def rank(matrix):
    return DomainMatrix.from_Matrix(matrix).rank()


def edge():
    x, y, z, w = sp.symbols("x y z w")
    slopes = [0, z+w, 2*z+5*w, 4*z+16*w]
    first = [-5*z+w, -13*z+23*w, 2*z+113*w, 8*z+95*w]
    second = [z+7*w, 7*z+43*w, 10*z+61*w, 4*z+25*w]
    middle = sp.diag(2*z-3*w, -z+2*w, 7*z+w, 5*z+13*w)
    for i, j in combinations(range(4), 2):
        entry = sp.cancel((first[j]*second[i]-first[i]*second[j])/(slopes[i]-slopes[j]))
        assert sp.Poly(entry, z, w).total_degree() == 1
        middle[i, j] = middle[j, i] = entry
    monad = sp.diag(*(x+a for a in slopes)).col_join(y*sp.eye(4)+middle)
    monad = monad.col_join(sp.Matrix([first, second]))
    omega = sp.zeros(10)
    for i, j in [(0, 4), (1, 5), (2, 6), (3, 7), (8, 9)]:
        omega[i, j], omega[j, i] = 1, -1
    assert (monad.T*omega*monad).applyfunc(sp.expand) == sp.zeros(4)
    coefficients = sp.Matrix.hstack(*(monad.diff(v) for v in (x, y, z, w)))
    assert rank(coefficients.T*omega*coefficients) == 10
    return sp.Matrix.hstack(*coefficients.nullspace())


def finite_projection(rel, prime=None):
    x, y, z = sp.symbols("x y z")
    plane = [x, y, z, 2*x+3*y+5*z]
    evaluation = sp.Matrix([[sum(plane[i]*rel[4*i+j, k] for i in range(4))
                             for j in range(4)] for k in range(6)])
    minors = [sp.Poly(evaluation[list(rows), :].det(method="domain-ge"), x, y, z)
              for rows in combinations(range(6), 4)]
    exponents = [(i, j, 4-i-j) for i in range(5) for j in range(5-i)]
    coefficients = sp.Matrix([[p.coeff_monomial(e) for p in minors] for e in exponents])
    domain_matrix = coefficients.to_DM()
    if prime is not None:
        domain_matrix = domain_matrix.convert_to(sp.GF(prime))
    assert domain_matrix.rank() == 15


def tangent_and_flip(rel):
    pairs = list(combinations(range(4), 2))
    sym = list(combinations_with_replacement(range(4), 2))

    def wedge(i, j, t, v, c):
        return int(i == t)*rel[4*j+v, c] - int(j == t)*rel[4*i+v, c]

    first = sp.Matrix([[wedge(i, j, t, v, c) for t in range(4) for c in range(6)]
                       for i, j in pairs for v in range(4)])
    second = sp.Matrix([[int(k == v)*wedge(i, j, t, h, c) + int(k == h)*wedge(i, j, t, v, c)
                         for t in range(4) for k in range(4) for c in range(6)]
                        for i, j in pairs for v, h in sym])
    return 60-rank(second), 24-rank(first)


def main():
    rel = edge()
    reverse = rel[[4*(i%4)+i//4 for i in range(16)], :]
    for presentation in (rel, reverse):
        finite_projection(presentation)
    assert tangent_and_flip(rel) == (1, 4)
    assert tangent_and_flip(reverse) == (0, 0)
    print("Both finite projections, rank-ten tangent, and directed tangent/flip coranks verified over Q.")


if __name__ == "__main__":
    main()
