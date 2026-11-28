"""A finite nonregular trisecant window and a generic obstruction to other middle roots."""

from itertools import combinations, product
from random import Random

import sympy as sp

from verify_opposite_root_triple_octics import monomials, octic_nullity
from verify_single_edge_instanton_extension import (
    P, coefficients_to_field, finite_projection, rank, standard_frobenius,
)


def null_columns(matrix):
    return matrix.to_DM().to_field().nullspace().to_Matrix().T


def middle_relations(transform):
    first = (sp.eye(5)[:, :4], sp.eye(5)[:, 1:])
    second = (sp.eye(6)[:, :5], sp.eye(6)[:, 1:])
    blocks = [v*transform*u for u in first for v in second]
    return sp.Matrix(16, 6, lambda a, b: blocks[a//4][b, a%4])


def finite_both(rels):
    # Maximal rank on a plane after good reduction certifies maximal rank over Q.
    small = sp.Matrix(coefficients_to_field(rels).tolist())
    finite_projection(small, prime=P)
    finite_projection(small[[4*(i%4)+i//4 for i in range(16)], :], prime=P)


def primitive_polynomial(rel):
    variables = sp.symbols("x:4")
    matrix = sp.Matrix(6, 4, lambda i, j: sum(variables[k]*rel[4*k+j, i] for k in range(4)))
    words = monomials(4)
    exponents = [tuple(word.count(i) for i in range(4)) for word in words]
    minors = [sp.Poly(matrix[list(rows), :].det(method="domain-ge"), *variables)
              for rows in combinations(range(6), 4)]
    quotient = null_columns(sp.Matrix([[p.coeff_monomial(e) for e in exponents] for p in minors])).T
    assert quotient.rows == 20

    def multiply(i):
        return quotient[:, [words.index(tuple(sorted((i,)+word))) for word in monomials(3)]]

    affine = multiply(0).inv()*multiply(1)
    polynomial = sp.Poly.from_list(affine.to_DM().charpoly(), sp.Symbol("z"), domain=sp.QQ)
    return polynomial.clear_denoms(convert=True)[1].primitive()[1]


def general_pair():
    transform = sp.Matrix([[-1, 5, -5, 1, -4], [1, 2, 1, 5, 4], [5, 0, -4, -1, 2],
                           [-3, -1, -1, 4, -5], [1, -2, -5, -2, 5]])
    rel = middle_relations(transform)
    finite_both(rel)
    finite_both(middle_relations(transform.inv()))
    assert octic_nullity(rel) == 0
    polynomial = primitive_polynomial(rel)
    for prime, degrees in [(31, [20]), (53, [1, 19]), (61, [1, 2, 17])]:
        reduction = polynomial.set_modulus(prime)
        assert reduction.degree() == 20 and reduction.gcd(reduction.diff()).degree() == 0
        assert sorted(f.degree() for f, _ in reduction.factor_list()[1]) == degrees


def common_quadric_transform():
    form = sp.zeros(5)
    form[0, 4] = form[4, 0] = 1
    form[1, 3] = form[3, 1] = 3
    form[2, 2] = -8
    rng = Random(576)
    transform = sp.eye(5)
    for _ in range(6):
        v = sp.Matrix([rng.randrange(-3, 4) for _ in range(5)])
        value = (v.T*form*v)[0]
        if value:
            transform = (sp.eye(5)-2*v*v.T*form/value)*transform
    transform = transform.T
    assert transform*form*transform.T == form
    return transform


def common_quadric():
    transform = common_quadric_transform()
    hankel = [sp.Matrix(5, 5, lambda i, j: int(i+j == h)) for h in range(9)]
    equations = sp.Matrix.hstack(*(m.reshape(25, 1) for m in hankel+[-transform.T*h*transform for h in hankel]))
    coefficients = null_columns(equations)[9:, :]
    assert coefficients.shape == (9, 4)
    forms = [sum((coefficients[j, k]*hankel[j] for j in range(9)), sp.zeros(5)) for k in range(4)]
    first = (sp.eye(5)[:, :4], sp.eye(5)[:, 1:])
    blocks = [[u.T*transform.T*g*v for g in forms] for u in first for v in first]
    tensor = sp.Matrix(16, 16, lambda a, b: blocks[a%4][b%4][a//4, b//4])
    standard_frobenius(tensor)
    # Every cyclic adjacent flattening has rank six over Q. Inspect both projections.
    indices = list(product(range(4), repeat=4))
    for shift in range(4):
        rotation = sp.zeros(16)
        for word in indices:
            a, b, c, d = word[shift:]+word[:shift]
            i, j, k, l = word
            rotation[4*a+b, 4*c+d] = tensor[4*i+j, 4*k+l]
        pivots = rotation.to_DM().rref()[1]
        assert len(pivots) == 6
        finite_both(rotation[:, list(pivots)])
    opposite = sp.Matrix(16, 16, lambda a, b: tensor[4*(a//4)+b//4, 4*(a%4)+b%4])
    assert rank(opposite) == 16

    transformed = [transform.T*g*transform for g in forms]
    phi, psi = [sp.Matrix(9, 4, lambda h, k: family[k][min(h, 4), max(h-4, 0)])
                for family in (forms, transformed)]
    indices = list(product(range(3), range(3), range(4)))
    rows = [[int(i == u)*phi[h+j, k] for i, j, k in indices] for u, h in product(range(3), range(7))]
    rows += [[int(j == v)*psi[h+i, k] for i, j, k in indices] for v, h in product(range(3), range(7))]
    solutions = null_columns(sp.Matrix(rows))
    assert solutions.shape == (36, 4)
    maps = [sp.Matrix(3, 3, lambda i, j: solutions[4*(3*i+j)+k, 0]) for k in range(4)]
    assert rank(sp.Matrix.hstack(*maps)) == rank(sp.Matrix.vstack(*maps)) == 3


if __name__ == "__main__":
    general_pair()
    common_quadric()
    print("Verified the S20 certificate, triple-octic exclusion, finite window, and (3,2,2,3) simple.")
