"""Nodal and cuspidal genus-one Bezout kernels, including the boundary Hankel determinant."""

from itertools import combinations
from random import Random

import sympy as sp
from flint import nmod_mat

from verify_opposite_root_triple_octics import monomials
from verify_single_edge_instanton_extension import P


def normal_form(kind, sections, t, z, w):
    a, b = sections
    aw, bw = a.subs(z, w), b.subs(z, w)
    if kind == "node":
        numerator = (z-t*w)*a*bw+(t*z-w)*aw*b
        degrees = (0, 1, 2, 3)
    else:
        numerator = (z-w)*(a*bw+aw*b)/t+z*w*(a*bw-aw*b)
        degrees = (0, 2, 3, 4)
    quotient, remainder = sp.div(numerator.expand(), z-w, z)
    assert remainder == 0
    polynomial = sp.Poly(quotient, z, w)
    matrix = sp.Matrix([[polynomial.coeff_monomial(z**i*w**j) for j in degrees] for i in degrees])
    basis = [1+z**4, z, z**2, z**3] if kind == "node" else [1, z**2, z**3, z**4]
    assert sp.expand(quotient-sum(matrix[i, j]*basis[i]*sp.sympify(basis[j]).subs(z, w)
                                   for i in range(4) for j in range(4))) == 0
    return matrix


def check(kind):
    z, w = sp.symbols("z w")
    rng = Random(1582)
    transform = sp.Matrix([[rng.randrange(-5, 6) for _ in range(4)] for _ in range(4)])
    assert transform.det()
    t = sp.Integer(3)
    if kind == "node":
        left, right = ([1+t*z**4, z, z**2, z**3], [1+z**4/t, z, z**2, z**3])
    else:
        left, right = ([1+t*z, z**2, z**3, z**4], [1-t*z, z**2, z**3, z**4])

    def evaluate(vector):
        a = sum(f*c for f, c in zip(left, vector))
        b = sum(f*c for f, c in zip(right, transform*vector))
        return normal_form(kind, (a, b), t, z, w), a, b

    ratios = []
    for _ in range(4):
        vector = sp.Matrix([rng.randrange(1, 8) for _ in range(4)])
        matrix, a, b = evaluate(vector)
        resultant = sp.resultant(a, b, z)
        if resultant:
            ratios.append(matrix.det()/resultant)
    assert len(ratios) >= 2 and ratios[0] and len(set(ratios)) == 1
    units = [sp.eye(4)[:, i] for i in range(4)]
    diagonals = [evaluate(e)[0] for e in units]
    pairs = monomials(2)
    components = [diagonals[i] if i == j else evaluate(units[i]+units[j])[0]-diagonals[i]-diagonals[j]
                  for i, j in pairs]
    change = sp.Matrix([[block[i, j] for block in components] for i, j in pairs])
    assert change.rank() == 10
    coords = sp.symbols("c:4")
    forms = change*sp.Matrix([coords[i]*coords[j] for i, j in pairs])
    n = sp.Matrix(4, 4, lambda i, j: forms[pairs.index(tuple(sorted((i, j))))])
    index = [tuple(word.count(i) for i in range(4)) for word in monomials(4)]
    rows = []
    for i in combinations(range(4), 2):
        for j in combinations(range(4), 2):
            poly = sp.Poly(n.extract(i, j).det(), *coords)
            row = []
            for m in index:
                a, b = poly.coeff_monomial(m).as_numer_denom()
                row.append(int(a)*pow(int(b), -1, P) % P)
            rows.append(row)
    rank = nmod_mat(rows, P).rank()
    assert rank > 15

    h = [sp.Integer(rng.randrange(-7, 8)) for _ in range(9)]
    if kind == "node":
        base = sp.Matrix(3, 3, lambda i, j: h[i+j+2])
        u, v = sp.Matrix(h[1:4]), sp.Matrix(h[5:8])
        coefficient = h[4]*base.det()-(u.T*base.adjugate()*v)[0]
        hankel = sp.Matrix(4, 4, lambda i, j: h[i+j+1])
        assert coefficient == -hankel.det()
    else:
        base = sp.Matrix(3, 3, lambda i, j: h[i+j+4])
        v = sp.Matrix(h[3:6])
        coefficient = -h[2]*base.det()+(v.T*base.adjugate()*v)[0]
        hankel = sp.Matrix(4, 4, lambda i, j: h[i+j+2])
        assert coefficient == -hankel.det()
    print(f"{kind}: resultant, Hankel boundary and {rank} independent quartic minors verified.")


if __name__ == "__main__":
    check("node")
    check("cusp")
