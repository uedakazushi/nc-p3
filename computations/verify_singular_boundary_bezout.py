"""Boundary matrices and residual polynomial eigenvectors for the two singular elliptic quartics."""

from itertools import combinations
from random import Random

import sympy as sp
from flint import nmod_mat

from verify_opposite_root_triple_octics import monomials
from verify_single_edge_instanton_extension import P


def divide(expression, z, w):
    quotient, remainder = sp.div(sp.expand(expression), z-w, z)
    assert remainder == 0
    return quotient


def boundary(a, b, kind, z, w):
    a0, a1, _, a3 = a
    az = sum(c*z**i for i, c in enumerate(a))
    bz = sum(c*z**i for i, c in enumerate(b))
    aw, bw = az.subs(z, w), bz.subs(z, w)
    bezout = divide(az*bw-aw*bz, z, w)
    if kind == "node":
        particular = (z+w)*az*aw+divide(a3*(z*az-w*aw)+a0*z*w*(z**4*aw-w**4*az), z, w)
        poly = sp.Poly(particular+z*w*bezout, z, w)
        basis = [1+z**4, z, z**2, z**3]
        degrees = [0, 1, 2, 3]
    else:
        particular = az*aw-a1*(z*az+w*aw)-a1*z*w*divide(az-aw, z, w)
        poly = sp.Poly(particular+(z*w)**2*bezout, z, w)
        basis = [1, z**2, z**3, z**4]
        degrees = [0, 2, 3, 4]
    matrix = sp.Matrix([[poly.coeff_monomial(z**i*w**j) for j in degrees] for i in degrees])
    assert matrix == matrix.T
    assert sp.expand(sum(matrix[i, j]*basis[i]*sp.sympify(basis[j]).subs(z, w)
                         for i in range(4) for j in range(4))-poly.as_expr()) == 0
    return matrix


def check(kind):
    z, w = sp.symbols("z w")
    a = sp.Matrix([2, -1, 1, 3])
    az = sum(c*z**i for i, c in enumerate(a))
    if kind == "node":
        u = a[3]/a[0]
        b = u*(az-a[0])/z-z*(az-a[3]*z**3)/u
    else:
        u = a[1]/a[0]
        b = (a[1]+(-u+u**2*z)*az)/z**2
    b = sp.Poly(sp.cancel(b+2*az), z)
    assert boundary(a, [b.nth(i) for i in range(4)], kind, z, w).rank() == 1

    rng = Random(584)
    transform = sp.Matrix(4, 4, lambda i, j: rng.randrange(-5, 6))
    units = [sp.eye(4)[:, i] for i in range(4)]
    diagonal = [boundary(v, transform*v, kind, z, w) for v in units]
    pairs = monomials(2)
    columns = [diagonal[i] if i == j else boundary(units[i]+units[j], transform*(units[i]+units[j]), kind, z, w)
               -diagonal[i]-diagonal[j] for i, j in pairs]
    change = sp.Matrix([[m[i, j] for m in columns] for i, j in pairs])
    assert change.rank() == 10
    x = sp.symbols("x:4")
    forms = change*sp.Matrix([x[i]*x[j] for i, j in pairs])
    matrix = sp.Matrix(4, 4, lambda i, j: forms[pairs.index(tuple(sorted((i, j))))])
    rows = []
    for i in combinations(range(4), 2):
        for j in combinations(range(4), 2):
            minor = sp.Poly(matrix.extract(i, j).det(), *x)
            rows.append([int(minor.coeff_monomial(sp.prod(x[k] for k in word))) % P for word in monomials(4)])
    assert nmod_mat(rows, P).rank() == 19


if __name__ == "__main__":
    check("node")
    check("cusp")
    print("Verified both boundary matrices, their eigenvector equations, and nineteen independent minors.")
