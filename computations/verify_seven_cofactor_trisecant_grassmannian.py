"""Bezoutian, Hankel, and common-quadric checks over Q in binary-form coordinates."""

from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def rank(matrix):
    return DomainMatrix.from_Matrix(matrix).rank()


def bezout(f, g, n, z, w):
    kernel = sp.div(f*g.subs(z, w) - f.subs(z, w)*g, z-w, z)[0]
    poly = sp.Poly(kernel, z, w)
    return sp.Matrix(n, n, lambda i, j: poly.coeff_monomial(z**i*w**j))


def verify_bezout():
    z, w = sp.symbols("z w")
    f = sp.prod(z-a for a in range(4))
    b = bezout(f, 1+z**4, 4, z, w)
    inverse = b.inv()
    assert all(inverse[i+1, j] == inverse[i, j+1] for i, j in product(range(3), repeat=2))
    h = 3+2*z-z**2+5*z**3
    col = sp.Matrix([h.coeff(z, i) for i in range(4)])
    assert bezout(z*h, h, 4, z, w) == col*col.T


def verify_overlap():
    T = sp.Matrix([[-1, 5, -5, 1, -4], [1, 2, 1, 5, 4], [5, 0, -4, -1, 2],
                   [-3, -1, -1, 4, -5], [1, -2, -5, -2, 5]])
    first = [sp.eye(5)[:, :4], sp.eye(5)[:, 1:]]
    second = [sp.eye(6)[:, :5], sp.eye(6)[:, 1:]]

    def relations(transform):
        return [v*transform*u for u in first for v in second]

    r, reverse = relations(T), relations(T.inv())
    r_next = [reverse[i] for i in (0, 2, 1, 3)]
    rows = []
    for x, middle, y in product(range(4), repeat=3):
        left = [r[middle][i, x] if y == j else 0 for i, j in product(range(6), range(4))]
        right = [r_next[middle][i, y] if x == j else 0 for i, j in product(range(6), range(4))]
        rows.append(left+right)
    degrees = sp.Matrix(rows)
    assert rank(degrees[:, :24]) == rank(degrees[:, 24:]) == 24
    assert rank(degrees) == 45
    hankel = [sp.Matrix(5, 5, lambda i, j: int(i+j == k)) for k in range(9)]
    both = hankel + [T.T*h*T for h in hankel]
    assert rank(sp.Matrix.hstack(*(m.reshape(25, 1) for m in both))) == 15


if __name__ == "__main__":
    verify_bezout()
    verify_overlap()
    print("Verified the Bezout identities and the generic cubic overlap by both exact rank calculations.")
