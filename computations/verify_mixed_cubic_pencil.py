"""The five base planes of the cubic Pluecker system for the split symmetric net."""
from itertools import combinations, product

import sympy as sp

from verify_prescribed_discriminant_seeds import bracket, cross

X = sp.symbols('x y z')
V = sp.Matrix(X)
INDICES = list(combinations(range(5), 3))
G = []
for i in range(3):
    m = sp.zeros(3)
    j, k = (i+1) % 3, (i+2) % 3
    m[j, k] = m[k, j] = 1
    G.append(m)
BASIS = G+[sp.diag(1, -1, 0), sp.diag(0, 1, -1)]
A = sp.Matrix([[sp.trace(cross(V)*bracket(g, m))/2 for m in BASIS] for g in G])


def coefficients(polynomial):
    poly = sp.Poly(polynomial, *X)
    return [poly.coeff_monomial(X[0]**i*X[1]**j*X[2]**(3-i-j))
            for i in range(4) for j in range(4-i)]


def verify_base(plane):
    pivots = plane.T.rref()[1]
    plane = plane*plane.extract(pivots, range(3)).inv()
    mixed = A*plane
    assert sp.expand(mixed.det()) == 0
    columns = []
    for row in range(5):
        if row in pivots:
            continue
        for col in range(3):
            variation = mixed.copy()
            variation[:, col] = A[:, row]
            columns.append(coefficients(variation.det()))
    assert sp.Matrix(columns).rank() == 6


def main():
    cubics = [A[:, indices].det() for indices in INDICES]
    assert sp.Matrix([coefficients(p) for p in cubics]).rank() == 6
    planes = [sp.eye(5)[:, :3]]
    verify_base(planes[0])
    for j, k in product((1, -1), repeat=2):
        v = sp.Matrix([1, j, k])
        assert all(p.subs(dict(zip(X, v))) == 0 for p in cubics)
        common_eigenline = sp.Matrix.hstack(*[m*v for m in BASIS])
        plane = sp.Matrix.hstack(*(cross(v)*common_eigenline).nullspace())
        assert plane.shape == (5, 3)
        verify_base(plane)
        planes.append(plane)
    assert sp.Matrix([[plane.extract(rows, range(3)).det() for rows in INDICES]
                      for plane in planes]).rank() == 4
    print("PASS: six independent cubics and five transverse zero planes in Gr(3,5)")


if __name__ == '__main__':
    main()
