"""Rational identities for the six-dimensional symmetrizer lifts in Moore coordinates."""
from itertools import combinations_with_replacement

import sympy as sp
from verify_prescribed_discriminant_seeds import axial, bracket


def moore(p, v):
    return sp.Matrix(3, 3, lambda k, j: p[(j-k) % 3]*v[(-k-j) % 3])


def check():
    a, b, c, r = sp.symbols('a b c r', nonzero=True)
    v = sp.symbols('x y z')
    p = [a, b, c]
    matrix = moore(p, v)
    hesse = (a**3+b**3+c**3)*sp.prod(v)-a*b*c*sum(x**3 for x in v)
    assert sp.expand(matrix.det()-hesse) == 0

    s = [-3*b*c*(a**3-b**3)+r*a/b, 3*a*a*(b**3-c**3)+r*c*c/(b*b), r*c/a]
    multipliers = [c/b, b/a, a/c]
    opposite = [u*m for u, m in zip(s, multipliers)]
    k = [(opposite[(i+1) % 3]-multipliers[i]*s[(i-1) % 3]).cancel() for i in range(3)]
    differences = [c*k[2]/b, b*k[0]/a, a*k[1]/c]
    assert sp.cancel(sum(differences)) == 0
    diagonal = [(2*differences[0]+differences[1])/3,
                (-differences[0]+differences[1])/3, (-differences[0]-2*differences[1])/3]
    q = [sp.diag(*diagonal), sp.zeros(3), sp.zeros(3)]
    for i in range(3):
        q[1][i, (i+1) % 3] = s[i]
        q[2][i, (i-1) % 3] = opposite[i]
    equations = sp.Matrix.hstack(*[matrix.diff(x) for x in v])
    for i, j in combinations_with_replacement(range(3), 2):
        h = sp.zeros(3)
        h[i, j] = h[j, i] = 1
        column = sp.Matrix.vstack(*[axial(h*u-u.T*h) for u in q])
        assert (equations*column).applyfunc(sp.cancel).is_zero_matrix

    doubled = [a*(b**3-c**3), b*(c**3-a**3), c*(a**3-b**3)]
    assert sp.expand(hesse.subs(dict(zip(v, doubled)))) == 0
    assert sp.expand(sum(hesse.diff(x).subs(dict(zip(v, p)))*y for x, y in zip(v, doubled))) == 0
    target = sum((v[i]*bracket(q[(i+1) % 3], q[(i+2) % 3]) for i in range(3)), sp.zeros(3))
    scale = sp.diag(*[s[i]*k[i]/doubled[2*i % 3] for i in range(3)])
    expected = scale*moore(doubled, v)[[0, 2, 1], :]
    assert (target-expected).applyfunc(sp.cancel).is_zero_matrix
    assert sp.cancel(scale.det()) != 0

    symmetry = []
    for i in range(9):
        h = sp.zeros(3)
        h[i] = 1
        symmetry.append(sp.Matrix.vstack(*[axial(h*matrix.diff(x)-(h*matrix.diff(x)).T) for x in v]))
    det = sp.Matrix.hstack(*symmetry).det(method='domain-ge')
    assert sp.expand(det-(a**3-b**3)*(a**3-c**3)*(b**3-c**3)) == 0


if __name__ == '__main__':
    check()
    print("PASS: Moore symmetrizer equations, half-period obstruction and commutator tangent-point identity")
