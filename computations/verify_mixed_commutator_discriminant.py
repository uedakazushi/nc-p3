"""Symbolic five-dimensional alternating form and cyclic reconstruction identities."""
from itertools import combinations

import sympy as sp


def comm(a, b):
    return a*b - b*a


def pf(a):
    return a[0, 1]*a[2, 3] - a[0, 2]*a[1, 3] + a[0, 3]*a[1, 2]


def verify():
    z = sp.symbols("x y z")
    x, y, t = z
    v = sp.Matrix(z)
    k = sp.Matrix([[0, -t, y], [t, 0, -x], [-y, x, 0]])
    base = [sp.diag(1, -1, 0), sp.diag(0, 1, -1)]
    for i, j in combinations(range(3), 2):
        element = sp.zeros(3)
        element[i, j] = element[j, i] = 1
        base.append(element)
    omega = sp.Matrix(5, 5, lambda i, j: (k*comm(base[i], base[j])).trace())
    vector = sp.zeros(3)
    for i, element in enumerate(base):
        indices = [j for j in range(5) if j != i]
        vector += (-1)**i * pf(omega.extract(indices, indices))*element
    expected = v*v.T - v.dot(v)*sp.eye(3)/3
    ratio = sp.cancel(vector[0, 1]/expected[0, 1])
    assert ratio.is_Rational and ratio != 0
    assert (vector-ratio*expected).applyfunc(sp.expand) == sp.zeros(3)
    entries = sp.symbols("a:6")
    mixed = sp.zeros(4)
    for (i, j), coefficient in zip(combinations(range(4), 2), entries):
        mixed[i, j], mixed[j, i] = coefficient, -coefficient
    assert sp.expand(mixed.extract([0, 1, 2], [0, 1, 3]).det()-mixed[0, 1]*pf(mixed)) == 0
    # Exercise the reconstruction at a nonsemisimple cyclic X.
    X = sp.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    Y = sp.Matrix([[1, 2, 0], [0, 0, 3], [1, 0, -1]])
    Z = sp.Matrix([[0, 3, 2], [4, 0, -1], [1, 0, 0]])
    a, b, c, d = sp.symbols("a b c d")
    C, D = comm(X, Y), comm(X, Z)
    corrected = comm(Y+a*X+b*X**2, Z+c*X+d*X**2)-comm(Y, Z)
    linear = a*D+b*(X*D+D*X)-c*C-d*(X*C+C*X)
    assert (corrected-linear).applyfunc(sp.expand) == sp.zeros(3)


if __name__ == "__main__":
    verify()
    print("PASS: mixed determinant, all five Pfaffians, and cyclic reconstruction at a Jordan block")
