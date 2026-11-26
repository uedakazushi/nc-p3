"""Exact identities for the five-space conic proof, over the rationals."""
from itertools import product

import sympy as sp


def exterior_restriction(left, right):
    rows = []
    for i, j in ((0, 1), (0, 2), (1, 2)):
        form = sp.zeros(3)
        form[i, j], form[j, i] = 1, -1
        restricted = left.T * form * right
        rows.append(list(restricted))
    return sp.Matrix(rows).rank()


def main():
    eye = sp.eye(3)
    plane = eye[:, :2]
    other_plane = eye[:, 1:]
    for left, right in ((eye, eye), (plane, eye), (eye, plane), (plane, other_plane)):
        assert exterior_restriction(left, right) == 3

    c1, c2, d1, d2, h, k12, k21 = sp.symbols("c1 c2 d1 d2 h k12 k21")
    c = sp.Matrix([c1, c2, h])
    d = sp.Matrix([d1, d2, h])
    rank_one = c * d.T
    residual = c1 * d2 - c2 * d1 - h * (k12 - k21)
    correction = rank_one[:2, :2] - h * sp.Matrix([[0, k12], [k21, 0]])
    assert (correction[0, 1] - correction[1, 0]).expand() == residual.expand()
    v1, v2, a, b = sp.symbols("v1 v2 a b")
    v = sp.Matrix([v1, v2, 0])
    tangent = v * d.T + c * v.T
    assert tangent[2, 2] == 0
    assert tangent[:2, 2] == h * v[:2, :]
    assert tangent[2, :2] == h * v[:2, :].T
    shear = sp.Matrix([-b, a]) * sp.Matrix([[c1, c2]])
    assert (a * c1 + b * c2 + shear[0, 1] - shear[1, 0]).expand() == 0
    for i, j in product(range(3), repeat=2):
        assert rank_one[i, 2] * rank_one[2, j] == rank_one[2, 2] * rank_one[i, j]
    print("PASS: quotient exterior forms, conic residual, tangent and symmetric shear")


if __name__ == "__main__":
    main()
