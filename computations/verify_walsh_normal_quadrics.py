"""Two generic rank-one obstructions and quaternion transport of missing quadratic words."""

from itertools import product

import sympy as sp

from flint_rational_function import rational_function_field
from verify_klein_center_lift import echelon, field_paths
from verify_rank_fifteen_genus_two_unit_line import tensor_windows


def windows(matrix):
    return tensor_windows({
        (i, i ^ b, i ^ b ^ e, i ^ e): value * (-1)**((i >> 1) * ((b ^ e) & 1))
        for i, b, e in product(range(4), range(3), range(3))
        if (value := matrix[b][e])
    })


def rank_one_obstructions():
    field = rational_function_field(("u", "v"))
    u, v = map(field.from_sympy, sp.symbols("u v"))

    def conic(t):
        return (1 + t*t, 1 - t*t, 2*t)
    for negative in (False, True):
        matrix = [[x*y*(-1)**int(negative and (b, e) == (0, 1))
                   for e, y in enumerate(conic(v))] for b, x in enumerate(conic(u))]
        relations = windows(matrix)
        for h, expected in enumerate((12, 8, 8, 8)):
            rows = []
            for pos, side in enumerate(relations[:3]):
                for color, relation in side:
                    for j in range(4):
                        outside = (j, j ^ color ^ h)
                        rows.append({outside[:pos] + letters + outside[pos:]: c
                                     for letters, c in relation.items()})
            assert 64 - len(echelon(rows)) == expected


def reversal_and_orientation():
    for i, j, k in product(range(4), repeat=3):
        b, e, h = i ^ j, j ^ k, i ^ j ^ k
        flip = (b & 1) * (e >> 1)
        first = (e >> 1) * (k & 1)
        second = (b >> 1) * (j & 1) + int(i == 3) + int(j == 3)
        assert (flip + first - (e >> 1) * (h & 1)) % 2 == 0
        assert (flip + second - (b & 1) * (h >> 1)) % 2 == 0
    matrix = [[0, 1, sp.I], [sp.I, 0, sp.I], [-1, 1, 0]]
    basis, image = field_paths(windows([[sp.sympify(v) for v in row] for row in matrix]), sp.S.One)
    for p in range(4):
        for a in range(4):
            left = [image(p, (i, i ^ 3, a)) for i in range(4)]
            right = [image(p, (a, i ^ 3, i)) for i in range(4)]
            columns = right + left
            linear = sp.Matrix(len(basis[p, 3, a ^ 3]), 8, lambda j, i: columns[i].get(j, 0))
            reduced, pivots = linear.rref()
            assert pivots == (0, 1, 2, 3)
            transport = reduced[:4, 4:]
            assert transport.det() == -1
            if a == 0 and p == 0:
                assert transport == sp.diag(1, 1, -1, 1)


def quaternion_equations():
    one = sp.eye(2)
    e1, e2 = sp.diag(1, -1), sp.Matrix([[0, 1], [1, 0]])
    basis = [one, e1, e2, e1*e2]
    q = sp.Matrix(2, 2, sp.symbols("a b c d"))
    qinv = q.adjugate()  # common (invertible) determinant is cleared
    nxt = e1 * qinv * e2
    for a, b in product(range(4), repeat=2):
        ea, eb = basis[a], basis[b]
        assert ea*eb == (-1)**((a >> 1)*(b & 1))*basis[a ^ b]
        arrow_sign = (-1)**(int(a != 0) + (a >> 1)*(b & 1) + (a & 1)*(b >> 1))
        kappa = (-1)**(int(b != 0) + (b >> 1))
        actual = e1 * ea * eb.adjugate() * ea * qinv * e2
        assert actual == arrow_sign * kappa * ea.det() * eb * nxt


if __name__ == "__main__":
    rank_one_obstructions()
    reversal_and_orientation()
    quaternion_equations()
    print("PASS: two rank-one quartic obstructions, orientation and all quadratic normal transports")
