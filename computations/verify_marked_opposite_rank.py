#!/usr/bin/env python3
"""Polynomial commutator identity and transversality at one marked zero of the symmetric window."""
from random import Random

import sympy as sp

from verify_finite_square_chart import frame
from verify_reversible_component import local_matrices


CYCLIC = ((1, 2), (2, 0), (0, 1))


def alternating_column(matrix):
    return sp.Matrix([matrix[i, j]-matrix[j, i] for i, j in CYCLIC])


def commutators(matrices, inverse_metric):
    return sp.Matrix.hstack(*(alternating_column(matrices[i]*inverse_metric*matrices[j])
                              for i, j in CYCLIC))


def verify_polynomial_identity():
    p, q, a, b, c = sp.symbols("p q a b c")
    r = -p-q
    weight = sp.diag(p, q, r)
    diagonal = sp.diag(q-r, r-p, p-q)
    net = []
    for scalar, (i, j) in zip((a, b, c), CYCLIC):
        matrix = scalar*diagonal
        matrix[i, j] = matrix[j, i] = 1
        net.append(matrix)
    d = commutators(net, sp.eye(3))
    moments = sp.Matrix(3, 3, lambda i, j: sp.trace(weight*net[i]*net[j]))
    assert (d*moments*d.T-d.det()*weight).applyfunc(sp.expand) == sp.zeros(3)


def verify_marked_open():
    rel, h, _, _, tangent, quadric, _ = local_matrices()
    x, y = frame()
    e = sp.Matrix.vstack(x[0, :], sp.eye(4)[1:, :])
    f = sp.Matrix.vstack(y[0, :], sp.eye(4)[1:, :])
    kron = sp.kronecker_product
    adapted = kron(e, f)*rel
    row = adapted[1:4, :]
    column = adapted[[4, 8, 12], :]
    jacobian = column.col_join(row).T
    inverse = (row*h*row.T).inv()
    matrices = [adapted[4*i+1:4*i+4, :]*h*row.T for i in range(1, 4)]
    assert commutators(matrices, inverse) == sp.diag(-sp.Rational(15, 32), 0, 0)

    def normal_block(direction):
        drelation = direction.reshape(6, 16).T*h.inv()
        dpoint = -jacobian.inv()*(kron(e[0, :], f[0, :])*drelation).T
        de, df = sp.zeros(4), sp.zeros(4)
        de[0, 1:], df[0, 1:] = dpoint[:3, :].T, dpoint[3:, :].T
        delta = (kron(de, f)+kron(e, df))*rel + kron(e, f)*drelation
        assert delta[0, :] == sp.zeros(1, 6)
        drow = delta[1:4, :]
        dinverse = -inverse*(drow*h*row.T + row*h*drow.T)*inverse
        dm = [delta[4*i+1:4*i+4, :]*h*row.T + adapted[4*i+1:4*i+4, :]*h*drow.T
              for i in range(1, 4)]
        full = sp.Matrix.hstack(*(alternating_column(dm[i]*inverse*matrices[j]
                                 + matrices[i]*dinverse*matrices[j] + matrices[i]*inverse*dm[j])
                                 for i, j in CYCLIC))
        return full[1:, 1:]

    random = Random(15)
    directions = [sp.Matrix([random.randrange(-2, 3) for _ in range(tangent.cols)]) for _ in range(2)]
    t = sp.symbols("t")
    line = directions[0] + t*directions[1]
    cone = (line.T*quadric*line)[0].expand()
    first, second = (normal_block(tangent*v) for v in directions)
    minor = (first+t*second).det().expand()
    assert sp.expand(15*cone+928*t*t+704*t+232) == 0
    assert sp.expand(256*minor-112*t*t-184*t+37) == 0
    assert sp.gcd(cone, minor) == 1


if __name__ == "__main__":
    verify_polynomial_identity()
    verify_marked_open()
    print("PASS marked commutator identity in five indeterminates and coprime tangent quadrics over Q")
