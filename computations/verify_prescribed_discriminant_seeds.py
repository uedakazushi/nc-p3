"""Exact finite-edge certificates at the three nondiagonalizable trace normals."""
from itertools import combinations

import sympy as sp


def bracket(x, y):
    return x*y-y*x


def cross(x):
    a, b, c = x
    return sp.Matrix([[0, -c, b], [c, 0, -a], [-b, a, 0]])


def axial(m):
    return sp.Matrix([m[2, 1], m[0, 2], m[1, 0]])


def self_adjoint(a, b, c, d, e):
    return sp.Matrix([[a, b, c], [d, -2*a, b], [e, d, a]])


def jacobian():
    q = [-cross(sp.eye(3)[:, i]) for i in range(3)]
    basis = []
    for i in range(8):
        x = sp.zeros(3)
        x[i] = 1
        x[2, 2] = -x.trace()
        basis.append(x)
    columns = []
    for pos in range(3):
        for x in basis:
            column = []
            for i, j in [(1, 2), (2, 0), (0, 1)]:
                v = sp.zeros(3)
                if pos == i:
                    v += bracket(x, q[j])
                if pos == j:
                    v += bracket(q[i], x)
                column.extend(list(v)[:8])
            columns.append(column)
    assert sp.Matrix(columns).T.det() == 6912


def certify(normal, rows):
    metric = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    g = [self_adjoint(*row) for row in rows]
    assert all((normal*x).trace() == 0 for x in g)
    cols = []
    for i, j in [(1, 2), (2, 0), (0, 1)]:
        m = sp.Matrix(3, 3, lambda a, c: (metric*bracket(g[a], g[c]))[i, j])
        assert m == -m.T
        cols.append(axial(m))
    target = sp.Matrix.hstack(*cols)
    assert target.det() != 0
    l = target.adjugate().T
    values = []
    for i in range(5):
        for j in range(5-i):
            v = sp.Matrix([i, j, 1])
            star = i+2*j+3
            m = sp.zeros(3, 1).row_join(star*sp.eye(3)+cross(l*v))
            m = m.col_join(v.row_join(sp.Matrix.hstack(*[-x*v for x in g])))
            values.append([int(m.extract(r, range(4)).det()) % 101
                           for r in combinations(range(6), 4)])
    degree_four = int(sp.Matrix(values).det(method="domain-ge")) % 101
    assert degree_four
    # Infinitesimal projective orthogonal stabilizer of the ordered span.
    gmat = sp.Matrix.hstack(*(x.reshape(9, 1) for x in g))
    left = sp.Matrix.hstack(*gmat.T.nullspace()).T
    skew = [metric*cross(sp.eye(3)[:, i]) for i in range(3)]
    stab = sp.Matrix.hstack(*(sp.Matrix.vstack(*[left*bracket(k, x).reshape(9, 1)
                                                for x in g]) for k in skew))
    assert stab.rank() == 3
    return target.det(), degree_four


if __name__ == "__main__":
    jacobian()
    seeds = [
        ((0, 0, 1, 0, 0), [(-2, -1, 2, -2, 0), (-2, -2, 1, -2, 0), (2, 0, 2, -2, 0)]),
        ((0, 1, 0, 0, 0), [(2, 2, 0, 0, 1), (-2, -1, -2, 0, -1), (2, -2, 1, 0, 0)]),
        ((1, 0, 1, 0, 0), [(2, 12, 6, -6, -12), (0, 0, 12, -12, 0), (-1, 12, 0, -6, 6)]),
    ]
    for normal, rows in seeds:
        print(normal, certify(self_adjoint(*normal), rows))
    print("PASS: all special-normal finite edges, stabilizers, and the skew commutator Jacobian")
