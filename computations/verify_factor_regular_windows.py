"""Exact relation/overlap ranks for transverse and common-radical regular windows."""
from itertools import product

from sympy import Matrix, QQ, diag, eye
from sympy.polys.matrices import DomainMatrix


def relations(radical, associative):
    basis = eye(4).columnspace()

    def ternary(x, y):
        if associative:
            return x.multiply_elementwise(y)
        return Matrix([x.dot(y), x[0]*y[1]+x[1]*y[0], x[0]*y[2]+x[2]*y[0]])

    def value(x, y):
        u, v = x[2:, :], y[2:, :]
        if radical:
            q = diag(1, 2)
            mixed = (x[0]*v+y[0]*u).col_join(x[1]*q*v+y[1]*q.inv()*u)
            height = x[0]*y[0]+x[1]*y[1]+u[0]*v[1]-u[1]*v[0]
            return mixed.col_join(Matrix([height, u.dot(v)]))
        if associative:
            # The first coordinate represents the unit (1,1,1).
            unit = Matrix([1, 1, 1])
        else:
            unit = Matrix([1, 0, 0])
        def vector(s, v):
            return s*unit+Matrix([0, v[0], v[1]])
        a = ternary(vector(x[0], u), vector(y[0], v))
        b = ternary(vector(x[1], u), vector(y[1], diag(1, 2)*v))
        return a.col_join(b)

    first = Matrix.hstack(*(value(x, y) for x, y in product(basis, repeat=2)))
    return [first, first[:, [4*j+i for i, j in product(range(4), repeat=2)]]]


def nonradical_units():
    # In each ternary cubic the two complementary polars have kernels r,e.
    # Both marked units are independent of r. A binary change separates J_U,J_V.
    g = diag(1, 2)
    def projected(unit_x, unit_y, a, b):
        first = unit_x*unit_y+a[1]*b[1]
        second = (unit_x*b[1]+a[1]*unit_y+a[0]*b[1]+a[1]*b[0])
        return Matrix([first, second])

    def value(x, y, correction):
        u, v = x[2:, :], y[2:, :]
        u2, v2 = g.inv().T*u, g*v
        height = (x[0]*v[0]+u[0]*y[0]+x[1]*v2[0]+u2[0]*y[1]
                  +correction*u[1]*v[1])
        return (projected(x[0], y[0], u, v).col_join(projected(x[1], y[1], u2, v2))
                .col_join(Matrix([height, u.dot(v)])))

    h = diag(eye(4), Matrix([[0, 1], [1, 0]]))
    basis = eye(4).columnspace()
    def matrix(c):
        return Matrix.hstack(*(value(x, y, c) for x, y in product(basis, repeat=2)))
    def defect(m, i, j, k, l):
        def pair(i, j, k, l):
            return (m[:, 4*i+j].T*h*m[:, 4*k+l])[0]
        return pair(i, j, k, l)-pair(i, l, k, j)
    v0 = defect(matrix(0), 2, 2, 3, 3)
    v1 = defect(matrix(1), 2, 2, 3, 3)
    first = matrix(v0/(v0-v1))
    assert all(defect(first, i, j, k, l) == 0 for i, j, k, l in product(range(4), repeat=4))
    assert first[:, 0] != first[:, 5] and first.rank() == 6
    return [first, first[:, [4*j+i for i, j in product(range(4), repeat=2)]]]


def dimension(matrices, n, start):
    rows = {}
    for pos in range(n-1):
        r = matrices[(start+pos) % 2]
        for prefix in range(4**pos):
            stride = 4**(n-pos-2)
            for tail, rel in product(range(stride), range(6)):
                row = {(prefix*16+i)*stride+tail: QQ.convert(r[rel, i])
                       for i in range(16) if r[rel, i]}
                rows[len(rows)] = row
    return 4**n-DomainMatrix(rows, (len(rows), 4**n), QQ).rank()


if __name__ == "__main__":
    for name, radical, associative in [('transverse', False, False),
                                      ('radical', True, False),
                                      ('associative', False, True)]:
        r = relations(radical, associative)
        dims = [[dimension(r, n, p) for n in (3, 4)] for p in range(2)]
        if not associative:
            assert dims == [[20, 35], [20, 35]]
        else:
            assert Matrix(4, 4, list(r[0][0, :])).rank() == 1
        print(name, dims)
    r = nonradical_units()
    dims = [[dimension(r, n, p) for n in (3, 4)] for p in range(2)]
    assert dims == [[20, 35], [20, 35]]
    print("root-line with two nonradical units", dims)
