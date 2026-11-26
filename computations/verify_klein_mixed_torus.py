"""Diagonal Walsh normalization and its AS cofactor on the homogeneous mixed torus."""
from itertools import product

import sympy as sp
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


def net(plus, minus):
    result = []
    for i in range(3):
        j, k = (i+1) % 3, (i+2) % 3
        m = sp.zeros(3)
        m[j, k], m[k, j] = plus[i], minus[i]
        result.append(m)
    return result


def main():
    u, a, b = sp.symbols('u a b', nonzero=True)
    v = [a, b, 1/(u**3*a*b)]
    r = [1/(u-e) for e in v]
    c = [1, 0, u*v[1]*r[1]/r[0]]
    c[1] = c[2]*u*v[0]*r[0]/r[2]
    h = -1/u
    q, p = net([u]*3, [1]*3), net([e*w for e, w in zip(v, r)], r)
    g = net(c, [h*c[(i+1) % 3] for i in range(3)])
    f = net([-r[(i-1) % 3]/c[(i-1) % 3] for i in range(3)], [-r[i]/(h*c[i]) for i in range(3)])
    def coefficient(i, j, k, l):
        if i == k == 0:
            return int(j == l)
        if j == l == 0:
            return int(i == k)
        if 0 in (i, j, k, l) and (i, j, k, l).count(0) >= 2:
            return 0
        if not i:
            return f[k-1][j-1, l-1]
        if not j:
            return p[l-1][i-1, k-1]
        if not k:
            return g[i-1][j-1, l-1]
        if not l:
            return q[j-1][i-1, k-1]
        return (g[i-1]*f[k-1])[j-1, l-1]+(q[j-1]*p[l-1])[i-1, k-1]

    for i, j, k, l in product(range(3), repeat=4):
        left = (g[i]*f[k]-f[k]*g[i])[j, l]
        right = (q[j]*p[l]-p[l]*q[j])[i, k]
        assert sp.cancel(left+right) == 0
    rows, values = [], []
    for i in range(4):
        for b0 in range(1, 4):
            for e in range(1, 4):
                indices = (i, i^b0, i^b0^e, i^e)
                row = [0]*25
                for pos, label in enumerate(indices):
                    row[4*pos+label] = 1
                row[16+3*(b0-1)+e-1] = -1
                rows.append(row)
                sign = (-1)**((i >> 1)*((b0^e) & 1))
                values.append(sp.factor(sp.Integer(sign)/coefficient(*indices)))
    matrix = sp.Matrix(rows)
    smith, left, right = [m.to_Matrix() for m in smith_normal_decomp(DomainMatrix.from_Matrix(matrix))]
    assert left*matrix*right == smith
    assert list(smith.diagonal()) == [1]*14+[2]*5+[4]*2+[0]*4
    def evaluate(row):
        return sp.factor(sp.prod(val**int(e) for e, val in zip(row, values)))
    assert all(evaluate(left.row(i)) == 1 for i in range(21, 36))
    wanted = sp.zeros(1, 25)
    for i, j, exponent in [(0, 0, 2), (1, 1, 2), (0, 1, -2), (1, 0, -2)]:
        wanted[16+3*i+j] = exponent
    certificate = sp.zeros(1, 36)
    for index, power in [(9, -1), (10, 1), (12, 1), (13, -1),
                         (18, -1), (19, 1), (21, 1), (22, -1)]:
        certificate[index] = power
    assert certificate*matrix == wanted
    ratio = evaluate(certificate)
    assert sp.factor(ratio-1-(u**6-1)*(a*b*u**4-1)/(u**4*(u-a)*(u-b))) == 0
    print('PASS: 81 mixed equations, 15 integral Walsh normalization relations and the u^6-1 cofactor')


if __name__ == '__main__':
    main()
