"""Exact ternary graph tests, including a nonreduced cubic, and normal-line transport."""

from itertools import combinations, product

import sympy as sp


def relation(n, terms):
    matrix = sp.zeros(n)
    for i, j, c in terms:
        matrix[i, j] += c
    return matrix


def overlap(relations):
    n = len(relations)
    indices = tuple(product(range(n), repeat=3))
    left = sp.Matrix([[int(i == a)*r[j, k] for a in range(n) for r in relations]
                      for i, j, k in indices])
    right = sp.Matrix([[r[i, j]*int(k == a) for r in relations for a in range(n)]
                       for i, j, k in indices])
    return [v for v in left.row_join(-right).nullspace()]


def graph_test():
    x = sp.Matrix(sp.symbols('x:3'))
    comm = [relation(3, [(j, i, 1), (i, j, -1)]) for i, j in combinations(range(3), 2)]
    jordan = [r.copy() for r in comm]
    jordan[0][0, 0] = -1
    squares = [relation(3, [(i, i, 1)]) for i in range(3)]
    for rels, expected, good in [(comm, 0, True), (jordan, x[0]**2*x[2], True),
                                 (squares, sp.prod(x), False)]:
        for transposed in (False, True):
            evaluated = sp.Matrix.vstack(*[x.T*(r.T if transposed else r) for r in rels])
            assert sp.expand(evaluated.det()**2-expected**2) == 0
            minors = [evaluated.extract(a, b).det()
                      for a in combinations(range(3), 2) for b in combinations(range(3), 2)]
            all_rank_two = True
            for var in x:
                gb = sp.groebner([m.subs(var, 1) for m in minors], *(v for v in x if v != var))
                all_rank_two &= gb.reduce(1+0*var)[1] == 0
            assert bool(all_rank_two) == good
        if good:
            w, = overlap(rels)
            assert sp.Matrix(3, 3, w[:9]).det() and sp.Matrix(3, 3, w[9:]).det()


def check_transport():
    rels = [relation(4, [(j, i, 1), (i, j, -1)]) for i, j in combinations(range(4), 2)]
    rels[0][0, 0] = -1
    basis = sp.eye(4)
    space = sp.Matrix.hstack(*[r.reshape(16, 1) for r in rels])
    for k in (0, 3):
        t = basis[:, k]
        left = sp.Matrix.hstack(*[(t*basis[:, j].T).reshape(16, 1) for j in range(4)])
        right = sp.Matrix.hstack(*[(basis[:, j]*t.T).reshape(16, 1) for j in range(4)])
        eq = left.row_join(space)
        assert eq.rank() == 10
        g = eq.gauss_jordan_solve(right)[0][:4, :]
        assert g.det() == 1 and g*t == t
        assert space.row_join(sp.kronecker_product(g, g)*space).rank() == 6
        expect = sp.eye(4)
        if k == 0:
            expect[0, 1] = 1
        assert g == expect


if __name__ == '__main__':
    graph_test()
    check_transport()
    print('PASS: plane graph, nonreduced cubic, rank-one exclusion, central and noncentral regular line transports')
