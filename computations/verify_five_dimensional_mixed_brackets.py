"""Both inverse sheets and the central one-period presents of the rank-five bracket."""
from itertools import combinations_with_replacement
import sympy as sp
from verify_prescribed_discriminant_seeds import bracket


def main():
    t = sp.symbols('t', nonzero=True)
    eye = sp.eye(3)
    K, D = [], []
    for i in range(3):
        j, k = (i+1) % 3, (i+2) % 3
        unit = eye[:, j]*eye[:, k].T
        K.append(unit-unit.T)
        D.append(unit+unit.T)
    q = [k+t*d for k, d in zip(K, D)]
    p = [k-t*d for k, d in zip(K, D)]
    g = [d+t*k for k, d in zip(K, D)]
    f = [-d+t*k for k, d in zip(K, D)]
    def same(a, b):
        return (a-b).applyfunc(sp.cancel) == sp.zeros(3)
    for i in range(3):
        for j in range(3):
            c = bracket(q[i], p[j])
            assert same(c, bracket(K[i]-D[i]/t, -t*t*K[j]-t*D[j]))
            assert same(bracket(g[i], f[j]), bracket(p[i], -q[j]))
            if i == j:
                assert same(c, -2*t*bracket(K[i], D[i]))
            else:
                assert same(c, (1-t*t)*bracket(K[i], K[j]))
                assert same(bracket(K[i], D[j]), -bracket(K[j], D[i]))
    # The two character-block maximal minors for the trace-space dimension.
    parts = [[eye[:, i]*eye[:, i].T-eye[:, 2]*eye[:, 2].T for i in range(2)]]
    parts += [[K[i], D[i]] for i in range(3)]
    c = [[bracket(x, y) for y in p] for x in q]
    for char, rows, cols, det in [(0, [5, 7], [0, 2], -4*(t*t-1)**2),
                                   (1, [0, 4, 8, 9], [0, 1, 2, 4], -768*t*t*(t*t-1)**2)]:
        triples = []
        for i in range(3):
            for h in parts[(i+1)^char]:
                triple = [sp.zeros(3) for _ in range(3)]
                triple[i] = h
                triples.append(triple)
        trace = sp.Matrix([[(s[i]*c[k][j]+s[k]*c[i][j]).trace() for s in triples]
                           for i, k in combinations_with_replacement(range(3), 2) for j in range(3)])
        assert sp.factor(trace[rows, cols].det()-det) == 0
    # After G=P,F=-Q the two relations at edge 01 identify the star as a central letter.
    for a in range(3):
        for b in range(3):
            for n in range(3):
                assert q[b][a, n] == -p[a][b, n]
    forward = sp.Matrix([[0, 1, -1], [-1, 0, 1], [1, -1, 0]])
    for sign in (1, -1):
        assert (sp.eye(3)-sign*forward**2).det() == (1+3*sign)**2
    print('PASS: both inverse choices, trace ranks, central presentations and endpoint tori')


if __name__ == '__main__':
    main()
