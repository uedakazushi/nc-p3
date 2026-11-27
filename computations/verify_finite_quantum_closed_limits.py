"""Exact quadratic obstruction for the polynomial chamber of finite point graphs."""
from itertools import combinations, product

import sympy as sp

K = range(4)
WEIGHTS = (0, 3, 5, 9)
PIVOTS = [4*i+j for i,j in combinations(K, 2)]


def words(degree):
    return [q for q in product(K, repeat=4) if sum(WEIGHTS[v] for v in q) == 17+degree]


def cut(tensor, shift):
    def entry(i, j):
        q = (i//4, i%4, j//4, j%4)
        return tensor.get(q[-1:]+q[:-1] if shift else q, 0)
    return sp.Matrix(16, 16, entry)


def main():
    volume = {q: sp.LeviCivita(*q) for q in words(0)}
    cuts = []
    for shift in (0,1):
        m = cut(volume, shift)
        left = sp.Matrix.hstack(*m.T.nullspace()).T
        right = sp.Matrix.hstack(*m.nullspace())
        cuts.append((shift, left, right, m[PIVOTS, PIVOTS].inv()))

    def linear(tensor):
        return sp.Matrix.vstack(*[(l*cut(tensor, sh)*r).reshape(100,1) for sh,l,r,_ in cuts])

    def coefficients(h):
        return sp.Matrix.hstack(*[linear({q: 1}) for q in words(h)])

    first = coefficients(1)
    assert first.shape == (200,14) and first.rank() == 11
    a,b,c = sp.symbols('a b c')
    delta = dict(zip(words(1), sp.Matrix.hstack(*first.nullspace())*sp.Matrix([a,b,c])))
    rhs = []
    for sh,l,r,inv in cuts:
        d = cut(delta, sh)
        rhs.append((l*d[:, PIVOTS]*inv*d[PIVOTS, :]*r).reshape(100,1))
    rhs = sp.Matrix.vstack(*rhs)
    second = coefficients(2)
    columns = list(second.rref()[1])
    assert second.shape == (200,12) and len(columns) == 7
    image = second[:,columns]
    rows = list(image.T.rref()[1])
    remainder = (rhs-image*image[rows,:].inv()*rhs[rows,:]).applyfunc(sp.factor)
    assert set(remainder) == {0, -4*a*b}

    m, d, inv = cut(volume,0), cut(delta,0), cuts[0][3]
    edge = (d[:,PIVOTS]-m[:,PIVOTS]*inv*d[PIVOTS,PIVOTS])*inv
    assert [edge[5,PIVOTS.index(2)], edge[12,PIVOTS.index(6)], edge[10,PIVOTS.index(3)]] == [c,-b,-a]
    print('The polynomial cut obstruction is 4BC. Both nonzero coefficients are necessary for a finite planar graph.')


if __name__ == '__main__':
    main()
