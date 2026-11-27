"""Transport, changes of marking, and quantum coefficients on the semisimple line."""
from itertools import combinations, product

import sympy as sp

from verify_split_inverse_windows import BITS, PAIRS, borders, reciprocal, transport

I = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Z = sp.diag(1, -1)
A = (X*Z, Z, X, I)
B = (Z, X*Z, I, X)
F = (X,)*4


def torus(a, b):
    return tuple(sp.diag(1, v) for v in (a, b, 1/a, 1/b))


def multiply(left, right):
    return tuple(m*n for m, n in zip(left, right))


def gauge(maps, basis):
    return tuple(basis[(i+2)%4].inv().T*m*basis[i].T for i, m in enumerate(maps))


def compare(left, right):
    for m, n in zip(left, right):
        assert all(sp.cancel(m[i]*n[j]-m[j]*n[i]) == 0 for i, j in combinations(range(4), 2))


def line_tensor(s, maps):
    p = {q:(sp.Integer(1) if q[0] == q[1] else s) if q[:2] == q[2:] else 0 for q in BITS}
    return reconstruct(p, maps)


def reconstruct(p, maps):
    targets = [transport(reciprocal(p, i), basis, PAIRS[i])
               for i, basis in enumerate(((I,)*4, maps))]
    scale = sp.cancel(targets[0][BITS[0]]/targets[1][BITS[0]])
    w = {}
    for q in BITS:
        assert sp.cancel(targets[0][q]-scale*targets[1][q]) == 0
        def word(pattern):
            return tuple(b+2*(j%2 ^ c) for j, (b, c) in enumerate(zip(q, pattern)))
        w[word((0,)*4)] = p[q]
        w[word((1,)*4)] = targets[0][q]
        for i, basis in enumerate(((I,)*4, maps)):
            for k, value in enumerate(borders(basis, PAIRS[i], i, q)):
                pat = tuple(int(j not in (i, i+1)) ^ k for j in range(4))
                w[word(pat)] = value*(scale if i == 1 and k == 0 else 1)
    return {q:v for q, v in w.items() if v != 0}


def face_values(w):
    tensor = {tuple(x ^ int(i >= 2) for i, x in enumerate(q)):v for q, v in w.items()}
    assert len(tensor) == 24 and all(len(set(q)) == 4 for q in tensor)
    faces = []
    for i in range(4):
        def value(q):
            return tensor[q[-i:]+q[:-i] if i else q]
        edges = {}
        for a, b in combinations(range(4), 2):
            c, d = [j for j in range(4) if j not in (a, b)]
            ratio = sp.cancel(-value((a,b,c,d))/value((b,a,c,d)))
            assert sp.cancel(ratio*value((b,a,d,c))+value((a,b,d,c))) == 0
            edges[a,b] = ratio
        faces.append(tuple(sp.factor(edges[a,b]*edges[b,c]/edges[a,c])
                           for a,b,c in combinations(range(4), 3)))
    assert len(set(faces)) == 1
    return faces[0]


def verify():
    a, b, c, d, s = sp.symbols('a b c d s', nonzero=True)
    compare(multiply(A, B), F)
    for e, f in product((0, 1), repeat=2):
        flip = multiply(A if e else (I,)*4, B if f else (I,)*4)
        def form(u, v):
            return multiply(torus(u, v), flip)
        maps = form(a, b)
        line_tensor(sp.I, maps)
        compare(gauge(maps, A), form(-1/a, -b))
        compare(gauge(maps, B), form(-a, -1/b))
        compare(gauge(maps, torus(c, d)), form(a*c**(2*(1-e)), b*d**(2*(1-f))))
    mixed = multiply(torus(a, sp.Integer(1)), A)
    rotate = tuple(mixed[(i-1)%4].inv() for i in range(4))
    swap = tuple((X if i%2 == 0 else I)*rotate[PAIRS[1][i]].T*(X if i%2 == 0 else I)
                 for i in range(4))
    compare(swap, multiply(torus(a, -sp.Integer(1)), A))
    line_tensor(s, (I,)*4)
    faces = face_values(line_tensor(s, multiply(torus(a, b), F)))
    print('Indexed quantum triangle values', faces)
    cmat = sp.Matrix([[a, b], [c, d]])
    line_tensor(sp.Integer(1), (cmat, X, cmat.inv().T, X))
    print('Both cuts, four flip types, their gauge actions and the mixed one-shift agree.')


if __name__ == '__main__':
    verify()
