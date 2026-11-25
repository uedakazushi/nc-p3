"""Clifford transport on a rank-three quadratic cone over the dual numbers."""

from itertools import product

import sympy as sp

I = sp.eye(2)
ZERO = sp.zeros(2)
NIL = sp.Matrix([[0, 1], [0, 0]])
DIAG = sp.diag(1, -1)
FLIP = sp.diag(1, 1, -1, -1)


def even(constant, linear):
    return constant.row_join(linear).col_join(ZERO.row_join(constant))


def odd(traceless, scalar):
    return even(traceless, scalar*I)*FLIP


def verify():
    vectors = [odd(NIL, 0), odd(NIL.T, 0), odd(DIAG, 0), odd(ZERO, 1)]
    paths = [sp.eye(4)]
    for dim in (4, 7, 8, 8):
        paths = [u*v for u, v in product(paths, vectors)]
        assert sp.Matrix.hstack(*[w.reshape(16, 1) for w in paths]).rank() == dim
    bad = even(ZERO, NIL)
    ideal = [u*bad*v for u, v in product(vectors, repeat=2)]
    assert sp.Matrix.hstack(*[w.reshape(16, 1) for w in ideal]).rank() == 4
    assert all(not any(w[:2, :2]) for w in ideal)
    cases = [(I, ZERO, 0), (I, NIL, 0), (I, DIAG, 0), (sp.diag(2, 3), ZERO, 0),
             (I+NIL, ZERO, 0), (DIAG, ZERO, 0), (DIAG, ZERO, 1)]
    for p, n, d in cases:
        a = even(p, (n+d*I)*p)
        next_a = even(p, (n-d*I)*p)
        assert not a[:2, 2:].trace() and not a.inv()[:2, 2:].trace()
        for u in vectors:
            image = a*u*next_a.inv()*FLIP
            assert image == even(image[:2, :2], image[0, 2]*I)
            assert not image[:2, :2].trace()
        assert bool(((p-p.trace()*I/2)*a[:2, 2:]).trace()) == bool(d)
    k = NIL+DIAG
    translate = even(I, k)
    for u in (DIAG, NIL, NIL.T):
        assert translate*odd(u, 0)*translate.inv() == odd(u, (k*u).trace())


if __name__ == '__main__':
    verify()
    print('PASS: cone Clifford dimensions, nilpotent exclusion, all normal forms and radical translations')
