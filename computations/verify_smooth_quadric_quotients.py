"""Exact spin block and normal quadratic transport checks."""

from itertools import product

import sympy as sp


def spin(u):
    return sp.Matrix.vstack(sp.Matrix.hstack(sp.zeros(2), u),
                            sp.Matrix.hstack(u.adjugate(), sp.zeros(2)))


def verify():
    units = [sp.eye(4)[:, i].reshape(2, 2) for i in range(4)]
    vectors = [spin(u) for u in units]
    paths = [sp.eye(4)]
    for expected in (4, 7, 8, 8):
        paths = [a*b for a, b in product(paths, vectors)]
        assert sp.Matrix.hstack(*[w.reshape(16, 1) for w in paths]).rank() == expected
    nil = sp.diag(sp.Matrix([[0, 1], [0, 0]]), sp.zeros(2))
    middle = [x*nil*y for x, y in product(vectors, repeat=2)]
    assert all(not any(w[:2, :]) for w in middle)
    assert sp.Matrix.hstack(*[w.reshape(16, 1) for w in middle]).rank() == 4
    r = sp.symbols('r', nonzero=True)
    forms = [(sp.diag(2, 3), sp.diag(2, 3), sp.diag(2, 3)),
             (sp.eye(2), sp.Matrix([[1, 1], [0, 1]]), sp.Matrix([[1, 1], [0, 1]])),
             (sp.diag(1, -1), sp.diag(r, -r), sp.diag(1/r, -1/r))]
    for p, q, qnext in forms:
        a = sp.diag(p, q)
        b = sp.diag(p, qnext)
        assert p.trace() == q.trace() == qnext.trace()
        assert p.inv().trace() == q.inv().trace()
        for v in vectors:
            image = a*v*b.inv()
            assert image == spin(image[:2, 2:])


if __name__ == '__main__':
    verify()
    print('PASS: spin dimensions, middle pure-block exclusion, semisimple, Jordan and alternating transports')
