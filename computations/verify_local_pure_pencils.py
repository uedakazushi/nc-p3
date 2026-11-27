"""Check the singular-direction lengths of the Kronecker plus regular blocks."""
from itertools import combinations
from functools import reduce

import sympy as sp


def regular(parts, t):
    return sp.diag(*(sp.Matrix(size, size, lambda i, j:
                         (t-root) if i == j else int(j == i+1))
                     for root, size in enumerate(parts)))


def verify():
    t = sp.symbols('t')
    for e in range(4):
        n = 3-e
        types = {0: [()], 1: [(1,)], 2: [(2,), (1, 1)],
                 3: [(3,), (1, 2), (1, 1, 1)]}[n]
        for parts in types:
            singular = sp.Matrix(e+1, e, lambda i, j: int(i == j) + t*int(i == j+1))
            block = regular(parts, t) if parts else sp.zeros(0)
            matrix = sp.diag(singular, block) if e else block
            minors = [matrix.extract(rows, range(3)).det() for rows in combinations(range(matrix.rows), 3)]
            gcd = sp.Poly(reduce(sp.gcd, minors), t).monic()
            assert gcd.degree() == n
            for root, size in enumerate(parts):
                assert matrix.subs(t, root).rank() == 2
                assert gcd.rem(sp.Poly((t-root)**size, t)).is_zero
    print("PASS: all finite pure center block sizes, divisor lengths and exceptional ranks")


if __name__ == '__main__':
    verify()
