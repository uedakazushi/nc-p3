"""Check the coordinate frame, the two core determinants and their pivots over Q."""
from itertools import product

import sympy as sp

# Terms of h indexed by powers of x1, t1, and by the odd coordinate.
H = ((1, 1, 0, 1), (1, 0, 1, -1), (0, 0, 2, 1))


def flattenings(B, P, Q, D):
    first, second = sp.zeros(16), sp.zeros(16)

    def add(i, v, j, u, value):
        first[4*i+v, 4*j+u] += value
        second[4*v+j, 4*u+i] += value

    for i, j in product(range(4), repeat=2):
        add(i, 3, j, 3, D[i, j])
    for i, j in product(range(3), repeat=2):
        for k, l, v, sign in H:
            add(i+l, 3, j+k, v, sign*P[i, j])
            add(i+k, v, j+l, 3, sign*Q[i, j])
    for i, j in product(range(2), repeat=2):
        for k, l, v, sign in H:
            for m, n, u, other in H:
                add(i+k+n, v, j+l+m, u, sign*other*B[i, j])
    return first, second


def verify():
    B = sp.Matrix(2, 2, sp.symbols('a:d'))
    # Zero-based row and column numbers in the two cyclic orderings.
    pivots = ((2, 4, 6, 8, 10, 12), (1, 2, 3, 8, 9, 10))
    for matrix, rows in zip(flattenings(B, sp.zeros(3), sp.zeros(3), sp.zeros(4)), pivots):
        det = sp.factor(matrix.extract(rows, rows).det())
        assert sp.simplify(det**2 - B.det()**6) == 0
        assert matrix.subs(dict(zip(B, [1, 2, 3, 7]))).rank() == 6

    # In the three odd coordinates count all choices of a pair relation
    # annihilating a moving edge and the two endpoint conditions.
    for even_first in (True, False):
        cols = list(product(range(4), range(3)))
        moving = [[int((v in (1, 2) and i+v-1 == degree) if even_first
                       else (v in (0, 1) and i+v == degree))
                   for i, v in cols] for degree in range(5)]
        fixed = ([(0, 0), (0, 1)] if even_first else [(3, 1), (3, 2)])
        rows = moving + [[int(c == p) for c in cols] for p in fixed]
        assert sp.Matrix(rows).rank() == 6
    assert sum((d+1)**2 for d in (1, 2, 2, 3)) == 38
    print("PASS: two core determinants, both relation spaces and 38 coefficients")


if __name__ == '__main__':
    verify()
