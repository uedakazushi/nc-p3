"""Recover the marked rank-ten window and its cubic ideal from four matrix triples."""
from itertools import product
from math import comb

from flint import fmpz_mat
from sympy import Matrix, eye

from verify_twenty_point_companion import gram_half


def main():
    gram = gram_half(1).T*gram_half(-1)

    def w(a, b, c, d):
        return gram[4*a+b, 4*c+d]

    def triple(position):
        arrays = []
        for i in range(1, 4):
            if position == 0:
                value = lambda a, b: w(0, a+1, i, b+1)
            elif position == 1:
                value = lambda a, b: w(a+1, 0, b+1, i)
            elif position == 2:
                value = lambda a, b: w(i, a+1, 0, b+1)
            else:
                value = lambda a, b: w(a+1, i, b+1, 0)
            arrays.append(Matrix(3, 3, value))
        return arrays

    f, p, g, q = (triple(i) for i in range(4))
    assert Matrix(3, 3, lambda b, d: w(0, b+1, 0, d+1)) == eye(3)
    assert Matrix(3, 3, lambda a, c: w(a+1, 0, c+1, 0)) == eye(3)
    assert all(m.trace() == 0 for t in (f, p, g, q) for m in t)
    for a, b, c, d in product(range(3), repeat=4):
        entry = w(a+1, b+1, c+1, d+1)
        assert entry == (g[a]*f[c])[b, d] + (q[b]*p[d])[a, c]
        assert entry == (f[c]*g[a])[b, d] + (p[d]*q[b])[a, c]
    # Obtain the actual eight cubic relations at each start by reducing a single mixed letter.
    relations = []
    for i in range(4):
        def rotate(a, b, c, d):
            indices = (a, b, c, d)
            return w(*(indices[-i:]+indices[:-i]))
        middle = [Matrix(3, 3, lambda a, c: rotate(a+1, b+1, c+1, 0)) for b in range(3)]
        rows = []
        for k in range(9):
            unit = Matrix(3, 3, lambda a, c: int(3*a+c == k))
            comm = [m*unit-unit*m for m in middle]
            rows.append([int(comm[b][a, c]) for a, b, c in product(range(3), repeat=3)])
        assert fmpz_mat(rows).rank() == 8
        relations.append(rows)
    for i, n in product(range(4), range(3, 6)):
        words = list(product(range(3), repeat=n))
        rows = []
        for offset in range(n-2):
            for outer in product(range(3), repeat=n-3):
                for rel in relations[(i+offset) % 4]:
                    rows.append([rel[9*w[offset]+3*w[offset+1]+w[offset+2]]
                                 if w[:offset]+w[offset+3:] == outer else 0 for w in words])
        assert 3**n-fmpz_mat(rows).rank() == comb(n+3, 3)-1
    print("PASS: two Schur complements, four eight-dimensional commutator spaces and cubic Hilbert through five")


if __name__ == "__main__":
    main()
