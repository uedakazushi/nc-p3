"""Exact restriction bounds, Artin reduction counts, and a bisymmetric four-slot evaluation."""
from itertools import combinations, combinations_with_replacement

import sympy as sp

from verify_finite_square_chart import frame, products


def main():
    for length in range(1,10):
        least = 1 + (length >= 2) + (length >= 4) + (length >= 8)
        assert min(length,least+1) * 2 > length
    quadratics = list(combinations_with_replacement(range(3),2))
    cubics = list(combinations_with_replacement(range(3),3))
    for pair in combinations(quadratics,2):
        assert sum(all(factor in pair for factor in combinations(word,2)) for word in cubics) <= 2
    u,v = frame()
    adjacent = sp.Matrix.hstack(*products(u,v).columnspace())
    four = products(adjacent,adjacent,True)
    trace, = four.T.nullspace()
    assert all(trace) and (adjacent.T*sp.diag(*trace)*adjacent).is_zero_matrix
    assert four.rank()==19
    for m in range(1,20):
        positions = list(range(0,20,2))+list(range(1,20,2))
        a,b = adjacent[positions[:m],:].rank(),adjacent[positions[m:],:].rank()
        assert a==m-10+b and 2*a>m
    print('Half-length restrictions, residual duality and the unique nondegenerate four-slot trace verified.')


if __name__ == '__main__':
    main()
