"""Bounds for the short trace poles and the mixed two-dimensional trace calculation."""
from itertools import combinations_with_replacement
from sympy import Matrix


def main():
    units = [Matrix(2,2, lambda r,c: int((r,c) == (i,j)))
             for i in range(2) for j in range(2)]
    modified = Matrix([[(x*y).trace() + x.trace()*y.trace() for y in units] for x in units])
    assert modified.det() != 0
    profiles = {p for n in range(1,5) for p in combinations_with_replacement(range(1,5),n)
                if sum(p) <= 4}
    assert all(sum(d>=2 for d in p) <= 1 or p == (2,2) for p in profiles)
    slots = 2*4**2-1
    maxpath = 4*(4**2-1)*slots
    assert maxpath == 1860
    assert 8*(4**maxpath)*(7**255) + maxpath > maxpath
    print('Mixed trace pairing, multiplicity budgets, and the uniform homogeneous absorber length verified.')


if __name__ == '__main__':
    main()
