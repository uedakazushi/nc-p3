"""Check the semipositive degree enumeration and the diagonal quotient dimensions."""
from itertools import combinations_with_replacement, product
from math import comb

import sympy as sp

from verify_conductor_pencil import coefficients, t


def distributions(n):
    for support in combinations_with_replacement(range(n), 4):
        yield tuple(support.count(i) for i in range(n))


def arc_signature(degrees):
    n = len(degrees)
    cyclic = degrees*2
    return tuple(int(total >= 2) + int(total >= 3)
                 for i in range(n) for length in range(1, n)
                 for total in [sum(cyclic[i:i+length])])


def check_arcs():
    for n in range(2, 9):
        groups = {}
        for d in distributions(n):
            groups.setdefault(arc_signature(d), []).append(d)
        survivors = []
        for group in groups.values():
            for i, d in enumerate(group):
                for e in group[:i]:
                    if all(x+y for x, y in zip(d, e)):
                        survivors.append(sorted(zip(d, e)))
        if n == 2:
            assert len(survivors) == 2
            assert all(s in ([(0,1),(4,3)], [(1,0),(3,4)]) for s in survivors)
        elif n == 3:
            assert survivors == [[(0,1),(1,0),(3,3)]]*3
        else:
            assert not survivors


def check_projection():
    b = [1,t,t**2,t**3]
    w = [1+2*t+4*t**2+7*t**6+3*t**7, t+4*t**4+t**7]
    quotient = sp.Matrix.hstack(*coefficients(w,7).T.nullspace()).T
    for conductor in (t, t**2):  # node at 0,infinity or cusp at 0
        reduced = [conductor*t**i for i in range(3)]
        matrix = quotient*coefficients([a*c for a,c in product(reduced,b)],7)
        full = coefficients([t**(i+j) for i,j in product(range(3),range(4))],5)
        assert matrix.rank() == full.rank() == matrix.col_join(full).rank() == 6
    assert 9*2**7 == 1152 > comb(20,3) == 1140


if __name__ == "__main__":
    check_arcs()
    check_projection()
    print("PASS: both binary projections, all semipositive arc signatures, degree 17 bound 1152 > 1140")
