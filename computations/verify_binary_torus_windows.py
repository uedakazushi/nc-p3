"""Cube layers for arbitrary binary torus weights and the exceptional Euler test."""
from itertools import combinations, combinations_with_replacement, product
from math import prod

import sympy as sp

from verify_torus_euler_profiles import PERMS


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cube_layers():
    shapes = set()
    for n in range(1, 5):
        vertices = list(product((0, 1), repeat=n))
        partial = [v for v in product((-1, 0, 1), repeat=n) if any(v) and (n < 4 or 0 in v)]
        for size in range(1, n+1):
            for points in combinations(vertices, size):
                differences = [tuple(a-b for a, b in zip(w, points[0])) for w in points[1:]]
                if any(all(dot(v, row) == 0 for row in differences) for v in partial):
                    continue
                matrix = sp.Matrix(differences) if differences else sp.zeros(0, n)
                if matrix.rank() != size-1:
                    continue
                assert n == size == 4
                normal, = matrix.nullspace()
                if min(normal) <= 0:
                    continue
                denom = sp.ilcm(*[x.q for x in normal])
                gaps = tuple(int(x*denom) for x in normal)
                shapes.add((tuple(sorted(gaps)), dot(gaps, points[0])))
    assert shapes == {((1,1,1,1), j) for j in (1,2,3)} | {((1,1,1,2), j) for j in (2,3)}


def check_larger_gap():
    gaps = (1, 1, 1, 2)
    # Invert all characters to reduce the top weight to 2.
    t, y = 11, 7
    allowed = {prod(1-y*t**j for j in powers) for powers in combinations_with_replacement(range(6),4)
               if sum(powers) == 8}
    for high in product(range(1,4), repeat=4):
        if dot(high, gaps) != 8:
            continue
        pairs = []
        for i in range(4):
            j = (i+1) % 4
            counts = [0]*4
            for a, b in product((0,1), repeat=2):
                counts[a*gaps[i]+b*gaps[j]] += (high[i] if a else 4-high[i])*(high[j] if b else 4-high[j])
            pairs.append(counts)
        def choices(i):
            for a in range(7):
                for b in range(7-a):
                    row = (a,b,6-a-b)
                    if all(v <= min(pairs[i][j], pairs[i+2][2-j]) for j,v in enumerate(row)):
                        yield row
        for left,right in product(choices(0), choices(1)):
            relations = (left,right,left[::-1],right[::-1])
            rows = []
            for i in range(4):
                prev = (i-1)%4
                terms = (1+y*t**2, -(4-high[i]+high[i]*t**gaps[i]),
                         sum(v*t**j for j,v in enumerate(relations[i])),
                         -((4-high[prev])*t**2+high[prev]*t**(2-gaps[prev])))
                rows.append([terms[(j-i)%4]*(y if j<i else 1) for j in range(4)])
            value = sum(sign*prod(rows[i][perm[i]] for i in range(4)) for perm,sign in PERMS)
            assert value not in allowed


if __name__ == '__main__':
    cube_layers()
    check_larger_gap()
    print('Only four equal binary gaps survive. The exceptional (1,1,1,2) slice fails the Euler identity.')
