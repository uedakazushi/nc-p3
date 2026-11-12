#!/usr/bin/env python3
"""Exact witness checker for period p on a strict four-window quadratic presentation.

This certifies a PROVIDED collection of four invertible 4x4 matrices g_i.
Failure does NOT prove non-periodicity: another witness may exist.
Only degree-2 relation-space preservation is checked, not AS regularity.
"""
from fractions import Fraction
import json

D = 4

def f(x):
    return Fraction(x)

def rank(vectors):
    a = [[f(x) for x in v] for v in vectors]
    if not a:
        return 0
    nrows, ncols = len(a), len(a[0])
    row = 0
    for col in range(ncols):
        pivot = next((i for i in range(row, nrows) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        d = a[row][col]
        a[row] = [x / d for x in a[row]]
        for i in range(nrows):
            if i != row and a[i][col]:
                v = a[i][col]
                a[i] = [x - v * y for x, y in zip(a[i], a[row])]
        row += 1
        if row == nrows:
            break
    return row

def identity():
    return [[int(i == j) for j in range(D)] for i in range(D)]

def change_relation(v, g, h):
    output = [Fraction(0) for _ in range(D * D)]
    for a in range(D):
        for b in range(D):
            if v[D * a + b]:
                for c in range(D):
                    for d in range(D):
                        output[D * c + d] += f(v[D * a + b]) * g[c][a] * h[d][b]
    return output

def relation(*terms):
    v = [Fraction(0)] * (D * D)
    for a, b, coefficient in terms:
        v[D * a + b] += Fraction(coefficient)
    return v

def family(lam):
    # Order X,Y,Z,W = 0,1,2,3. Six relations for each initial vertex.
    w = []
    for i in range(4):
        eps = (-1)**i
        w.append([
            relation((0,1,1),(1,0,1),(2,2,-1),(3,3,-1)),
            relation((0,2,1),(2,0,1),(3,1,-eps*lam)),
            relation((0,3,1),(3,0,1),(2,2,-1)),
            relation((1,2,1),(2,1,1),(3,3,-1)),
            relation((1,3,1),(3,1,1)),
            relation((2,3,1),(3,2,1)),
        ])
    return w

def commutative():
    return [[relation((a,b,1),(b,a,-1))
             for a in range(D) for b in range(a+1,D)]
            for _ in range(4)]

def has_period_witness(r, g, p):
    assert len(r) == 4 and len(g) == 4
    for i in range(4):
        if rank(g[i]) != D:
            return False
        source, target = r[i], r[(i + p) % 4]
        if rank(source) != 6 or rank(target) != 6:
            raise ValueError("Not a 6-dimensional quadratic relation subspace")
        transformed = [change_relation(v, g[i], g[(i+1)%4]) for v in source]
        if rank(transformed + target) != 6:
            return False
    return True

def verify():
    g = [identity() for _ in range(4)]
    r = family(Fraction(1))
    assert has_period_witness(r, g, 2)
    assert not has_period_witness(r, g, 1)  # ONLY identity witness fails
    assert has_period_witness(commutative(), g, 1)
    assert family(Fraction(-1)) == [r[(i+1)%4] for i in range(4)]
    return {
      "result": "PASS", "arithmetic": "exact Fraction",
      "four_vertex_relation_ranks": [rank(x) for x in r],
      "A(1)_period2_identity_witness": True,
      "A(1)_period1_identity_witness": False,
      "commutative_period1_identity_witness": True,
      "shift_sign_relation_identity": True,
      "limitation": "failure of a specific witness never proves non-periodicity; no Koszul or AS test"
    }

if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
