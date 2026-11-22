"""Exact mixed spans for matched branches and attachments in pairs of additive quartics."""
from itertools import product

import sympy as sp

t = sp.Symbol("t")


def sections(degrees):
    """Common value, sum of first derivatives zero, independent higher terms."""
    n = len(degrees)
    basis = [[sp.S.One]*n]
    for i in range(n-1):
        row = [sp.S.Zero]*n
        row[i], row[-1] = t, -t
        basis.append(row)
    for i, degree in enumerate(degrees):
        for power in range(2, degree+1):
            row = [sp.S.Zero]*n
            row[i] = t**power
            basis.append(row)
    assert len(basis) == 4
    return basis


def rank(degrees, other_degrees, components):
    left, right = sections(degrees), sections(other_degrees)
    rows = []
    for i, x, j, y in components:
        polys = [sp.Poly(a[i].subs(t,x)*b[j].subs(t,y),t)
                 for a,b in product(left,right)]
        max_degree = max(p.degree() for p in polys)
        rows.extend([p.nth(k) for p in polys] for k in range(max_degree+1))
    return sp.Matrix(rows).rank()


def main():
    zero = sp.S.Zero
    pair = (3,1)
    graph = [(0,t,0,t), (1,t,1,t)]
    assert rank(pair,pair,graph) == 8
    assert rank(pair,pair,[graph[0], (1,t,1,2*t)]) == 9
    # One shared cubic and two loose lines, attached at the singular center.
    central = [graph[0], (1,t,1,zero), (1,zero,1,t)]
    assert rank(pair,pair,central) == 8
    # Moving the vertical attachment away from the center breaks the rank bound.
    assert rank(pair,pair,[*central[:2], (1,1,1,t)]) == 9
    # Disjoint centers, with no matched branch or a matched line or cubic.
    cross = [(i,t,0,1) for i in range(2)] + [(0,1,j,t) for j in range(2)]
    assert rank(pair,pair,cross) == 7
    for matched, expected in [(0,9), (1,8)]:
        i = matched
        components = [(i,t,i,t-1), (1-i,t,i,-1), (i,1,1-i,t)]
        assert rank(pair,pair,components) == expected
    degrees = (2,1,1)
    assert rank(degrees,degrees,[(i,t,i,t) for i in range(3)]) == 8
    print("PASS: additive graph spans 7,8,9, distinct centers, attachments and shared tangents")


if __name__ == "__main__":
    main()
