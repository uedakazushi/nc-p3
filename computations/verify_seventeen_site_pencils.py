"""Infinitesimal conic pencil and confluent square identities for 17-site graphs."""
from itertools import combinations

import sympy as sp


def check_contact():
    a, b, c, d, x = sp.symbols("A B C D x")
    modulus = x**2-b*x-a
    z = c+d*x
    assert sp.rem(x*z, modulus, x) == a*d+(b*d+c)*x
    reduced = sp.rem(z**2, modulus, x).subs(c, -d*b).expand()
    assert reduced.coeff(x) == -b*d**2
    for s in range(1, 14):
        for f in range(s+1):
            for va in range(s+1):
                for vb in range(s+1):
                    if f+va >= s and 2*f+vb >= s:
                        assert 2*f+min(va, vb) >= s


def check_jets():
    # Exact linear identities used when the missing square is the edge 01.
    # D_i has coefficient a_j at single j and -a_j at edge ij for j!=i.
    a = sp.symbols("a:4", nonzero=True)
    edges = list(combinations(range(4), 2))
    J = sp.Matrix([[a[j]+(sum(a) if i == j else 0) for j in range(4)]
                   + [-(a[v] if i == u else a[u] if i == v else 0) for u, v in edges]
                   for i in range(4)])
    missed = 4+edges.index((0, 1))
    cross = 4+edges.index((0, 2))
    assert J[1, missed] != 0 and J[1, cross] == 0
    assert J[2, missed] == 0
    assert [J[0, cross], J[2, cross], J[2, 0]] == [-a[2], -a[0], a[0]]
    assert max(11-k*(4-k) for k in (1, 2, 3)) == 8
    assert (15-2+1)//2+1 > 7
    assert 2*(4+3) == 14 and 4+3+1 > 7


if __name__ == "__main__":
    check_contact()
    check_jets()
    print("PASS: two-root contact equations, truncated valuations, and missing-edge first jets")
