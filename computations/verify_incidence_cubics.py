#!/usr/bin/env python3
"""Exact finite reduced six-bilinear-equation control for the cubic restriction."""
from itertools import combinations, combinations_with_replacement

import sympy as sp


def main():
    # Equations l_i(x) m_i(y)=0; the two lists evaluate cubic polynomials at six nodes.
    nodes = (0, 1, 2, 3, 4, 5)
    other_nodes = (6, 8, 10, 12, 14, 16)
    triples = list(combinations(range(6), 3))
    monomials = list(combinations_with_replacement(range(4), 3))
    t = sp.symbols("t")

    def point(values):
        return sp.Poly(sp.prod(t-v for v in values), t).all_coeffs()[::-1]

    left, right = [], []
    for split in triples:
        x = point(nodes[i] for i in split)
        y = point(other_nodes[i] for i in range(6) if i not in split)
        assert all(sum(a*n**d for d, a in enumerate(x)) *
                   sum(b*m**d for d, b in enumerate(y)) == 0
                   for n, m in zip(nodes, other_nodes))
        for vector, rows in ((x, left), (y, right)):
            rows.append([vector[i]*vector[j]*vector[k] for i, j, k in monomials])
    assert sp.Matrix(left).rank() == sp.Matrix(right).rank() == 20
    print("PASS 20 distinct intersections and both cubic restriction ranks 20 over Q")


if __name__ == "__main__":
    main()
