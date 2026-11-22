#!/usr/bin/env python3
"""Symbolic checks of the five-coefficient normal form and an exact four-edge quotient."""
from itertools import combinations, product

import sympy as sp


def main():
    a, b, c, d, n = coordinates = sp.symbols("a b c d n")
    u, v = sp.symbols("u v")
    small = sp.Matrix([[a, b], [c, d]])
    matrix = small.row_join(n*sp.eye(2)).col_join(sp.zeros(2).row_join(small.adjugate()))
    quadric = a*d-b*c
    assert sp.expand(matrix.det()-quadric**2) == 0
    triples = list(combinations(range(4), 3))
    minors = [matrix.extract(rows, cols).det() for rows, cols in product(triples, repeat=2)]
    assert sp.Poly(sp.gcd_list(minors), *coordinates).monic().as_expr() == quadric
    front = matrix*sp.Matrix([v, -u, 0, 0])
    back = sp.Matrix([[0, 0, u, v]])*matrix
    assert sp.expand(front[0]-back[3]) == sp.expand(front[1]+back[2]) == 0
    assert front[2:] == [0, 0] and list(back)[:2] == [0, 0]

    # A rank-five tensor for the quotient check; no regularity is assumed.
    def units(*positions):
        result = sp.zeros(4)
        for i, j in positions:
            result[i, j] = 1
        return result

    opposite = [units((0, 0), (1, 1)), units((3, 2)),
                units((1, 0), (2, 1)), units((0, 3)), units((3, 3))]
    coeff = [matrix.diff(t) for t in coordinates]
    def entry(word):
        i, j, k, l = word
        return sum(A[k, i]*B[l, j] for A, B in zip(coeff, opposite))
    words = list(product(range(4), repeat=4))
    paired = sp.Matrix(16, 16, lambda i, j: entry((i//4, j//4, i%4, j%4)))
    assert paired.rank() == 5
    # Keep the first source covector, its matched opposite, and complements to both supports.
    keep = ({0}, {2, 3}, {3}, {3})
    for i in range(4):
        j = (i+1) % 4
        assert not any(entry(word) for word in words if word[i] in keep[i] and word[j] in keep[j])
    assert sp.prod(len(slot) for slot in keep) == 2
    print("PASS determinant, cubic gcd, matched binary planes and four zero edge images")


if __name__ == "__main__":
    main()
