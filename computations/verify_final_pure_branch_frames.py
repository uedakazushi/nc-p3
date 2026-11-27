"""Linear incidence on the two residual cyclic frames, over the rationals."""
from itertools import product

import sympy as sp

# The three-coordinate parts of L, M are (a,b,0), (0,a,b).
LINE = [{(0, 0)}, {(0, 1)}, {(1, 2)}, {(2, 2)}, {(0, 2), (1, 1)}]


def kernel(edges, alphabets):
    words = list(product(*alphabets))
    rows = []
    for i, groups in enumerate(edges):
        for group in groups:
            for other in product(alphabets[(i+2) % 4], alphabets[(i+3) % 4]):
                row = [int((w[i], w[(i+1) % 4]) in group
                           and (w[(i+2) % 4], w[(i+3) % 4]) == other)
                       for w in words]
                if any(row):
                    rows.append(row)
    return words, sp.Matrix(rows).nullspace()


def four():
    words, basis = kernel([LINE]*4, [range(3)]*4)
    assert len(basis) == 1
    core = dict(zip(words, basis[0]))
    pairs = [(1, 0), (1, 1), (2, 0), (2, 1)]
    block = sp.Matrix([[core[r+c] for c in pairs] for r in pairs])
    assert block.det() and block.inv()[0, 0] == block.inv()[3, 3] == 0
    assert all(core[w] == core[w[1:]+w[:1]] for w in words)
    words, lifts = kernel([LINE]*4, [[3]]+[range(3)]*3)
    insertions = [sp.Matrix([core[(i,)+w[1:]] for w in words]) for i in range(3)]
    extra = sp.Matrix([int(w == (3, 2, 1, 0)) for w in words])
    assert len(lifts) == 4 and sp.Matrix.hstack(*insertions, extra).rank() == 4
    assert sp.Matrix.hstack(*lifts, *insertions, extra).rank() == 4


def three():
    start = [{(0, 0)}, {(0, 1)}] + [
        {(i, j+1) for i, j in product(range(4), range(2)) if i+j == n}
        for n in range(5)]
    end = [{(1, 3)}, {(2, 3)}] + [
        {(j, i) for i, j in product(range(4), range(2)) if i+j == n}
        for n in range(5)]
    words, basis = kernel([start, LINE, LINE, end], [range(4)]+[range(3)]*3)
    assert len(basis) == 2
    # In the two kernel vectors the four binary cubic coefficients differ by one shift.
    first, second = [dict(zip(words, vector)) for vector in basis]
    for i, a, b, c in words:
        assert second[i, a, b, c] == (first[i-1, a, b, c] if i else 0)
        if i == 3:
            assert first[i, a, b, c] == 0


if __name__ == '__main__':
    four()
    three()
    print("PASS: cyclic core, complementary inverse entries, single-letter lifts and common linear root")
