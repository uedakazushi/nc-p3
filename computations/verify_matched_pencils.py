#!/usr/bin/env python3
"""Exact ranks of both adjacent restrictions on every positive-dimensional plane fiber."""
from itertools import product

import sympy as sp


def pencil(indices):
    # Drop unused output rows. Each column is multiplication by s or by t.
    return tuple(sp.diag(*(sp.eye(e+1)[:, k:k+e] for e in indices)) for k in (0, 1))


def restriction_rank(left, right):
    return (sp.kronecker_product(left[0], right[0]) +
            sp.kronecker_product(left[1], right[1])).rank()


def quadric_cut():
    multiply = pencil((3,))
    contract = tuple(m.T for m in multiply)
    central = sp.Matrix([[2, 1, 0], [0, 3, 1], [1, 0, 1]])
    partner = sp.Matrix([[1, 2, 3], [0, 2, 1], [1, 0, 4]])
    assert central.det() and partner.det()
    # Independent changes of the two parameter bases in the second ruling.
    mult_other = (multiply[0]+multiply[1], multiply[0]-2*multiply[1])
    contr_other = (2*contract[0]+contract[1], contract[0]+contract[1])
    cut = sp.zeros(16)
    for i, j in product(range(2), repeat=2):
        left = multiply[j]*central*contract[i]
        right = mult_other[j]*partner*contr_other[i]
        cut += sp.kronecker_product(left, right)
    outer = sum((sp.kronecker_product(a,b) for a,b in zip(multiply, mult_other)), sp.zeros(16,9))
    inner = sum((sp.kronecker_product(a,b) for a,b in zip(contract, contr_other)), sp.zeros(9,16))
    assert cut == outer*sp.kronecker_product(central, partner)*inner
    assert outer.rank() == inner.rank() == cut.rank() == 9


def main():
    rows = [(3,), (2,), (1, 1)]
    columns = [(1,), (2,), (1, 1), (3,)]
    expected = [[(3, 4), (6, 8), (6, 8), (9, 12)],
                [(2, 3), (4, 6), (4, 6), (6, 8)],
                [(2, 4), (4, 6), (4, 8), (6, 8)]]
    for i, j in product(range(3), range(4)):
        left, right = pencil(rows[i]), pencil(columns[j])
        ranks = []
        admissible = True
        for pair in (right, tuple(a.T for a in right)):
            rank = restriction_rank(left, pair)
            ranks.append(rank)
            proper_support = min(left[0].rows, pair[0].rows) < 4
            admissible &= rank <= 6 - int(proper_support)
        assert tuple(ranks) == expected[i][j]
        assert admissible == (j == 0)
    quadric_cut()
    print("PASS all 12 pencil comparisons and independent Segre cut rank nine")


if __name__ == "__main__":
    main()
