"""Exact balanced weights and exposing cocharacters for every sparse regular doublet support."""
from itertools import combinations, product
from pathlib import Path

import sympy as sp

from verify_klein_atlas_codes import WORDS, normalizer
from verify_klein_single_missing_supports import SURVIVORS, decode

K = range(4)
ORDER = [(i, i^h, k, k^h) for h, i, k in product(K, repeat=3)]
REGULAR = (50,51,52,54,55,56,57,58,59,62,63,64,65,66,76,78,82,85,88,89,90,
           94,95,96,112,115,117,118,124,129,134,139,142,143,144,145,146,147)


def incidence(support):
    return sp.Matrix([[int(q[i] == j) for q in support] for i, j in product(K, repeat=2)])


def balanced(support):
    mat = incidence(support)
    rows = mat[list(mat.T.rref()[1]), :]
    inverse = rows.T*(rows*rows.T).inv()
    ones = sp.ones(rows.rows, 1)
    weight = inverse*ones
    if min(weight) <= 0:
        guess = weight.applyfunc(lambda x: max(x, sp.Rational(1, 32)))
        weight = guess+inverse*(ones-rows*guess)
    assert min(weight) > 0 and mat*weight == sp.ones(16, 1)


def terms(mask):
    return frozenset(q for i, q in enumerate(ORDER) if mask >> i & 1)


def relabelings(support):
    return {frozenset(r for q, r in zip(WORDS, images) if q in support)
            for images in normalizer()}


def sparse():
    masks = [int(line, 16) for line in Path(__file__).with_name(
        'klein_doublet_parallel_masks.txt').read_text().splitlines()]
    # Each array is the ordered list of four traceless diagonal cocharacters.
    to_50 = (
        ((51,62,63), (0,2,-1,-1, -1,-1,0,2, 2,0,-1,-1, -1,-1,0,2)),
        ((54,55,64,65,76,78,112,115,117), (1,1,3,-5, -5,3,1,1, 1,1,3,-5, 3,-5,1,1)),
        ((56,), (2,0,-1,-1, -1,-1,2,0, 2,0,-1,-1, -1,-1,0,2)),
        ((57,), (0,2,-1,-1, -1,-1,0,2, 0,2,-1,-1, -1,-1,2,0)),
        ((58,), (1,-1,0,0, 0,0,1,-1, 1,-1,0,0, 0,0,-1,1)),
    )
    to_52 = (((59,66), (1,1,-1,-1, -1,-1,1,1, 1,1,-1,-1, -1,-1,1,1)),)
    handled = set()
    for destination, certificates in ((50, to_50), (52, to_52)):
        target = terms(masks[destination])
        balanced(sorted(target))
        orbit = relabelings(target)
        for indices, weights in certificates:
            assert all(sum(weights[4*i:4*i+4]) == 0 for i in K)
            for n in indices:
                support = terms(masks[n])
                degrees = {q: sum(weights[4*i+j] for i, j in enumerate(q)) for q in support}
                assert min(degrees.values()) == 0 and max(degrees.values()) > 0
                assert frozenset(q for q,d in degrees.items() if d == 0) in orbit
                handled.add(n)
        for start, color in product((0, 1), K):
            def present(i, j):
                q = (i, i^color, j, j^color)
                return int((q[-start:]+q[:-start] if start else q) in target)
            mat = sp.Matrix(4, 4, present)
            if color < 2:
                assert any(abs(mat.extract(I,J).det()) == 1
                           for I,J in product(combinations(K,2),repeat=2))
            else:
                assert any(mat)
    assert len(handled) == 17
    for n in set(REGULAR)-handled:
        balanced(sorted(terms(masks[n])))
    for mask in (0x0ff066663c3caa55, 0x0ff066667e7eaff5):
        balanced(sorted(terms(mask)))
    balanced(sorted(decode(SURVIVORS[-1])))


if __name__ == '__main__':
    sparse()
    dense = tuple(WORDS)
    for pair, expected_rank in (((1,1),12), ((1,2),13)):
        thin = [q for q in dense if q[0]^q[1] not in (0,pair[0]) or q[1]^q[2] not in (0,pair[1])]
        assert len(thin) == 48 and incidence(thin).rank() == expected_rank
        assert incidence(thin)*sp.ones(48,1) == 12*sp.ones(16,1)
    nonconstant = [q for q in dense if len(set(q)) != 1]
    assert incidence(nonconstant).rank() == incidence(dense).rank() == 13
    assert incidence(nonconstant)*sp.ones(60,1) == 15*sp.ones(16,1)
    print('All sparse supports checked. Seventeen nonclosed doublet supports limit to 50 or 52.')
