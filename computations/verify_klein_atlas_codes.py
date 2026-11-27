"""Saturated torus codes and finite normalizer for index-preserving Klein markings."""
from functools import lru_cache
from itertools import permutations, product
from math import prod

import sympy as sp

from verify_negative_component_quotient import lattice, words

WORDS = tuple(words)


@lru_cache(None)
def characters(support):
    return tuple(tuple(map(int, row)) for row in lattice(support).tolist())


def torus_code(coefficients):
    support = tuple(sorted(q for q, value in coefficients.items() if value != 0))
    values = [sp.sympify(coefficients[q]) for q in support]
    return support, tuple(prod(v**n for v, n in zip(values, row) if n)
                          for row in characters(support))


def normalizer():
    for linear in permutations((1, 2, 3)):
        linear = (0, *linear)
        for a, b, c in product(range(4), repeat=3):
            offsets = a, b, c, a ^ b ^ c
            yield tuple(tuple(linear[i] ^ d for i, d in zip(q, offsets)) for q in WORDS)


def equivalent_marked(left, right):
    target = torus_code(right)
    for image in normalizer():
        transformed = {r: left[q] for q, r in zip(WORDS, image) if left.get(q, 0)}
        if tuple(sorted(transformed)) == target[0] and torus_code(transformed) == target:
            return True
    return False


def finite_relabeling_check():
    # Three affine permutations determine the fourth by total degree zero.
    # Restricting the equation to one free index fixes their common linear part.
    affine = list(permutations(range(4)))
    maps = set()
    for first, second, third in product(affine, repeat=3):
        fourth = tuple(first[i] ^ second[0] ^ third[0] for i in range(4))
        if any(first[0] ^ second[i] ^ third[0] != fourth[i] or
               first[0] ^ second[0] ^ third[i] != fourth[i] for i in range(4)):
            continue
        maps.add(tuple(tuple(f[i] for f, i in zip((first, second, third, fourth), q))
                       for q in WORDS))
    assert len(maps) == 384 and maps == set(normalizer())


if __name__ == '__main__':
    finite_relabeling_check()
    assert len(characters(WORDS)) == 51
    sparse = {q: 1+i for i, q in enumerate(WORDS) if (3*i+1)%7 < 3}
    scales = ((2,3,1,7), (5,1,3,2), (11,3,2,1), (3,1,1,2))
    scaled = {q: value*prod(row[i] for row, i in zip(scales, q)) for q, value in sparse.items()}
    assert torus_code(sparse) == torus_code(scaled)
    image = list(normalizer())[143]
    relabeled = {r: scaled[q] for q, r in zip(WORDS, image) if q in scaled}
    assert equivalent_marked(sparse, relabeled)
    changed = dict(relabeled)
    changed[next(iter(changed))] *= 37
    assert not equivalent_marked(sparse, changed)
    print('384 simultaneous affine relabelings, saturated full-support rank 51, sparse comparisons passed.')
