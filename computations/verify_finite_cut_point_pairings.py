"""Evaluate all cut point pairings on the twenty lifted marks, with no random sampling."""
from collections import Counter
from itertools import product

import numpy as np

from verify_finite_three_slot_evaluations import pivot_orders
from verify_twenty_mixed_bracket_invariants import (
    MOD, P, configurations, lift, lift_points, tensor,
)


def pair_matrix(w, left, opposite):
    return np.tensordot(np.tensordot(left, w, axes=1), opposite, axes=(1, 0)) % MOD


def main():
    w = tensor(lift())
    middle, first, third = configurations()
    marks = np.array([lift_points(w, mark) for mark in zip(first, middle, third, middle)])
    generic = Counter()
    special = Counter()

    for a, b in product(range(20), repeat=2):
        ranks = []
        for i in range(4):
            order = [(i+j) % 4 for j in range(4)]
            matrix = pair_matrix(w.transpose(order), marks[a, i], marks[b, (i+2) % 4])
            pivots = pivot_orders(matrix)
            ranks.append(pivots.count(0))
            assert len(pivots) == (3 if a == b else 2)
        assert len(set(ranks)) == 1
        generic[a == b] += 1
        special[a == b, ranks[0]] += 1
    assert generic == {True: 20, False: 380}
    print(f'PASS: four phases, 20 self-pairings and 380 ordered off-diagonal pairings; mod {P}: {special}')


if __name__ == '__main__':
    main()
