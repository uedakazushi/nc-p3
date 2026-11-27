"""Check support bounds and the three-letter evaluations on finite-window fixtures."""
import numpy as np

from verify_finite_square_chart import frame, products
from verify_twenty_mixed_bracket_invariants import (
    DIGITS, MOD, P, configurations, lift, lift_points, tensor,
)


def restriction_bound(length):
    return sum(length >= threshold for threshold in (1, 2, 4, 8))


def pivot_orders(matrix):
    """Smith pivot valuations below the working precision, using only unit reductions."""
    matrix = matrix.copy() % MOD
    modulus = MOD
    order = 0
    orders = []
    while min(matrix.shape) and order < DIGITS:
        pivots = np.argwhere(matrix % P != 0)
        if not len(pivots):
            matrix //= P
            modulus //= P
            order += 1
            continue
        i, j = pivots[0]
        matrix[[0, i]] = matrix[[i, 0]]
        matrix[:, [0, j]] = matrix[:, [j, 0]]
        first = matrix[0] * pow(int(matrix[0, 0]), -1, modulus) % modulus
        matrix = (matrix[1:, 1:] - matrix[1:, :1]*first[None, 1:]) % modulus
        orders.append(order)
    return orders


def main():
    for length in range(1, 6):
        assert 2*restriction_bound(length) > length
    for length in range(1, 8):
        assert length-restriction_bound(length)+1 <= 5
    assert max(a+b for a in range(1, 21) for b in range(1, 21)
               if restriction_bound(a)+restriction_bound(b) <= 4) == 8
    u, v = frame()
    assert products(products(u, u, True), v).to_DM().rank() == 16
    w = tensor(lift())
    middle, first, third = configurations()
    points = np.array([lift_points(w, mark) for mark in zip(first, middle, third, middle)])
    for slot in range(4):
        a, b, c = [points[:, (slot+k) % 4, :] for k in range(3)]
        triple = (a[:, :, None, None]*b[:, None, :, None]*c[:, None, None, :]).reshape(20, 64)
        orders = pivot_orders(triple)
        expected = [0]*20 if slot % 2 == 0 else [0]*16 + [1]*4
        assert orders == expected
        print(f'slot {slot}: rank 20 with pivot valuations {orders}')
    print('PASS: support inequalities, symmetric evaluation and lifted evaluation matrices')


if __name__ == '__main__':
    main()
