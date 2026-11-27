"""Linear completions of the nineteen-point marked edge at all twenty split marks."""
from itertools import product

import numpy as np
from flint import nmod_mat

from verify_mixed_bracket_reconstruction import G, Q, P, brackets, units
from verify_twenty_mixed_bracket_invariants import configurations, marked_tensor, tensor


def centered(array):
    identity = np.eye(3, dtype=np.int64)
    return (array-np.array([a.trace()*pow(3,-1,P)*identity for a in array])) % P


def companion_rank(t):
    G0 = centered(t[1:, 1:, 0, 1:])
    F0 = centered(t[0, 1:, 1:, 1:].transpose(1, 0, 2))
    Q0 = centered(t[1:, 1:, 1:, 0].transpose(1, 0, 2))
    P0 = centered(t[1:, 0, 1:, 1:].transpose(2, 0, 1))
    def equations(g, p):
        return (brackets(g, F0)+brackets(Q0, p).transpose(2, 3, 0, 1)).ravel() % P
    assert not np.any(equations(G0, P0)) and np.any(G0) and np.any(P0)
    zero = np.zeros_like(G0)
    columns = []
    for h, basis in product(range(3), units(P)):
        slot = zero.copy()
        slot[h] = basis
        columns.extend((equations(slot, zero), equations(zero, slot)))
    assert nmod_mat(np.array(columns).tolist(), P).rank() == 47


if __name__ == '__main__':
    W = tensor((G, -G, Q, Q))
    odd, first, third = configurations()
    for track in zip(first, odd, third, odd):
        t = (marked_tensor(W, np.array(track, dtype=object)) % P).astype(np.int64)
        for i in range(4):
            companion_rank(t.transpose(tuple((i+j)%4 for j in range(4))))
    print('PASS: all 80 adjacent nineteen-point data have a one-dimensional companion kernel')
