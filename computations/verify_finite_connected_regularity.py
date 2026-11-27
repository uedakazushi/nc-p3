"""Exact weight layers for the finite additive cuts and the exceptional torus support."""
from itertools import product

import sympy as sp

from verify_finite_additive_windows import jordan, tensor_sum


def layers(first, second):
    def weights(parts):
        return [2*j-d+1 for d in parts for j in range(d)]
    basis = list(product(weights(first), weights(second)))
    total = tensor_sum(jordan(first), jordan(second))
    return total, {h: [i for i, pair in enumerate(basis) if sum(pair) == h]
                   for h in set(map(sum, basis))}


def main():
    total, weights = layers((3, 1), (2, 1, 1))
    # After deleting zz, D^2 lands entirely in the highest left letter.
    assert (total**2)[:, 1:][[i for i in range(16) if i//4 != 2], :].is_zero_matrix
    for first, second, bottom, length, needed in [
            ((4,), (3, 1), -3, 4, 2), ((4,), (4,), -4, 5, 1),
            ((3, 1), (3, 1), -2, 3, 2)]:
        total, weights = layers(first, second)
        cols = weights[bottom]
        lift = sp.eye(16)[:, cols]
        allowed = lift*sp.Matrix.hstack(*((total**length)*lift).nullspace())
        assert allowed[:4, :].rank() == needed
        assert (total**(length-1)*allowed).rank() == allowed.cols
        assert all((total**length)[:, weights[h]].is_zero_matrix for h in weights if h > bottom)
    # Any submodule of L4 x L4 of type L5 + L1 maps in degree zero
    # to alternating tensors, including when the scalar block drops rank.
    n = jordan((4,))
    total, weights = layers((4,), (4,))
    flip = sp.Matrix(16, 16, lambda i,j: int(i == 4*(j%4)+j//4))
    for bottom, length in [(-4, 5), (0, 1)]:
        lift = sp.eye(16)[:, weights[bottom]]
        start = lift*sp.Matrix.hstack(*((total**length)*lift).nullspace())
        assert start.cols == 1
        assert ((flip+sp.eye(16))*start).is_zero_matrix
    a, b = [0, 1, 2, 3], [0, 1, 3, 4]
    word = (3, 0, 0, 3)
    sums = [v[word[i]] + w[word[(i+1)%4]] for i,(v,w) in enumerate(((a,b),(b,a),(a,b),(b,a)))]
    assert sums == [3, 0, 4, 7]
    print('All mixed and uniform Jordan layer constraints and the conflicting mixed torus coefficient verified.')


if __name__ == '__main__':
    main()
