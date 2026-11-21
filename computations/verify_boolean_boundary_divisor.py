#!/usr/bin/env python3
"""Three-cut tangent, directed symmetrizers, and transverse determinant at the Boolean boundary."""
import argparse
import json
from itertools import product
from pathlib import Path

import sympy as s
from frobenius_tensor import nullspace_sparse
from tensor_symmetrization import swap_factors, symmetrizer_space
from verify_nonkoszul_boundary import cubic_frame, tensor_at
from verify_pure_anticommutative_slice import cut_equations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    w = tensor_at(0)
    affine_dimension = len(nullspace_sparse(cut_equations(w), 256))
    assert affine_dimension == 70
    ranks = []
    for target in range(4):
        moved = swap_factors(w, 0, target)
        row = []
        for source in range(4):
            if source == target:
                row.append(None)
                continue
            slot = target if source == 0 else source
            generator, = symmetrizer_space(moved, slot)
            row.append(generator.rank())
        ranks.append(row)
    assert ranks == [[None, 4, 1, 4], [4, None, 1, 4], [3, 3, None, 3], [4, 4, 1, None]]

    t, monomials, frame = cubic_frame()
    octad = [v+(1,) for v in product(range(2), repeat=3)]
    weights = [s.cancel((-1)**sum(v[:3])*t/(t+sum(v[:3]))) for v in octad]
    def trace(word):
        return s.factor(sum(a*s.prod(v[k] for k in word) for v, a in zip(octad, weights)))
    free = [(0,1,2), (1,3,3), (2,3,3), (3,3,3)]
    H = s.Matrix(4, 4, lambda i, j: trace((i,)+free[j]))
    all_cubics = s.Matrix(4, 20, lambda i, j: trace((i,)+monomials[j]))
    assert (H*frame.T-all_cubics).applyfunc(s.cancel) == s.zeros(4, 20)
    assert H.subs(t, 0) == s.diag(0, 0, 0, 1)
    determinant = s.factor(H.det())
    assert s.cancel(determinant+2*t**3/((t+1)**3*(t+2)**3*(t+3))) == 0
    reverse = (t*H.inv()).applyfunc(s.cancel).subs(t, 0)
    assert reverse.rank() == 3 and reverse*H.subs(t, 0) == s.zeros(4)
    result = {"result": "PASS", "arithmetic": "exact rational and QQ(t)",
              "three_cut_projective_tangent_dimension": affine_dimension-1,
              "directed_symmetrizer_ranks_target_by_source": ranks,
              "all_off_diagonal_symmetrizer_spaces_dimension": 1,
              "normalized_determinant": str(determinant), "order_at_zero": 3,
              "regularized_inverse_rank_at_zero": 3}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print("PASS smooth-boundary tangent, all 12 oriented lines, transverse order three")


if __name__ == "__main__":
    main()
