#!/usr/bin/env python3
"""Exact tangent spaces at a symmetric complete intersection with no square in its pencil."""
import argparse
import json
from fractions import Fraction as F
from itertools import combinations_with_replacement, product
from pathlib import Path

from frobenius_tensor import nullspace_sparse
from single_negative_slice import opposite_rank_six
from tensor_symmetrization import encode, swap_factors
from verify_pure_anticommutative_slice import cut_equations
from verify_quartic_obstruction import exact_rank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    co = [1, 2, 3, 5]
    quartics = list(combinations_with_replacement(range(4), 4))
    indices = {word: i for i, word in enumerate(quartics)}
    quadrics = [[(i, i, 1), ((i+1)%4, (i+2)%4, c)] for i, c in enumerate(co)]
    rows = []
    for quadric in quadrics:
        for tail in combinations_with_replacement(range(4), 2):
            rows.append({indices[tuple(sorted((a, b)+tail))]: F(v) for a, b, v in quadric})
    trace, = nullspace_sparse(rows, 35)
    w = {encode(word): trace[indices[tuple(sorted(word))]] for word in product(range(4), repeat=4)
         if indices[tuple(sorted(word))] in trace}
    first, second, _ = opposite_rank_six(w)
    P = [first[4*i:4*i+4, :] for i in range(4)]
    Q = [second[:, 4*i:4*i+4] for i in range(4)]

    # One linearized equation for each matrix entry of each commutator/quadric.
    equations = quadrics + [[(i, j, 1), (j, i, -1)] for i in range(4) for j in range(i)]
    tangent = []
    for equation in equations:
        for r, t in product(range(4), repeat=2):
            row = {}
            for a, b, scale in equation:
                for k in range(6):
                    p = 24*a+6*r+k
                    q = 96+24*b+4*k+t
                    row[p] = row.get(p, 0)+scale*F(Q[b][k, t])
                    row[q] = row.get(q, 0)+scale*F(P[a][r, k])
            tangent.append(row)
    rank = exact_rank(tangent)
    assert rank == 125 and 192-rank == 4**2+6**2+4**2-1

    kernel = nullspace_sparse(cut_equations(w), 256)
    assert len(kernel) == 70
    symmetrized = []
    for t in kernel:
        swapped = swap_factors(t, 0, 1)
        symmetrized.append({k: t.get(k, 0)+swapped.get(k, 0) for k in t.keys() | swapped.keys()})
    assert exact_rank(symmetrized) == 55
    result = {"result": "PASS", "arithmetic": "exact rational",
              "quadrics": "x_i^2+c_i*x_(i+1)*x_(i+2), cyclic indices", "c": co,
              "quartic_trace_dimension": 1, "module_linearization_rank": rank,
              "fixed_ring_module_tangent": 67, "three_cut_affine_tangent": 70,
              "swap01_fixed_affine_tangent": 55, "three_cut_projective_tangent": 69}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print("PASS fixed-CI module tangent 67, swap-fixed tangent 55, projective tangent 69")


if __name__ == "__main__":
    main()
