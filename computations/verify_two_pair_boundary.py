#!/usr/bin/env python3
"""Exact two-pair degeneration with disjoint binary planes and a nonzero smoothing derivative."""
import argparse
import json
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

import sympy as s
from frobenius_tensor import flatten, reconstructed_windows, orthogonal_windows, nullspace_sparse
from tensor_symmetrization import encode, swap_factors, symmetrizer_space
from verify_pure_anticommutative_slice import cut_equations
from verify_quartic_obstruction import exact_rank, insertion_columns, relation_matrix


def fixture():
    alpha, epsilon = (1,0,2,1,3), (2,1,3,0,5)
    beta = ((2,3,5,7), (11,13,17,19))
    delta = ((2,1,0,3), (2,4,1,-1))
    pairs = [(0,0), (0,1), (1,1)]
    def binary_block(quartic):
        return s.Matrix(3, 3, lambda i, j: quartic[sum(pairs[i]+pairs[j])])
    A, E = binary_block(alpha), binary_block(epsilon)
    B = s.Matrix(3, 4, lambda i, j: beta[j%2][sum(pairs[i])+j//2])
    D = s.Matrix(4, 3, lambda i, j: delta[i//2][i%2+sum(pairs[j])])
    left, right = B.T*A.inv()*B, D*E.inv()*D.T
    derivative = right[0, 3]-right[1, 2]
    scale = -(left[0, 3]-left[1, 2])/derivative
    assert scale == -s.Rational(23, 162)
    delta = (delta[0], tuple(scale*v for v in delta[1]))
    D[2:, :] = scale*D[2:, :]
    C = left+D*E.inv()*D.T
    assert C[0, 3] == C[1, 2]
    def coefficient(a, b, c, d):
        first_r = [i for i in (a,b) if i < 2]
        last_r = [i for i in (c,d) if i < 2]
        first_s = [i-2 for i in (a,b) if i >= 2]
        last_s = [i-2 for i in (c,d) if i >= 2]
        if first_s and last_r:
            return 0
        r, t = first_r+last_r, first_s+last_s
        if len(r) == 4:
            return alpha[sum(r)]
        if len(r) == 3:
            return beta[t[0]][sum(r)]
        if len(t) == 4:
            return epsilon[sum(t)]
        if len(t) == 3:
            return delta[r[0]][sum(t)]
        return C[2*r[0]+t[0], 2*r[1]+t[1]]
    w = {}
    for word in product(range(4), repeat=4):
        value = coefficient(*word)
        if value:
            w[encode(word)] = F(value)
    return w, derivative


def no_pure_mixed_relations(window):
    relations = relation_matrix(window)
    x, y = s.symbols("x y")
    for side in ("A", "D"):
        excluded = [4*a+b for a, b in product(range(4), repeat=2)
                    if (a >= 2 if side == "A" else b < 2)]
        kernel = relations[excluded, :].nullspace()
        assert len(kernel) == 3
        mats = [(relations*v).reshape(4,4) for v in kernel]
        images = ([x*m[0,:].T+y*m[1,:].T for m in mats] if side == "A"
                  else [x*m[:,2]+y*m[:,3] for m in mats])
        pencil = s.Matrix.hstack(*images)
        minors = [pencil[list(rows), :].det() for rows in combinations(range(4), 3)]
        assert s.Poly(s.gcd_list(minors), x, y).total_degree() == 0


def verify_order(w):
    windows = reconstructed_windows(w)
    dual = orthogonal_windows(windows)
    for i in range(4):
        assert [len(nullspace_sparse(insertion_columns(dual, i, d), 4**d)) for d in (3,4)] == [4,1]
    return [256-exact_rank(insertion_columns(windows, i, 4)) for i in range(4)]


def opposite_plane_quotients():
    """Check both projections on two binary quartics supported on complementary planes."""
    quartics = ((1, 0, 2, 1, 3), (2, 1, 3, 0, 5))
    tensor = {}
    for plane, quartic in enumerate(quartics):
        for word in product(range(2), repeat=4):
            tensor[encode(tuple(2*plane+i for i in word))] = F(quartic[sum(word)])
    windows = [relation_matrix(window) for window in reconstructed_windows(tensor)]
    first, second = s.eye(4)[:2, :], s.eye(4)[2:, :]

    def project(arrows):
        return [s.kronecker_product(arrows[i], arrows[(i+1) % 4])*window
                for i, window in enumerate(windows)]

    edges = project([second, first, first, second])
    edge_ranks = [edge.rank() for edge in edges]
    assert edge_ranks == [0, 3, 0, 3]

    triples = project([first, s.eye(4), first, second])
    triple_ranks = [edge.rank() for edge in triples]
    assert triple_ranks == [3, 3, 0, 0]
    left = s.kronecker_product(triples[0], s.eye(2))
    right = s.kronecker_product(s.eye(2), triples[1])
    triple_dimension = left.rows - left.row_join(right).rank()
    assert triple_dimension == 8
    return {"first_edge_ranks": edge_ranks, "second_edge_ranks": triple_ranks,
            "second_length_four_block_dimension": 2*triple_dimension}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    w, derivative = fixture()
    assert [flatten(w, (i,)) for i in range(4)] == [4]*4
    assert [flatten(w, (0,j)) for j in (1,2,3)] == [6]*3
    assert swap_factors(w,0,1) == w and swap_factors(w,2,3) == w
    ranks = []
    for i in range(4):
        moved = swap_factors(w,0,i)
        row = []
        for j in range(4):
            if i == j:
                row.append(None)
            else:
                h, = symmetrizer_space(moved, i if j == 0 else j)
                row.append(h.rank())
        ranks.append(row)
    assert all(ranks[i][j] == (4 if i//2 == j//2 else 2)
               for i, j in product(range(4), repeat=2) if i != j)
    assert len(nullspace_sparse(cut_equations(w), 256)) == 70
    no_pure_mixed_relations(reconstructed_windows(w)[1])
    adjacent = verify_order(w)
    opposite = verify_order(swap_factors(w, 1, 2))
    assert adjacent == [41,42,41,42]
    result = {"result": "PASS", "arithmetic": "exact rational",
              "mixed_Schur_obstruction_derivative": str(derivative),
              "coefficient_scale": "-23/162", "all_singleton_ranks": 4, "all_pair_ranks": 6,
              "directed_symmetrizer_ranks_target_by_source": ranks,
              "projective_three_cut_tangent_dimension": 69,
              "finite_Frobenius_both_orders": True, "mixed_window_has_no_pure_tensor": True,
              "quartic_dimensions_adjacent_pairs": adjacent, "quartic_dimensions_opposite_pairs": opposite,
              "complementary_plane_projection_control": opposite_plane_quotients()}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print("PASS two-pair fixture, both cyclic orders, smooth tangent, Schur derivative and plane projections")


if __name__ == "__main__":
    main()
