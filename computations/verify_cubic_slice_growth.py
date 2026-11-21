#!/usr/bin/env python3
"""Compute the repeated three-letter quotient from four relation spaces and a covector."""
import argparse
import json
from itertools import combinations, product
from pathlib import Path

import sympy as s
from frobenius_tensor import reconstructed_windows
from general_cut_tor import word_frame
from verify_indexed_anticommutative_ore import window
from verify_nonkoszul_boundary import tensor_at
from verify_quartic_obstruction import relation_matrix


def quotient_map(columns):
    functionals = columns.T.nullspace()
    return s.Matrix.vstack(*(v.T for v in functionals)) if functionals else s.zeros(0, columns.rows)


def block_dimension(windows, marker, covector):
    ell = s.Matrix([covector])
    if ell.shape != (1, 4) or not any(covector):
        raise ValueError("the cut needs a nonzero four-dimensional covector")
    right, first, second, left = [relation_matrix(windows[(marker+i)%4]) for i in range(4)]
    identity = s.eye(4)
    X = quotient_map(s.kronecker_product(ell, identity)*right)
    Z = quotient_map(s.kronecker_product(identity, ell)*left)
    if X.rows == 0 or Z.rows == 0:
        return {"outer_dimensions": [X.rows, Z.rows], "block_dimension": 0}
    r0 = s.kronecker_product(X, identity)*first
    r1 = s.kronecker_product(identity, Z)*second
    insertions = s.kronecker_product(r0, s.eye(Z.rows)).row_join(s.kronecker_product(s.eye(X.rows), r1))
    rank = insertions.rank()
    return {"outer_dimensions": [X.rows, Z.rows], "projected_relations": [r0.rank(), r1.rank()],
            "insertion_rank": rank, "block_dimension": insertions.rows-rank}


def marker_matrix(windows, marker, ell):
    """The linear 20 by 120 matrix, evaluated at a given covector."""
    start = (marker+1) % 4
    dual, _ = word_frame(windows, start, 3)
    _, tails = word_frame(windows, (start+1) % 4, 2)
    _, heads = word_frame(windows, start, 2)
    columns = []
    for r in windows[marker]:
        contracted = [sum(ell[a]*r.get((a,b), 0) for a in range(4)) for b in range(4)]
        for word in tails:
            columns.append([sum(v*row.get(16*b+word, 0) for b, v in enumerate(contracted)) for row in dual])
    for r in windows[(marker-1) % 4]:
        contracted = [sum(ell[b]*r.get((a,b), 0) for b in range(4)) for a in range(4)]
        for word in heads:
            columns.append([sum(v*row.get(4*word+a, 0) for a, v in enumerate(contracted)) for row in dual])
    return s.Matrix(columns).T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    p = [0, 0, 0, 1]
    bad = block_dimension(reconstructed_windows(tensor_at(0)), 2, p)
    assert bad == {"outer_dimensions": [3, 3], "projected_relations": [5, 5],
                   "insertion_rank": 26, "block_dimension": 10}
    matrix = marker_matrix(reconstructed_windows(tensor_at(0)), 2, p)
    assert matrix.shape == (20, 120) and matrix.rank() == 10
    commutators = [{(a,b): 1, (b,a): -1} for a, b in combinations(range(4), 2)]
    ore = [window([0]*9, [1,0,0], i) for i in range(4)]
    records = []
    for marker, ell in product(range(4), (p, [1, 0, 0, 0], [1, 2, 3, 5])):
        polynomial = block_dimension([commutators]*4, marker, ell)["block_dimension"]
        signed = block_dimension(ore, marker, ell)["block_dimension"]
        assert polynomial == 1 and signed <= 1
        assert 20-marker_matrix([commutators]*4, marker, ell).rank() == polynomial
        assert 20-marker_matrix(ore, marker, ell).rank() == signed
        records.append({"marker": marker, "covector": ell, "polynomial": polynomial, "signed_Ore": signed})
    result = {"result": "PASS", "arithmetic": "exact rational",
              "Boolean_t0_marker2_pure_cube": bad, "Boolean_marker_matrix_rank": 10, "regular_controls": records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print("PASS Boolean block dimension 10; all 24 regular control block dimensions at most one")


if __name__ == "__main__":
    main()
