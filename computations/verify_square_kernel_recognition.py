#!/usr/bin/env python3
"""Recover the squarefree normal form after independent rational basis changes."""
import argparse
import json
from itertools import combinations, permutations, product
from pathlib import Path

import sympy as s
from frobenius_tensor import reconstructed_windows, orthogonal_windows, nullspace_sparse
from single_negative_slice import opposite_rank_six
from tensor_symmetrization import act_factor, encode
from verify_quartic_obstruction import insertion_columns


def verify(changes):
    standard = {encode(word): 1 for word in permutations(range(4))}
    w = standard
    for slot, change in enumerate(changes):
        w = act_factor(w, change, slot)
    first, second, _ = opposite_rank_six(w)
    axes = changes[0].inv()
    P, Q = [], []
    for axis in range(4):
        P.append(sum((axes[axis, a]*first[4*a:4*a+4, :] for a in range(4)), s.zeros(4, 6)))
        Q.append(sum((axes[axis, a]*second[:, 4*a:4*a+4] for a in range(4)), s.zeros(6, 4)))
    u = [q.nullspace()[0] for q in Q]
    assert all(p.rank() == q.rank() == 3 and p*q == s.zeros(4) for p, q in zip(P, Q))
    assert s.Matrix.hstack(*u).det()
    for i, j in combinations(range(4), 2):
        assert P[i]*Q[j] == P[j]*Q[i] and (P[i]*Q[j]).rank() == 2
    for j in range(1, 4):
        lhs, rhs = Q[0]*u[j], Q[j]*u[0]
        row = next(k for k in range(6) if lhs[k])
        u[j] = u[j]*(rhs[row]/lhs[row])
    pairs = [Q[i]*u[j] for i, j in combinations(range(4), 2)]
    assert s.Matrix.hstack(*pairs).det()
    for i, j in combinations(range(4), 2):
        assert Q[i]*u[j] == Q[j]*u[i]
    basis = []
    for omitted in range(4):
        i, j, k = sorted(set(range(4))-{omitted})
        basis.append(P[i]*Q[j]*u[k])
    frame = s.Matrix.hstack(*basis)
    assert frame.det()
    for i, j, k in product(range(4), repeat=3):
        expected = s.Matrix([standard.get(encode((i, j, a, k)), 0) for a in range(4)])
        assert P[i]*Q[j]*u[k] == frame*expected

    # Independently reconstruct the fourth slot from the two primal windows.
    relations = reconstructed_windows(w)[:2]
    cubic = nullspace_sparse(insertion_columns(orthogonal_windows(relations), 0, 3), 64)
    inclusion = s.Matrix(64, len(cubic), lambda i, j: cubic[j].get(i, 0))
    actual = s.Matrix(64, 4, lambda i, j: w.get(4*i+j, 0))
    coordinate_change, parameters = inclusion.gauss_jordan_solve(actual)
    assert parameters.rows == 0 and coordinate_change.det()
    assert inclusion*coordinate_change == actual
    return {"input_factor_changes": [g.tolist() for g in changes],
            "kernel_frame_det": str(s.Matrix.hstack(*u).det()),
            "pair_frame_det": str(s.Matrix.hstack(*pairs).det()),
            "triple_frame_det": str(frame.det()),
            "two_window_fourth_slot_det": str(coordinate_change.det())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    g = s.Matrix([[1,2,0,1],[0,1,3,0],[0,0,1,2],[0,0,0,1]])
    h = s.Matrix([[1,0,0,0],[1,2,0,0],[2,3,1,0],[1,1,1,3]])
    fixtures = [verify([s.eye(4)]*4), verify([g, g, h, g.T*h])]
    result = {"result": "PASS", "arithmetic": "exact rational",
              "checks": "all six pair ranks and normalized exchanges, full tensor, two-window cubic inclusion",
              "fixtures": fixtures}
    if args.output:
        args.output.write_text(json.dumps(result, default=int, indent=2)+"\n")
    print("PASS square-axis normalization and two-window reconstruction in both bases")


if __name__ == "__main__":
    main()
