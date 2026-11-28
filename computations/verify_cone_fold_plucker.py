"""Verify descent of the symmetric fold to the binary Hankel and Pluecker systems."""
from itertools import combinations, combinations_with_replacement

import numpy as np
from flint import nmod_mat

from verify_quadric_cone_middle import PRIME, binary_coordinates, cone_edges, field, fixed_tangent


def main():
    pairs = list(combinations(range(4), 2))
    symmetric = list(combinations_with_replacement(range(4), 2))
    equations = []
    for i,j in pairs:
        for k in range(3):
            row = np.zeros(24, dtype=int)
            row[6*i+j+k] = 1
            row[6*j+i+k] = -1
            equations.append(row)
    assert field(np.array(equations)).rank() == 15
    _, edge = cone_edges()
    _, transfer = binary_coordinates(edge)
    brackets = np.zeros((9, 6), dtype=int)
    for k,(i,j) in enumerate(pairs):
        for h in range(6):
            brackets[i+h,k] += int(transfer[h,j])
            brackets[j+h,k] -= int(transfer[h,i])
    assert field(brackets).rank() == 5
    dual, dim = field(brackets).transpose().nullspace()
    assert dim == 4
    lam = np.array([[int(dual[i,j]) for i in range(9)] for j in range(4)], dtype=int)
    cut = []
    for i,j in symmetric:
        row = [lam[:,i+j+k] for k in range(3)]
        row.append(sum(int(transfer[h,j])*lam[:,i+h] for h in range(6)) % PRIME)
        cut.append(np.concatenate(row).tolist())
    assert nmod_mat(cut, PRIME).rank() == 10
    fixed_tangent()
    print("PASS: Hankel conditions leave nine coefficients, sample bracket rank five and opposite cut rank ten; rank-six tangent verified")


if __name__ == "__main__":
    main()
