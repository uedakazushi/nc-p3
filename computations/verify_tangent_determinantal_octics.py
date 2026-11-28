"""Compute the bordered octic of infinitesimal square roots of a rational instanton edge."""

from itertools import combinations, combinations_with_replacement, permutations
from math import prod

import numpy as np
import sympy as sp
from flint import nmod_mat, nmod_mpoly_ctx

from verify_finite_square_tangent_instantons import edge
from verify_normal_configuration_octics import smooth_plane_form
from verify_opposite_root_triple_octics import monomials, quotient_projection
from verify_semisimple_period_obstruction import P
from verify_single_edge_instanton_extension import coefficients_to_field


def tangent_matrices(rel):
    quartics, projection = quotient_projection(rel, 4)
    wedge = list(combinations(range(4), 2))
    pairs = monomials(2)
    rows = []
    for (i, j), (k, l) in combinations_with_replacement(wedge, 2):
        block = np.zeros((20, 100), dtype=np.int64)
        # Evaluate (x_i e_j - x_j e_i)^t H(x) (x_k e_l - x_l e_k).
        for a, b, c, d, sign in [(i, k, j, l, 1), (i, l, j, k, -1),
                                 (j, k, i, l, -1), (j, l, i, k, 1)]:
            index = pairs.index(tuple(sorted((c, d))))
            for t, word in enumerate(pairs):
                block[:, 10*index+t] += sign*projection[:, quartics[tuple(sorted((a, b)+word))]]
        rows.extend((block % P).tolist())
    return nmod_mat(rows, P).nullspace()


def bordered_section(basis, column, coords):
    pairs = monomials(2)
    entries = [[0 for _ in range(5)] for _ in range(5)]
    for index, (i, j) in enumerate(pairs):
        value = sum(int(basis[10*index+t, column])*coords[k]*coords[l] for t, (k, l) in enumerate(pairs))
        entries[i][j] = entries[j][i] = value
    for i in range(4):
        entries[4][i] = entries[i][4] = coords[i]
    return sum((-1)**sum(p[i] > p[j] for i, j in combinations(range(5), 2))
               * prod(entries[i][p[i]] for i in range(5)) for p in permutations(range(5)))


def main():
    rel = sp.Matrix(coefficients_to_field(edge()).tolist())
    basis, size = tangent_matrices(rel)
    assert size == 17  # Sixteen infinitesimal basis changes and one Grassmann tangent.
    ctx = nmod_mpoly_ctx.get(("x", "y", "z"), P)
    x, y, z = ctx.gens()
    coords = (x, y, z, x+2*y+5*z)
    section = next(f for j in range(size) if (f := bordered_section(basis, j, coords)))
    smooth_plane_form(section, P)
    print("Verified one tangent modulo basis changes and a smooth plane section of its bordered octic.")


if __name__ == "__main__":
    main()
